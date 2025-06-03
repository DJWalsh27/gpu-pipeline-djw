import os
import argparse
import logging
import subprocess
import time
import csv
import psutil
import shutil
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pipeline_functions import run_parabricks_fq2bam, run_parabricks_haplotypecaller

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

FASTQ_DIR = "./fastq_files"
ALIGNMENT_DIR = "./alignment"
VCF_DIR = "./vcf_files"
LOG_DIR = "./logs"
QUALITY_DIR = "./quality_reports"
REFERENCE = "./reference_genome/reference.fa"

def check_conda_env():
    current_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    if current_env != "bioinfo_env":
        logging.error("This script must be run in the 'bioinfo_env' conda environment.")
        exit(1)

def count_reads(fastq_file):
    try:
        cmd = ["gunzip", "-c", fastq_file] if fastq_file.endswith(".gz") else ["cat", fastq_file]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        num_lines = sum(1 for _ in proc.stdout)
        proc.stdout.close()
        proc.wait()
        return num_lines // 4
    except Exception as e:
        logging.error(f"Error counting reads in {fastq_file}: {e}")
        return None

def count_aligned_reads(bam_file):
    try:
        result = subprocess.run(
            ['conda', 'run', '-n', 'bioinfo_env', 'samtools', 'flagstat', bam_file],
            capture_output=True, text=True, check=True
        )
        total_reads = properly_paired = None
        for line in result.stdout.splitlines():
            if "in total" in line:
                total_reads = int(line.split()[0])
            if "properly paired" in line:
                properly_paired = int(line.split()[0])
        if total_reads and properly_paired:
            return total_reads, properly_paired
        logging.error(f"Could not parse flagstat output for {bam_file}")
        return None, None
    except Exception as e:
        logging.error(f"Error reading alignment stats from {bam_file}: {e}")
        return None, None

def count_variants(vcf_file):
    try:
        result = subprocess.run(
            ['conda', 'run', '-n', 'bioinfo_env', 'bcftools', 'stats', vcf_file],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            if "number of records:" in line:
                return int(line.split("number of records:")[1].strip().split()[0])
        return None
    except Exception as e:
        logging.error(f"Error counting variants in {vcf_file}: {e}")
        return None

def get_cpu_ram_usage():
    cpu_percent = psutil.cpu_percent(interval=1)
    ram_usage = psutil.virtual_memory().used / (1024 * 1024)
    return cpu_percent, ram_usage

def get_gpu_usage_monitor(duration):
    peak_util = peak_mem = 0
    start_time = time.time()
    while time.time() - start_time < duration:
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, check=True
            )
            for line in result.stdout.strip().splitlines():
                util, mem = map(int, line.split(','))
                peak_util = max(peak_util, util)
                peak_mem = max(peak_mem, mem)
            time.sleep(1)
        except Exception as e:
            logging.error(f"Error monitoring GPU usage: {e}")
            return None, None
    return peak_util, peak_mem

def get_disk_io(sample_id):
    if shutil.which("iostat") is None:
        logging.error("iostat not found. Install sysstat for disk I/O stats.")
        return None, None
    try:
        result = subprocess.run(["iostat", "-d", "-k", "1", "30"], capture_output=True, text=True, check=True)
        lines = result.stdout.strip().splitlines()
        peak_read = peak_write = 0
        for line in lines[3:]:
            if line.strip() and not line.startswith("Device"):
                parts = line.split()
                read_speed = float(parts[2]) / 1024
                write_speed = float(parts[3]) / 1024
                peak_read = max(peak_read, read_speed)
                peak_write = max(peak_write, write_speed)
        return peak_read, peak_write
    except Exception as e:
        logging.error(f"Error getting disk I/O stats for {sample_id}: {e}")
        return None, None

def parse_parabricks_output(output_str):
    metrics = {}
    for key, pattern in [
        ("bwa_time", r"bwalib run finished in ([\d.]+) seconds"),
        ("sorting_time", r"Sorting and Marking: ([\d.]+) seconds"),
        ("bqsr_time", r"BQSR and writing final BAM:\s*([\d.]+) seconds"),
        ("time_reading", r"Time spent reading:\s*([\d.]+)\s*seconds"),
        ("monitor_time", r"Time spent monitoring.*?:\s*([\d.]+)"),
        ("regions_processed", r"Regions-Processed\s+\d+\s+\d+\s+(\d+)\s+\d+"),
        ("regions_per_minute", r"Regions/Minute\s+(\d+)")
    ]:
        m = re.search(pattern, output_str)
        if m:
            metrics[key] = float(m.group(1))
    rate_block = re.search(r"Rate stats.*?:\s*min rate:\s*([\d.]+).*?max rate:\s*([\d.]+).*?avg rate:\s*([\d.]+)", output_str, re.DOTALL)
    if rate_block:
        metrics["min_rate"] = float(rate_block.group(1))
        metrics["max_rate"] = float(rate_block.group(2))
        metrics["avg_rate"] = float(rate_block.group(3))
    if "regions_processed" not in metrics:
        last_progress = re.findall(r"Regions-Processed\s+\d+\s+\d+\s+(\d+)\s+\d+", output_str)
        if last_progress:
            metrics["regions_processed"] = float(last_progress[-1])
    if "regions_per_minute" not in metrics:
        last_rate = re.findall(r"Regions/Minute\s+(\d+)", output_str)
        if last_rate:
            metrics["regions_per_minute"] = float(last_rate[-1])
    return metrics

def parse_quality_reports(sample_id):
    quality_dir = os.path.join(QUALITY_DIR, sample_id)
    metrics = {}
    try:
        with open(os.path.join(quality_dir, "insert_size.txt"), 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "MEAN_INSERT_SIZE" in line and i + 1 < len(lines):
                    metrics["mean_insert_size"] = float(lines[i + 1].split()[0])

        with open(os.path.join(quality_dir, "gcbias_summary.txt"), 'r') as f:
            for line in f:
                if line.startswith("All Reads"):
                    parts = line.split()
                    metrics["gc_bias"] = float(parts[8])  # GC_NC_40_59

        with open(os.path.join(quality_dir, "mean_quality_by_cycle.txt"), 'r') as f:
            qualities = []
            header_passed = False
            for line in f:
                if line.strip() and not line.startswith("#"):
                    if not header_passed:
                        header_passed = True
                        continue
                    parts = line.split()
                    if len(parts) > 1 and parts[1].replace('.', '').isdigit():
                        qualities.append(float(parts[1]))
            metrics["mean_quality"] = sum(qualities) / len(qualities) if qualities else 0
    except Exception as e:
        logging.error(f"Error parsing quality reports for {sample_id}: {e}")
    return metrics

def process_sample(sample_id, num_threads, num_gpus, low_memory):
    check_conda_env()
    os.makedirs(LOG_DIR, exist_ok=True)
    fastq1 = os.path.join(FASTQ_DIR, sample_id, f"{sample_id}_1.fastq.gz")
    fastq2 = os.path.join(FASTQ_DIR, sample_id, f"{sample_id}_2.fastq.gz")
    if not os.path.exists(fastq1) or not os.path.exists(fastq2):
        logging.error(f"FASTQ files missing for {sample_id}. Skipping.")
        return
    output_bam = os.path.join(ALIGNMENT_DIR, f"{sample_id}_aligned.bam")
    output_vcf = os.path.join(VCF_DIR, f"{sample_id}.vcf")
    logging.info(f"Processing {sample_id} with {num_threads} threads, {num_gpus} GPUs, low_memory={low_memory}.")

    total_reads = count_reads(fastq1)
    if total_reads is None:
        return
    bases_processed = total_reads * 150

    sample_tmp_dir = os.path.join("tmp", sample_id)
    os.makedirs(sample_tmp_dir, exist_ok=True)
    import pipeline_functions
    pipeline_functions.TMP_DIR = os.path.abspath(sample_tmp_dir)

    with ThreadPoolExecutor(max_workers=1) as executor:
        gpu_future = executor.submit(get_gpu_usage_monitor, 180)
        start_time_fq2bam = time.time()
        fq2bam_output, success = run_parabricks_fq2bam(
            fastq1, fastq2, REFERENCE, output_bam, QUALITY_DIR, sample_id, num_gpus, num_threads, low_memory
        )
        fq2bam_time = time.time() - start_time_fq2bam
        if not success:
            return
        fq2bam_metrics = parse_parabricks_output(fq2bam_output)

        start_time_hc = time.time()
        hc_output, success_hc = run_parabricks_haplotypecaller(output_bam, REFERENCE, output_vcf, num_gpus)
        haplotypecaller_time = time.time() - start_time_hc
        if not success_hc:
            return
        hc_metrics = parse_parabricks_output(hc_output)

    total_reads_aligned, properly_paired_aligned = count_aligned_reads(output_bam)
    if total_reads_aligned is None or properly_paired_aligned is None:
        return
    alignment_rate = (properly_paired_aligned / total_reads_aligned) * 100 if total_reads_aligned > 0 else 0
    num_variants = count_variants(output_vcf)
    cpu_usage, ram_usage = get_cpu_ram_usage()
    gpu_util, gpu_mem = gpu_future.result() if gpu_future else (0, 0)
    disk_read_speed, disk_write_speed = get_disk_io(sample_id)
    quality_metrics = parse_quality_reports(sample_id)

    temp_metrics_file = os.path.join(LOG_DIR, f"{sample_id}_metrics.csv")
    with open(temp_metrics_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = [
            "Sample", "Reads Processed", "Bases Processed", "fq2bam Time (sec)", "HaplotypeCaller Time (sec)",
            "BWA Time (sec)", "Sorting Time (sec)", "BQSR Time (sec)", "Alignment Rate (%)", "Variants Called",
            "Regions Processed", "Regions/Minute", "CPU Usage (%)", "RAM Usage (MB)", "Peak GPU Usage (%)",
            "Peak GPU Memory (MB)", "Disk Read Speed (MB/s)", "Disk Write Speed (MB/s)", "Time Spent Reading (sec)",
            "Min Rate (bases/GPU/min)", "Max Rate (bases/GPU/min)", "Avg Rate (bases/GPU/min)", "Monitor Time (sec)",
            "Mean Insert Size", "GC Bias", "Mean Quality Score"
        ]
        writer.writerow(header)
        row = [
            sample_id, total_reads, bases_processed, fq2bam_time, haplotypecaller_time,
            fq2bam_metrics.get("bwa_time", "NA"), fq2bam_metrics.get("sorting_time", "NA"),
            fq2bam_metrics.get("bqsr_time", "NA"), alignment_rate, num_variants,
            hc_metrics.get("regions_processed", "NA"), hc_metrics.get("regions_per_minute", "NA"),
            cpu_usage, ram_usage, gpu_util, gpu_mem, disk_read_speed, disk_write_speed,
            fq2bam_metrics.get("time_reading", "NA"), fq2bam_metrics.get("min_rate", "NA"),
            fq2bam_metrics.get("max_rate", "NA"), fq2bam_metrics.get("avg_rate", "NA"),
            fq2bam_metrics.get("monitor_time", "NA"), quality_metrics.get("mean_insert_size", "NA"),
            quality_metrics.get("gc_bias", "NA"), quality_metrics.get("mean_quality", "NA")
        ]
        writer.writerow(row)

def process_samples(sample_file, concurrent, threads, gpus, low_memory):
    with open(sample_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    with ThreadPoolExecutor(max_workers=concurrent) as executor:
        futures = {executor.submit(process_sample, sample, threads, gpus, low_memory): sample for sample in samples}
        for future in as_completed(futures):
            sample = futures[future]
            try:
                future.result()
            except Exception as e:
                logging.error(f"Error processing sample {sample}: {e}")

def merge_metrics(threads, gpus, concurrent):
    merged_filename = f"processing_metrics_t{threads}g{gpus}c{concurrent}.csv"
    merged_rows = []
    header_written = False
    for filename in os.listdir(LOG_DIR):
        if filename.endswith("_metrics.csv"):
            filepath = os.path.join(LOG_DIR, filename)
            with open(filepath, 'r') as csvfile:
                reader = csv.reader(csvfile)
                rows = list(reader)
                if rows:
                    if not header_written:
                        merged_rows.append(rows[0])
                        header_written = True
                    merged_rows.extend(rows[1:])
    with open(merged_filename, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerows(merged_rows)
    logging.info(f"Merged metrics saved to {merged_filename}")
    return merged_filename

def write_realclock_metric(merged_filename, realclock_time):
    base, ext = os.path.splitext(merged_filename)
    realclock_filename = f"{base}_realclock{ext}"
    with open(realclock_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Real Clock Duration (sec)"])
        writer.writerow([realclock_time])
    logging.info(f"Real clock metric saved to {realclock_filename}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--gpus', type=int, default=1)
    parser.add_argument('--samples', default='samples.txt')
    parser.add_argument('--concurrent', type=int, default=4)
    parser.add_argument('--low-memory', action='store_true', help="Enable low-memory mode for fq2bam")
    args = parser.parse_args()

    start_time = time.time()
    process_samples(args.samples, args.concurrent, args.threads, args.gpus, args.low_memory)
    merged_filename = merge_metrics(args.threads, args.gpus, args.concurrent)
    realclock_time = time.time() - start_time
    write_realclock_metric(merged_filename, realclock_time)
    logging.info("Processing, merging, and real clock measurement completed.")

if __name__ == "__main__":
    main()
