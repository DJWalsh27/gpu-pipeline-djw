import os
import subprocess
import logging
import requests

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

BASE_URL = "http://storage.googleapis.com/caendr-site-private-bucket/bam/c_elegans"
TMP_DIR = os.path.join(os.getcwd(), "tmp")
os.makedirs(TMP_DIR, exist_ok=True)

def download_bam(sample_id, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    bam_file = f"{sample_id}.bam"
    bam_path = os.path.join(output_dir, bam_file)
    if os.path.exists(bam_path):
        logging.info(f"{bam_file} already exists in {output_dir}. Skipping download.")
        return bam_path
    logging.info(f"Downloading {bam_file} to {output_dir}...")
    url = f"{BASE_URL}/{bam_file}"
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        temp_path = f"{bam_path}.incomplete"
        with open(temp_path, 'wb') as f_out:
            for chunk in response.iter_content(chunk_size=8192):
                f_out.write(chunk)
        os.rename(temp_path, bam_path)
        logging.info(f"Downloaded {bam_file} successfully.")
        return bam_path
    except Exception as e:
        logging.error(f"Failed to download {bam_file}: {e}")
        return None

def run_bam2fq(input_bam, output_dir, n_threads=2):
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.basename(input_bam).replace('.bam', '')
    fastq1 = os.path.join(output_dir, f"{base_name}_1.fastq.gz")
    fastq2 = os.path.join(output_dir, f"{base_name}_2.fastq.gz")
    uid, gid = os.getuid(), os.getgid()
    cmd = [
        'docker', 'run', '--rm', '--gpus', 'all', '--user', f'{uid}:{gid}',
        '-v', f"{os.path.abspath(os.path.dirname(input_bam))}:/data/input",
        '-v', f"{os.path.abspath(output_dir)}:/data/output",
        '-v', f"{os.path.abspath(TMP_DIR)}:/data/tmp",
        'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1', 'pbrun', 'haplotypecaller',
        '--in-bam', f"/data/input/{os.path.basename(input_bam)}",
        '--out-prefix', f"/data/output/{base_name}",
        '--out-suffixF', '_1.fastq.gz', '--out-suffixF2', '_2.fastq.gz',
        '--num-threads', str(n_threads), '--tmp-dir', '/data/tmp'
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logging.info(f"Converted BAM to FASTQ for {os.path.basename(input_bam)} with {n_threads} threads:\n{result.stdout}")
        return fastq1, fastq2
    except subprocess.CalledProcessError as e:
        logging.error(f"bam2fq failed for {input_bam}: {e}\n{e.stderr}")
        return None, None

def check_and_download_reference(reference_genome_dir, reference_genome_file):
    if not os.path.exists(reference_genome_dir):
        os.makedirs(reference_genome_dir)
    if not os.path.isfile(reference_genome_file):
        logging.info(f"{reference_genome_file} not found. Downloading...")
        ref_url = "https://ftp.ensembl.org/pub/release-113/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
        gz_file = reference_genome_file + '.gz'
        try:
            subprocess.run(['wget', '-O', gz_file, ref_url], check=True)
            subprocess.run(['gunzip', '-f', gz_file], check=True)
            logging.info("Reference genome downloaded and unzipped.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to download or unzip reference genome: {e}")
            return False
    return True

def check_and_create_bwa_index(reference_genome_file):
    index_files = [f"{reference_genome_file}.{ext}" for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']]
    if all(os.path.exists(f) for f in index_files):
        logging.info("BWA index files already exist.")
        return True
    try:
        subprocess.run(['bwa', 'index', reference_genome_file], check=True)
        logging.info("BWA index files created successfully.")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error creating BWA index files: {e}")
        return False

def run_parabricks_fq2bam(fastq1, fastq2, reference, output_bam, metrics_base_dir, sample_id, num_gpus=1, num_threads=8, low_memory=False):
    out_bam_dir = os.path.dirname(output_bam)
    os.makedirs(out_bam_dir, exist_ok=True)
    sample_metrics_dir = os.path.join(metrics_base_dir, sample_id)
    os.makedirs(sample_metrics_dir, exist_ok=True)
    uid, gid = os.getuid(), os.getgid()
    log_file = os.path.join(TMP_DIR, f"{sample_id}_fq2bam.log")
    cmd = [
        'docker', 'run', '--rm', '--gpus', 'all', '--user', f'{uid}:{gid}',
        '-v', f'{os.path.dirname(fastq1)}:/data/fastq',
        '-v', f'{os.path.dirname(reference)}:/data/reference',
        '-v', f'{out_bam_dir}:/data/output',
        '-v', f'{sample_metrics_dir}:/data/metrics',
        '-v', f"{os.path.abspath(TMP_DIR)}:/data/tmp",
        'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1', 'pbrun', 'fq2bam',
        '--ref', f'/data/reference/{os.path.basename(reference)}',
        '--in-fq', f'/data/fastq/{os.path.basename(fastq1)}'
    ]
    if fastq2:
        cmd.append(f'/data/fastq/{os.path.basename(fastq2)}')
    cmd.extend([
        '--out-bam', f'/data/output/{os.path.basename(output_bam)}',
        '--num-gpus', str(num_gpus), '--bwa-cpu-thread-pool', str(num_threads),
        '--out-qc-metrics-dir', '/data/metrics', '--tmp-dir', '/data/tmp',
        '--bwa-options', '-K 1000000',  # Add this line
        '--verbose', '--logfile', f'/data/tmp/{os.path.basename(log_file)}'
    ])
    if low_memory:
        cmd.append('--low-memory')
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(log_file, 'r') as f:
            full_output = result.stdout + f.read()
        logging.info(f"Parabricks fq2bam output for {sample_id}:\n{full_output}")
        return full_output, True
    except subprocess.CalledProcessError as e:
        logging.error(f"Parabricks fq2bam failed for {sample_id}: {e}\n{e.stderr}")
        return e.stderr, False

def run_parabricks_haplotypecaller(input_bam, reference, output_vcf, num_gpus=1):
    out_vcf_dir = os.path.dirname(output_vcf)
    os.makedirs(out_vcf_dir, exist_ok=True)
    uid, gid = os.getuid(), os.getgid()
    log_file = os.path.join(TMP_DIR, f"{os.path.basename(input_bam)}_hc.log")
    cmd = [
        'docker', 'run', '--rm', '--gpus', 'all', '--user', f'{uid}:{gid}',
        '-v', f'{os.path.dirname(input_bam)}:/data/input',
        '-v', f'{os.path.dirname(reference)}:/data/reference',
        '-v', f'{out_vcf_dir}:/data/output',
        '-v', f"{os.path.abspath(TMP_DIR)}:/data/tmp",
        'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1', 'pbrun', 'haplotypecaller',
        '--ref', f'/data/reference/{os.path.basename(reference)}',
        '--in-bam', f'/data/input/{os.path.basename(input_bam)}',
        '--out-variants', f'/data/output/{os.path.basename(output_vcf)}',
        '--num-gpus', str(num_gpus),
        '--tmp-dir', '/data/tmp', '--verbose', 
        '--logfile', f'/data/tmp/{os.path.basename(log_file)}'
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(log_file, 'r') as f:
            full_output = result.stdout + f.read()
        logging.info(f"Parabricks haplotypecaller output:\n{full_output}")
        return full_output, True
    except subprocess.CalledProcessError as e:
        logging.error(f"Parabricks haplotypecaller failed: {e}\n{e.stderr}")
        return e.stderr, False
