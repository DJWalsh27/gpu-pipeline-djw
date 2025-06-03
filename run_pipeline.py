import os
import argparse
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pipeline_functions import download_bam, run_bam2fq

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

# Define key directories
RAW_BAM_DIR = "./raw_bams"
FASTQ_DIR = "./fastq_files"
LOG_DIR = "./logs"

def check_conda_environment(env_name, python_version):
    """Ensure the specified Conda environment exists."""
    try:
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True, check=True)
        if env_name in result.stdout:
            logging.info(f"Conda environment '{env_name}' already exists.")
        else:
            logging.info(f"Creating Conda environment '{env_name}'...")
            subprocess.run(['conda', 'create', '-y', '-n', env_name, f'python={python_version}'], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error while checking or creating Conda environment: {e}")
        exit(1)

def prepare_fastq(sample_id, num_threads, single_end):
    """Download BAM and convert it to FASTQ."""
    os.makedirs(LOG_DIR, exist_ok=True)
    log_file = os.path.join(LOG_DIR, f"{sample_id}.log")

    fastq_dir = os.path.join(FASTQ_DIR, sample_id)
    os.makedirs(fastq_dir, exist_ok=True)

    bam_path = download_bam(sample_id, RAW_BAM_DIR)
    if not bam_path:
        logging.error(f"Failed to download BAM for {sample_id}. Exiting.")
        return

    logging.info(f"Converting BAM to FASTQ for sample {sample_id}. Logs: {log_file}")
    
    # Handle single-end and paired-end correctly
    if single_end:
        fastq1, _ = run_bam2fq(bam_path, fastq_dir, num_threads)
        fastq2 = None
    else:
        fastq1, fastq2 = run_bam2fq(bam_path, fastq_dir, num_threads)

    if not fastq1 or (not single_end and not fastq2):
        logging.error(f"Failed to generate FASTQ files for {sample_id}.")
        return

    logging.info(f"FASTQ generation completed for {sample_id}.")

def process_samples(sample_file, concurrent_samples, num_threads, single_end):
    """Process all samples from the provided file in parallel."""
    with open(sample_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]

    with ThreadPoolExecutor(max_workers=concurrent_samples) as executor:
        future_to_sample = {
            executor.submit(prepare_fastq, sample, num_threads, single_end): sample for sample in samples
        }
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                future.result()
            except Exception as e:
                logging.error(f"Error processing sample {sample}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Prepare FASTQ files for genomic analysis.")
    parser.add_argument('--env', default='bioinfo_env', help="Conda environment name.")
    parser.add_argument('--threads', type=int, default=4, help="Number of CPU threads to use for bam2fq.")
    parser.add_argument('--gpus', type=int, default=1, help="Number of GPUs to use (not used in bam2fq).")
    parser.add_argument('--samples', default='samples.txt', help="File with list of sample IDs.")
    parser.add_argument('--concurrent', type=int, default=2, help="Number of concurrent samples to process.")
    parser.add_argument('--single-end', action='store_true', help="Process samples as single-end reads.")
    args = parser.parse_args()

    for d in [RAW_BAM_DIR, FASTQ_DIR, LOG_DIR]:
        os.makedirs(d, exist_ok=True)

    check_conda_environment(args.env, '3.9')

    # Install dependencies inside the conda environment
    install_cmd = ["conda", "run", "-n", args.env, "python", "install_dependencies.py"]
    try:
        logging.info("Installing dependencies...")
        subprocess.run(install_cmd, check=True)
        logging.info("Dependencies installed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to install dependencies: {e}")
        exit(1)

    process_samples(args.samples, args.concurrent, args.threads, args.single_end)
    logging.info("FASTQ preparation completed.")

if __name__ == "__main__":
    main()
