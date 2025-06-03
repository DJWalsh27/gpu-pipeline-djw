import argparse
import logging
import os
from pipeline_functions import (
    download_bam,
    run_bam2fq,
    run_parabricks_fq2bam,
    run_parabricks_haplotypecaller
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

RAW_BAM_DIR = "./raw_bams"
FASTQ_DIR = "./fastq_files"
ALIGNMENT_DIR = "./alignment"
VCF_DIR = "./vcf_files"
METRICS_DIR = "./quality_reports"

# Reference file assumed to be set up by install_dependencies.py
REFERENCE = "./reference_genome/reference.fa"

def main(sample_id, num_threads, num_gpus, single_end=False):
    """Main pipeline execution for a single sample."""

    # 1. Ensure key directories exist
    os.makedirs(RAW_BAM_DIR, exist_ok=True)
    os.makedirs(FASTQ_DIR, exist_ok=True)
    os.makedirs(ALIGNMENT_DIR, exist_ok=True)
    os.makedirs(VCF_DIR, exist_ok=True)
    os.makedirs(METRICS_DIR, exist_ok=True)

    # 2. Download the BAM file
    bam_path = download_bam(sample_id, RAW_BAM_DIR)
    if not bam_path:
        logging.error(f"Failed to download BAM for sample {sample_id}. Exiting.")
        return

    # 3. Convert BAM to FASTQ
    fastq_dir = os.path.join(FASTQ_DIR, sample_id)
    os.makedirs(fastq_dir, exist_ok=True)

    if single_end:
        fastq1, _ = run_bam2fq(bam_path, fastq_dir)
        if not fastq1:
            logging.error(f"Failed to generate single-end FASTQ for {sample_id}. Exiting.")
            return
        fastq2 = None
    else:
        fastq1, fastq2 = run_bam2fq(bam_path, fastq_dir)
        if not fastq1 or not fastq2:
            logging.error(f"Failed to generate paired-end FASTQ for {sample_id}. Exiting.")
            return

    # 4. Define output paths
    output_bam = os.path.join(ALIGNMENT_DIR, f"{sample_id}_aligned.bam")
    output_vcf = os.path.join(VCF_DIR, f"{sample_id}.vcf")

    # 5. Run FASTQ to BAM and variant calling
    if run_parabricks_fq2bam(
        fastq1, 
        fastq2, 
        REFERENCE, 
        output_bam, 
        METRICS_DIR, 
        sample_id, 
        num_gpus, 
        num_threads
    ):
        if run_parabricks_haplotypecaller(output_bam, REFERENCE, output_vcf, num_gpus):
            logging.info(f"VCF generation completed for sample {sample_id}.")

            # 6. Cleanup
            if os.path.exists(bam_path):
                os.remove(bam_path)
                logging.info(f"Deleted raw BAM: {bam_path}")
            if os.path.exists(output_bam):
                os.remove(output_bam)
                logging.info(f"Deleted pipeline-generated BAM: {output_bam}")
        else:
            logging.error(f"VCF generation failed for sample {sample_id}.")
    else:
        logging.error(f"Alignment failed for sample {sample_id}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline Script for Single Sample Processing")
    parser.add_argument('--sample', required=True, help="Sample ID to process")
    parser.add_argument('--threads', required=True, type=int, help="Number of threads for processing")
    parser.add_argument('--gpus', required=True, type=int, help="Number of GPUs to use")
    parser.add_argument('--single-end', action='store_true', help="Perform single-end alignment instead of paired-end")
    args = parser.parse_args()

    main(args.sample, args.threads, args.gpus, single_end=args.single_end)
