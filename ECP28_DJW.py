import os
import subprocess
import argparse
import pandas as pd

def check_file_exists(filepath):
    """Check if a required file exists, exit if not found."""
    if not os.path.isfile(filepath):
        print(f"Error: Required file {filepath} not found.")
        exit(1)

def extract_ad_values(vcf_file, output_file):
    """Extract AD values from VCF using vcftools."""
    cmd = [
        "vcftools", 
        "--gzvcf", vcf_file,
        "--extract-FORMAT-info", "AD",
        "--out", output_file
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Allelic Depth values extracted successfully to {output_file}.AD.FORMAT")
    except subprocess.CalledProcessError as e:
        print(f"Error running vcftools: {e}")
        exit(1)

def format_ad_output(raw_file, formatted_output):
    """Format the extracted AD file into a clean tab-separated output."""
    try:
        df = pd.read_csv(raw_file, sep='\t', header=0)
        df.columns = [col.strip() for col in df.columns]  # Clean column names
        
        # Save as a formatted TSV
        df.to_csv(formatted_output, sep='\t', index=False)
        print(f"Formatted AD values saved to {formatted_output}")
    except Exception as e:
        print(f"Error processing AD values: {e}")
        exit(1)

def process_samples(sample_file, vcf_dir, output_dir):
    """Process AD extraction for all samples listed in the sample file."""
    check_file_exists(sample_file)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    ad_reports_dir = os.path.join(output_dir, "AD_reports")
    if not os.path.isdir(ad_reports_dir):
        os.makedirs(ad_reports_dir, exist_ok=True)
    
    with open(sample_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    
    for sample in samples:
        vcf_file = os.path.join(vcf_dir, f"{sample}.vcf.gz")
        output_file = os.path.join(ad_reports_dir, f"{sample}_AD_values")
        
        if os.path.isfile(vcf_file):
            extract_ad_values(vcf_file, output_file)
            format_ad_output(f"{output_file}.AD.FORMAT", f"{output_file}_formatted.tsv")
        else:
            print(f"Warning: VCF file for {sample} not found. Skipping.")

def main():
    parser = argparse.ArgumentParser(description="Extract Allelic Depth (AD) values from VCF for multiple samples.")
    parser.add_argument("--samples", required=True, help="File listing sample names")
    parser.add_argument("--vcf_dir", required=True, help="Directory containing VCF files")
    parser.add_argument("--output_dir", required=True, help="Directory to save extracted AD values")
    args = parser.parse_args()
    
    # Ensure directories exist
    if not os.path.isdir(args.vcf_dir):
        print(f"Error: VCF directory {args.vcf_dir} not found.")
        exit(1)
    
    # Process all samples
    process_samples(args.samples, args.vcf_dir, args.output_dir)
    
if __name__ == "__main__":
    main()
