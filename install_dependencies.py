import subprocess
import sys
import logging
import time
import os
import json
import requests
from pipeline_functions import (
    check_and_download_reference,
    check_and_create_bwa_index
)

# Configure logging with timestamps
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

PARABRICKS_IMAGE = os.getenv('PARABRICKS_IMAGE', 'nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1')

def check_internet_connection():
    """Check for internet connectivity before proceeding with downloads."""
    try:
        requests.get("https://8.8.8.8", timeout=3)
        logging.info("Internet connection confirmed.")
    except requests.RequestException:
        logging.error("No internet connection. Please check your connection and retry.")
        sys.exit(1)

def configure_conda_channels():
    """Configure Conda channels for package installation with error handling."""
    channels = ['defaults', 'bioconda', 'conda-forge']
    try:
        for channel in channels:
            subprocess.run(['conda', 'config', '--add', 'channels', channel], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to configure Conda channel {channel}: {e}")
        sys.exit(1)

def install_package(package_name, channel=None):
    """Install a package using Conda if not already installed with the required version."""
    logging.info(f"Checking if {package_name} is installed.")
    try:
        if "==" in package_name:
            base_package_name, required_version = package_name.split("==")
        else:
            base_package_name = package_name
            required_version = None

        result = subprocess.run(
            ['conda', 'list', base_package_name, '--json'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        installed_packages = json.loads(result.stdout)

        for pkg in installed_packages:
            # Match name and version if specified
            if pkg['name'] == base_package_name and (not required_version or pkg['version'] == required_version):
                logging.info(f"{package_name} is already installed with the required version.")
                return

        logging.info(f"{package_name} not installed with required version. Installing via Conda...")
        install_cmd = ['conda', 'install', '-y']
        if channel:
            install_cmd += ['-c', channel]
        install_cmd.append(package_name)

        start_time = time.time()
        subprocess.run(install_cmd, check=True)
        elapsed_time = time.time() - start_time
        logging.info(f"{package_name} installed in {elapsed_time / 60:.2f} minutes.")

    except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
        logging.error(f"Failed to install {package_name} via Conda: {e}")
        sys.exit(1)

def install_packages_sequentially(packages):
    """Install multiple packages sequentially."""
    failed_packages = []
    for pkg in packages:
        try:
            if isinstance(pkg, tuple):
                install_package(pkg[0], pkg[1])
            else:
                install_package(pkg)
        except Exception as e:
            logging.error(f"Failed to install {pkg}: {e}")
            failed_packages.append(pkg)
    if failed_packages:
        logging.warning(f"The following packages failed to install: {', '.join(failed_packages)}")

def check_docker_installed():
    """Check if Docker is installed and accessible."""
    try:
        subprocess.run(['docker', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info("Docker is installed and accessible.")
    except FileNotFoundError:
        logging.error("Docker is not installed. Please install Docker before proceeding.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error checking Docker installation: {e}")
        sys.exit(1)

def check_nvidia_docker():
    """Check if NVIDIA Docker is set up for GPU usage."""
    try:
        subprocess.run(
            ['docker', 'run', '--rm', '--gpus', 'all', 'nvidia/cuda:12.0.1-devel-ubi8', 'nvidia-smi'], check=True
        )
        logging.info("Docker GPU setup verified successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Docker GPU verification failed: {e}")
        sys.exit(1)

def pull_parabricks_image():
    """Pull the Clara Parabricks Docker image if not already available."""
    logging.info(f"Checking for Clara Parabricks Docker image: {PARABRICKS_IMAGE}")
    try:
        result = subprocess.run(['docker', 'images', '-q', PARABRICKS_IMAGE], stdout=subprocess.PIPE, text=True)
        if result.stdout.strip():
            logging.info(f"Clara Parabricks Docker image {PARABRICKS_IMAGE} already exists.")
            return
        logging.info(f"Pulling Clara Parabricks Docker image: {PARABRICKS_IMAGE}")
        subprocess.run(['docker', 'pull', PARABRICKS_IMAGE], check=True)
        logging.info("Clara Parabricks Docker image pulled successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to pull Clara Parabricks Docker image: {e}")
        sys.exit(1)

def setup_reference_genome(base_dir):
    """Set up reference genome and ensure BWA index is created."""
    reference_genome_dir = os.path.join(base_dir, 'reference_genome')
    reference_genome_file = os.path.join(reference_genome_dir, 'reference.fa')
    os.makedirs(reference_genome_dir, exist_ok=True)

    if not check_and_download_reference(reference_genome_dir, reference_genome_file):
        logging.error("Failed to prepare reference genome.")
        sys.exit(1)

    if not check_and_create_bwa_index(reference_genome_file):
        logging.error("Failed to create BWA index for the reference genome.")
        sys.exit(1)

def main():
    total_start_time = time.time()
    check_internet_connection()  # Ensure connectivity

    configure_conda_channels()

    required_packages = [
        'samtools==1.16.1', 'bwa==0.7.18',
        'vcftools==0.1.16', 'wget==1.21.4', 'unzip==6.0',
        'plink2==2.00a5.12', 'pandas==2.2.3', 'openssl==3.3.2',
        'pigz==2.3.4',
        ('bcftools==1.16', 'bioconda'),
    ]
    install_packages_sequentially(required_packages)  # Sequential package installation

    check_docker_installed()
    pull_parabricks_image()
    check_nvidia_docker()

    base_dir = os.getcwd()
    setup_reference_genome(base_dir)

    total_end_time = time.time()
    total_time_elapsed = total_end_time - total_start_time
    logging.info(f"Setup complete in {total_time_elapsed // 60:.0f} minutes and {total_time_elapsed % 60:.0f} seconds.")

if __name__ == "__main__":
    main()
