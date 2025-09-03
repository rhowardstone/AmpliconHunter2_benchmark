#!/bin/bash

# Exit on error
set -e

# Clone the benchmarking repository
echo "Cloning AmpliconHunter2_benchmark repository..."
git clone https://github.com/rhowardstone/AmpliconHunter2_benchmark.git
cd AmpliconHunter2_benchmark

# Install AmpliconHunterv1.1
cd src/
git clone https://github.com/rhowardstone/AmpliconHunter.git
cd AmpliconHunter
pip install -e .
cd ..

# Download the file from Google Drive
echo "Downloading AmpliconHunter2_test-genomes.tar from Google Drive..."

# Check if gdown is available
if command -v gdown &> /dev/null; then
    echo "gdown found, downloading automatically..."
    gdown --id 1Nt7MjwfL3pIa5Axa3I2z3xoO2Ait_ito -O AmpliconHunter2_test-genomes.tar
else
    echo "gdown not found. Please either:"
    echo "  1. Install gdown with: pip install gdown"
    echo "  2. Download manually from: https://drive.google.com/file/d/1Nt7MjwfL3pIa5Axa3I2z3xoO2Ait_ito/view?usp=drive_link"
    echo ""
    echo "Press Enter when AmpliconHunter2_test-genomes.tar is in the current directory..."
    read
fi

# Check if file exists
if [ ! -f "AmpliconHunter2_test-genomes.tar" ]; then
    echo "Error: AmpliconHunter2_test-genomes.tar not found!"
    exit 1
fi

# Extract the main tar file
echo "Extracting AmpliconHunter2_test-genomes.tar..."
tar -xf AmpliconHunter2_test-genomes.tar

# Create temporary directory for all extracted genomes
mkdir -p extracted_genomes

# Find and extract all tar.xz files
echo "Extracting all tar.xz files..."
for archive in *.tar.xz; do
    if [ -f "$archive" ]; then
        echo "Extracting $archive..."
        tar -xJf "$archive" -C extracted_genomes/
    fi
done

# Create the target directory structure
echo "Creating directory structure..."
mkdir -p data/fasta/G204800
mkdir -p data/fasta/G102400
mkdir -p data/fasta/G051200
mkdir -p data/fasta/G025600
mkdir -p data/fasta/G012800
mkdir -p data/fasta/G006400

# Move exactly 204,800 genome files to G204800
echo "Moving 204,800 genome files to G204800..."
find extracted_genomes -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) | head -n 204800 | while read file; do
    mv "$file" data/fasta/G204800/
done

# Get list of files in G204800 for progressive halving
cd data/fasta/G204800
files=(*)
total_files=${#files[@]}

# Copy first 102,400 files to G102400 (half of 204,800)
echo "Copying 102,400 files to G102400..."
for ((i=0; i<102400 && i<$total_files; i++)); do
    cp "${files[$i]}" ../G102400/
done

# Copy first 51,200 files to G051200 (half of 102,400)
echo "Copying 51,200 files to G051200..."
for ((i=0; i<51200 && i<$total_files; i++)); do
    cp "${files[$i]}" ../G051200/
done

# Copy first 25,600 files to G025600 (half of 51,200)
echo "Copying 25,600 files to G025600..."
for ((i=0; i<25600 && i<$total_files; i++)); do
    cp "${files[$i]}" ../G025600/
done

# Copy first 12,800 files to G012800 (half of 25,600)
echo "Copying 12,800 files to G012800..."
for ((i=0; i<12800 && i<$total_files; i++)); do
    cp "${files[$i]}" ../G012800/
done

# Copy first 6,400 files to G006400 (half of 12,800)
echo "Copying 6,400 files to G006400..."
for ((i=0; i<6400 && i<$total_files; i++)); do
    cp "${files[$i]}" ../G006400/
done

# Return to base directory
cd ../../../

# Cleanup - commented out as requested
# echo "Cleaning up temporary files and extra genomes..."
# rm -rf extracted_genomes/
# rm -f *.tar.xz
# rm -f AmpliconHunter2_test-genomes.tar

# Print summary
echo "Setup complete!"
echo "Directory structure created:"
echo "  data/fasta/G204800/ - 204,800 files"
echo "  data/fasta/G102400/ - 102,400 files"
echo "  data/fasta/G051200/ - 51,200 files"
echo "  data/fasta/G025600/ - 25,600 files"
echo "  data/fasta/G012800/ - 12,800 files"
echo "  data/fasta/G006400/ - 6,400 files"

# Verify counts
echo ""
echo "Verifying file counts:"
for dir in data/fasta/G*; do
    count=$(find "$dir" -type f | wc -l)
    echo "  $dir: $count files"
done

echo ""
echo "If there are no errors above, you may now run:"
echo "  bash run_benchmarks.sh results"
echo ""
