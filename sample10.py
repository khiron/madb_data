import os
import random
import shutil

# Define the source and destination directories
source_dir = '/home/richard/source/ensembl/primates10_113/output'
destination_dir = '/home/richard/source/ensembl/primates10_113/output/sample10'

# Create the destination directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# Set the random seed for reproducibility
random.seed(42)

# Get a list of all FASTA files in the source directory
fasta_files = [f for f in os.listdir(source_dir) if f.endswith('.fa')]

# Randomly sample 10 FASTA files
sampled_files = random.sample(fasta_files, 10)

# Move the sampled files to the destination directory
for file in sampled_files:
    shutil.move(os.path.join(source_dir, file), os.path.join(destination_dir, file))

print(f"Moved {len(sampled_files)} files to {destination_dir}")