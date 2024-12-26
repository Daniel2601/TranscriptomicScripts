# Load necessary libraries
library(Rsamtools)
library(ShortRead)
library(R.utils)

# Define paths to your files
fastq_files <- c("sample1.fastq", "sample2.fastq") # Replace with your FASTQ files
reference_genome <- "reference.fasta"              # Replace with your reference genome
output_dir <- "output_bam"                         # Directory to store BAM files

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to check if a command-line tool is available
check_tool <- function(tool) {
  cmd <- paste("which", tool)
  path <- system(cmd, intern = TRUE)
  if (length(path) == 0) {
    stop(paste(tool, "is not installed or not in your PATH. Please install it before proceeding."))
  }
}

# Check for required tools
check_tool("bwa")      # Replace with "bowtie2" if you prefer Bowtie2
check_tool("samtools")

# Index the reference genome (only need to do this once)
index_reference <- function(ref_genome) {
  # Check if BWA index files exist
  required_files <- paste0(ref_genome, c(".bwt", ".pac", ".ann", ".amb", ".sa"))
  if (!all(file.exists(required_files))) {
    message("Indexing the reference genome with BWA...")
    cmd <- paste("bwa index", ref_genome)
    system(cmd)
  } else {
    message("Reference genome already indexed.")
  }
}

index_reference(reference_genome)

# Function to align FASTQ to reference and convert to BAM
fastq_to_bam <- function(fastq, ref_genome, output_dir) {
  # Define output file names
  sam_file <- file.path(output_dir, paste0(basename(fastq), ".sam"))
  bam_file <- file.path(output_dir, paste0(basename(fastq), ".bam"))
  sorted_bam <- file.path(output_dir, paste0(basename(fastq), ".sorted.bam"))
  
  # Align FASTQ to reference using BWA MEM
  message(paste("Aligning", fastq, "to", ref_genome, "using BWA MEM..."))
  align_cmd <- paste("bwa mem", ref_genome, fastq, ">", sam_file)
  system(align_cmd)
  
  # Convert SAM to BAM
  message(paste("Converting SAM to BAM for", sam_file))
  samtools_view_cmd <- paste("samtools view -bS", sam_file, "-o", bam_file)
  system(samtools_view_cmd)
  
  # Sort BAM file
  message(paste("Sorting BAM file", bam_file))
  samtools_sort_cmd <- paste("samtools sort", bam_file, "-o", sorted_bam)
  system(samtools_sort_cmd)
  
  # Index the sorted BAM file
  message(paste("Indexing BAM file", sorted_bam))
  samtools_index_cmd <- paste("samtools index", sorted_bam)
  system(samtools_index_cmd)
  
  # Optionally, remove intermediate SAM and unsorted BAM files
  file.remove(sam_file)
  file.remove(bam_file)
  
  message(paste("BAM file created and sorted:", sorted_bam))
  
  return(sorted_bam)
}

# Process each FASTQ file
bam_files <- lapply(fastq_files, function(fastq) {
  fastq_to_bam(fastq, reference_genome, output_dir)
})

# Optional: Verify BAM files
message("Verifying BAM files...")
for (bam in bam_files) {
  if (file.exists(bam)) {
    message(paste("BAM file exists:", bam))
  } else {
    warning(paste("BAM file not found:", bam))
  }
}

message("Conversion from FASTQ to BAM completed successfully.")
