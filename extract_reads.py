import argparse
import gzip
from Bio import SeqIO
from bloom_filter import BloomFilter
import multiprocessing as mp

# Function to load k-mers into a Bloom Filter
def load_kmers(kmer_file, expected_size=2304212, error_rate=0.001):
    bloom = BloomFilter(max_elements=expected_size, error_rate=error_rate)

    with open(kmer_file, 'r') as f:
        for line in f:
            kmer, _ = line.strip().split()  # Ignore frequency, just load kmer
            bloom.add(kmer)
    return bloom

# Function to check if a read contains any k-mer from the bloom filter
def contains_kmer(sequence, bloom_filter, kmer_size):
    for i in range(len(sequence) - kmer_size + 1):
        kmer = sequence[i:i+kmer_size]
        if kmer in bloom_filter:
            return True
    return False

# Function to process paired-end reads and extract matching ones
def process_reads(bloom_filter, kmer_size, r1_file, r2_file, output_r1, output_r2, output_orphans, start, end):
    with gzip.open(r1_file, 'rt') as r1_handle, gzip.open(r2_file, 'rt') as r2_handle, \
         gzip.open(output_r1, 'wt') as out_r1, gzip.open(output_r2, 'wt') as out_r2, \
         gzip.open(output_orphans, 'wt') as out_orphans:

        r1_records = SeqIO.parse(r1_handle, "fastq")
        r2_records = SeqIO.parse(r2_handle, "fastq")

        for i, (r1_record, r2_record) in enumerate(zip(r1_records, r2_records)):
            if i < start:
                continue
            if i >= end:
                break

            r1_contains_kmer = contains_kmer(str(r1_record.seq), bloom_filter, kmer_size)
            r2_contains_kmer = contains_kmer(str(r2_record.seq), bloom_filter, kmer_size)

            # If both reads contain k-mers, write them to output files
            if r1_contains_kmer and r2_contains_kmer:
                SeqIO.write(r1_record, out_r1, "fastq")
                SeqIO.write(r2_record, out_r2, "fastq")
            # If only one read contains a k-mer, write it to the orphan file
            elif r1_contains_kmer or r2_contains_kmer:
                if r1_contains_kmer:
                    SeqIO.write(r1_record, out_orphans, "fastq")
                if r2_contains_kmer:
                    SeqIO.write(r2_record, out_orphans, "fastq")

# Parallel processing function
def process_in_parallel(bloom_filter, kmer_size, r1_file, r2_file, output_r1, output_r2, output_orphans, num_workers=4):
    total_reads = sum(1 for _ in SeqIO.parse(gzip.open(r1_file, 'rt'), "fastq"))

    # Determine chunk size based on the total reads and number of workers
    chunk_size = total_reads // num_workers

    # Prepare arguments for multiprocessing
    jobs = []
    for i in range(num_workers):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < num_workers - 1 else total_reads
        jobs.append((bloom_filter, kmer_size, r1_file, r2_file, output_r1, output_r2, output_orphans, start, end))

    # Use multiprocessing to parallelize the work
    with mp.Pool(num_workers) as pool:
        pool.starmap(process_reads, jobs)

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Extract reads containing k-mers from paired-end FASTQ files.")
    parser.add_argument("--kmer_file", required=True, help="Path to the k-mer file (text format)")
    parser.add_argument("--r1_file", required=True, help="Path to the paired-end R1 FASTQ file (gzipped or plain)")
    parser.add_argument("--r2_file", required=True, help="Path to the paired-end R2 FASTQ file (gzipped or plain)")
    parser.add_argument("--output_r1", required=True, help="Output filtered R1 FASTQ file")
    parser.add_argument("--output_r2", required=True, help="Output filtered R2 FASTQ file")
    parser.add_argument("--output_orphans", required=True, help="Output orphan reads FASTQ file")
    parser.add_argument("--kmer_size", type=int, default=31, help="Size of the k-mers (default is 31)")
    parser.add_argument("--num_workers", type=int, default=4, help="Number of workers for parallel processing")
    return parser.parse_args()

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # Load k-mers into the Bloom Filter
    bloom_filter = load_kmers(args.kmer_file)

    # Process reads and extract matching pairs using parallel processing
    process_in_parallel(bloom_filter, args.kmer_size, args.r1_file, args.r2_file, args.output_r1, args.output_r2, args.output_orphans, args.num_workers)
