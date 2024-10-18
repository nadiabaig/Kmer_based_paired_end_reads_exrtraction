# KmerReadExtractor

**KmerReadExtractor** is a bioinformatics tool designed to extract reads containing specific k-mers from paired-end sequencing data. It leverages a Bloom Filter for efficient k-mer membership testing and supports parallel processing to handle large datasets quickly.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Command-Line Arguments](#command-line-arguments)
- [Dependencies](#dependencies)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Features

- **K-mer based paired end reads extraction**: Efficiently extracts reads containing specified k-mers from paired-end FASTQ files.
- **Bloom Filter**: Utilizes a Bloom Filter for rapid k-mer membership checking, minimizing memory usage.
- **Parallel Processing**: Processes large datasets in parallel, improving performance.
- **Flexible Input**: Supports gzipped and plain FASTQ files.
- **Customizable**: Allows users to specify k-mer sizes and the number of parallel nodes.

## Installation

To run **KmerExtractor**, ensure you have Python installed (version 3.6 or higher). You can install the required libraries using `pip`:

```bash
pip install biopython bloom-filter
```

## Usage

To run **KmerExtractor**, use the following command format:

```bash
python extract_reads.py --kmer_file <kmer_file> --r1_file <reads_R1.fastq.gz> --r2_file <reads_R2.fastq.gz> --output_r1 <output_R1.fastq.gz> --output_r2 <output_R2.fastq.gz> --output_orphans <orphans.fastq.gz>  --kmer_size <kmer_size> --num_workers <num_nodes>
```

### Command-Line Arguments

| Argument       | Description                                                             |
|----------------|-------------------------------------------------------------------------|
| `--kmer_file`  | Path to the k-mer file (text format)                                   |
| `--r1_file`    | Path to the paired-end R1 FASTQ file (gzipped or plain)               |
| `--r2_file`    | Path to the paired-end R2 FASTQ file (gzipped or plain)               |
| `--output_r1`  | Output filtered R1 FASTQ file                                           |
| `--output_r2`  | Output filtered R2 FASTQ file                                           |
| `--output_orphans`  | Output orphan reads                                           |
| `--kmer_size`  | Size of the k-mers (default is 31)                                     |
| `--num_workers`| Number of workers for parallel processing (default is 4)                |

## Dependencies

- **Python 3.6 or higher**
- **Biopython**: For reading and writing FASTQ files.
- **Bloom-Filter**: For efficient k-mer membership checking.

## Examples

Here are some examples of how to run **KmerExtractor**:

### Basic Example

```bash
python extract_reads.py --kmer_file kmers.txt --r1_file reads_R1.fastq.gz --r2_file reads_R2.fastq.gz --output_r1 output_R1.fastq.gz --output_r2 output_R2.fastq.gz --output_orphans <orphans.fastq.gz>  --kmer_size 31 --num_workers 4
```

### Adjusting K-mer Size and Workers

```bash
python extract_reads.py --kmer_file kmers.txt --r1_file reads_R1.fastq.gz --r2_file reads_R2.fastq.gz --output_r1 output_R1.fastq.gz --output_r2 output_R2.fastq.gz --output_orphans <orphans.fastq.gz>  --kmer_size 21 --num_workers 8
```

## Contribution

Contributions are welcome! If you have suggestions for improvements or features, feel free to fork the repository and submit a pull request.

1. Fork the repository.
2. Create your feature branch (`git checkout -b feature/YourFeature`).
3. Commit your changes (`git commit -m 'Add some feature'`).
4. Push to the branch (`git push origin feature/YourFeature`).
5. Open a pull request.


---

Feel free to modify this guideline according to your preferences or specific requirements for the GitHub repository!
