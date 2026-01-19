# Demux_Illu_ONT
Script to demultiplexing Illumina NGS libraries from Nanopore sequensing fastq files

*This script is designed to demultiplex FASTQ files (including .gz) by S5 and N7 indices using the base.xlsx and Index.xlsx files*


# **It allows you to:**
* Distribute reads across samples.
* Save the identification scheme in CSV and Excel.
* Build a summary of the number of reads per sample.
* Work with the streaming approach (low memory consumption, batch processing).
* Manage the maximum number of rows exported to Excel.


# **Dependencies**
* Requires Python 3.8+ and the following libraries:
* pip install biopython pandas openpyxl


# **Example run:**
```
python demultiplex_v0.15.py \
-i input.fastq.gz \
-b base.xlsx \
-x Index.xlsx \
-o results \
-t 16 \
--batch 50000 \
--max-excel-rows 1000000
```


# **Command Line Arguments Description Default**
```
-i, --input FASTQ or FASTQ.GZ input file is required
-b, --base base.xlsx file with indexes (Name, Forward, Reverse) is required
-x, --index Index.xlsx file with sample table (Sample, S5, N7) is required
-o, --output Results folder is required
-t, --threads Number of threads 4
--batch Batch size for parallel processing 100000
--max-excel-rows Maximum number of rows to save in Excel 1,000,000
```


# **Results**
After execution, the following will be created in the results/ folder:
```
reads_scheme.csv — the full distribution scheme (unlimited)
reads_scheme.xlsx — the scheme, but only the first --max-excel-rows rows (for quick viewing)
sample_summary.csv — sample summary with read counts
sample_summary.xlsx — the same summary in Excel
SampleName.fastq.gz — individual FASTQs for each sample
unmultiplexed.fastq.gz — reads without identification
demultiplex.log — script execution log
```


# **Operation Logic**
- Loading the base.xlsx and Index.xlsx tables.
- Building S5 and N7 index mapping dictionaries.
- Reading FASTQs as a stream (SeqIO.parse), without loading them entirely into memory.
- Processing reads in batches (--batch) using multiprocessing.


# **For each read:**
- Determine the scheme (Forward/Reverse).
- Search for S5/N7 indices at the beginning and end of the read.
- Check for a match (S5, N7) in Index.xlsx. Write the result to CSV, FASTQ for a sample, or unmultiplexed.


# **After completion:**
- CSV > Excel (row limit).
- The sample_summary summary is generated.

Thus, version 0.15 is stable and scalable for files >10 GB.

