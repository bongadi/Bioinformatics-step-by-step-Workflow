# Bioinformatics-Workflow

### Comprehensive Bioinformatics Workflow Explanation

#### 1. Downloading Data from NCBI
- **Downloading SRA Files**:
  - **Tool**: SRA Toolkit
  - **Installation**:
    ```bash
    conda install -c bioconda sra-tools
    ```
  - **Usage**: To download sequencing data from NCBI's Sequence Read Archive (SRA), use the `fasterq-dump` command:
    - Single-end reads:
      ```bash
      fasterq-dump SRR_ID
      ```
    - Paired-end reads:
      ```bash
      fasterq-dump --split-files SRR_ID
      ```
    The `fasterq-dump` tool is faster than `fastq-dump` and outputs reads in the FASTQ format.

- **Downloading non-SRA Files from NCBI**:
  - NCBI hosts many types of sequence data in formats like FASTA or FASTQ that aren’t part of the SRA. These can be directly downloaded using command-line tools such as `wget` or `curl`, or simply using the download links provided by NCBI.
  - Example using `wget`:
    ```bash
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/filename.fasta
    ```

#### 2. Quality Control (QC)
- **Why Perform QC?**: Before proceeding to downstream analysis, it is critical to assess the quality of your sequencing reads. QC helps you identify issues such as low-quality reads, adapter contamination, over-represented sequences, or poor sequencing quality, which can negatively affect assembly or alignment.
- **Tool**: FastQC
- **Installation**:
  ```bash
  conda install -c bioconda fastqc
  ```
- **Usage**:
  - For single or paired-end reads:
    ```bash
    fastqc input.fastq
    ```
  FastQC will generate an HTML report with detailed graphs and statistics on the quality of your reads, allowing you to decide if filtering and trimming are necessary.

#### 3. Filtering and Trimming Reads
- **Why Filter and Trim?**:
  - **Filtering**: Removes low-quality reads or reads that do not meet length requirements. This step is crucial to improve the overall accuracy and efficiency of downstream processes such as assembly or alignment.
  - **Trimming**: Involves removing low-quality bases or adapter sequences that might be present at the beginning or end of your reads. These unwanted regions can introduce errors in assembly and alignment.
  
- **Tool (Illumina Reads)**: fastp (a fast, comprehensive tool for both trimming and filtering)
  - **Installation**:
    ```bash
    conda install -c bioconda fastp
    ```
  - **Usage**:
    - Single-end reads:
      ```bash
      fastp -i input.fastq -o output_trimmed.fastq -q 20 -l 50
      ```
    - Paired-end reads:
      ```bash
      fastp -i input_R1.fastq -I input_R2.fastq -o output_R1_trimmed.fastq -O output_R2_trimmed.fastq -q 20 -l 50
      ```

- **Tool (Long Reads)**: filtlong (for nanopore and PacBio long-read data)
  - **Installation**:
    ```bash
    conda install -c bioconda filtlong
    ```
  - **Usage**:
    ```bash
    filtlong --min_length 1000 input.fastq > filtered_output.fastq
    ```

#### 4. Read Assembly
- **Why Perform Assembly?**: Once reads are cleaned and filtered, genome assembly is performed to reconstruct the entire genome or transcriptome from the sequencing reads. This process is especially crucial for de novo assemblies, where a reference genome is not available.
- **Tool (Illumina Short Reads)**: SPAdes
  - **Installation**:
    ```bash
    conda install -c bioconda spades
    ```
  - **Usage**:
    - Single-end:
      ```bash
      spades.py -s input.fastq -o output_directory
      ```
    - Paired-end:
      ```bash
      spades.py -1 input_R1.fastq -2 input_R2.fastq -o output_directory
      ```

- **Tool (Long Reads)**: Flye (for assembling long reads)
  - **Installation**:
    ```bash
    conda install -c bioconda flye
    ```
  - **Usage**:
    ```bash
    flye --nano-raw input.fastq --out-dir output_directory
    ```

- **Tool (Hybrid Assembly)**: Unicycler (combining short and long reads)
  - **Installation**:
    ```bash
    conda install -c bioconda unicycler
    ```
  - **Usage**:
    ```bash
    unicycler -1 input_R1.fastq -2 input_R2.fastq -l long_reads.fastq -o output_directory
    ```

#### 5. Annotation
- **Why Perform Annotation?**: Genome annotation helps identify the biological features (e.g., genes, regulatory regions) within the assembled genome. Annotation assigns functional information to the genomic sequence, making it useful for further biological interpretation.
- **Tool**: Prokka (bacterial genome annotation)
  - **Installation**:
    ```bash
    conda install -c bioconda prokka
    ```
  - **Usage**:
    ```bash
    prokka --outdir output_dir --prefix annotation genome.fasta
    ```

#### 6. Read Alignment
- **Why Perform Alignment?**: Alignment is the process of mapping sequencing reads back to a reference genome to identify sequence similarities, variations, or to quantify expression levels in RNA-Seq experiments.
- **Tool**: BWA (Burrows-Wheeler Aligner for aligning reads to reference genomes)
  - **Installation**:
    ```bash
    conda install -c bioconda bwa
    ```
  - **Usage**:
    ```bash
    bwa mem reference.fasta input.fastq > output.sam
    ```

#### 7. Phylogenetic Analysis
- **Why Perform Phylogenetic Analysis?**: This step is critical for understanding evolutionary relationships between sequences by constructing a phylogenetic tree. It provides insights into how sequences are related or diverged over time.
- **Tool**: IQ-TREE (phylogenetic tree building)
  - **Installation**:
    ```bash
    conda install -c bioconda iqtree
    ```
  - **Usage**:
    ```bash
    iqtree -s alignment.fasta -m GTR+G -bb 1000 -nt AUTO
    ```

#### 8. Visualization and Interpretation
- **Why Visualize Results?**: Visualization is essential for interpreting the quality of data, the assembly, and the functional elements within the genome. Tools like IGV help visualize alignments and annotations, providing a graphical representation for better understanding.

- **Tool**: IGV (Integrative Genomics Viewer for viewing reads aligned to reference genomes)
  - **Installation**:
    ```bash
    conda install -c bioconda igv
    ```

---

### Summary:
This bioinformatics workflow guides you through downloading data, performing quality control, filtering, trimming, assembling genomes, annotating them, and performing phylogenetic analysis. Each step has been carefully designed to ensure high-quality data for downstream analysis, and tools like FastQC, fastp, SPAdes, Flye, and Prokka help optimize the pipeline. All tools are easily installed using Conda, making the workflow efficient and reproducible.
