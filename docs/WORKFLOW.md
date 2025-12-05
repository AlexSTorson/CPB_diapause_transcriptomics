# Detailed Analysis Workflow

This document provides step-by-step instructions for reproducing the analyses presented in Torson et al. (202X).

## Table of Contents
1. [Data Acquisition](#data-acquisition)
2. [Read Processing](#read-processing)
3. [Transcriptome Assembly](#transcriptome-assembly)
4. [Annotation](#annotation)
5. [Differential Expression Analysis](#differential-expression-analysis)
6. [Network Analysis](#network-analysis)
7. [Functional Enrichment](#functional-enrichment)

---

## Data Acquisition

### Raw Sequencing Data

Download raw FASTQ files from NCBI SRA:

```bash
# Install SRA Toolkit if needed
# conda install -c bioconda sra-tools

# Download all libraries
for SRX in {SRX6428284..SRX6428305}; do
    prefetch $SRX
    fasterq-dump --split-files --outdir raw_data/ $SRX
done
```

### Reference Genome

```bash
# Download L. decemlineata RefSeq genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/GCF_000500325.1_Ldec_2.0/GCF_000500325.1_Ldec_2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/GCF_000500325.1_Ldec_2.0/GCF_000500325.1_Ldec_2.0_genomic.gff.gz

# Decompress
gunzip *.gz
```

---

## Read Processing

### Quality Trimming with fastp

```bash
# Run fastp on all samples (Phred ≥ 15)
for sample in raw_data/*_1.fastq; do
    base=$(basename $sample _1.fastq)
    fastp -i raw_data/${base}_1.fastq \
          -I raw_data/${base}_2.fastq \
          -o trimmed_data/${base}_1_trimmed.fastq \
          -O trimmed_data/${base}_2_trimmed.fastq \
          --qualified_quality_phred 15 \
          --html reports/${base}_fastp.html \
          --json reports/${base}_fastp.json
done
```

### Quality Assessment

```bash
# Pre-trimming QC
fastqc -o qc_reports/pre_trim/ raw_data/*.fastq

# Post-trimming QC
fastqc -o qc_reports/post_trim/ trimmed_data/*.fastq

# Aggregate reports
multiqc qc_reports/ -o qc_reports/multiqc/
```

---

## Read Mapping and Assembly

### Build HISAT2 Index

```bash
# Build index from reference genome
hisat2-build GCF_000500325.1_Ldec_2.0_genomic.fna ldec_index
```

### Map Reads with HISAT2

```bash
# Map each sample
for sample in trimmed_data/*_1_trimmed.fastq; do
    base=$(basename $sample _1_trimmed.fastq)
    hisat2 -x ldec_index \
           -1 trimmed_data/${base}_1_trimmed.fastq \
           -2 trimmed_data/${base}_2_trimmed.fastq \
           -S alignments/${base}.sam \
           --threads 8
    
    # Convert to BAM and sort
    samtools view -bS alignments/${base}.sam | samtools sort -o alignments/${base}.sorted.bam
    samtools index alignments/${base}.sorted.bam
done
```

---

## Transcriptome Assembly

### Individual Assembly with StringTie

```bash
# Assemble each sample individually
for bam in alignments/*.sorted.bam; do
    base=$(basename $bam .sorted.bam)
    stringtie $bam \
              -G GCF_000500325.1_Ldec_2.0_genomic.gff \
              -o assemblies/${base}.gtf \
              -l ${base}
done
```

### Merge Assemblies with TACO

```bash
# Create assembly list
ls assemblies/*.gtf > assembly_list.txt

# Merge with TACO
taco_run -o merged_assembly/ \
         --gtf-expr-attr TPM \
         assembly_list.txt

# This produces merged_assembly/assembly.gtf
```

### Extract Count Data

```bash
# Use StringTie's prepDE.py script
python prepDE.py -i assembly_list.txt \
                 -g gene_count_matrix.csv \
                 -t transcript_count_matrix.csv
```

---

## Annotation

### BLASTx for Novel Transcripts

```bash
# Extract novel transcript sequences (not in RefSeq)
# This requires custom filtering based on StringTie class codes

# Run BLASTx against NCBI nr (Arthropoda)
blastx -query novel_transcripts.fasta \
       -db nr \
       -taxids 6656 \
       -evalue 1e-3 \
       -outfmt 5 \
       -out blast_results.xml \
       -num_threads 16
```

### Functional Annotation with Blast2GO

Use OmicsBox GUI or command-line tools:
1. Import BLASTx results
2. Map GO terms
3. Annotate sequences
4. Export annotation table

### KEGG Annotation with KAAS

1. Submit novel transcript sequences to KAAS web server
2. Select "bi-directional best hit" method
3. Use Eukaryotic + Cole