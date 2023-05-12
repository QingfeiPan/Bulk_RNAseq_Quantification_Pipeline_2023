# Bulk RNA-seq Quantification Pipeline 2023

## Overview

Bulk RNA sequencing (RNA-Seq) is a highly sensitive and accurate tool for meansuring expression across the transcriptome. In addition to the transcriptome quantification, RNA-Seq also allows researchers to detect new splicing junctions (e.g. TOPHAP/TOPHAP2-regtools), novel transcripts (e.g. Cufflinks), gene fusion (e.g. STAR-Fusion, Arriba), single nucleotide variants (e.g. STAR-GATK), and other features. **This pipeline is for transcriptome quantification purpose only.**

The current bulk RNA-Seq quantification methods can be grouped into two categories (as summarized in the table below). To ensure the accuracy of quantification, we employ one signature method from each of these two categories for cross-validation: **1)** **Bowtie2/STAR-RSEM**, the most highly cited alignment-based method which shows the highest accuracy in most benchmarks; **2)** **Salmon**, one wicked-fast and highly-accurate alignment-free method which is recently further enhanced by integrating selective alignment and decoy sequences.

In addition, we also introduce **3) STAR-HTSeq**, another alignment-based method recomended by GDC. 

|            | Alignment-based methods                                      | Alignment-free methods                                       |
| ---------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Definition | Methods that quanfity from the **alignments to either transcritpome or genome** (i.e. BAM/SAM files). The coordinates of mapped reads are provided. | Methdos that quantify from **read-transcript matches**. The coordinates of mapped reads are **NOT** provided. |
| Principle  | "Seed and extend" for aligners.<br />Quantifiers vary in rules and weighting methods to count reads. | "*k-mer*s" based indexing;<br />Multiple models to identify read-transcript matches,<br /> e.g. SEMEs for  Salmon, T-DBG for Kallisto |
| Examples   | **Aligner**: Bowtie/Bowtie2, STAR, BWA, HISAT2, TopHat2, et. al.<br />**Quantifier**: RSEM, HTSeq, featureCounts, IsoEM, Cufflinks et. al. | Salmon, Kallisto, Sailfish, Fleximer, RNA-Skim, RapMap, et. al. |
| Accuracy   | High                                                         | a little bit lower or equal                                  |
| Speed      | Slow, a few hours for a typical run                          | Super-fast, a few minutes for a type run                     |



Here is an overview of the pipelines:



## Preprocessing

The RNA-Seq data you start with could **very in format** (e.g., FASTQ, BAM, FASTA) and usually **contain noisy sequences** (e.g., adapters leftovers, poor quality bases and other contaminations). So, we need to pre-process these data to generate the standard-in-format, clean-in-sequence FASTQ files which can be directly proceed to quantification analysis.

### 1. Prepare raw FASTQ files

Usually, you have two FASTQ files (R1, R2) for each paired-end sequencing sample, or one single FASTQ file for each single-end sequencing sample. If so, you are good to move forward.

However, if this is not the case, you will need to generate the raw FASTQ files by yourself:

* If you start with multiple FASTQ files for each mate, usually generated in different lanes, you need to merge them into one:

  ```bash
  ## For single-end sequencing
  cat sample1_L001.fq.gz sample1_L002.fq.gz sample1_L003.fq.gz sample1_L004.fq.gz > sample1.raw.fq.gz
  
  ## For paired-end sequencing
  cat sample1_L001_R1.fq.gz sample1_L002_R1.fq.gz sample1_L003_R1.fq.gz sample1_L004_R1.fq.gz > sample1_R1.raw.fq.gz
  cat sample1_L001_R2.fq.gz sample1_L002_R2.fq.gz sample1_L003_R2.fq.gz sample1_L004_R2.fq.gz > sample1_R2.raw.fq.gz
  ```

* If you start with BAM files, usually collected from other sources, you need to convert them into FASTQ files:

  ```bash
  bedtools=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/bedtools
  
  ## For single-end sequencing
  $bedtools bamtofastq -i sample1.bam -fq sample1.raw.fq
  
  ## For paired-end sequencing
  $bedtools bamtofastq -i sample1.bam -fq sample1_R1.raw.fq -fq2 sample1_R2.raw.fq
  ```

  > NOTE:
  >
  > The BAM files of paired-end sequencing, **MUST BE SORTED BY NAME**. To sort the BAM files, please use:
  >
  > ```bash
  > samtools=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/samtools
  > $samtools sort -n -o sample1.sorted sample1.bam
  > ```

### 2. Quality control

Quality control of raw FASTQ files provides not only the metrics to evaluate the quality of FASTQ files, but also key informations used in subsequent analysis, like Encoding of Phred Scores, Sequencing Length and Adapter Content.

```bash
fastqc=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/fastqc

## For single-end sequencing
$fastqc sample1.raw.fq.gz -o output_dir

## For paired-end sequencing
$fastqc sample_R1.raw.fq.gz -o output_dir
$fastqc sample_R2.raw.fq.gz -o output_dir
```

This analysis will generate a .html file with all metrics integerated. For more details of this report, please check out here: https://rtsf.natsci.msu.edu/sites/_rtsf/assets/File/FastQC_TutorialAndFAQ_080717.pdf

 From this report, you can figure out:

* **Encoding of Phred Scores**: Phred+33 is denoted as Illumina 1.9/Sanger, while Phred+64 encoding as illumina 1.5 or lower. (Find more details here: https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)
* **Sequenc Length**: The most common values are 46/45, 76/75, 101/100 or 151/150.
* **Adapter Content**: Illumina Universal Adapter (AGATCGGAAGAG), Illumina Small RNA 3' Adapter (TGGAATTCTCGG), Illumina Small RNA 5' Adapter (GATCGTCGGACT), Nextera Transposase Sequence (CTGTCTCTTATA) and SOLID Small RNA Adapter(CGCCTTGGCCGT).

### 3. Adapter Trimming

Adapter trimming analysis trims not only the **adapter sequences**, but also the **sequences of unknown or low-quality bases**. It also discards the reads of **too-short length**. So, even though no significant adapter content was found in quality control analysis, it is still highly recommended to perform this analysis to  remove the low-quality reads.

In this pipeline, Cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html) is used for adapter trimming. Another round of quality control by FastQC is recommended to confirmed the adapter triming works well.

```bash
cutadapter=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/cutadapt

## For single-end sequencing
$cutadapt -a AGATCGGAAGAG --trim-n --max-n=0.5 --quality-base=33 -q 30 -m 30 --cores=8 -o sample1.clean.fq.gz sample1.raw.fq.gz
$fastqc sample1.clean.fq.gz -o output_dir

## For paired-end sequencing
$cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --trim-n --max-n=0.5 --quality-base=33 -q 30,20 -m 30 --cores=8 -o sample1_R1.clean.fq.gz -p sample1_R2.clean.fq.gz sample1_R1.raw.fq.gz sample1_R2.raw.fq.gz
$fastqc sample1_R1.clean.fq.gz -o output_dir
$fastqc sample1_R2.clean.fq.gz -o output_dir
```

> Key options:
>
> * -a/-A: used to specify the 3' adapters. Remove them if no significant sequences was found in quality control analysis.
> * --trim-n/--max-n: used to handle the unknown bases in reads. --trim-n removes flanking N bases from each reads, and --max-n discards reads containing higher proportion of unknown reads than n.
> * --quality-base: used to specify the encoding of Phred scores. By default, 33 is used as the cutoff. This should always be the case with the data generated in recent years. Some old FASTQ files generated years ago could be encoded as ASCII (phred quality + 64), and --quality-base=64 must be used.
> * -m: used to specify the minimum length of qualified reads. After the trimming of adapter sequences, unknown bases and low-quanlity bases, reads shorter than m will be discarded. -m=30 works well for must cases, and -m=20 is recommend for those with raw sequences of <= 50nt.

## Quantification by Salmon

Salmon is a bulk RNA-Seq quantifier with **wicked-fast speed** and **comparable accuracy** (https://salmon.readthedocs.io/en/latest/salmon.html). It provides two working modes:

* **Mapping-based mode**: this is the feature mode that makes Salmon famous. Samlon employes a SEME (super maximal exact match)-based chaining algorithm to find and score potential mapping loci, making it super fast (because no alignment is needed) and comprably accurate. From version 1.0.0 and higher, Salmon introduced the **Selective Alignment** and **Decoy Sequences** to further improve its quantification accuracy (for more details: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8).
* **Alignment-based mode**: Salmon is not an aligner, so this mode, even though called alignment-based, does take FASTQ files. Instead, you can simpley provide: 1) BAM/SAM files of alignments generated by other aligners, e.g., STAR, Bowtie; and 2) FASTA file of reference transcriptome. This mode doesn't require indexing.

In this pipeline, we suggest to use the mapping-based mode ONLY.

```bash
salmon=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/salmon

## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/Salmon/index_decoy
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/Salmon/index_decoy
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/Salmon/index_decoy
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/Salmon/index_decoy

## transcript-to-gene mapping files
tr2gene_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/gencode.v43.primary_assembly.annotation.gtf
tr2gene_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/gencode.v43lift37.annotation.gtf
tr2gene_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/gencode.vM32.primary_assembly.annotation.gtf
tr2gene_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/gencode.vM25.primary_assembly.annotation.gtf

## Quantification by Salmon
# For single-end sequencing
$salmon quant -i $hg38v43_index -l A -p 8 -g $hg38v43_tr2gene -r sample1.clean.fq.gz --validateMappings -o sample1_quantSalmon

# For paired-end sequencing
$salmon quant -i $hg38v43_index -l A -p 8 -g $hg38v43_tr2gene -1 sample1_R1.clean.fq.gz -2 sample1_R2.clean.fq.gz --validateMappings -o sample1_quantSalmon
```

> Key options:
>
> Only a few options need to be specified by the users (as shown above). Salmon could figure out the others by itself.
>
> * -l: Library type. Salmon employs a three-letter string to denote the library type (https://salmon.readthedocs.io/en/latest/library_type.html#fraglibtype): 1) the relative orientation of two matched mates, including I = inward, O = outward and M = matching; 2) the protocol is stranded or unstranded, including S = stranded and U = unstranded; 3) the strand from which the read originates in a strand-specific protocol, including F = read 1 (or single-end read) from the forward stand and R =  read 1 (or single-end read) from the reverse stand. As a result, there are 9 library types: 3 for single-end library (IU, MU, OU) and 6 for paired-end library (ISF, ISR, MSF, MSR, OSR, OSF). If you don't which is the one used in your case, then just use "A". Salmon will learn it from the real data.
> * -g: 

## Quantification by RSEM

RSEM (RNA-Seq by Expectation Maximization) is an accurate and user-friendly software tool for quantifying transcript and gene abundances from RNA-seq data. 

* RSEM is a quantifier, not an aligner. RSEM can directly take the FASTQ files, but it does align the reads by itself. Instead, it empolys the Bowtie2 (by default) or STAR/HISAT2 (optional) for read alignment.
* RSEM doesn't rely on the existence of a reference genome. It quantifies from the alignments to reference transcriptome.
* RSEM is famouse for its ability to effectively use ambiguously-mapping reads. This is the main reason for its high accuracy in quantification.

In this pipeline, we provide codes for both Bowtie2-RSEM and STAR-RSEM pipelines:

```bash
rsem=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/rsem-calculate-expression

## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/RSEM
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/RSEM
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/RSEM
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/RSEM

## Quantification by Bowtie2-RSEM
$rsem --num-threads 8 \
--bowtie2 --bowtie2-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin \
--bowtie2-sensitivity-level sensitive \
--strandedness reverse --phred33-quals \
--output-genome-bam --sort-bam-by-coordinate \
--paired-end sample1_R1.clean.fq.gz sample1_R2.clean.fq.gz \
$index_hg38v43/index_bowtie2 sample1_quantRSEM

## Quantification by STAR-RSEM
$rsem -num--threads 8 \
--star --star-path /research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin \
--star-gzipped-read-file \
--strandedness reverse --phred33-quals \
--output-genome-bam --sort-bam-by-coordinate \
--paired-end sample1_R1.clean.fq.gz sample1_R2.clean.fq.gz \
$index_hg38v43/index_star sample1_quantRSEM
```

> Key options:
>
> * --strandedness: defines the strandedness of the RNA-Seq reads. This is similar to the library type of Salmon. There are three options: none, forward and reverse.
> * --phred33-quals: encoding method of Phred score.
> * --star-gzipped-read-file: use it if the FASTQ files are in .fq.gz format.

## Quantification by HTSeq

In addition to Salmon and RSEM, we also provide the STAR-HTSeq pipeline for RNA-Seq quantification, because:

* STAR-HTSeq strategy is recommended by GDC.
* STAR

```bash
star=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/STAR
htseq=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/conda_env/bulkRNAseq-2023/bin/htseq-count

## Indexing directories
index_hg38v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg38/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_hg19v43=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/hg19/gencode.release43/bulkRNAseq/STAR/index_overhang100
index_mm39vM32=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm39/gencode.releaseM32/bulkRNAseq/STAR/index_overhang100
index_mm10vM25=/research_jude/rgs01_jude/groups/yu3grp/projects/software_JY/yu3grp/yulab_databases/references/mm10/gencode.releaseM25/bulkRNAseq/STAR/index_overhang100

## Alignment by STAR
$star --genomeDir $index_hg38oh100 --readFilesIn sample1_R1.clean.fq.gz sample1_R2.clean.fq.gz --readFilesCommand zcat \
--runThreadN 8 --twopassMode Basic --quantMode TranscriptomeSAM \
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif \
--outSAMattributes All --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:SampleName SM:SampleName LB:Illumina PL:Illumina PU:Illumina --outFileNamePrefix $outdir/sample/

## Quantification by HTSeq
$htseq -f bam -r pos -s reverse -a 10 -t exon -i gene_id -m intersection-nonempty --nonunique none --secondary-alignments score --supplementary-alignments score -n 8 $outdir/sample/Aligned.out.bam $gtf_hg38 > $outdir/sample/htseq_counts.txt

```





