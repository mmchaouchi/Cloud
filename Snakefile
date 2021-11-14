
################################################
## The main entry point of Snakemake workflow ##
# Run FASTQC on ATACseq samples
# Author: Mohamed Malek CHAOUCHI
# AFFILIATION: Clermont Auvergne University
# CONTACT: mohamed_malek.chaouchi@uca.fr
################################################
configfile: "config/config.yaml",
import glob
##########################################################
#### Target rule with final output of the workflow
#########################################################
# Collecting runs in 'ech' variable
ech, num=glob_wildcards("data/mydatalocal/atacseq/subset/{samplename}_{number}.fastq.gz")
rule all:
    input:
        expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
        expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
        expand("results/trim/{samplename}_1_trim_paired.fastq.gz", samplename=ech), 
        expand("results/trim/{samplename}_2_trim_paired.fastq.gz", samplename=ech),
        expand("results/fastqc_post/{sample}_1_trim_paired_fastqc.zip", sample =ech),
        expand("results/fastqc_post/{sample}_2_trim_paired_fastqc.html", sample=ech),
        expand("results/bowtie2/{samplename}.bam", samplename=ech),
        expand("results/bowtie2_cleaned/{samplename}_cleaned.bam", samplename=ech), 
        "results/npz/results.npz",
        expand("results/deeptools_cor/heatmap_SpearmanCorr_readCounts.png"),
        expand("results/deeptools_cor/SpearmanCorr_readCounts.tab"),
        expand("results/deeptools_cov/coverage.png"),
        expand("results/deeptools_cov/coverage.tab"),
        expand("results/macs2/{sample}_peaks.narrowPeak", sample =ech),
        expand("results/macs2/{sample}_summits.bed", sample =ech),
        expand("results/macs2/{sample}_model.r", sample =ech),
        expand("results/macs2/{sample}_peaks.xls", sample =ech)
rule unzip:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("tmp/{sample}.fastq")
    shell:
        """
        mkdir -p tmp
        gunzip -c {input} > {output}
        """
rule fastqc_init:
    input:
        "tmp/{sample}.fastq"
    output:
        html="results/fastqc_init/{sample}_fastqc.html",
        zip="results/fastqc_init/{sample}_fastqc.zip"
    conda:
        "config/env/qc.yaml"
    threads: 2
    shell:
        """
        mkdir -p results/fastqc_init
        fastqc {input} -o "results/fastqc_init" -t {threads}
        """
rule trimming:
    input:
        r1="data/mydatalocal/atacseq/subset/{samplename}_1.fastq.gz",
        r2="data/mydatalocal/atacseq/subset/{samplename}_2.fastq.gz"
    output:
        fwd_P="results/trim/{samplename}_1_trim_paired.fastq.gz",
        rvr_P="results/trim/{samplename}_2_trim_paired.fastq.gz",
        fwd_U="results/trim/{samplename}_1_trim_unpaired.fastq.gz",
        rvr_U="results/trim/{samplename}_2_trim_unpaired.fastq.gz"
    conda:
        "config/env/trim.yaml"
    threads: 2
    shell:
        """
        mkdir -p results/trim
        trimmomatic PE -threads {threads} \
        -trimlog results/trim/trim.log -summary results/trim/stats \
        {input.r1} {input.r2} \
        {output.fwd_P} \
        {output.fwd_U} \
        {output.rvr_P} \
        {output.rvr_U} ILLUMINACLIP:data/mydatalocal/atacseq/subset/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:33
        """
rule fastqc_post:
    input:
        "results/trim/{sample}.fastq.gz"
    output:
        html="results/fastqc_post/{sample}_fastqc.html",
        zip="results/fastqc_post/{sample}_fastqc.zip"
    conda:
        "config/env/qc.yaml"
    threads: 2
    shell:
        """
        mkdir -p results/fastqc_post
        fastqc {input} -o "results/fastqc_post" -t {threads}
        """
rule bowtie2:
    input:
        r1="results/trim/{samplename}_1_trim_paired.fastq.gz",
        r2="results/trim/{samplename}_2_trim_paired.fastq.gz"
    output:
        aln="results/bowtie2/{samplename}.bam"
    conda:
        "config/env/bowtie2.yaml"
    threads: 2
    shell:
        """
        bowtie2  --very-sensitive -p 1 -k 10  -x data/mydatalocal/bowtie2/all -1 {input.r1}  -2 {input.r2} | samtools view -q 2 -bS  -  |  samtools sort - -o {output.aln}
        samtools index -b {output.aln}
        """

rule picard:
    input:
        aln_a="results/bowtie2/{samplename}.bam"        
    output:
        aln_clean="results/bowtie2_cleaned/{samplename}_cleaned.bam",
        aln_clean_txt="results/bowtie2_cleaned/{samplename}_cleaned.txt"
    conda:
        "config/env/picard.yaml"
    threads: 2
    shell:
        """
        picard MarkDuplicates \
        I={input.aln_a} \
        O={output.aln_clean} \
        M={output.aln_clean_txt} \
        REMOVE_DUPLICATES=true
        samtools index -b {output.aln_clean}
        """
rule deeptools_npz:
    input:
        expand("results/bowtie2_cleaned/{samplename3}_cleaned.bam", samplename3=ech)
    output:
        "results/npz/results.npz"
    conda:
        "config/env/deeptools.yaml"
    threads: 2
    shell:
        """
        multiBamSummary bins -b {input} \
        -o {output:q}
        """

rule deeptools_cor:
    input:
        "results/npz/results.npz"
    output:
        png="results/deeptools_cor/heatmap_SpearmanCorr_readCounts.png",
        matrix="results/deeptools_cor/SpearmanCorr_readCounts.tab"
    conda:
        "config/env/deeptools.yaml"
    threads: 2
    shell:
        """
        plotCorrelation \
        -in {input}  \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.png}   \
        --outFileCorMatrix {output.matrix}
        """

rule deeptools_cov:
    input:
       expand("results/bowtie2_cleaned/{samplename}_cleaned.bam", samplename=ech)
    output:
        png="results/deeptools_cov/coverage.png",
        matrix="results/deeptools_cov/coverage.tab"
    conda:
        "config/env/deeptools.yaml"
    threads: 2
    shell:
        """
        plotCoverage --bamfiles {input} \
        --plotFile {output.png} \
        --plotTitle "Coverage" \
        --outRawCounts {output.matrix} \
        --ignoreDuplicates
        """

rule macs2:
    input:
        "results/bowtie2_cleaned/{sample}_cleaned.bam"
    output:
        "results/macs2/{sample}_model.r", "results/macs2/{sample}_peaks.narrowPeak","results/macs2/{sample}_peaks.xls", "results/macs2/{sample}_summits.bed"
    conda:
        "config/env/macs2.yaml"
    threads: 2
    shell:
        """
        macs2 callpeak -t {input}  \
        -f BAM \
        -n {wildcards.sample} \
        --outdir results/macs2
        """
