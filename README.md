# Cloud_Snakemake

Ps : Could not push the git repository due to commits on large files that gihub could not take, instead, i downloaded the necassary files and created a new git repository localy, the cloud VM is still on and could be requested on demand.

####################################################

ATAC-seq workflow.

Mohamed Malek CHAOUCHI

Université Clermont auvergne, 2021

####################################################

### Dependencies, in order of use:
		- FastQC
		- Trimmomotic                             
		- FastQC
		- Bowtie2                  
		- Picard(MarkDuplicates)                 
		- deepTools(multiBamSummary)
		- deepTools(plotCorrelation) 
		- deepTools(plotCoverage)        
		- MACS2(callpeak)

### Launching workflow:
- Execute Snakefile in snakemake conda env.

### Directory structure

		├── .gitignore
		├── README.md
		├──Snakefile
		├── config
		│   ├── config.yaml
		│   └── env
		│        ├── qc.yaml
		│        ├── trim.yaml
		│        ├── bowtie2.yaml
		│        ├── deeptools.yaml
		│        ├── picard.yaml
		│        └── macs2.yaml
		├── data
		│   └── mydatalocal
		│         ├── bowtie2
		│         └── atacseq
		│                └──subset
		└── results
		
- Datasets go into         data/mydatalocal/atacseq/subset
- Bowtie2 indexes go into  data/mydatalocal/bowtie2
- All output will be in    results
- Snakemake config is in   config
- Conda environment config files are in config/env
