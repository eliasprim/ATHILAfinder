# ATHILAfinder: Tool to detect ATHILA LTR retrotransposons in plant genomes

Current annotation tools are generally designed to identify all classes of LTR retrotransposons. However, pipelines specifically tailored to specific transposable elements (TEs) remain scarce. Our results demonstrate that TE-specific tools can substantially enhance the accuracy of TE identification and improve subsequent genome analyses. 

ATHILAfinder is a specialized tool developed for large-scale yet precise discovery and preliminary characterization of intact ATHILA elements and their solo LTRs in plant genomes. The program begins by searching for lineage-specific seed motifs corresponding to the LTR–internal junctions of ATHILA elements. It then reconstructs complete elements around these motifs using additional filtering criteria to ensure high-quality elements. A complementary homology-based step enables the identification of both intact elements and solo LTRs. In addition, ATHILAfinder estimates the insertion age of ATHILA elements and screens their internal regions with Hidden Markov Models (HMMs) representing LTR retrotransposon genes from the Pfam database. The tool generates multiple outputs suitable for downstream analyses, including phylogenetic and structural investigations of ATHILA elements. Validation of ATHILAfinder was performed using two published telomere-to-telomere *Arabidopsis thaliana* assemblies (including the Col-CEN genome available in this GitHub repository), the *Arabidopsis lyrata* reference genome, and several other *Brassicaceae* species such as *Brassica nigra*, *Draba nivalis*, *Erysimum cheiranthoides*, and *Megadenia pygmaea*.

## Update

ATHILAfinder has been modified to detect and characterise elements from any LTR retrotransposon lineage in plants. The default PBS and PPT signatures included are specific to ATHILA elements from the *Brassicaceae* family. Users wishing to analyze other lineages must provide the corresponding lineage-specific PBS and PPT signatures.

## Requirements

- Linux OS or macOS
- Installation of conda 


## Installation and execution

0. Clone the repository
```shell
git clone https://github.com/eliasprim/ATHILAfinder

cd ATHILAfinder
```

1. Create conda environment
```shell
conda env create --file athilafinder.yml
```

2. Activate conda environment

```shell
conda activate athilafinder
```

## Usage

1. Modify ATHILAfinder parameters

A detailed description of all available parameters is provided in the ATHILAfinder_runscript.sh file. Users can modify these parameters by selecting the appropriate command-line options for their operating system (Linux or macOS).  

2. Execute ATHILAfinder_runscript.sh

```shell
bash ATHILAfinder_runscript.sh
```

3. Get the results

Upon successful completion of an ATHILAfinder run, two output directories are generated: the main and the side outputs. The main output directory contains all files required for downstream analyses, whereas the side output directory includes supplementary files useful for validation and in-depth characterization of ATHILA elements.

## Demo data

*Arabidopsis thaliana* Col-CEN_v1.2.fasta is provided for a toy run.

## Main dependencies (all can be install using the yml file)

BEDTools: Quinlan,A.R. and Hall,I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26, 841–842.

BLAST: Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of molecular biology, 215(3), 403-410.

EMBOSS: Rice, P., Longden, I. and Bleasby, A. (2000) EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics 16, 276–277. 

HMMER: http://hmmer.org/

Vmatch: Kurtz,S. The Vmatch large scale sequence analysis software A Manual.

## Citing ATHILAfinder

The accompanying manuscript will be available soon. In the meanwhile, please refer to this repository.
