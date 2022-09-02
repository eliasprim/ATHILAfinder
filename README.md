```diff
+ATHILAfinder
```
<br />

### Overview:

The ATHILA clade of Ty3 LTR retrotransposons is found across plants and occupies ~2.7% of the *Arabidopsis thaliana* genome. A recent study showed that ATHILAs are the only transposable elements that have invaded the core centromeric arrays in *Arabidopsis thaliana*, disrupting their genetic and epigenetic organization. To better assess the significance of ATHILA elements in the function and evolution of centromeres in *Arabidopsis thaliana* and other species, a computational pipeline was needed that could identify ATHILAs with high efficiency and sensitivity. Current annotation tools may lack this precision as they are broadly designed to identify all types of LTR retrotransposons. 


Here, we present ATHILAfinder, an expertly built tool for the large-scale, yet accurate, discovery and preliminary analysis of intact ATHILA elements and their solo LTRs in plant genomes. ATHILAfinder uses as primary seeds motifs that are specific to the LTR-internal junctions of the ATHILA sequence. It then builds intact elements around these motifs using additional filters to return high quality elements. A homology-based step is included to identify intact elements and solo LTRs. ATHILAfinder calculates the insertion age of the ATHILA elements and scans the internal sequence with Hidden Markov Models of LTR retrotransposon genes retrieved by Pfam. Validation of ATHILAfinder was based on two recently published telomere-to-telomere *Arabidopsis thaliana* assemblies (Col-0 is included in this GitHub project) and the reference genome of *Arabidopsis lyrata*. 


Phylogenetic and structural analyses show that these are indeed true ATHILA sequences. These results were compared with the output of the widely-used EDTA pipeline. EDTA identified only ~30% of intact ATHILAs. The availability of pipelines tailored for specific TEs is limited, but our results suggest that they can significantly improve TE identification and downstream genome analysis.


The various outputs of ATHILAfinder allow several different analyses of ATHILA elements, like phylogenetic and structural analyses. 


ATHILAfinder was made for the model organism *Arabidopsis thaliana* and its closest related species *Arabidopsis lyrata*. However, with some homology-based steps related to PBS and PPT signatures of ATHILA elements, it can be applied to any plant genome.   


<br />
<br />

## **The ATHILAfinder dependancies are the following:** 

### Software: 

```diff
! Vmatch (2.3.1), BEDTools (v2.27.1), blastall (2.2.26), EMBOSS (6.6.0.0), HMMER (3.3)
```

### Programming Languages and (packages): 

```diff
- R (dplyr, tidyr, data.table, plyranges, stringr, tibble, Biostrings, seqinr)

- Python (sys, regex, from itertools import groupby, from Bio import SeqIO, from Bio.SeqRecord import SeqRecord)

- Perl

- Awk
```

<br />
<br />

### **More info about the execution and the parameters can be found in the run_ATHILAfinder_parameters.sh file** 
