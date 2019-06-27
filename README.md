# BenchmarkBPprediction
---

> **Cite as:** Assessment of branch point prediction tools to predict physiological branch points and their alteration by variants.
*Raphaël LEMAN, Hélène Tubeuf, Sabine Raad, Isabelle Tournier, Céline Derambure, Raphaël Lanos, Pascaline Gaildrat, Gaia Castelain, Julie Abdat, Audrey Kilian, Stéphanie Baert-Desurmont, Angelina Legros, Nicolas Goardon, Céline Quesnelle, Agathe Ricou, Laurent Castera, Dominique Vaur, Gérald Le Gac, Chandran Ka, Yann Fichou, Françoise Bonnet-Dorion, Nicolas Sevenet, Marine Guillaud-Bataille, Nadia Boutry-Kryza, Inès Schultz, Virginie Caux-Moncoutier, Maria Rossing, Logan C. Walker, Amanda B. Spurdle, Claude Houdayer, Alexandra Martins, Sophie Krieger*

This repository contains the data and scripts used for this study.
Three sets of data were used: the Ensembl data, the RNAseq data and variant data. The scripts used to compare bioinformatics tools (HSF, SVM-BPfinder, BPP, Branchpointer, LaBranchoR and RNABPS) are in R language.

## Installation and Usage

To run these scripts, the following dependencies are needed:
+ R (v3.0 or later)
```Bash
sudo apt-get update
sudo apt-get install r-base r-base-dev
```
+ package 'ROCR'
```R
install.packages('ROCR')
```
+ package 'gplots'
```R
install.packages('gplots')
```
+ package 'ggplot2'
```R
install.packages('ggplot2')
```
+ package 'VennDiagram'
```R
install.packages('VennDiagram')
```

Download the repository and the Ensembl data
```Bash
git clone https://github.com/raphaelleman/BenchmarkBPprediction
cd ./BenchmarkBPprediction
```

Download [Branch point data from Ensembl dataset](https://mega.nz/#!C753GIKL!tspYhFUaNC6ZgUxwCnVV9bvFctRGjU-4Vtos_lZc1C4 "tittle"). Then put it in the 'data' folder.

```Bash
#Usage Example
Rscript ./scripts/studyRNAseqBP.r
```

## SCRIPTS

The script [**studyEnsemblBPscore.r**](https://github.com/raphaelleman/BenchmarkBPprediction/blob/master/scripts/studyEnsemblBPscore.r "tittle") is deicated to the Ensembl data and permits to perform:
* ROC analysis
* VENN diagram

The script [**studyRNAseqBP.r**](https://github.com/raphaelleman/BenchmarkBPprediction/blob/master/scripts/studyRNAseqBP.r "tittle") is deicated to the RNAseq data and permits to perform:
* ROC analysis
* VENN diagram
* Boxplot expression according to presence or not of Branch point
* Correlation expression and branch point score

The script [**studyVariantBPscore.r**](https://github.com/raphaelleman/BenchmarkBPprediction/blob/master/scripts/studyVariantsBPscore.r "tittle") is deicated to the variant data and permits to perform:
* Histogram of variant repartition
* ROC analysis
* Histogram of relatif repartition
* Score combination by regresion logistic

## DATA

Hereafter, the scheme of each data set use in this study.

### Ensembl data

The set of natural 3'ss and control AG

| Column names | Description |
|------------:|:--------:|
| Chr | Chromosone of 3'ss |
| strand | Strand of transcript |
| sstype | "Acc", acceptor splice site |
| pos | Position of 3'ss (hg19) |
| UsedSite | 1- alternative 3'ss, 0- control AG |
| MES | MES score |
| SSF | SSF-like score |
| ESR | ESRseq score |
| Branch_class | Score above optimal threshold (Branchpointer) |
| Branch_pos | Position of predicted BP by Branchpointer |
| Branch_score | Score of predicted BP by Branchpointer |
| zscBPP | Score of predicted BP by BPP |
| zscSVM | Score of predicted BP by SVM-BPfinder |
| scoreLB | Score of predicted BP by LaBranchoR |
| scoreRNABPS | Score of predicted BP by RNABPS |

### RNAseq data

The set of alternative 3'ss indentified from RNAseq data

| Column names | Description |
|------------:|:--------:|
| chr | Chromosone of 3'ss |
| Start_hg19 | Start of the junction (End if strand -) |
| Pos_hg19 | Position of 3'ss (hg19) |
| Strand | Strand of transcript |
| Transcript | Transcript name (RefSeq) |
| Gene | Gene symbol |
| UsedSite | 1- alternative 3'ss, 0- control AG |
| Expression_Average_(%) | Average expression in % ("." for control AG) |
| Expression_Average_(read) | Average expression in read count ("." for control AG) |
| nbSamp | Nb samples supporting the alternative 3'ss ("." for control AG) |
| MES | MES score |
| SSF | SSF-like score |
| bp_pos | Position of predicted BP by BPP |
| bp_zsc | Score of predicted BP by BPP |
| BPP_class | Score above optimal threshold (BPP) |
| ss_dist | Position of predicted BP by SVM-BPfinder |
| svm_scr | Score of predicted BP by SVM-BPfinder |
| SVM_class | Score above optimal threshold (SVM-BPfinder) |
| branchpoint_dist | Position of predicted BP by Branchpointer |
| branchpoint_prob | Score of predicted BP by Branchpointer |
| branchpoint_class | Score above optimal threshold (Branchpointer) |
| LB_pos | Position of predicted BP by LaBranchoR |
| LB_score | Score of predicted BP by LaBranchoR |
| LB_class | Score above optimal threshold (LaBranchoR) |
| RNABPS_pos | Position of predicted BP by RNABPS |
| RNABPS_score | Score of predicted BP by RNABPS |
| RNABPS_class | Score above optimal threshold (RNABPS) |
| Cum_class | Overlapping of predictions |

### Variant data

The collection of variant with their *in vitro* RNA studies.

| Column names | Description |
|------------:|:--------:|
| ID | Name of variant |
| chr | Chromosome |
| strand | Strand of transcript (+: forward, -: reverse) |
| gene | Gene Name |
| transcript | Transcript ID (RefSeq) |
| intron | Intron number |
| cNomen | Transcriptomic coordinate |
| gNomen | Genomic coordinate (hg19) |
| distSS | Distance between variant and Acceptor site |
| varType | Type of mutation |
| nb_analyse | Number of study for the variant |
| Assay.Method | Materiel and Method used |
| Authora | The contributor |
| published | If the RNA studies were published or not |
| Result | Overall impact of variant on splicing |
| class_effect | Class of splicing alteration (0: no impact, 1: splicing alteration) |
| MES_wt | MES score of acceptor splice site (wildType) |
| MES_mut | MES score of acceptor splice site (mutated) |
| delta_MES | Delta MES score |
| SSF_wt | SSF-like score of acceptor splice site (wildType) |
| SSF_mut | SSF-like score of acceptor splice site (mutated) |
| delta_SSF | Delta SSF-like score |
| to_3prime | Distance between variant and Acceptor site |
| BP_num_REF | Number of branch point found by branchpointer (wildType) |
| BP_num_ALT | Number of branch point found by branchpointer (mutated) |
| deleted_n | Number of deleted branch points (branchpointer) |
| created_n | Number of created branch points (branchpointer) |
| max_prob_REF | Branchpointer score, maximal value (wildType) |
| max_prob_ALT | Branchpointer score, maximal value (mutated) |
| max_U2_REF | U2 score (acceptor site, wildType) |
| max_U2_ALT | U2 score (acceptor site, mutated) |
| PosPB_Branch | Position of branch point with maximal score |
| ProbPBarea_Branch | Interval of 4-mer motif of branch point |
| MutInPBarea_Branch | If variant was in the 4-mer of branch point |
| Delta_prob | Delta Branchpointer score |
| bps_WT | Sequence of wildType Branch point with maximal score |
| bps_MUT | Sequence of mutated Branch point with maximal score |
| bp_pos_WT | Relative position of branch point to the acceptor site (wildType) |
| bp_pos_MUT | Relative position of branch point to the acceptor site (mutated) |
| zsc_WT | BPP score (wildType) |
| zsc_MUT | BPP score (mutated) |
| PosPB_BPP | Position of branch point with maximal score |
| ProbPBarea_BPP | Interval of 4-mer motif of branch point |
| MutInPBarea_BPP | If variant was in the 4-mer of branch point |
| Delta_BPP | Delta score of BPP |
| ss_dist_WT | Relative position of branch point to the acceptor site (wildType) |
| ss_dist_MUT | Relative position of branch point to the acceptor site (mutated) |
| bp_seq_WT | Sequence of wildType Branch point with maximal score |
| bp_seq_MUT | Sequence of mutated Branch point with maximal score |
| svm_scr_WT | SVM-BPfinder score (wildType) |
| svm_scr_MUT | SVM-BPfinder score (mutated) |
| PosPB_SVM | Position of branch point with maximal score |
| ProbPBarea_SVM | Interval of 4-mer motif of branch point |
| MutInPBarea_SVM | If variant was in the 4-mer of branch point |
| Delta_SVM | Delta score of SVM-BPfinder |
| score_wt | HSF score (wildType) |
| score_mut | HSF score (mutated) |
| delta_HSF | Delta score of HSF |
| interpret | Overall prediction of HSF |
| LB_pos_WT | Relative position of branch point to the acceptor site (wildType) |
| LB_pos_MUT | Relative position of branch point to the acceptor site (mutated) |
| LB_score_WT | LaBranchoR score (wildType) |
| LB_score_MUT | LaBranchoR score (mutated) |
| PosPB_LB | Position of branch point with maximal score |
| ProbPBarea_LB | Interval of 4-mer motif of branch point |
| MutInPBarea_LB | If variant was in the 4-mer of branch point |
| Delta_LB | Delta score of LaBranchor |
| RBPS_posA_WT | Relative position of branch point to the acceptor site (wildType) |
| RBPS_posA_MUT | Relative position of branch point to the acceptor site (mutated) |
| RBPS_score_WT | RNABPS score (wildType) |
| RBPS_score_MUT | RNABPS score (mutated) |
| PosPB_RBPS | Position of branch point with maximal score |
| ProbPBarea_RBPS | Interval of 4-mer motif of branch point |
| MutInPBarea_RBPS | If variant was in the 4-mer of branch point |
| Delta_RBPS | Delta score of RNABPS |
