# Skin Microbiome Project Example
Trial QIIME2 Pipeline in collected data
## Publication
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5710430/#note-DOI170016-1-s
## Collect data using SRAtoolkit
1. trace ncbi: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP016977&o=acc_s%3Aa
2. Download data
## Import data to QIIME2
## Visualization quality control
## Denoising using DADA2
## Generate Silva database
### Generate 16S V3-V4 amplicon reference
```bash
qiime rescript get-silva-data --p-version '138' --p-target 'SSURef_NR99' --p-include-species-labels --o-silva-sequences silva-138-ssu-nr99-seqs.qza --o-silva-taxonomy silva-138-ssu-nr99-tax.qza
qiime rescript cull-seqs --i-sequences silva-138-ssu-nr99-seqs.qza --o-clean-sequences silva-138-ssu-nr99-seqs-cleaned.qza
qiime taxa filter-seqs --i-sequences silva-138-ssu-nr99-seqs-cleaned.qza --i-taxonomy silva-138-ssu-nr99-tax.qza --p-exclude 'd__Eukaryota' --p-mode 'contains' --o-filtered-sequences silva138_noEuk_seqs.qza
qiime rescript filter-seqs-length-by-taxon --i-sequences silva138_noEuk_seqs.qza --i-taxonomy silva-138-ssu-nr99-tax.qza --p-labels Archaea Bacteria --p-min-lens 900 1200 --o-filtered-seqs silva138_noEuk_AB_seqs.qza --o-discarded-seqs silva138_Euk_seqs_discard.qza
qiime rescript dereplicate --i-sequences silva138_noEuk_AB_seqs.qza --i-taxa silva-138-ssu-nr99-tax.qza --p-threads 12 --o-dereplicated-sequences silva138_noEuk_AB_seqs_uniq.qza --o-dereplicated-taxa silva138_noEuk_AB_tax_uniq.qza
qiime feature-classifier extract-reads --i-sequences silva138_noEuk_AB_seqs_uniq.qza --p-f-primer ACTCCTAYGGGRBGCASCAG --p-r-primer AGCGTGGACTACNNGGGTATCTAAT --p-n-jobs 12 --o-reads silva138_AB_V3-V4seqs.qza
```
## Taxonomy classifier
## Alpha/Beta diversity measurement
