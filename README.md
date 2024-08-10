# Skin Microbiome Project Example Tutorial
Trial QIIME2 Pipeline in collected data
## Publication
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5710430/#note-DOI170016-1-s

Biopsy specimens were analyzed at the Department of Microbiology and Infection Control, Statens Serum Institut, Copenhagen, Denmark. DNA was extracted using a kit (QIAamp DNA Mini Kit; Qiagen) according to the manufacturer’s instructions for tissues. For each batch of DNA extraction, a “negative” control was included containing buffers but no sample material for downstream analysis. DNA was amplified using a 2-step polymerase chain reaction using custom 341F/806R primers targeting the V3-V4 16S regions, as well as 3 primer sets targeting the hypervariable regions V3-V4 of the 18SrDNA gene, and amplicons were sequenced on a desktop sequencer (MiSeq; Illumina, Inc) using the v2 reagent kit. For details concerning primer design and library preparation, see the eAppendix in the Supplement. Sequence data are available at the European Nucleotide Archive (accession number PRJEB15266).
## Collect data using SRAtoolkit
1. Trace ncbi: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=ERP016977&o=acc_s%3Aa
2. Download data
3. Install SRAtoolkit
  ```bash
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --show channels
  ```
    Install SRAtoolkit:
  ```bash
    conda create -n sratool sra-tools
    conda env list 
  ```
4. Download data:
   
  Get this file in Trace link
  ```bash
    SraAccList.txt
  ```
  After that run this command:
    
  ```bash
    prefetch --option-file SraAccList.txt
    find . -name '*.sra' -print0 | xargs -0 mv -t . 
    find . -type d -empty -delete
    ls *.sra | parallel -j0 fastq-dump --split-files --origfmt {}
    mkdir fastq
    mv *.fastq fastq
    cd fastq
    gzip fastq
    mkdir sra
    mv *.sra sra
  ```
## Import Fastq to QIIME2
  We perform importing data to QIIME2 following manifest protocol:
  
  Step 1: Generate manifest file
```bash
  echo -e 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath' > manifest.tsv
  for FOR in reads/*_1*gz;
  do ID=$(basename $FOR | cut -f1 -d_);
  REV=${FOR/_1/_2};
  echo -e "${ID}\t${PWD}/${FOR}\t${PWD}/${REV}";
  done >>manifest.tsv
```
  Step 2: Import to QIIME2
```bash
  qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path reads.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
## Visualization quality control
```bash
  qiime demux summarize \
  --i-data reads.qza \
  --o-visualization reads.qzv
```
## Denoising using DADA2
```bash
  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads.qza \
  --p-trim-left-f 1 \
  --p-trim-left-r 1 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```
Visualization
  ```bash
  qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file sample-metadata.tsv
  qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv
  qiime metadata tabulate \
    --m-input-file denoising-stats.qza \
    --o-visualization denoising-stats.qzv
  ```
## Generate Silva database
### Generate 16S V3-V4 amplicon reference
```bash
qiime rescript get-silva-data --p-version '138' --p-target 'SSURef_NR99' --p-include-species-labels --o-silva-sequences silva-138-ssu-nr99-seqs.qza --o-silva-taxonomy silva-138-ssu-nr99-tax.qza
qiime rescript cull-seqs --i-sequences silva-138-ssu-nr99-seqs.qza --o-clean-sequences silva-138-ssu-nr99-seqs-cleaned.qza
qiime taxa filter-seqs --i-sequences silva-138-ssu-nr99-seqs-cleaned.qza --i-taxonomy silva-138-ssu-nr99-tax.qza --p-exclude 'd__Eukaryota' --p-mode 'contains' --o-filtered-sequences silva138_noEuk_seqs.qza
qiime rescript filter-seqs-length-by-taxon --i-sequences silva138_noEuk_seqs.qza --i-taxonomy silva-138-ssu-nr99-tax.qza --p-labels Archaea Bacteria --p-min-lens 900 1200 --o-filtered-seqs silva138_noEuk_AB_seqs.qza --o-discarded-seqs silva138_Euk_seqs_discard.qza
qiime rescript dereplicate --i-sequences silva138_noEuk_AB_seqs.qza --i-taxa silva-138-ssu-nr99-tax.qza --p-threads 12 --o-dereplicated-sequences silva138_noEuk_AB_seqs_uniq.qza --o-dereplicated-taxa silva138_noEuk_AB_tax_uniq.qza
qiime feature-classifier extract-reads --i-sequences silva138_noEuk_AB_seqs_uniq.qza --p-f-primer ACTCCTAYGGGRBGCASCAG --p-r-primer AGCGTGGACTACNNGGGTATCTAAT --p-n-jobs 12 --o-reads silva138_AB_V3-V4seqs.qza
qiime rescript dereplicate --i-sequences silva138_AB_V3-V4seqs.qza --i-taxa silva138_noEuk_AB_tax_uniq.qza --o-dereplicated-sequences silva138_AB_V3-V4seqs_uniq.qza --o-dereplicated-taxa silva138_AB_V3-V4taxa_uniq.qza
qiime rescript evaluate-fit-classifier --i-sequences silva138_AB_V3-V4seqs_uniq.qza --i-taxonomy silva138_AB_V3-V4taxa_uniq.qza --o-classifier silva138_AB_V3-V4_classifier.qza --o-observed-taxonomy silva138_AB_V3-V4_predicted_taxonomy.qza --o-evaluation silva138_AB_V3-V4_classifier_eval.qzv --p-n-jobs 0
```
## Taxonomy classifier
## Alpha/Beta diversity measurement
