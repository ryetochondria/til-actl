
# **RNAseq pipeline log**

 2018 -

## PIPELINE PREPARATION

### Backup raw data

1. Fastq data backup to Google cloud storage

### Reference files generated using the GENCODE v26 build (GRCh38)

1. Build STAR index
2. Build RSEM index
3. Build genome reference fasta
4. Build gene annotation file, includes ERCC genes

**Constructing reference genome files:** Reference files were built from GRCh38.p10 with annotations GENCODE v26 (release date March 2017). HLA genes, alt contigs/scaffolds/loci, and decoys were excluded. ERCC spike-in genome was included (ERCC92).

Original source: (GRCh38.p10)[https://www.gencodegenes.org/human/release_26.html]

**Genome Reference**
* Hommo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta

**Gene Annotations**
* gencode.v26.GRCh38.ERCC.genes.gtf

**STAR index**
* Built in STAR for 2x75bp read-length
* bash scripts/index_build_star.sh
* STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v26_oh75

**RSEM reference**
* bash scripts/index_build_rsem.sh
* rsem_reference_GRCh38_gencode26_ercc