#!/bin/bash

# Cria a pasta data se não existir
mkdir -p data
cd data

echo "🔽 Baixando arquivos de alinhamento e genoma de referência..."

# CRAM e CRAI
curl -C - -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA06994/exome_alignment/NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram
curl -C - -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA06994/exome_alignment/NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai

# BED das regiões exônicas
curl -C - -O https://www.twistbioscience.com/sites/default/files/resources/2022-12/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed

# Genoma de referência e index
curl -C - -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
curl -C - -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

echo "🔽 Baixando arquivos para VerifyBamID..."

# Arquivos do VerifyBamID
curl -C - -O https://github.com/Griffan/VerifyBamID/raw/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.UD
curl -C - -O https://github.com/Griffan/VerifyBamID/raw/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.V
curl -C - -O https://github.com/Griffan/VerifyBamID/raw/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.bed
curl -C - -O https://github.com/Griffan/VerifyBamID/raw/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.mu

echo "✅ Download concluído com sucesso!"
