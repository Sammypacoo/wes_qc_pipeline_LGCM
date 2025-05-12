#!/bin/bash

# Caminho do diretório onde estão os arquivos
DATA_DIR="./data"

# Arquivos e seus MD5 esperados
declare -A FILES_MD5
FILES_MD5["NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram"]="3d8d8dc27d85ceaf0daefa493b8bd660"
FILES_MD5["NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai"]="15a6576f46f51c37299fc004ed47fcd9"
FILES_MD5["hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed"]="c3a7cea67f992e0412db4b596730d276"

echo "🔍 Verificando integridade dos arquivos via MD5..."
echo

for file in "${!FILES_MD5[@]}"; do
    filepath="$DATA_DIR/$file"
    expected_md5="${FILES_MD5[$file]}"
    
    if [ ! -f "$filepath" ]; then
        echo "❌ Arquivo não encontrado: $filepath"
        continue
    fi

    calculated_md5=$(md5sum "$filepath" | awk '{print $1}')

    if [ "$calculated_md5" = "$expected_md5" ]; then
        echo "✅ $file: OK"
    else
        echo "❌ $file: MD5 inválido"
        echo "   Esperado : $expected_md5"
        echo "   Encontrado: $calculated_md5"
    fi
done
