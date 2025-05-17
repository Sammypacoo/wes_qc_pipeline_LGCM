
# Identificação 

Estou participando do **Desafio Técnico – Bioinformata**, proposto pelo **Laboratory of Genetics and Molecular Cardiology**.

- **Nome completo:** Samantha Lucciola Guedes Paço
- **Celular:** (21) 99929-3774
- **E-mail institucional:** samanthalucciola@gmail.com
- **GitHub:** [https://github.com/Sammypacoo](https://github.com/Sammypacoo)
- **LinkedIn:** [https://www.linkedin.com/in/samanthapaco/](https://www.linkedin.com/in/samanthapaco/)


#  Dados de Entrada
script ./scripts/download_data.sh cria uma pasta ./data dentro do diretorio do github e baixa os arquivos abaixo :

Os arquivos utilizados neste desafio foram organizados na pasta `data/`. Eles se dividem em duas categorias: arquivos de alinhamento e anotação genômica, e arquivos de referência para análise de contaminação com o software **VerifyBamID**.

##  Arquivos de Alinhamento e Referência Genômica
Eles foram baixados no site abaixo 
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA06994/exome_alignment/

- [CRAM do exoma - NA06994 (1000 Genomes Project)](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA06994/exome_alignment/NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram)  
- [CRAI (índice do CRAM)](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA06994/exome_alignment/NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai)  
- [BED das regiões exônicas (Twist Bioscience hg38 v2.0.2)](https://www.twistbioscience.com/sites/default/files/resources/2022-12/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed)  
- [Genoma de referência - GRCh38 full + decoy + HLA](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa)  
- [Index do genoma de referência (.fai)](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai)  

## Arquivos para Análise de Contaminação (VerifyBamID)

Esses arquivos devem ser baixados e salvos na mesma pasta `data/` para a execução do **VerifyBamID**:

- [1000g.phase3.10k.b38.exome.vcf.gz.dat.UD](https://github.com/Griffan/VerifyBamID/blob/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.UD)  
- [1000g.phase3.10k.b38.exome.vcf.gz.dat.V](https://github.com/Griffan/VerifyBamID/blob/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.V)  
- [1000g.phase3.10k.b38.exome.vcf.gz.dat.bed](https://github.com/Griffan/VerifyBamID/raw/refs/heads/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.bed)  
- [1000g.phase3.10k.b38.exome.vcf.gz.dat.mu](https://github.com/Griffan/VerifyBamID/blob/master/resource/exome/1000g.phase3.10k.b38.exome.vcf.gz.dat.mu)  


# Organizacão do Snakemaker  
##  Rule 1:  Fazer os arquivos serem executáveis
  ```bash
  chmod +x scripts/*.sh
  ```
##  Rule 2: Verificação Hashes MD5 dos arquivos cram,crai e bed baixados

### Hashes MD5 esperados:

  - Arquivo cram`NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram`  
    `3d8d8dc27d85ceaf0daefa493b8bd660`

  - Arquivo crai `NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai`  
    `15a6576f46f51c37299fc004ed47fcd9`

  - Arquivo BED (Twist Exome):  
    `c3a7cea67f992e0412db4b596730d276`

  Foi utilizado a pipeline presente no script md5_check.sh, que aplicou o método  abaixo e depois foi comparado os valores MD5 calculados e esperados

  ```bash
  md5sum arquivo.txt
  ```
## Rule 3: Conversão do arquivo cram em bam e também a criacao do arquivo index bai
Recorri ao samtools pra conversão do arquivo bam, os arquivos output bam são encaminhados para pasta ./results

  ```bash
  samtools view -T {input.ref} -b -o {output.bam} {input.cram}
  samtools index  {output.bam}
  ```


##  Rule 4: Cálculo de cobertura 
Recorri ao mosdept para o cálculo da cobertura

  ```bash
  mosdepth \
      --by {input.bed} \
      --thresholds 10,30 \
      --threads 4 \
         ./results/cobertura \
      {input.bam}
  ```
 O mosdept calcula a cobertura (número de leituras alinhadas) em regiões do genoma a partir de arquivos BAM/CRAM. É comumente usado para:
  - Avaliar a qualidade de sequenciamento

  - Verificar se regiões de interesse estão suficientemente cobertas

  - Preparar dados para análise de CNV ou filtragem por cobertura

### Arquivos de saída do mosdepth foram encaminhados para pasta ./results
A ferramenta [`mosdepth`](https://github.com/brentp/mosdepth) foi utilizada com a opção `--by`, aplicando um arquivo BED para calcular a profundidade de cobertura apenas nas regiões especificadas, nesse projeto o exoma.

| Arquivo                             | Descrição                                                                                                     |
|-----------------------------------|---------------------------------------------------------------------------------------------------------------|
| `{prefix}.mosdepth.global.dist.txt` | Distribuição da profundidade por base em todo o arquivo BAM, independente das regiões do BED.                   |
| `{prefix}.mosdepth.regions.dist.txt` | Distribuição da profundidade por base das regiões do BED. 
| `{prefix}.mosdepth.summary.txt`      | Estatísticas de cobertura agregadas por região do BED (média, mediana, min, max, etc.).                         |
| `{prefix}.regions.bed.gz`              | Cobertura média (ou mediana, se opção usada) por região do arquivo BED, em formato BED compactado.             |
| `{prefix}.thresholds.bed.gz`           | Contagem de bases em cada região do BED com cobertura ≥ aos thresholds configurados (ex: ≥10x, ≥30x).          |                     |
| `{prefix}.per-base.bed.gz`              | **Não gerado quando se usa `--by`.**  


Iremos utilizar os aquivos {prefix}.regions.bed.gz e {prefix}.thresholds.bed.gz 

## Rule 5: Relatório da Cobertura

Este passo executa o script `scripts/cobertura_report.py` para gerar o relatório de cobertura a partir dos arquivos `{prefix}.regions.bed.gz` e `{prefix}.thresholds.bed.gz`.

### O que é feito:

- **Profundidade média:** média de vezes que cada base foi sequenciada.
- **Cobertura ≥10x e ≥30x:** indicam regiões com cobertura adequada para detectar variantes comuns e raras, respectivamente.
- **Gráficos gerados:**
  - Histograma em escala logarítmica para visualizar melhor a distribuição da cobertura.
  - Heatmap mostrando a cobertura ao longo dos cromossomos.

## Rule 6: Inferência do Sexo Genético
Este passo executa o script `scripts/infer_sex.py`

A inferência foi feita com base em três pilares:

1. **Razões de cobertura:**  
   Calculou-se a cobertura média dos cromossomos **chr1** (autossômico), **chrX** e **chrY**.  
   Razões X/1 próximas de 2 indicam perfil feminino (XX); razões Y/1 acima de 0.1 indicam masculino (XY).

2. **Visualização (boxplot):**  
   Gráfico da cobertura nos cromossomos chr1, chrX e chrY em escala logarítmica, evidenciando diferenças na cobertura entre os sexos.

3. **Cobertura do gene SRY:**  
   Presença de cobertura no gene SRY, exclusivo do cromossomo Y, confirma perfil masculino; ausência sugere perfil feminino.


## Rule 7: Estimativa de Contaminação com VerifyBamID2

Esta etapa utiliza o software **VerifyBamID2** para estimar a **contaminação do DNA** na amostra analisada.

### Objetivo:
Detectar **misturas indesejadas de DNA** (ex. contaminação cruzada entre amostras), o que pode afetar a qualidade das análises genômicas.

### Entradas principais:
- `BAM`: arquivo de alinhamento da amostra.
- `Reference`: genoma de referência.
- Arquivos auxiliares com PCA e variância baseados no painel 1000 Genomes (pré-processados).

FREEMIX é uma estimativa da proporção de DNA contaminante presente em uma amostra de sequenciamento.

Ele é calculado pelo programa VerifyBamID2, que compara as leituras do seu arquivo BAM/CRAM com um painel de referência de variantes esperadas (como o 1000 Genomes).
### Saída:
- Arquivo `.selfSM` contendo a **fração estimada de contaminação**.

## Rule 8: Relatório da Contaminação
### Interpretação rápida:
- **FREEMIX < 0.02** → sem contaminação significativa.
- **FREEMIX > 0.05** → atenção: contaminação pode comprometer resultados.


# ✅ Pré-requisitos

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) instalado
- Git 

---

# 🛠️ Rodar o programa

1. Clonar o repositório:

```bash
git clone https://github.com/seu_usuario/wes_qc_pipeline_LGCM.git
cd wes_qc_pipeline_LGCM.git
```

2. Criar o ambiente Conda com Snakemake
```bash
conda env create -f snakemake_base_env.yaml
conda activate snakemake_base
```

3. Baixar os arquivos de entrada na pasta data
```bash
chmod +x ./scripts/download_data.sh
./scripts/download_data.sh
```

4. Chamar a Pipeline
```bash
snakemake --use-conda --conda-frontend mamba --cores 4
```