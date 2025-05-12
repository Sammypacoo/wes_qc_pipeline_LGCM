### Dados de entrada

Os arquivos necessários foram baixados manualmente dos seguintes links:

- [Arquivo CRAM](http://ftp.1000genomes.ebi.ac.uk/.../NA06994.exome.cram)
- [Arquivo CRAI](http://ftp.1000genomes.ebi.ac.uk/.../NA06994.exome.cram.crai)
- [BED das regiões exônicas](https://www.twistbioscience.com/...)


### Hashes MD5 esperados:

- `NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram`  
  `3d8d8dc27d85ceaf0daefa493b8bd660`

- `NA06994.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram.crai`  
  `15a6576f46f51c37299fc004ed47fcd9`

- Arquivo BED (Twist Exome):  
  `c3a7cea67f992e0412db4b596730d276`

### Verificação 


## ✅ Pré-requisitos

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) instalado
- Git 

---

## 🛠️ Rodar o programa

1. Clonar o repositório:

```bash
git clone https://github.com/seu_usuario/wes_qc_pipeline_LGCM.git
cd seu_pipeline
```

2. Criar o ambiente Conda com Snakemake
```bash
conda env create -f snakemake_base_env.yaml
conda activate snakemake_base
```

3. Chamr a Pipeline
```bash
snakemake --use-conda --conda-frontend mamba --cores 4
```