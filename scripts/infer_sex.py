import gzip
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import natsort
import sys

# Usando input/output do Snakemake
bed_path = snakemake.input.regions
report_path = snakemake.output.sex_report
plot_path = report_path.replace(".md", ".png")

# Leitura e parsing
data = []
with gzip.open(bed_path, 'rt') as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if not row or row[0].startswith("#") or len(row) < 4:
            continue
        try:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            coverage = float(row[-1])
            data.append((chrom, start, end, coverage))
        except ValueError:
            continue

df = pd.DataFrame(data, columns=["chrom", "start", "end", "coverage"])
df["chrom"] = pd.Categorical(df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True)
df = df.sort_values(["chrom", "start"])

mean_cov = df.groupby("chrom")["coverage"].mean()
cov_chr1 = mean_cov.get("chr1", 1)
cov_chrX = mean_cov.get("chrX", 0)
cov_chrY = mean_cov.get("chrY", 0)
ratio_X = cov_chrX / cov_chr1
ratio_Y = cov_chrY / cov_chr1

sexo = "Masculino (XY)" if ratio_Y > 0.1 else "Feminino (XX)"
interpretacao_cobertura = (
    "A razão X/1 próxima de 2.0 e a ausência de cobertura no cromossomo Y indicam um perfil feminino (XX). "
    "No boxplot, observa-se que o cromossomo X apresenta uma mediana de cobertura similar ou superior ao cromossomo 1, "
    "enquanto o cromossomo Y apresenta cobertura desprezível, reforçando o perfil XX."
    if sexo == "Feminino (XX)"
    else
    "A presença de cobertura no cromossomo Y e a razão X/1 próxima de 1.0 indicam um perfil masculino (XY). "
    "No boxplot, a cobertura do cromossomo X está próxima da de chr1, e o cromossomo Y apresenta cobertura clara, "
    "reforçando a presença do cromossomo Y e o perfil XY."
)

# Gera gráfico
df_filtered = df[df['chrom'].isin(['chr1', 'chrX', 'chrY'])].copy()
df_filtered['chrom'] = pd.Categorical(df_filtered['chrom'], categories=['chr1', 'chrX', 'chrY'], ordered=True)

plt.figure(figsize=(8, 5))
sns.boxplot(data=df_filtered, x='chrom', y='coverage', palette='pastel')
plt.yscale("log")
plt.ylabel("Cobertura (log scale)")
plt.title("Distribuição da cobertura por cromossomo")
plt.grid(True, axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()
plt.savefig("/home/venus/mar/alunos/slucciola/teste_LGCM/relatorio_sex.png")
plt.close()

# Abrir e ler arquivo BED com pandas, assumindo que há cabeçalho e colunas suficientes
with gzip.open("cobertura.regions.bed.gz", 'rt') as f:
    df = pd.read_csv(f, sep='\t', header=None, comment='#')

# Ver as primeiras colunas para entender a estrutura
print(df.head())

# Se você encontrou o nome "SRY" na coluna 3 (índice 3), por exemplo:
sry_df = df[df[3].astype(str).str.contains("SRY", case=False, na=False)]
sry_df=sry_df.iloc[0,-1]
sry_df=round(sry_df,2)
## Imprimir o resultado do SRY
# Cobertura do gene SRY (exemplo)

comentario_sry = (
    "A ausência de cobertura nessa região é compatível com indivíduos do sexo feminino."
    if sry_df < 1
    else "A presença de cobertura nessa região é compatível com indivíduos do sexo masculino."
)

# Escreve no .md
with open(report_path, "w") as f:
    f.write("## Inferência de Sexo Genético\n\n")
    f.write(f"- Cobertura média:\n")
    f.write(f"A partir da análise da cobertura média por cromossomo, foi realizada a inferência do sexo genético com base na razão de cobertura dos cromossomos sexuais (X e Y) em relação ao cromossomo autossômico chr1.\n")
    f.write(f"  - chr1: {cov_chr1:.2f}\n")
    f.write(f"  - chrX: {cov_chrX:.2f}\n")
    f.write(f"  - chrY: {cov_chrY:.2f}\n\n")
    f.write(f"- Razão X/1: **{ratio_X:.2f}**, Razão Y/1: **{ratio_Y:.2f}**\n\n")
    f.write(f"  - Interpretação: {interpretacao_cobertura} \n")
    f.write(f"![Boxplot de cobertura]({'/home/venus/mar/alunos/slucciola/teste_LGCM/relatorio_sex.png'.split('/')[-1]})\n")
    f.write(f"- Verificação adicional: **Gene SRY** \n")
    f.write(f", marcador exclusivo do cromossomo Y, foi utilizado como validação adicional da inferência.\n")
    f.write(f"  - Cobertura média sobre SRY: **{sry_df}x**.\n")
    f.write(f" {comentario_sry}\n\n")
    f.write(f"**Inferência: {sexo}**\n\n")