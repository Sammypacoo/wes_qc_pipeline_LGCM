import gzip
import csv
import polars as pl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
import os

conteudo = """

###  Cobertura com escala logarítmica

A escala logarítmica **comprime valores altos** e **expande os valores baixos**. Isso ajuda a visualizar melhor **regiões de baixa cobertura**, mas pode dificultar a interpretação exata dos números sem prestar atenção ao eixo.

No gráfico acima, o eixo X está em `log10` (base 10), o que significa:

- `10^0` = **1x**
- `10^1` = **10x**
- `10^1.5` ≈ **32x**
- `10^2` = **100x**

---

### Como comparar visualmente

- Identifique onde estão os valores de **10x e 30x** no eixo logarítmico:
  - `10x` → `log10(10)` = **1**
  - `30x` → `log10(30)` ≈ **1.48**

- Trace linhas imaginárias verticais nesses pontos para observar quantas regiões estão abaixo ou acima desses limiares.

- **Interpretação prática do gráfico:**
  - Se houver um **grande número de regiões à esquerda de 1 (10x)**, há **muita cobertura baixa**, o que pode comprometer a análise.
  - Se o **pico do histograma está à direita de 1.5 (30x)**, significa que a **maioria das regiões tem cobertura adequada ou ótima**.
"""




def mean_coverage(file):
    data = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            try:
                coverage = float(fields[-1])
            except ValueError:
                continue
            data.append((fields[0], fields[1], fields[2], coverage))
    df = pl.DataFrame(data, schema=["chrom", "start", "end", "coverage"])
    mean_cov = df.select(pl.col("coverage").mean()).item()
    return round(mean_cov, 2)

def thresholds_percent(file):
    data = []
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                co_10x = float(fields[-2])
                co_30x = float(fields[-1])
                data.append((chrom, start, end, co_10x, co_30x))
            except ValueError:
                continue
    df = pl.DataFrame(data, schema=["chrom", "start", "end", "co_10x", "co_30x"])
    df = df.with_columns((pl.col("end") - pl.col("start")).alias("length"))
    pct_10x = round((df["co_10x"].sum() / df["length"].sum()) * 100, 2)
    pct_30x = round((df["co_30x"].sum() / df["length"].sum()) * 100, 2)
    return pct_10x, pct_30x

def plot_histogram(df, output_path):
    sns.histplot(df["coverage"], log_scale=True)
    plt.xlabel("Cobertura")
    plt.ylabel("Número de regiões")
    plt.title("Histograma da cobertura")
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_heatmap(df, output_path):
    df = df.sort_values(["chrom", "start"])
    df["region_bin"] = df.groupby("chrom").cumcount()
    pivot = df.pivot(index="chrom", columns="region_bin", values="coverage").fillna(0)
    plt.figure(figsize=(16, 6))
    sns.heatmap(pivot, cmap="YlGnBu", cbar_kws={"label": "Cobertura (x)"}, norm=LogNorm())
    plt.title("Heatmap da Profundidade de Cobertura por Região")
    plt.xlabel("Regiões ordenadas (bins)")
    plt.ylabel("Cromossomos")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def main():
    regions = snakemake.input.regions
    thresholds = snakemake.input.thresholds
    hist_path = "results/hist_cobertura.png"
    heatmap_path = "results/heatmap_cobertura.png"
    output = snakemake.output.cobertura_report

    # Coleta dados para os gráficos
    data = []
    with gzip.open(regions, 'rt') as f:
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

    # Gera gráficos
    plot_histogram(df, hist_path)
    plot_heatmap(df, heatmap_path)

    # Cálculo de métricas
    mean_cov = mean_coverage(regions)
    pct_10x, pct_30x = thresholds_percent(thresholds)
   
    '''
    with open(output, "w") as f:
       

        f.write("##  Interpretação das Métricas\n")
        f.write("- A profundidade média representa a média de vezes que cada base foi sequenciada. Valores acima de 30x são geralmente considerados confiáveis para análise.\n")
        f.write("- A porcentagem de cobertura ≥10x e ≥30x indica o quanto do ou exoma alvo foi suficientemente sequenciado.\n")
        f.write("  - ≥10x: adequado para detecção de variantes comuns.\n")
        f.write("  - ≥30x: adequado para detecção precisa de variantes raras.\n\n")
        
        f.write(f"{conteudo}\n\n")
        
        f.write("# Relatório de Cobertura dos resultados analisados \n\n")
        f.write(f"**Profundidade média:** {mean_cov}x\n\n")
        f.write(f"**Cobertura ≥10x:** {pct_10x}%\n")
        f.write(f"**Cobertura ≥30x:** {pct_30x}%\n\n")

        f.write("# Histograma da Cobertura\n")
        f.write("Este gráfico mostra a distribuição da profundidade de cobertura por região. A escala logarítmica permite visualizar melhor as regiões com baixa e alta cobertura.\n\n")
        f.write(f'<img src="{os.path.basename(hist_path)}" alt="Histograma de cobertura" width="600">\n\n')

        f.write("# Heatmap de Cobertura\n")
        f.write("O heatmap exibe a cobertura ao longo dos cromossomos. Cada linha representa um cromossomo e cada coluna representa uma região ordenada. Valores mais claros indicam maior cobertura.\n\n")
        f.write(f'<img src="{os.path.basename(heatmap_path)}" alt="Heatmap de cobertura" width="1500">\n\n')
      '''
     
    with open(output, "w") as f:
            f.write("## 📌 Interpretação das Métricas\n")
            f.write("- A **profundidade média** representa a média de vezes que cada base foi sequenciada.\n")
            f.write("  - ✅ Valores acima de 30x são geralmente considerados confiáveis para análise.\n")
            f.write("- A **porcentagem de cobertura ≥10x e ≥30x** indica o quanto do exoma alvo foi suficientemente sequenciado:\n")
            f.write("  - 🧪 ≥10x: adequado para detecção de variantes comuns.\n")
            f.write("  - 🧬 ≥30x: adequado para detecção precisa de variantes raras.\n\n")


            f.write("# ✅ **Resumo dos Resultados de Cobertura**\n\n")
            f.write(f"- 🔬 **Profundidade média:** **{mean_cov}x**\n")
            f.write(f"- 🧪 **Cobertura ≥10x:** **{pct_10x}%** (boa para variantes comuns)\n")
            f.write(f"- 🧬 **Cobertura ≥30x:** **{pct_30x}%** (essencial para variantes raras)\n\n")

            f.write("> ⚠️ Regiões com cobertura abaixo de 10x podem comprometer a detecção de variantes.\n")
            f.write("> ✅ A maioria das regiões com ≥30x indica uma boa confiabilidade dos dados.\n\n")

            f.write("## 📊 Visualização dos Resultados\n\n")
            f.write("---\n\n")

            f.write("### 📉 Histograma da Cobertura (escala logarítmica)\n")
            f.write("Este gráfico mostra a distribuição da profundidade de cobertura por região. A escala logarítmica permite visualizar melhor as regiões com baixa e alta cobertura.\n\n")
            f.write("**Entendendo a escala logarítmica (base 10):**\n")
            f.write("- 10^0 = 1x\n")
            f.write("- 10^1 = 10x\n")
            f.write("- 10^1.5 ≈ 32x\n")
            f.write("- 10^2 = 100x\n\n")
            f.write("🧠 Para interpretar:\n")
            f.write("- 🔸 10x → log10(10) = 1\n")
            f.write("- 🔸 30x → log10(30) ≈ 1.48\n")
            f.write("- Se muitas regiões estão à esquerda de 1 → baixa cobertura.\n")
            f.write("- Se o pico está à direita de 1.5 → cobertura adequada.\n\n")
            f.write(f'<img src="{os.path.basename(hist_path)}" alt="Histograma de cobertura" width="600">\n\n')

            f.write("# 🔥 Heatmap de Cobertura\n")
            f.write("O heatmap exibe a cobertura ao longo dos cromossomos. Cada linha representa um cromossomo e cada coluna representa uma região ordenada.\n")
            f.write("🟦 **Regiões mais escuras representam melhor cobertura**, enquanto 🟥 **regiões mais claras indicam baixa cobertura ou ausência de dados.** Isso ajuda a identificar áreas mal sequenciadas que podem comprometer a análise.\n\n")
            f.write(f'<img src="{os.path.basename(heatmap_path)}" alt="Heatmap de cobertura" width="1200">\n\n')



if __name__ == "__main__":
    main()
