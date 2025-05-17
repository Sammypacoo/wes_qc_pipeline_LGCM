import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Caminho para o arquivo .selfSM
selfsm_path = snakemake.input.selfsm
report_path = snakemake.output.contaminacao_report

# Ler o arquivo
df = pd.read_csv(selfsm_path, sep='\t')

# Pegar a primeira amostra (ou loop depois)
sample = df.iloc[0]
seq_id = sample["#SEQ_ID"]
freemix = sample["FREEMIX"]
avg_dp = sample["AVG_DP"]
snps = sample["#SNPS"]
reads = sample["#READS"]

# Interpretar FREEMIX
if freemix < 0.005:
    status = "✅ Sem evidência de contaminação"
elif freemix < 0.03:
    status = "⚠️ Potencial leve contaminação"
else:
    status = "❌ Contaminação significativa"


# Criar conteúdo do relatório em Markdown
md_content = f"""
# Relatório de Contaminação – VerifyBamID2
---

## 📊 Tabela de Referência – Interpretação do FREEMIX

| Intervalo FREEMIX | Interpretação                     |
|------------------:|-----------------------------------|
| `< 0.005`         | ✅ Sem evidência de contaminação   |
| `0.005 – 0.03`    | ⚠️ Leve contaminação               |
| `> 0.03`          | ❌ Alta contaminação               |

---
**Amostra:** `{seq_id}`  
**Total de SNPs:** {snps}  
**Número de Leituras:** {reads}  
**Profundidade Média:** {avg_dp:.2f}  
**FREEMIX (estimativa de contaminação):** `{freemix:.6f}`  
**Status de contaminação:** {status}

"""

# Salvar como arquivo .md
with open("report_path", "w") as f:
    f.write(md_content)

print("✅ Relatório Markdown salvo como 'report_verifybamid2.md'")
