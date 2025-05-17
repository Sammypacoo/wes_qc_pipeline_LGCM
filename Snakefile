# Snakefile
configfile: "config.yaml"

rule all:
    input:
        "logs/scripts_ready.done",
        "logs/md5_check.done",
        "data/sample.bam",
        "results/coverage_stats.txt",
        "results/sex_inference.txt",
        "results/contamination.txt"

# Rule 1:  make scripts executable
rule make_scripts_executable:
    input:
        expand("scripts/{script}", script=["md5_check.sh","add_percent_column.sh"])
    output:
        touch("logs/scripts_ready.done")
    shell:
        """
        chmod +x scripts/*.sh
        touch {output}
        """

# Rule 2: Verificação dos arquivos via md5sum
rule check_md5:
    input:
        cram=config["cram"],
        crai=config["crai"],
        bed=config["bed"]
    output:
        done="logs/md5_check.done"
    log:
        "logs/md5_check.log"
    shell:
        """
        bash scripts/md5_check.sh > {log} 2>&1 && touch {output.done}
        """
# Rule 3: Conversão de CRAM para BAM
rule convert_cram:
    input:
        cram=config["cram"],
        crai=config["crai"],
        ref=config["reference"]
    output:
       bam= config["bam"],
       bai= config["bai"]   
    conda:
        "envs/env.yml"
    shell:
        """
        samtools view -T {input.ref} -b -o {output.bam} {input.cram}
        samtools index  {output.bam}

        """

# Rule 4: Cálculo de cobertura com mosdepth
rule calculate_coverage:
    input:
        bam= config["bam"],
        bed=config["bed"]
    output:
        "results/cobertura.thresholds.bed.gz",
        "results/cobertura.regions.bed.gz"
    conda:
        "envs/env.yml"
    shell:
        """
        mosdepth \
            --by {input.bed} \
            --thresholds 10,30 \
            --threads 4 \
             cobertura \
            {input.bam}
        """

# Add uma coluna com a porcentagem de 10x e 30x
rule add_percent_column:
    input:
        cobertura_10_30="results/cobertura.thresholds.bed.gz"
    output:
        "results/cobertura.thresholds_percentual.tsv"

    shell:
        """
        bash scripts/add_percent_column.sh 
        """

rule cobertura_report:
    input:
        "results/cobertura.thresholds.bed.gz",
        "results/cobertura.regions.bed.gz"
    output:
        "results/report.md"
    conda:
        "envs/report_env.yaml"
    shell:
        "python scripts/generate_report.py"

# Inferência de sexo com script Python
rule infer_sex:
    input:
        regions="results/cobertura.regions.bed.gz",
    output:
        report="results/relatorio_sex.md"
    conda:
        "envs/env.yml"
    script:
        "scripts/infer_sex.py"

# Estimativa de contaminação com verifyBamID2
rule contamination:
    input:
        bam="data/sample.bam",
        vcf=config["vcf"]
    output:
        "results/contamination.txt"
    shell:
        """
        verifybamid2 --bam {input.bam} --vcf {input.vcf} --out results/contam && \
        mv results/contam.selfSM {output}
        """
