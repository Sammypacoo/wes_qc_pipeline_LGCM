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

# make scrripts executable
rule make_scripts_executable:
    input:
        expand("scripts/{script}", script=["md5_check.sh"])
    output:
        touch("logs/scripts_ready.done")
    shell:
        """
        chmod +x scripts/*.sh
        touch {output}
        """

# Verificação dos arquivos via md5sum
rule check_md5:
    input:
        "data/md5_checksums.txt",
        "data/{cram}",
        "data/{crai}",
        "data/{bed}"
    output:
        "logs/md5_check.done"
    shell:
         """
        bash scripts/md5_check.sh && touch {output}
        """
# Conversão de CRAM para BAM
rule convert_cram:
    input:
        cram="data/{cram}",
        crai="data/{crai}",
        ref=config["reference"]
    output:
        "data/sample.bam"
    shell:
        """
        samtools view -T {input.ref} -b -o {output} {input.cram}
        """

# Cálculo de cobertura com bedtools
rule calculate_coverage:
    input:
        bam="data/sample.bam",
        bed="data/{bed}"
    output:
        "results/coverage_stats.txt"
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.bam} > {output}
        """

# Inferência de sexo com script Python
rule infer_sex:
    input:
        bam="data/sample.bam"
    output:
        "results/sex_inference.txt"
    shell:
        """
        python scripts/infer_sex.py {input.bam} > {output}
        """

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
