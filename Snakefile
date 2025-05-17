# Snakefile
configfile: "config.yaml"

rule all:
    input:
        cobertura_report=config["cobertura_report"],
        sex_report=config["sex_report"],
        selfsm=config["selfsm"]

# Rule 1:  make scripts executable
rule make_scripts_executable:
    input:
        expand("scripts/{script}", script=["md5_check.sh","cobertura_report.py","infer_sex.py"])
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
        scripts/md5_check.sh > {log} 2>&1 && touch {output.done}
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
    log:
        "logs/convert_cram.log" 
    conda:
        "envs/env.yaml"
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
        regions= config["regions"],
        thresholds=config["thresholds"]
    conda:
        "envs/env.yaml"
    log:
        "logs/calculate_coverage.log"
    shell:
        """
        mosdepth \
            --by {input.bed} \
            --thresholds 10,30 \
            --threads 4 \
               ./results/cobertura \
            {input.bam}
        """
# Rule 5: Relatório da cobertura
rule cobertura_report:
    input:
        regions= config["regions"],
        thresholds=config["thresholds"]
    output:
        cobertura_report=config["cobertura_report"]
    conda:
        "envs/env.yaml"
    shell:
        "python scripts/cobertura_report.py"

# Rule 6 :Inferência de sexo com script Python
rule infer_sex:
    input:
        regions= config["regions"],
    output:
        sex_report=config["sex_report"]
    conda:
        "envs/env.yaml"
    script:
        "scripts/infer_sex.py"

# Rule 7 :Estimativa de contaminação com verifyBamID2
rule contamination:
    input:
        bam=config["bam"],
        ref=config["reference"],
        vb_UD=config["vb_UD"],
        vb_bed=config["vb_bed"],
        vb_mu=config["vb_mu"],
        vb_V=config["vb_V"]
    output:
        selfsm=config["selfsm"]

    log:
        "logs/verifybamid2.log"

    shell:
        """
        verifybamid2 \
        --BamFile {input.bam} \
        --Reference  {input.ref}\
        --SVDPrefix /data/1000g.phase3.10k.b38.exome.vcf.gz.dat \
        --Output /results/NA06994_verifybamid2_result

        """

# Rule 8: Relatório da contaminação
rule contaminacao_report:
    input:
        selfsm=config["selfsm"]
    output:
        contaminacao_report=config["contaminacao_report"]
    conda:
        "envs/env.yaml"
    shell:
        "python scripts/contaminacao_report.py"



