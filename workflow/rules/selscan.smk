# Copyright 2025 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


wildcard_constraints:
    method="ihs|nsl"

rule download_selscan:
    output:
        selscan="resources/tools/selscan/selscan-2.0.3",
        norm="resources/tools/selscan/norm",
    log:
        "logs/download/download_selscan.log",
    shell:
        """
        rm -rf selscan-2.0.3* >> {log} 2>&1
        wget -c https://github.com/szpiech/selscan/archive/refs/tags/v2.0.3.tar.gz >> {log} 2>&1
        tar -xvzf v2.0.3.tar.gz >> {log} 2>&1
        mv selscan-2.0.3/bin/linux/selscan-2.0.3 {output.selscan} >> {log} 2>&1
        mv selscan-2.0.3/bin/linux/norm {output.norm} >> {log} 2>&1
        rm -rf v2.0.3.tar.gz selscan-2.0.3/ >> {log} 2>&1
        """

rule prepare_vcf_for_selscan:
    input:
        vcf=config["vcf_file"]
    output:
        filtered_vcf="results/selscan_input/{dataset}/chr{i}.filtered.vcf.gz",
        index="results/selscan_input/{dataset}/chr{i}.filtered.vcf.gz.tbi"
    params:
        maf=0.05,
        dataset=config["dataset"]
    log:
        "logs/selscan/{dataset}/chr{i}.prepare_vcf_for_selscan.log"
    shell:
        """
        bcftools view {input.vcf} \
            --min-alleles 2 \
            --max-alleles 2 \
            --types snps \
            --min-af {params.maf}:minor \
            -Oz -o {output.filtered_vcf} 2> {log}

        tabix -p vcf {output.filtered_vcf} 2>> {log}
        """


rule extract_snp_pos:
    input:
        vcf=rules.prepare_vcf_for_selscan.output.filtered_vcf
    output:
        map="results/selscan_input/{dataset}/chr{i}.biallelic.snps.map"
    log:
        "logs/selscan/{dataset}/chr{i}.extract_snp_pos.log"
    shell:
        """
        bcftools query -f "%CHROM\\t%CHROM:%POS:%REF:%ALT\\t%POS\\t%POS\\n" {input.vcf} > {output.map} 2> {log}
        """


rule estimate_selscan_scores:
    input:
        selscan=rules.download_selscan.output.selscan,
        vcf=rules.prepare_vcf_for_selscan.output.filtered_vcf,
        map=rules.extract_snp_pos.output.map
    output:
        out="results/selscan/{dataset}/{method}_{maf}/chr{i}.{method}.out",
        formatted_out="results/selscan/{dataset}/{method}_{maf}/chr{i}.{method}.formatted.out"
    params:
        output_prefix="results/selscan/{dataset}/{method}_{maf}/chr{i}",
    resources:
        cpus=8
    log:
        "logs/selscan/{dataset}/chr{i}.{method}_{maf}.estimate_selscan_scores.log"
    shell:
        """
        {input.selscan} --vcf {input.vcf} --map {input.map} \
            --{wildcards.method} \
            --out {params.output_prefix} \
            --threads {resources.cpus} \
            --maf {wildcards.maf} 2> {log}

        awk -v chr={wildcards.i} 'BEGIN{{OFS="\\t"}} NR==1{{print; next}} {{print chr,$2,$3,$4,$5,$6}}' \
            {output.out} > {output.formatted_out} 2>> {log}
        """


rule normalize_selscan_scores:
    input:
        norm=rules.download_selscan.output.norm,
        scores=rules.estimate_selscan_scores.output.formatted_out
    output:
        normalized="results/selscan/{dataset}/{method}_{maf}/chr{i}.{winsize}bp.{method}_{maf}.normalized",
        windows="results/selscan/{dataset}/{method}_{maf}/chr{i}.{winsize}bp.{method}_{maf}.windows",
        log_file="results/selscan/{dataset}/{method}_{maf}/chr{i}.{winsize}bp.{method}_{maf}.norm.log"
    params:
        winsize=lambda wildcards: int(wildcards.winsize),
        winsize_kb=lambda wildcards: int(wildcards.winsize) // 1000
    log:
        "logs/selscan/{dataset}/chr{i}.{winsize}bp.{method}_{maf}.normalize_selscan_scores.log"
    shell:
        """
        {input.norm} --files {input.scores} \
            --log {output.log_file} \
            --{wildcards.method} \
            --bp-win \
            --winsize {params.winsize} 2>> {log}

        mv {input.scores}.100bins.norm {output.normalized} 2>> {log}
        mv {input.scores}.100bins.norm.{params.winsize_kb}kb.windows {output.windows} 2>> {log}
        """


rule format_selscan_manhattan:
    input:
        scores=rules.normalize_selscan_scores.output.normalized
    output:
        formatted="results/selscan/{dataset}/chr{i}.{winsize}bp.{method}_{maf}.manhattan.txt"
    params:
        chrom="{i}",
        winsize=lambda wildcards: int(wildcards.winsize),
        method="{method}"
    log:
        "logs/selscan/{dataset}/chr{i}.{winsize}bp.{method}_{maf}.format_selscan_manhattan.log"
    shell:
        """
        awk -v chr="{params.chrom}" -v winsize={params.winsize} -v method={params.method} 'BEGIN{{OFS="\\t"}}
             NR==1{{print "SNP", "CHR", "BP", "window_start", "window_end", "normalized_"method; next}}
             {{
                 pos = $2
                 window_start = int(pos / winsize) * winsize
                 window_end = window_start + winsize
                 print $1":"$2, chr, pos, window_start, window_end, $7
             }}' \
        {input.scores} > {output.formatted} 2> {log}
        """

rule plot_selscan:
    input:
        scores=rules.format_selscan_manhattan.output.formatted
    output:
        candidates="results/plots/selscan/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.top_candidates.txt",
        plot="results/plots/selscan/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.manhattan.png"
    params:
        score_column=lambda wildcards: f"normalized_{wildcards.method}",
        use_absolute="TRUE",
        cutoff=lambda wildcards: float(wildcards.cutoff),
        width=640,
        height=240,
        color1="#56B4E9",
        color2="#F0E442"
    log:
        "logs/plots/selscan/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.plot_selscan.log"
    script:
        "../scripts/manhattan.R"


rule pyplot_selscan_windows:
    input:
        scores=rules.normalize_selscan_scores.output.windows
    output:
        plot="results/plots/selscan/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.pyplot_windows.png"
    params:
        chrom="{i}",
        title=lambda wildcards: f"Normalized {wildcards.method.upper()} Windows - Chr {wildcards.i} (window size: {int(wildcards.winsize)/1000:.0f}kb)"
    log:
        "logs/plots/selscan/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.pyplot_windows.log"
    script:
        "../scripts/pyplot_windows.py"


rule create_candidate_bed_windows:
    input:
        scores=rules.normalize_selscan_scores.output.windows
    output:
        bed="results/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.window_selection_regions.bed"
    params:
        chrom="{i}",
        cutoff=lambda wildcards: float(wildcards.cutoff)
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.create_candidate_bed_windows.log"
    script:
        "../scripts/create_selscan_windowed_bed.py"


rule extract_single_selscan_region:
    input:
        scores=rules.normalize_selscan_scores.output.windows
    output:
        bed="results/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.single_region.bed"
    params:
        chrom="{i}",
        format="selscan",
        cutoff=lambda wildcards: float(wildcards.cutoff)
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.extract_single_region.log"
    script:
        "../scripts/extract_strongest_region.py"



rule merge_selscan_bed:
    input:
        bed=rules.create_candidate_bed_windows.output.bed
    output:
        merged="results/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.window_selection_regions_merged.bed"
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{winsize}bp_{cutoff}.{method}_{maf}.merge_selscan_bed.log"
    shell:
        """
        bedtools sort -i {input.bed} | bedtools merge > {output.merged} 2> {log}
        """
