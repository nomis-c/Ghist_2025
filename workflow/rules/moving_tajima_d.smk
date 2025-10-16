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


rule moving_tajima_d:
    input:
        vcf=config["vcf_file"]
    output:
        scores="results/tajima_d/{dataset}/chr{i}.{window}snps_{step}.moving_tajima_d.scores.txt"
    params:
        window_size=lambda wildcards: int(wildcards.window),
        step_size_ratio=lambda wildcards: float(wildcards.step),
        dataset=config["dataset"]
    log:
        "logs/tajima_d/{dataset}/chr{i}.{window}snps_{step}.moving_tajima_d.log"
    script:
        "../scripts/tajima_d.py"


rule format_moving_tajima_d_manhattan:
    input:
        scores=rules.moving_tajima_d.output.scores
    output:
        formatted="results/tajima_d/{dataset}/chr{i}.{window}snps_{step}.moving_tajima_d.manhattan.txt"
    params:
        chrom="{i}"
    log:
        "logs/tajima_d/{dataset}/chr{i}.{window}snps_{step}.moving_tajima_d.format_manhattan.log"
    shell:
        """
        awk -v chr="{params.chrom}" 'BEGIN{{OFS="\\t"}}
             NR==1{{print "SNP", "CHR", "BP", "window_start", "window_end", "tajima_d"}}
             NR>1 && $3<0 {{print chr":"$1, chr, $1, $1, $2, $3}}' \
        {input.scores} > {output.formatted} 2> {log}
        """

rule plot_moving_tajima_d_manhattan:
    input:
        scores=rules.format_moving_tajima_d_manhattan.output.formatted
    output:
        candidates="results/plots/moving_tajima_d/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.top_candidates.txt",
        plot="results/plots/moving_tajima_d/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.manhattan.png"
    params:
        score_column="tajima_d",
        cutoff=lambda wildcards: float(wildcards.cutoff),
        use_absolute="FALSE",
        width=640,
        height=240,
        color1="#56B4E9",
        color2="#F0E442"
    log:
        "logs/plots/moving_tajima_d/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.manhattan.log"
    script:
        "../scripts/manhattan.R"


rule pyplot_moving_tajima_d:
    input:
        scores=rules.format_moving_tajima_d_manhattan.output.formatted
    output:
        plot="results/plots/moving_tajima_d/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.pyplot.png"
    params:
        score_column="tajima_d",
        chrom="{i}",
        title=lambda wildcards: f"Moving Tajima's D - Chr {wildcards.i} (window size: {wildcards.window} SNPs, step size: {int(int(wildcards.window)*float(wildcards.step))} SNPs)"
    log:
        "logs/plots/moving_tajima_d/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.pyplot.log"
    script:
        "../scripts/pyplot.py"


rule create_moving_tajima_d_bed:
    input:
        scores=rules.moving_tajima_d.output.scores
    output:
        bed="results/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.selection_regions.bed"
    params:
        cutoff=lambda wildcards: float(wildcards.cutoff),
        chrom="{i}"
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.create_bed.log"
    script:
        "../scripts/create_tajima_d_bed.py"


rule extract_single_moving_tajima_region:
    input:
        scores=rules.moving_tajima_d.output.scores
    output:
        bed="results/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.single_region.bed"
    params:
        chrom="{i}",
        format="moving_tajima",
        cutoff=lambda wildcards: float(wildcards.cutoff)
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.extract_single.log"
    script:
        "../scripts/extract_strongest_region.py"


rule merge_moving_tajima_d_bed:
    input:
        bed=rules.create_moving_tajima_d_bed.output.bed
    output:
        merged="results/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.selection_regions_merged.bed"
    log:
        "logs/selection_candidates/{dataset}/chr{i}.{window}snps_{step}_{cutoff}.moving_tajima_d.merge_bed.log"
    shell:
        """
        bedtools sort -i {input.bed} | bedtools merge > {output.merged} 2> {log}
        """
