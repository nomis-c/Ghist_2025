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


import pandas as pd

# Get inputs
input_file = snakemake.input.scores
output_bed = snakemake.output.bed
chrom = snakemake.params.chrom
file_format = snakemake.params.format

if file_format == 'windowed_tajima':
    # Read Tajima's D scores
    df = pd.read_csv(input_file, sep='\t')
    df['abs_score'] = df['tajima_d'].abs()
    start_col = 'window_start'
    end_col = 'window_end'
    strongest_idx = df['tajima_d'].idxmin()

elif file_format == 'moving_tajima':
    df = pd.read_csv(input_file, sep='\t')
    df['abs_score'] = df['tajima_d'].abs()
    start_col = 'window_start'
    end_col = 'window_end'
    strongest_idx = df['tajima_d'].idxmin()


elif file_format == 'selscan':
    # Read selscan windows
    df = pd.read_csv(
        input_file,
        sep='\s+',
        header=None,
        names=['start', 'end', 'n_snps', 'frac_extreme', 'percentile', 'max_score']
    )
    # Filter out invalid windows
    df['max_score_num'] = pd.to_numeric(df['max_score'], errors='coerce')
    df = df[(df['percentile'] != -1) & (df['max_score_num'].notna())].copy()
    df['abs_score'] = df['max_score_num'].abs()
    strongest_idx = df['abs_score'].idxmax()
    start_col = 'start'
    end_col = 'end'


# Get the region with highest absolute score
strongest = df.loc[strongest_idx]

# Write single-line BED
with open(output_bed, 'w') as f:
    f.write(f"{chrom}\t{int(strongest[start_col])}\t{int(strongest[end_col])}\n")
