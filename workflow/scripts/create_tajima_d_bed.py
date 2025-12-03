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

# Snakemake inputs
scores_file = snakemake.input.scores
output_bed = snakemake.output.bed
cutoff = float(snakemake.params.cutoff)
chrom = snakemake.params.chrom

# Read scores
df = pd.read_csv(scores_file, sep='\t')

df_sorted = df.sort_values('tajima_d', ascending=True)

# Get top X% most negative candidates
n_top = max(1, int(len(df_sorted) * cutoff))
top_candidates = df_sorted.head(n_top)

# Create BED format
bed_data = []
for idx, row in top_candidates.iterrows():
    bed_data.append([
        chrom,
        int(row['window_start']),
        int(row['window_end']),
    ])

# Sort by genomic position
bed_df = pd.DataFrame(bed_data, columns=['chr', 'start', 'end'])
bed_df = bed_df.sort_values(['chr', 'start'])

# Write BED file
bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
