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

# Read and filter windows
windows = pd.read_csv(
    snakemake.input.scores,
    sep='\s+',
    header=None,
    names=['start', 'end', 'n_snps', 'frac_extreme', 'percentile', 'max_score']
)

# Convert max_score to numeric, filter out NA and unreliable windows
windows['max_score_num'] = pd.to_numeric(windows['max_score'], errors='coerce')
valid = windows[(windows['percentile'] != -1) & (windows['max_score_num'].notna())].copy()

# Get absolute max_score for sorting
valid['abs_max_score'] = valid['max_score_num'].abs()

# Sort by absolute max_score
valid_sorted = valid.sort_values('abs_max_score', ascending=False)

n_top = max(1, int(len(valid_sorted) * float(snakemake.params.cutoff)))
top = valid_sorted.head(n_top)

# Create and write BED
bed = pd.DataFrame({
    'chr': snakemake.params.chrom,
    'start': top['start'].astype(int),
    'end': top['end'].astype(int),
}).sort_values(['chr', 'start'])

bed.to_csv(snakemake.output.bed, sep='\t', header=False, index=False)
