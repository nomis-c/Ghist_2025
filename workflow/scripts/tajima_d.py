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


import allel
import numpy as np
import pandas as pd

# Snakemake inputs
vcf_file = snakemake.input.vcf
output_scores = snakemake.output.scores
window_size = snakemake.params.window_size
step_size_ratio = snakemake.params.step_size_ratio
step_size = int(step_size_ratio * window_size)

# Read VCF
callset = allel.read_vcf(vcf_file, fields=['variants/POS', 'calldata/GT'])
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']

# Calculate allele counts
ac = gt.count_alleles()

# Calculate Tajima's D
d = allel.moving_tajima_d(ac, size=window_size, step=step_size)

window_start = allel.moving_statistic(pos, np.min, size=window_size, step=step_size)
window_end = allel.moving_statistic(pos, np.max, size=window_size, step=step_size)

# Remove NaN values
valid_mask = ~np.isnan(d)
window_start = window_start[valid_mask]
window_end = window_end[valid_mask]
d = d[valid_mask]

# Save results
results_df = pd.DataFrame({
    'window_start': window_start.astype(int),
    'window_end': window_end.astype(int),
    'tajima_d': d
})
results_df.to_csv(output_scores, sep='\t', index=False)
