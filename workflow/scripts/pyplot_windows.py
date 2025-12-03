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
import matplotlib.pyplot as plt

# Read windows file
windows = pd.read_csv(
    snakemake.input.scores,
    sep=r'\s+',
    header=None,
    names=['start', 'end', 'n_snps', 'frac_extreme', 'percentile', 'max_score']
)

# Filter valid windows only (those norm considers reliable)
windows['max_score_num'] = pd.to_numeric(windows['max_score'], errors='coerce')
valid = windows[(windows['percentile'] != -1) & (windows['max_score_num'].notna())].copy()

# Plot using window start position
plt.figure(figsize=(12, 4))
plt.plot(valid['start'], valid['max_score_num'], linewidth=0.5, alpha=0.8)
plt.title(snakemake.params.title)
plt.xlabel('Chromosome Position')
plt.ylabel('Max Score per Window')
plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
