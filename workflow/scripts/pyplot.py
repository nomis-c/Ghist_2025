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

# Load data
data = pd.read_csv(snakemake.input.scores, sep='\t')
data = data.sort_values('BP')

# Plot
plt.figure(figsize=(12, 4))
plt.plot(data['BP'], data[snakemake.params.score_column], linewidth=0.5, alpha=0.8)
plt.title(snakemake.params.title)
plt.xlabel('Chromosome Position')
plt.ylabel(snakemake.params.score_column)
plt.savefig(snakemake.output.plot, dpi=300, bbox_inches='tight')
plt.close()
