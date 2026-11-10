import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# load data
df = pd.read_csv("r2_maf_imputed.txt", sep="\t")


# some variants are not imputed, remove
df = df.dropna(subset=["R2", "MAF"])
df = df[df["R2"] > 0]

# ---- PLOT R2 vs MAF
# create MAF bins, and compute the R2 mean for each bin
bins = np.logspace(np.log10(1e-4), np.log10(0.5), 150)
df["MAF_bin"] = pd.cut(df["MAF"], bins=bins, include_lowest=True, right=True)

# compute mean R2 per bin
# get bin centers
binned = df.groupby("MAF_bin", observed=True)["R2"].mean().reset_index()
binned["MAF_center"] = binned["MAF_bin"].apply(lambda x: x.mid)

# plot
figsize=(7,5)
plt.plot(binned["MAF_center"], binned["R2"])

# custom ticks spacing (as in the paper https://pmc.ncbi.nlm.nih.gov/articles/PMC7973508/#Abs1)
plt.xscale("log")   
plt.xlim(1e-4, 5e-1)
ticks = [1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3,1e-2, 5e-2, 1e-1, 5e-1]
# Then update the xticks label
ax = plt.gca()
ax.set_xticks(ticks)
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:.0e}"))
plt.xticks(ticks, rotation=45)

plt.ylim(0, 1.1)

plt.xlabel("Binned MAF", fontsize=16)
plt.ylabel("Mean R2", fontsize=16)
plt.title("R2 vs MAF", fontsize=16)
plt.show()
plt.grid()

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

plt.savefig("R2_vs_MAF_imputation.png", dpi=120)
plt.close()

# ------ PLOT FREQUENCE R2
figsize=(7,5)

plt.hist(df["R2"], bins=120)

plt.xlabel("R2", fontsize=16)
plt.ylabel("Frequency", fontsize=16)
plt.show()
plt.grid()

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.tight_layout()

plt.savefig("R2_freq_imputation.png", dpi=120)
plt.close()

