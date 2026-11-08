import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

# The parser will automatically set the boundaries around european samples
# of the input dataset, based on the reference panel 1000 genomes
# These boundaries will be relaxed with margin of 15% fraction wider than
# the delimitation within the reference panel
 
parser = argparse.ArgumentParser(
    description="Plot PC1/PC2 and define the EUR boundaries from the 1000G reference panel")
parser.add_argument("--marginX", type=float, default=0.3,
                    help="Margin at x-axis around the EUR samples, as a fraction of the EUR range: "
                         "0.3 = box 30%% wider on each side (default: 0.3)")
parser.add_argument("--marginY", type=float, default=0.3,
                    help="Margin at y-axis around the EUR samples, as a fraction of the EUR range: "
                         "0.30 = box 30%% wider on each side (default: 0.3)")

parser.add_argument("--out_bounds", default="EUR_boundaries.txt")
args = parser.parse_args()

def eur_boundaries(df, marginX, marginY):
    """Return (xmin, xmax, ymin, ymax) of the EUR reference samples,
    extended on each side by a fraction (margin) of the EUR range."""
    eur = df[df["race"] == "EUR"]
    if eur.empty:
        sys.exit("ERROR: no EUR samples found in the reference panel")
 
    xmin, xmax = eur["PC1"].min(), eur["PC1"].max()
    ymin, ymax = eur["PC2"].min(), eur["PC2"].max()
 
    mx = marginX * (xmax - xmin)
    my = marginY * (ymax - ymin)
 
    return xmin - mx, xmax + mx, ymin - my, ymax + my


# load data
pca = pd.read_table("MDS_merge3_projected.sscore",sep="\t")
race = pd.read_table("racefile.txt", sep=r"\s+", header=0)


# keep only need col
pca.columns = pca.columns.str.replace("#", "")
pca = pca[["FID", "IID", "PC1_AVG", "PC2_AVG"]]
pca.columns = ["FID", "IID", "PC1", "PC2"]

mergedata = pd.merge(pca, race, on=["IID","FID"], how="inner")
mergedata = mergedata[mergedata["race"].isin(["EUR", "ASN", "AMR", "AFR", "OWN"])]

#print(mergedata.head())

# EUR boundaries
xmin, xmax, ymin, ymax = eur_boundaries(mergedata, args.marginX, args.marginY)
 
with open(args.out_bounds, "w") as f:
    f.write("xmin xmax ymin ymax\n")
    f.write(f"{xmin} {xmax} {ymin} {ymax}\n")


plt.figure(figsize=(8,8))

ax = sns.scatterplot(
    data=mergedata,
    x="PC1",
    y="PC2",
    hue="race",
    s=40,
    palette={
        "EUR": "green",
        "ASN": "red",
        "AMR": "orange",
        "AFR": "blue",
        "OWN": "black"
    }
)

# Change only the legend label
handles, labels = ax.get_legend_handles_labels()
labels = ["WTCCC" if label == "OWN" else label for label in labels]
ax.legend(handles, labels, fontsize=22)

ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin,
                       fill=False, linestyle="--", edgecolor="grey", linewidth=1.5))

plt.xlabel("PC1", fontsize=22) #16
plt.ylabel("PC2", fontsize=22)
plt.title("PCA", fontsize=22)

plt.xticks(fontsize=20) #14
plt.yticks(fontsize=20)


#plt.legend(fontsize=22) #13
plt.tight_layout()
plt.savefig("PCA_plot.png", dpi=150)

