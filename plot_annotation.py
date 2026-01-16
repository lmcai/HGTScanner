import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
import math, os, argparse


parser = argparse.ArgumentParser(description='Plot HGTScanner annotation and save to pdf.')
parser.add_argument('-f', metavar='tsv_file', help='summary tsv file', required=True)
parser.add_argument('-o', metavar='output', help='putput file name', required=True)
parser.add_argument('-l', metavar='wrapping_len(kb)', help='genome wrapping length in kb', required=True)
args = parser.parse_args()

sum_file = args.f
output_pdf = args.o
wrap_kb = int(args.l)

wrap_bp = wrap_kb * 1000

df = pd.read_csv(
    sum_file,
    sep="\t",
    header=None,
    names=["seqid", "start", "end", "type"]
)

color_map = {
    "VGT": "skyblue",             # vertical gene transfer / baseline
    "mtpt": "forestgreen",        # mitochondria-derived
    "moss HGT": "#FF7F0E",        # orange
    "algae HGT": "#E377C2",       # pink-magenta
    "seed HGT": "#D62728",        # red
    "Hachettea HGT": "#9467BD",   # purple
    "Podocarpaceae HGT": "#2CA02C", # dark green
    "gymnosperm HGT": "#8C564B"   # brown
}

palette = plt.colormaps["tab10"].colors
chroms = df["seqid"].unique()

with PdfPages(output_pdf) as pdf:
    for chrom in chroms:
        chrom_df = df[df["seqid"] == chrom]
        max_pos = chrom_df["end"].max()
        num_rows = math.ceil(max_pos / wrap_bp)
        fig, ax = plt.subplots(figsize=(14, 1.5 * num_rows + 1))
        # Assign colors to unseen annotation types
        all_types = chrom_df["type"].unique()
        for i, t in enumerate(all_types):
            if t not in color_map:
                color_map[t] = palette[i % len(palette)]
        for i in range(num_rows):
            row_start = i * wrap_bp
            row_end = (i + 1) * wrap_bp
            y_pos = num_rows - 1 - i  # top → bottom
            # Background row
            ax.add_patch(
                Rectangle(
                    (0, y_pos),
                    wrap_bp,
                    0.5,
                    color="#f0f0f0",
                    zorder=0
                )
            )
            # Row label
            ax.text(
                -wrap_bp * 0.01,
                y_pos + 0.25,
                f"{row_start // 1000} kb",
                ha="right",
                va="center",
                fontsize=9
            )
            # Features overlapping this row
            mask = (
                (chrom_df["start"] < row_end) &
                (chrom_df["end"] > row_start)
            )
            for _, feat in chrom_df[mask].iterrows():
                f_start = max(feat["start"], row_start)
                f_end = min(feat["end"], row_end)
                x_pos = f_start - row_start
                width = f_end - f_start
                ax.add_patch(
                    Rectangle(
                        (x_pos, y_pos),
                        width,
                        0.5,
                        color=color_map.get(feat["type"], "grey"),
                        edgecolor="black",
                        linewidth=0.4,
                        label=feat["type"]
                    )
                )
        ax.set_xlim(0, wrap_bp)
        ax.set_ylim(-0.5, num_rows)
        ax.set_yticks([])
        ax.set_xlabel("Position (bp)")
        ax.set_title(f"Genome Map: {chrom} (wrap = {wrap_kb} kb)")
        # Deduplicate legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(
            by_label.values(),
            by_label.keys(),
            loc="upper left",
            bbox_to_anchor=(1, 1),
            frameon=False
        )
        plt.tight_layout()
        pdf.savefig(fig)   # <-- ADD PAGE
        plt.close(fig)

print(f"Color-coded genome annotation PDF written to: {output_pdf}")
