import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
import math, os, argparse


parser = argparse.ArgumentParser(description='Plot HGTScanner annotation and save to pdf.')
parser.add_argument('-tsv', metavar='tsv_file', help='summary tsv file')
parser.add_argument('-bed', metavar='bed_file', help='annotation bed file')
parser.add_argument('-o', metavar='output_file', help='putput file name', required=True)
parser.add_argument('-l', metavar='wrapping_len(kb)', help='genome wrapping length in kb', required=True)
args = parser.parse_args()

output_pdf = args.o
wrap_kb = int(args.l)

wrap_bp = wrap_kb * 1000

if args.tsv:
	sum_file = args.tsv
	df = pd.read_csv(
    	sum_file,
    	sep="\t"
    )
elif args.bed:
	sum_file = args.bed
	df = pd.read_csv(
    	sum_file,
    	sep="\t",
    	header=None,
    	names=["Query", "Start", "End", "Classification"]
	)
else:
	print('############################################################\n\
#ERROR: Annotation file in tsv or bed format is missing!\n\
Usage:\n\
python HGTScanner.py -m mt -q <query sequence> -o <output prefix> -f <family> [optional] -mt_add_seq <reference fasta> -e <e value> -b <bed file for masking>')
	sys.exit()

color_map = {
	"inconclusive": "#D3D3D3", # light gray
	"ancestral mt transfer (high confidence)": "#C9A400", # mustard / dark yellow
	"ancestral mt transfer (putative)": "#F2E085", # pale yellow
	"high confidence alien MTPT": "#C51B8A",  # magenta
	"native MTPT": "#228B22", # forest green
    "VGT": "#0B3C5D", #navy blue
    "High confidence HGT": "#D62728", #strong red
    "Putative HGT": "#FF7F0E", #orange
}

palette = plt.colormaps["tab10"].colors
chroms = df["Query"].unique()

with PdfPages(output_pdf) as pdf:
    for chrom in chroms:
        chrom_df = df[df["Query"] == chrom]
        max_pos = chrom_df["End"].max()
        num_rows = math.ceil(max_pos / wrap_bp)
        fig, ax = plt.subplots(figsize=(14, 1.5 * num_rows + 1))
        # Assign colors to unseen annotation types
        all_types = chrom_df["Classification"].unique()
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
                (chrom_df["Start"] < row_end) &
                (chrom_df["End"] > row_start)
            )
            for _, feat in chrom_df[mask].iterrows():
                f_start = max(feat["Start"], row_start)
                f_end = min(feat["End"], row_end)
                x_pos = f_start - row_start
                width = f_end - f_start
                ax.add_patch(
                    Rectangle(
                        (x_pos, y_pos),
                        width,
                        0.5,
                        color=color_map.get(feat["Classification"], "grey"),
                        edgecolor="black",
                        linewidth=0.4,
                        label=feat["Classification"]
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
