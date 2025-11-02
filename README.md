# <img src="/images/logo.jpg" width="100" height="100"> 		HGTScanner		

**HGTScanner** is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Liming Cai and Rachel Cohen. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

**Contact**: Reach out to Liming Cai (cail@ufl.edu) for questions or open an [Issue](https://github.com/lmcai/HGTScanner/issues).

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner#v-general-guidelines-for-genome-skimming-data-collection)

[I. Prerequisites and installation](https://github.com/lmcai/HGTScanner#i-prerequisites-and-installation)
	
[II. Quick start](https://github.com/lmcai/HGTScanner#ii-quick-start)

[III. Complete tutorial for HGT detection in Aphyllon](https://github.com/lmcai/HGTScanner#iii-complete-tutorial-for-hgt-detection-in-aphyllon)

## I. General guidelines for experimental design

Examples of what HGTScanner can do :smiley: :

- Identify the genomic location and donor lineage (to family and genus level) of HGT from an **known or unknown** green plant mitochondrion to a parasitic plant mitochondrion, in both gene and **intergenic** regions.
 
- Identify the genomic location and donor lineage of HGT from an unknown fungus mitochondrion to a plant mitochondrion.

- Identify ancient, **shared** HGT among closely related HGT receiver species using phylogeny.
   
- Identify the genomic location, donor lineage (to genus level) of mitochondrial **plastid** DNAs in a group of plant species. No prior knowledge of the plastid DNA donor is known.

What is very challenging and **discoraged use** of HGTScanner :imp::

- Make claims regarding donor lineage at the species level

- Make claims of HGT between two closely-related sister lineages

- Apply HGTScanner to nuclear genomes (the primary concern is the all-by-all pairwise BLAST is not optimized for giga-sized genomes)

## I. Prerequisites and installation

**Installation**

To install HGTScanner, simply download it use `git`:
```
git clone https://github.com/lmcai/HGTScanner.git
```
Or from the source package:
```
#for v1.1
wget [html]
tar xzvf hgtscanner_v1.1.tar.gz
```
Then install dependencies via conda
```
cd HGTScanner
conda env create -f environment.yml
conda activate hgtscanner
```

To update your local version for any future releases, `cd` into the `PhyloHerb` directory then type
```
git fetch --prune origin
git reset --hard origin
git clean -f -d
```
