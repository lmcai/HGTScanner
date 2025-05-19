# HGTScanner		<img src="/images/logo.png" width="100" height="100">

**What is HGTScanner**: HGTScanner is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Cai, Liming, Cohen, Rachel. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner#v-general-guidelines-for-genome-skimming-data-collection)

[I. Prerequisites and installation](https://github.com/lmcai/HGTScanner#i-prerequisites-and-installation)
	
[II. Quick start](https://github.com/lmcai/HGTScanner#ii-quick-start)

[III. Complete tutorial for HGT detection in Aphyllon](https://github.com/lmcai/HGTScanner#iii-complete-tutorial-for-hgt-detection-in-aphyllon)

## I. General guidelines for experimental design

What HGTScanner **CAN** do:

- Identify the genomic location and donor lineage (to family and genus level) of HGT from an unknown green plant mitochondrion to a parasitic plant mitochondrion.
 
 	- Identify the genomic location and donor lineage of HGT from an unknown fungus mitochondrion to a mycoheterotrophic plant mitochondrion.
  
  	- Identify ancient, shared HGT among closely related HGT receiver species using phylogeny.
   
 	- Identify the genomic location, donor lineage (to genus level) of mitochondrial plastid DNAs in a group of plant species. No prior knowledge of the plastid DNA donor is known.

What is **very challenging** and discoraged use of HGTScanner:
	- Make claims regarding donor lineage at the species level
 	- 

## I. Prerequisites and installation

To process large datasets (>20 sp), high performance cluster is recommended. Mac and PC may suffer from insufficient memory during the assembly, alignment, or phylogenetic reconstruction. If you have difficulties installing PhyloHerb, please contact Liming Cai (lmcai@utexas.edu) or open an [Issue](https://github.com/lmcai/PhyloHerb/issues).

