# <img src="/images/logo.jpg" width="100" height="100"> 		HGTScanner		

**HGTScanner** is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Liming Cai and Rachel Cohen. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

**Contact**: Reach out to Liming Cai (cail@ufl.edu) for questions or open an [Issue](https://github.com/lmcai/HGTScanner/issues).

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner#v-general-guidelines-for-genome-skimming-data-collection)

[II. Prerequisites and installation](https://github.com/lmcai/HGTScanner#i-prerequisites-and-installation)
	
[III. Quick start](https://github.com/lmcai/HGTScanner#ii-quick-start)

[IV. Complete tutorial for HGT detection in Aphyllon](https://github.com/lmcai/HGTScanner#iii-complete-tutorial-for-hgt-detection-in-aphyllon)

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

## II. Prerequisites and installation

**Installation**

1. To install HGTScanner, simply download it use `git`:
```
git clone https://github.com/lmcai/HGTScanner.git
```
Or from the source package:
```
#for v1.1
wget [html]
tar xzvf hgtscanner_v1.1.tar.gz
```
2. Then install dependencies via conda, which will create a conda environment called `hgtscanner`
```
cd HGTScanner
conda env create -f environment.yml
```
The environment can be activated by
```
conda activate hgtscanner
```

3. To update your local version for any future releases, `cd` into the `PhyloHerb` directory then type
```
git fetch --prune origin
git reset --hard origin
git clean -f -d
```

## II. Quick start

### Usage:

```
HGTscanner.py		-q query -o output_prefix -f family [-mtpt] 
					[-pt_fix_id id_file] [-pt_add_seq fatsa] [-pt_add_id id_file] 
					[-hit number of hits] [-mt_add_seq fasta] [-b bed_file] [-e evalue]
```

### Options:
```
options:
  -h, --help           show this help message and exit
  -q query             Fasta file of the query genome
  -o output_prefix     Output prefix
  -f family            Family of the query for HGT classification
  -mtpt                Invoking the MTPT mode
  -pt_fix_id id_file   A file of user-selected GenBank accession numbers for MTPT detection.
  -pt_add_seq fatsa    A fasta file containing plastid references for MTPT detection.
  -pt_add_id id_file   A file user-selected GenBank accession numbers for MTPT detection.
  -hit number		   Number of best blast hits to be included.
  -mt_add_seq fasta    A fasta file containing mitochondrial references for mt HGT detection.
  -b bed_file          A bed file for regions to be masked
  -e evalue            BLAST evalue threshold
```

### 1. MTPT detection

*Input:* A fasta-formatted assembly of the query organelle genome. To identify MTPT, use the following command:

```
#Use the default 3722-genera plastid Viridiplantae database
python HGTscanner.py -mtpt -q [query.fas]  -o [output_prefix] -f [query_species_family]
```
The list of species included in our built-in plastid Viridiplantae database can be found [here](/database/pt_Viridiplantae_taxonomy.tsv). One representative species per genus (totalling 3722) has been selected from the entire NCBI plastid reference genome database. Users can add more species from the **built in** database by:
```
#Add species to the default database
python HGTscanner.py -mtpt -q [query.fas]  -o [output_prefix] -f [query_species_family] -pt_add_id [more_genbank_id.txt]
#Use a user-selected Genbank accessions as BLAST database (instead of adding to the default db)
python HGTscanner.py -mtpt -q [query.fas]  -o [output_prefix] -f [query_species_family] -pt_fix_id [genbank_id.txt]
```

*Output:* 

The following files will be generated:
```
#Most important
*.mtpt.sum.tsv: A summary spreadsheet of the location and classification of MTPT regions.
#Supporting evidence
*_HGTscanner_supporting_files: A folder containing the alignments and raw sequences of candidate MTPT regions.
*.mtpt.blast: BLAST result between the query and the plastid database.
*.mtpt.bed: A bed file of BLAST hit between the query and the database, with taxon information added. The last column is a unique ID for this BLAST hit.
*.mtpt.merged.bed: A bed file generated by merging BLAST hits from `*.mtpt.bed`.
```
### 2. HGT detection in mitochondrial genomes

*Input:* A fasta-formatted assembly of the query organelle genome. To identify mito HGT, use the following command. Note that the optional bed file for masking gene coding and MTPT regions is optional, but highly recommended to avoid excessive BLAST hits in these loci:
```
python HGTscanner.py -q [query.fas]  -o [output_prefix] -f [query_species_family] -b [bed_file_for_masking]
```

*Output:* 

The following files will be generated:

### 3. Adding custom sequences

Users can also add their own plastid (via `-pt_add_seq`) or mitochondrial (via `-mt_add_seq`) sequences to provide more comprehensive taxon sampling. The headers of these sequences need to be formatted as "family|species|[additional_header_string]". **No space (' ') should be included in the header.** This is necessary to use the taxonomy information for HGT detection. Unformatted sequences will not be used in the downstream analysis.
```
#Example of user-added fasta
>Fabaceae|Indigofera_tinctoria|my_seq_ID1
ATCGATCGATCG
>Geraniaceae|Geranium maculatum|my_seq_ID2
ATCGATCGATCG
...
```

For MTPT:
```
#Add custom sequences
python HGTscanner.py -mtpt -q [query.fas]  -o [output_prefix] -f [query_species_family] -pt_add_seq [fasta_file] 
#Or combining many sources
python HGTscanner.py -mtpt -q [query.fas]  -o [output_prefix] -f [query_species_family] -pt_add_id [more_genbank_id.txt] -pt_add_seq [fasta_file]
```
For mt HGT:
```
python HGTscanner.py -q [query.fas]  -o [output_prefix] -f [query_species_family] -b [bed_file_for_masking] -mt_add_seq [mt.fas]
```
