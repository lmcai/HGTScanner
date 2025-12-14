# <img src="/images/logo.jpg" width="100" height="100"> 		HGTScanner		

**HGTScanner** is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Liming Cai and Rachel Cohen. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

**Contact**: Reach out to Liming Cai (cail@ufl.edu) for questions or open an [Issue](https://github.com/lmcai/HGTScanner/issues).

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner/blob/main/README.md#i-general-guidelines-for-experimental-design)

[II. Installation](https://github.com/lmcai/HGTScanner#ii-installation)
	
[III. Quick start](https://github.com/lmcai/HGTScanner#iii-quick-start)

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

## II. Installation

HGTScanner has been test on macOS and Linux.

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

3. Download Viridiplantae database

4. To update your local version for any future releases, `cd` into the `PhyloHerb` directory then type
```
git fetch --prune origin
git reset --hard origin
git clean -f -d
```

## III. Quick start

### Usage:

```
HGTscanner.py [-h] -m mode -o output_prefix -taxon taxon_file [-q query] [-pt_fix_id id_file]
              [-pt_add_seq fatsa] [-pt_add_id id_file] [-hit integer] [-mt_add_seq fasta] [-wd dir]
              [-b bed_file] [-e evalue] [-nofasttree]
```

### Options:
```
options:
options:
  -h, --help          show this help message and exit
  -m mode             Choose from: mtpt, mtpt_eval, mt, mt_eval
  -o output_prefix    Output prefix
  -taxon taxon_file   A file containing the family of the query and a list of its close relatives for HGT classification
  -q query            Fasta file of the query genome
  -pt_fix_id id_file  A file of user-selected GenBank accession numbers for MTPT detection.
  -pt_add_seq fatsa   A fasta file containing plastid references for MTPT detection.
  -pt_add_id id_file  A file user-selected GenBank accession numbers for MTPT detection.
  -hit integer        Number of best blast hits to be included.
  -mt_add_seq fasta   A fasta file containing mitochondrial references for mt HGT detection.
  -wd dir             Path to working dir where *.sum.tsv and *_HGTscanner_supporting_files are located.
  -b bed_file         A bed file for regions to be masked
  -e evalue           BLAST evalue threshold
  -nofasttree         No FastTree phylogeny inference
```

### 1. MTPT detection

*Input:* A fasta-formatted assembly of the query organelle genome. To identify MTPT, use the following command:

```
#Use the default 3722-genera plastid Viridiplantae database
python HGTscanner.py -m mtpt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file]
```
The list of species included in our built-in plastid Viridiplantae database can be found [here](/database/pt_Viridiplantae_taxonomy.tsv). One representative species per genus (totalling 3722) has been selected from the entire NCBI plastid reference genome database. Users can add more species from the **built-in** database by:
```
#Add species to the default database
python HGTscanner.py -m mtpt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file] -pt_add_id [more_genbank_id.txt]
#Use a user-selected Genbank accessions as BLAST database (instead of adding to the default db)
python HGTscanner.py -m mtpt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file] -pt_fix_id [genbank_id.txt]
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
python HGTscanner.py -m mt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file] [optional] -b [bed_file_for_masking]
```

*Output:* 

The following files will be generated:

### 3. Adding custom sequences

Users can also add their own plastid (via `-pt_add_seq`) or mitochondrial (via `-mt_add_seq`) sequences to provide more comprehensive taxon sampling. The headers of these sequences need to be formatted as "family|binomial_species|[additional_header_string]". **No space (' ') should be included in the header.** This is necessary to use the taxonomy information for HGT detection. Unformatted sequences will not be used in the downstream analysis.

Example of user-added fasta
```
>Fabaceae|Indigofera_tinctoria|my_seq_ID1
ATCGATCGATCG
>Geraniaceae|Geranium_maculatum|my_seq_ID2
ATCGATCGATCG
...
```
Once the sequences are correctly prepared, the users can add them to the database to identify MTPT or mt HGT.

For MTPT:
```
#Add custom sequences
python HGTscanner.py -m mtpt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file]] -pt_add_seq [fasta_file] 
#Add custom sequences and add more sequences from the built-in database
python HGTscanner.py -m mtpt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file] -pt_add_id [more_genbank_id.txt] -pt_add_seq [fasta_file]
```
For mt HGT:
```
python HGTscanner.py -m mt -q [query.fas]  -o [output_prefix] -taxon [taxonomy_file] -mt_add_seq [mt.fas] [optional] -b [bed_file_for_masking] 
```
