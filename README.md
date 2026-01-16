# <img src="/images/logo.jpg" width="100" height="100"> 		HGTScanner		

**HGTScanner** is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Liming Cai and Rachel Cohen. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

**Contact**: Reach out to Liming Cai (cail@ufl.edu) for questions or open an [Issue](https://github.com/lmcai/HGTScanner/issues).

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner/blob/main/README.md#i-general-guidelines-for-experimental-design)

[II. Installation](https://github.com/lmcai/HGTScanner#ii-installation)
	
[III. Quick start](https://github.com/lmcai/HGTScanner#iii-quick-start)

[IV. Complete tutorial for HGT detection in Aeginetia](https://github.com/lmcai/HGTScanner#iv-complete-tutorial-for-mtpt-and-hgt-detection-in-aeginetia)

## I. General guidelines for experimental design

Examples of what HGTScanner can do :smiley: :

- Identify the genomic location and donor lineage (to family and, less confidently, genus level) of HGT from an **known or unknown** green plant mitochondrion to a parasitic plant mitochondrion, in both gene and **intergenic** regions.
 
- Identify the genomic location and donor lineage of HGT from an unknown fungus mitochondrion to a plant mitochondrion. With user-supplied fungus sequence database.

- Identify ancient, **shared** HGT among closely related species using phylogeny.
   
- Identify the genomic location, donor lineage (to family/genus level) of mitochondrial **plastid** DNAs in a group of plant species. No prior knowledge of the plastid DNA donor is known.

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
              [-b bed_file] [-e evalue] [-notree]
```

### Options:
```
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
  -notree         No FastTree phylogeny inference
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
*_HGTscanner_supporting_files: A folder containing the alignments and raw sequences of candidate MTPT regions. Sequences are named as 'Family|Genus_species|ID_start1_end1+start2_end2+...'
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
```
#Most important
*.hgt.sum.tsv: A summary spreadsheet of the location and classification of HGT regions.
#Supporting evidence
*_HGTscanner_supporting_files: A folder containing the alignments and raw sequences of candidate HGT regions. Sequences are named as 'Family|Genus_species|ID_start1_end1+start2_end2+...'
*.mt.blast: BLAST result between the query and the mitochondrial database.
*.mt.bed: A bed file of BLAST hit between the query and the database, with taxon information added. The last column is a unique ID for this BLAST hit.
*.merged.bed: A bed file generated by merging BLAST hits from `*.mtpt.bed`.
*.alnmap.bed: A bed file for the start and end location of the query sequence for each locus.
```

### 3. Adding custom sequences

Users can also add their own plastid (via `-pt_add_seq`) or mitochondrial (via `-mt_add_seq`) sequences to provide more comprehensive taxon sampling. The headers of these sequences need to be formatted as "Family|Genus_species|additional_seq_ID". **No space (' ') should be included in the header.** This is necessary to use the taxonomy information for HGT detection. Unformatted sequences will not be used in the downstream analysis.

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

### 4. Use IQ-TREE/RAxML instead of FastTree
Users can choose IQ-TREE/RAxML instead of FastTree for phylogenetic inference to improve accuracy. This also allows for the distribution of phylogenetic inference to high performance computer clusters. 

**Important**: The resulting phylogeny must be named in the format of `[prefix].mtpt.[id].aln.fas.treefile` for MTPT and `[prefix].hgt.[id].aln.fas.treefile` for HGT (e.g., Atr.hgt.1.aln.fas.treefile). 

These resulting trees should be placed in the same folder as the alignments in `[output_prefix]_HGTscanner_supporting_files`. Then the following command can be used to assess HGT:
```
#for MTPT
python HGTscanner.py -m mtpt_eval -o [output_prefix] -taxon [taxonomy_file] -wd [working_dir]
#for HGT
python HGTscanner.py -m mt_eval -o [output_prefix] -taxon [taxonomy_file] -wd [working_dir]
```
Note that the `-wd` argument should be the name of the directory containing both *.mtpt.sum.tsv and *_HGTscanner_supporting_files:
```
working_dir/
├── sp_HGTscanner_supporting_files/
│   ├── sp.mtpt.1.fas
│   └── sp.mtpt.1.aln.fas
│   └── sp.mtpt.1.aln.fas.treefile
│   ...
└── sp.mtpt.sum.tsv
```

### 5. Visualization
The built-in script `plot_annotation.py` allows users to visualize annotation results generated by HGTScanner. The input annotation file can be either the `*.sum.tsv` output from HGTScanner or a user-curated `*.bed` file.

To improve visualization of long genomes, users must specify a wrapping length parameter in kb (e.g., `-l 100`). The genome will be soft-wrapped at this length and displayed across multiple horizontal rows, while preserving the original genomic coordinates. This approach increases plot resolution and readability without fragmenting annotations.

a) For direct output from HGTScanner:
```
python plot_annotation.py -tsv ./example/Ain.mtpt.sum.tsv -l 100 -o Ain.mtpt.pdf
```
This result in plot like this:


b) For user-curated bed file `Ain.ann.bed`, which has the following format:
```
Ain_1	4409	5021	native MTPT
Ain_1	36973	37263	HGT
Ain_1	37676	38793	ancestral mt transfer (high confidence)
Ain_1	84935	85061	VGT
Ain_1	97432	97980	gene
...
```
The following command can be used:
```
python plot_annotation.py -tsv Ain.ann.bed -l 100 -o Ain.mtpt.pdf
```


## IV. Complete tutorial for MTPT and HGT detection in Aeginetia

Here, we demonstrate how to identify MTPT and HGT in Aeginetia indica, which is an Orobanchaceae parasitic plant specialized in grasses. All input and output files can be found in the `/example`.

### 1. Input preparation for MTPT

a. The query mitochondrial genome assembly `Aeginetia_indica.fas` (GenBank ID: NC_069194) in FASTA format.

b. List of ingroup families `taxon.txt`, which contained all families of Lamiales. The query family is labeled as `query` and others as `ingroup`. The the family names and labels are seperated by a tab.

```
Orobanchaceae	query
Lamiaceae	ingroup
Scrophulariaceae	ingroup
Oleaceae	ingroup
Phrymaceae	ingroup
Acanthaceae	ingroup
Verbenaceae	ingroup
Plantaginaceae	ingroup
Gesneriaceae	ingroup
Bignoniaceae	ingroup
```

c. Custom plastid sequence `pt.fasta`. The header should follow the format `Family|Genus_species|ID`:

```
>Orobanchaceae|Aeginetia_indica|Aeginetia_indica
ATATCATTATGATAAAATTGGTAAATTAATGCTGTTATGATGAAATTGGTAGACATGTTGCTTTTAGACAGCAATATTAA
...
>Orobanchaceae|Conopholis_alpina|LM013
ATATCATTATGATAAAATTGGTAAATTAATGCTGTTATGATGAAATTGGTAGACATGTTGCTTTTAGACAGCAATATTAA
...
```

### 2. Run HGTScanner for MTPT
Run the following commands:
```
python HGTscanner.py -m mtpt -q Aeginetia_indica.fas  -o Ain -taxon taxon.txt -pt_add_seq pt.fasta
```

### 3. Output interpretation

The following output files were generated

a. Summary spreadsheet `Ain.mtpt.sum.tsv` with the following columns:
- Locus_ID: numbered ID of the locus
- Target_scaffold: sequence header of the query FASTA
- Start: Start position of the loci 
- End: End position of the loci
- Classification: classification of MTPT based on the criteria presented in Cai and Cohen (2026).
- Support: Branch support for the placement of the query sequence. May range 0-1 for FastTree and 0-100 for IQ-TREE and RAxML.
- Sister_family: List of families in the immediate sister of the query, seperated by ',' 
- Sister_genus: List of genera in the immediate sister of the query, seperated by ',' 
- Sister_species: List of species in the immediate sister of the query, seperated by ',' 

b. `Ain.mtpt.blast`: the original BLASTN result

c. `Ain.mtpt.bed`: a bed file of BLASTN hits ordered by position and with added taxon information. The 11 columns contain the following information: 
- (1) query sequence ID
- (2) query start
- (3) query end
- (4) hit sequence ID
- (5) hit start
- (6) hit end
- (7) bit score
- (8) e-value
- (9) hit species
- (10) hit family
- (11) unique ID for the BLAST record

d. `Ain.mtpt.merged.bed`: a bed file for BLAST record IDs contained in each locus. The {i}th row represents the homologous sequence locations for {i}th locus in the summary spreadsheet. The list of numbers in the fourth column represents the unique BLAST record ID in `Ain.mtpt.bed`.

e. `Ain_HGTscanner_supporting_files`: a folder containing the raw sequence `.fas`, alignment `.aln.fas`, and phylogeny `.treefile` for each locus.

After examining the `Ain.mtpt.sum.tsv`, locus #1, 15, 26, 28 were classified as **native MTPT**, locus #14 and 18 are **high confidence alien MTPTs** from Poaceae. These loci should be masked in the next step. Others are either mitochondrial transfers or inconclusive.

### 4. Input preparation for mt HGT

a. The query mitochondrial genome assembly `Aeginetia_indica.fas` 

b. List of ingroup families `taxon.txt`

c. A bed file to mask gene regions and MTPTs `Ain.mask.bed`
```
Ain_1	4409	5021	mtpt
Ain_1	84935	85061	mtpt
Ain_1	86090	86359	mtpt
Ain_1	67	525 rpl16
Ain_1	39819	40223	nad3
...
```

### 5. Running mt HGT identification

Use the following command to identify homology, scaffold alignment, build phylogeny, and evaluate HGT:
```
python HGTscanner.py -m mt -q Aeginetia_indica.fas  -o Ain -taxon taxon.txt -b Ain.mask.bed
```
