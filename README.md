# <img src="/images/logo.jpg" width="100" height="100"> 		HGTScanner		

**HGTScanner** is a python package to identify the genomic location and evolutionary history of horizontal gene transfers (HGT) in bacteria-like organellar genomes

**License**: MIT License

**Citation**: Liming Cai and Rachel Cohen. HGTScanner: Detecting horizontal gene transfer at fine scale using whole organelle genomes.

**Contact**: Reach out to Liming Cai (cail@ufl.edu) for questions or open an [Issue](https://github.com/lmcai/HGTScanner/issues).

## Quick link
[I. General guidelines for experimental design](https://github.com/lmcai/HGTScanner/blob/main/README.md#i-general-guidelines-for-experimental-design)

[II. Installation](https://github.com/lmcai/HGTScanner#ii-installation)
	
[III. Quick start](https://github.com/lmcai/HGTScanner#iii-quick-start)

[IV. Complete tutorial for HGT detection in Aphyllon](https://github.com/lmcai/HGTScanner#iv-complete-tutorial-for-hgt-detection-in-aeginetia)

## I. General guidelines for experimental design

Examples of what HGTScanner can do :smiley: :

- Identify the genomic location and donor lineage (to family and genus level) of HGT from an **known or unknown** green plant mitochondrion to a parasitic plant mitochondrion, in both gene and **intergenic** regions.
 
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

## IV. Complete tutorial for MTPT and HGT detection in Aeginetia

Here, we demonstrate how to identify MTPT and HGT in Aeginetia indica, which is an Orobanchaceae parasitic plant specialized in grasses. All input and output files can be found in the `example`.

### 1. Input preparation for MTPT

a. The query mitochondrial genome assembly `Aeginetia_indica.fas` (GenBank ID: NC_069194) in FASTA format.

b. List of ingroup families `taxon.txt`, which contained all families of Lamiales. The query family is labeled as `query` and others as `ingroup`. The label and the family names are seperated by a tab.

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
>Orobanchaceae|Conopholis_alpina|LM013
ATATCATTATGATAAAATTGGTAAATTAATGCTGTTATGATGAAATTGGTAGACATGTTGCTTTTAGACAGCAATATTAA
```

### 2. Run HGTScanner for MTPT
Run the following commands:
```
python HGTscanner.py -m mtpt -q Aeginetia_indica.fas  -o Ain -taxon taxon.txt -pt_add_seq pt.fasta
```

### 3. Output interpretation

a. Summary spreadsheet: Ain.mtpt.sum.tsv

| Locus_ID	| Target_scaffold	| Start	| End	| Phylo_classification	| Support	| Sister_family	| Sister_genus	| Sister_species |
| 1	| Ain_1	| 4409	| 5021	| native MTPT	| 0.862	| Orobanchaceae	| Xylanche,Christisonia,Gleadovia,Schwalbea,Pedicularis,Centranthera,Aphyllon,Siphonostegia,Phtheirospermum,Harveya,Phacellanthus,Lathraea,Aeginetia,Pterygiella	| Orobanchaceae|LM037|LM037_10669_11244,Orobanchaceae|Castilleja_paramensis|NC_031805.1_141671_142246,Orobanchaceae|Castilleja_paramensis|Castilleja_paramensis_141671_142246,Orobanchaceae|LM001|LM001_138201_138776,Orobanchaceae|LM105|LM105_141668_142243,Orobanchaceae|Triphysaria_versicolor|NC_053793.1_141314_141889,Orobanchaceae|Pedicularis_chinensis|Pedicularis_chinensis_136262_136837,Orobanchaceae|LM040|LM040_140977_141552,Orobanchaceae|Phtheirospermum_japonicum|NC_053792.1_142148_142723,Orobanchaceae|Xylanche_himalaica|NC_068836.1_60244_60839,Orobanchaceae|Aeginetia_indica|Aeginetia_indica_17251_17792,Orobanchaceae|Christisonia_kwangtungensis|Christisonia_kwangtungensis_44660_44975,Orobanchaceae|Pedicularis_alaschanica|NC_080929.1_136003_136578,Orobanchaceae|LM162|LM162_128435_129010,Orobanchaceae|Melampyrum_koreanum|NC_057523.1_132698_133279,Orobanchaceae|Gleadovia_mupinensis|NC_086656.1_104068_104632,Orobanchaceae|Phacellanthus_tubiflorus|NC_068243.1_53992_54566,Orobanchaceae|Cistanche_sinensis|Cistanche_sinensis_20869_21397,Orobanchaceae|LM081|LM081_22108_22680,Orobanchaceae|Cistanche_tubulosa|Cistanche_tubulosa_20016_20577,Orobanchaceae|LM089|LM089_12913_13467,Orobanchaceae|LM099|LM099_66850_67425,Orobanchaceae|Aphyllon_californicum|NC_025651.1_109761_110336,Orobanchaceae|LM016|LM016_21663_22238,Orobanchaceae|LM098|LM098_10502_11077,Orobanchaceae|LM014|LM014_100642_101218,Orobanchaceae|LM020|LM020_22342_22922,Orobanchaceae|Aphyllon_purpureum|Aphyllon_purpureum_13430_14005,Orobanchaceae|Aphyllon_fasciculatum|NC039679_54369_54935,Orobanchaceae|LM033|LM033_141452_142027,Orobanchaceae|LM121|LM121_141964_142539,Orobanchaceae|LM136|LM136_141643_142218,Orobanchaceae|Lindenbergia_indica|NC_080204.1_121609_122185,Orobanchaceae|LM085|LM085_139232_139808,Orobanchaceae|Lathraea_squamaria|Lathraea_squamaria_139299_139874,Orobanchaceae|Lathraea_japonica|NC_086659.1_10667_11232,Orobanchaceae|Pterygiella_cylindrica|NC_063636.1_143550_144125,Orobanchaceae|Brandisia_cauliflora|NC_073570.1_142850_143425,Orobanchaceae|LM122|LM122_143312_143887,Orobanchaceae|Cymbaria_daurica|NC_064104.1_140623_141198,Orobanchaceae|LM161|LM161_25738_26313,Orobanchaceae|Siphonostegia_chinensis|NC_046038.1_138409_138984,Orobanchaceae|Schwalbea_americana|NC_023115.1_150473_151048,Orobanchaceae|LM051|LM051_119310_119885,Orobanchaceae|LM109|LM109_120167_120742,Orobanchaceae|LM143|LM143_137832_138407,Orobanchaceae|Striga_asiatica|NC_067570.1_148229_148528,Orobanchaceae|LM093|LM093_20680_21026,Orobanchaceae|Centranthera_grandiflora|NC_059747.1_136935_137510,Orobanchaceae|LM117|LM117_137010_137585,Orobanchaceae|Harveya_capensis|Harveya_capensis_10610_11179,Orobanchaceae|LM030|LM030_130664_131239,Orobanchaceae|Euphrasia_regelii|NC_045041.1_141798_142373
| 8	| Ain_1	| 36973	| 37263	| ancestral mt transfer	| NA	| NA	| NA	| NA |
