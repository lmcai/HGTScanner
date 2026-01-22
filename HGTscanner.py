import sys
ascii_art = r"""
  _  _  ___ _____ ___                           
 | || |/ __|_   _/ __|__ __ _ _ _  _ _  ___ _ _ 
 | __ | (_ | | | \__ / _/ _` | ' \| ' \/ -_| '_|
 |_||_|\___| |_| |___\__\__,_|_||_|_||_\___|_|  
                                                
"""
print(ascii_art)
print('############################################################\n\
HGTScanner v1.1\n\
A Python tool for genome-wise detection of horizontal gene transfers\n')

try:
    # Attempt to import all required (non-standard) modules
    from Bio import SeqIO
    from ete3 import Tree
    import ete3
    import pybedtools
    from numpy import median
    import numpy as np
    from statistics import mode
    import statistics as stats
    import pandas as pd
    #from scipy.spatial import cKDTree
    from hmmlearn.hmm import PoissonHMM
    #set np random seed to get a consistent outcome
    np.random.seed(27)
    
    # Import standard modules after the checks are passed
    import subprocess
    import os, argparse
    import datetime
    import warnings
    warnings.filterwarnings("ignore")
except ImportError as e:
    # If any module is missing, print a warning to the screen and exit
    print("------------------------------------------------------------")
    print("ERROR: Missing Required Python Module")
    print(f"The module '{e.name}' could not be imported.")
    print("\nPlease ensure all required dependencies are installed.")
    print("\nRequired modules include: biopython (Bio), ete3, pybedtools, numpy, pandas.")
    print("------------------------------------------------------------")
    sys.exit(1) # Exit the script with an error code


parser = argparse.ArgumentParser(description='HGTscanner identifies HGT blocks in organellar genomes.')
parser.add_argument('-m', metavar='mode', help='Choose from: mtpt, mtpt_eval, mt, mt_eval', required=True)
parser.add_argument('-o', metavar='output_prefix', help='Output prefix', required=True)
parser.add_argument('-taxon', metavar='taxon_file', help='A file containing the family of the query and a list of its close relatives for HGT classification', required=True)
parser.add_argument('-q', metavar='query', help='Fasta file of the query genome')
parser.add_argument('-pt_fix_id', metavar='id_file', help='A file of user-selected GenBank accession numbers for MTPT detection.')
parser.add_argument('-pt_add_seq', metavar='fatsa', help='A fasta file containing plastid references for MTPT detection.')
parser.add_argument('-pt_add_id', metavar='id_file', help='A file user-selected GenBank accession numbers for MTPT detection.')
parser.add_argument('-hit', metavar='integer', help='Number of best blast hits to be included.')
parser.add_argument('-mt_add_seq', metavar='fasta', help='A fasta file containing mitochondrial references for mt HGT detection.')
parser.add_argument('-wd', metavar='dir', help='Path to working dir where *.sum.tsv and *_HGTscanner_supporting_files are located.')
parser.add_argument('-b', metavar='bed_file', help='A bed file for regions to be masked')
parser.add_argument('-e', metavar='evalue', help='BLAST evalue threshold')
parser.add_argument('-notree', action='store_true', help='No FastTree phylogeny inference')

####################################
#I. pass argument values, check required arguments and input file format
################################

args = parser.parse_args()
script_path = os.path.abspath(sys.argv[0])
script_path = os.path.dirname(script_path)

def id2seq(ids,output_file):
	recs=SeqIO.parse(script_path+'/database/Viridiplantae_pt_aug2025.genome.fas','fasta')
	out=open(output_file,'a')
	for rec in recs:
		id=rec.id
		if id.split('.')[0] in ids:d=SeqIO.write(rec,out,'fasta')
		out.close()

if args.m =='mtpt':
	#mtpt mode
	#assemble the plastid dataset for MTPT
	try:
		sp=args.o
		query=args.q
		print(str(datetime.datetime.now())+'\tDetecting MTPT in '+sp+' using the query sequence '+query)
		print(str(datetime.datetime.now())+'\tChecking taxonomy file...')
		taxon = open(args.taxon).readlines()
		ingroup=[]
		fam=''
		for l in taxon:
			try:
				ingroup.append(l.split()[0])
				if l.split()[1].lower() == 'query':fam=l.split()[0]
			except IndexError:
				print(str(datetime.datetime.now())+'\tMalformatted taxonomy file. Please double check. Exit...')
				sys.exit()
		if fam=='':
			sys.exit(str(datetime.datetime.now())+'\tQuery family not set. Exit...')
		print(f"{datetime.datetime.now()}\tTaxonomy file looks OK")
		print(f"{datetime.datetime.now()}\tThe query belongs to family: {fam}; The following are close relatives: {', '.join(ingroup)}")
		if args.pt_fix_id:
			#Use custom list of pt db
			id_list=open(args.pt_fix_id)
			id2seq([i.strip() for i in id_list],sp+'.pt_db.fas')
			print(str(datetime.datetime.now())+'\tUsing custom list of plastid reference: '+args.pt_fix_id)
		else:
			#build db on top of default  3722-sp pt db
			add_seq=0
			if args.pt_add_seq:
				#add custom plastid sequences to database
				pt_reference = args.pt_add_seq
				print(str(datetime.datetime.now())+'\tAdd custom plastid fasta '+pt_reference+' in addition to the NCBI Viridiplantae plastid database')
				print(str(datetime.datetime.now())+'\tChecking file formats for '+pt_reference)
				#check header of the input seq
				headers=open(pt_reference).readlines()
				headers=[i[1:] for i in headers if i.startswith('>')]
				if all(i.count("|") == 2 for i in headers):
					S='cat '+pt_reference+' '+script_path+'/database/Viridiplantae_pt_aug2025.representative.fas >'+sp+'.pt_db.fas'
					os.system(S)
					add_seq=1
				else:
					sys.exit(str(datetime.datetime.now())+'\tMalformatted custom fasta file: '+pt_reference+'. All headers should be >FAMILY|SPECIES|ID. Exit...')
			if args.pt_add_id:
				add_seq=1
				print(str(datetime.datetime.now())+'\tAdd custom list of plastid based on '+args.pt_add_id+' in addition to the NCBI Viridiplantae plastid database')
				id2seq([i.strip() for i in id_list],sp+'.pt_db.fas')
				S='cat '+script_path+'/database/Viridiplantae_pt_aug2025.representative.fas >>'+sp+'.pt_db.fas'
			if add_seq==0:
				print(str(datetime.datetime.now())+'\tNo custom plastid reference provided. Using default...')
				S='cp '+script_path+'/database/Viridiplantae_pt_aug2025.representative.fas '+sp+'.pt_db.fas'
				os.system(S)				
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -m mtpt -q <query sequence> -o <output prefix> -taxon <file> [optional] -e <e value> -pt_add_seq <fasta> -pt_add_id <list of id file> -pt_fix_id <list of id file>')
	except Exception as e:
		print(f"Error: {type(e).__name__} - {e}")
		sys.exit()
elif args.m == 'mtpt_eval':
	if args.wd is None:parser.error("the following argument is required for 'mtpt_eval': -wd")
	try:
		sp=args.o
		print(str(datetime.datetime.now())+'\tClassify MTPT for '+sp)
		print(str(datetime.datetime.now())+'\tChecking taxonomy file...')
		taxon = open(args.taxon).readlines()
		ingroup=[]
		fam=''
		for l in taxon:
			try:
				ingroup.append(l.split()[0])
				if l.split()[1].lower() == 'query':fam=l.split()[0]
			except IndexError:
				print(str(datetime.datetime.now())+'\tMalformatted taxonomy file. Please double check. Exit...')
				sys.exit()
		if fam=='':
			sys.exit(str(datetime.datetime.now())+'\tQuery family not set. Exit...')
		print(f"{datetime.datetime.now()}\tTaxonomy file looks OK")
		print(f"{datetime.datetime.now()}\tThe query belongs to family: {fam}; The following are close relatives: {', '.join(ingroup)}")
		sum_path=args.wd
		if not (os.path.isfile(f"{sum_path}/{sp}.mtpt.sum.tsv") and os.path.isdir(f"{sum_path}/{sp}_HGTscanner_supporting_files")):
			sys.exit(str(datetime.datetime.now())+f"\tNo {sp}.mtpt.sum.tsv or {sp}_HGTscanner_supporting_files found in {sum_path}. Exit...")
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -m mtpt_eval -o <output prefix> -wd <dir> -taxon <file>')
	except Exception as e:
		print(f"Error: {type(e).__name__} - {e}")
		sys.exit()
elif args.m == 'mt_eval':
	if args.wd is None:parser.error("the following argument is required for 'mt_eval': -wd")
	try:
		sp=args.o
		print(str(datetime.datetime.now())+'\tClassify HGT for '+sp)
		print(str(datetime.datetime.now())+'\tChecking taxonomy file...')
		taxon = open(args.taxon).readlines()
		ingroup=[]
		fam=''
		for l in taxon:
			try:
				ingroup.append(l.split()[0])
				if l.split()[1].lower() == 'query':fam=l.split()[0]
			except IndexError:
				print(str(datetime.datetime.now())+'\tMalformatted taxonomy file. Please double check. Exit...')
				sys.exit()
		if fam=='':
			sys.exit(str(datetime.datetime.now())+'\tQuery family not set. Exit...')
		print(f"{datetime.datetime.now()}\tTaxonomy file looks OK")
		print(f"{datetime.datetime.now()}\tThe query belongs to family: {fam}; The following are close relatives: {', '.join(ingroup)}")
		sum_path=args.wd
		if not (os.path.isfile(f"{sum_path}/{sp}.alnmap.bed") and os.path.isdir(f"{sum_path}/{sp}_HGTscanner_supporting_files")):
			sys.exit(str(datetime.datetime.now())+f"\tNo {sp}.alnmap.bed or {sp}_HGTscanner_supporting_files found in {sum_path}. Exit...")
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -m mt_eval -o <output prefix> -wd <dir> -taxon <file>')
	except Exception as e:
		print(f"Error: {type(e).__name__} - {e}")
		sys.exit()
elif args.m == 'mt':
	#default mode for mtHGT
	try:
		query=args.q
		sp=args.o
		print(str(datetime.datetime.now())+'\tDetecting HGT in '+sp+' using the query sequence '+query)
		print(str(datetime.datetime.now())+'\tChecking taxonomy file...')
		taxon = open(args.taxon).readlines()
		ingroup=[]
		fam=''
		for l in taxon:
			try:
				ingroup.append(l.split()[0])
				if l.split()[1].lower() == 'query':fam=l.split()[0]
			except IndexError:
				print(str(datetime.datetime.now())+'\tMalformatted taxonomy file. Please double check. Exit...')
				sys.exit()
		if fam=='':
			sys.exit(str(datetime.datetime.now())+'\tQuery family not set. Exit...')
		print(f"{datetime.datetime.now()}\tTaxonomy file looks OK")
		print(f"{datetime.datetime.now()}\tThe query belongs to family: {fam}; The following are close relatives: {', '.join(ingroup)}")
		if args.mt_add_seq:
			#add custom references to database
			mt_reference = args.mt_add_seq
			print(str(datetime.datetime.now())+'\tAdded custom mitochondrial reference '+mt_reference+' to the NCBI Viridiplantae mitochondrial database')
			print(str(datetime.datetime.now())+'\tChecking file formats for '+mt_reference)
			#check header of the input seq
			headers=open(mt_reference).readlines()
			headers=[i[1:] for i in headers if i.startswith('>')]
			if all(i.count("|") == 2 for i in headers):
				S='cat '+mt_reference+' '+script_path+'/database/Viridiplantae_mt_aug2025.genome.fas >'+sp+'.mt_db.fas'
				os.system(S)
			else:
				sys.exit(str(datetime.datetime.now())+'\tMalformatted custom fasta file: '+pt_reference+'. All headers should be >FAMILY|SPECIES|ID. Exit...')
		else:
			S='cp '+script_path+'/database/Viridiplantae_mt_aug2025.genome.fas '+sp+'.mt_db.fas'
			os.system(S)
			print(str(datetime.datetime.now())+'\tNo custom mitochondrial reference provided. Will use the built-in Viridiplantae mitochondrial database')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -m mt -q <query sequence> -o <output prefix> -f <family> [optional] -mt_add_seq <reference fasta> -e <e value> -b <bed file for masking>')
	except Exception as e:
		print(f"Error: {type(e).__name__} - {e}")
		sys.exit()
else:
	print(f"Error!!! Choose one of the following running modes: mtpt, mtpt_eval, mt, mt_eval. Quitting...")
	sys.exit()
	
################################
#II. Define functions
################################

#MASK gene/MTPT regions
def mask_fasta_with_bed(fasta_file, bed_file, output_file):
	sequences = SeqIO.index(fasta_file, "fasta")
	bed = open(bed_file).readlines()
	recs={}
	for k in sequences.keys():
		recs[k]=str(sequences[k].seq)
	for l in bed:
   		sequence_id = l.split()[0]
   		start = int(l.split()[1])
   		end = int(l.split()[2])
   		sequence = recs[sequence_id]
   		masked_sequence = sequence[:(start-1)] + "N" * (end - start + 1) + sequence[end:]
   		recs[sequence_id] = masked_sequence
	out=open(output_file,'w')
	for k in recs.keys():
		d=out.write('>'+k+'\n'+recs[k]+'\n')

def id2bed(ids,bed_file):
	id_dict={}
	for l in bed_file:
		id_dict[l.split()[-1]]=l
	filtered_bed=[id_dict[j] for j in ids]
	return(filtered_bed)


#takes in bed file of aligned region, order them, and filter for ones selected to output
#e.g.
#Pch_1	133169	133918	BA000042.1	329847	330607	Nicotiana_tabacum
#Pch_1	133907	134389	BA000042.1	330960	331502	Nicotiana_tabacum
#Pch_1	134385	136416	BA000042.1	331627	333652	Nicotiana_tabacum
#Pch_1	133169	133912	CP129458.1	165029	164277	Quercus_variabilis
#Pch_1	133169	133299	CP129458.1	220967	221096	Quercus_variabilis #this row will be removed
#Pch_1	133907	135460	CP129458.1	163913	162355	Quercus_variabilis
def aln_scaffolder(bedtxt,split_block=False):
	rawtable = [line.split('\t') for line in bedtxt]
	sortedtable = sorted(rawtable, key=lambda x: (x[3], int(x[1])))
	filteredtable=[]
	cursp=''
	end_pos=0
	for l in sortedtable:
		if l[3]!=cursp:
			#a new sp
			filteredtable.append(l)
			cursp=l[3]
			end_pos=int(l[2])
		else:#same sp
			if int(l[1])>end_pos-50:
				#allow for 50 bp overlap in the reference
				filteredtable.append(l)
				end_pos=int(l[2])
	consolidatedtable={}
	consolidated_block_count={}
	for l in filteredtable:
		try:
			consolidatedtable[l[3]].append(l)
			consolidated_block_count[l[3]]=consolidated_block_count[l[3]]+1
		except KeyError:
			consolidatedtable[l[3]]=[l]
			consolidated_block_count[l[3]]=1
	mode_consolidated_block_count=mode([consolidated_block_count[k] for k in consolidated_block_count.keys()])
	if split_block==True and mode_consolidated_block_count>1:
		#break this scaffold into smaller chunks
		for k in consolidated_block_count.keys():
			outputtable=['X']
			if consolidated_block_count[k]==mode_consolidated_block_count:
				outputtable.append(['\t'.join(i[:3]) for i in consolidatedtable[k]])
				break
	else:
		outputtable=[]
		for k in consolidatedtable.keys():
			refpos=[l[1]+'-'+l[2] for l in consolidatedtable[k]]
			targetpos=[l[4]+'-'+l[5] for l in consolidatedtable[k]]
			#outputtable.append(consolidatedtable[k][0][0]+'\t'+';'.join(refpos)+'\t'+consolidatedtable[k][0][3]+'\t'+';'.join(targetpos)+'\t'+'\t'.join(consolidatedtable[k][0][6:]))
			outputtable.append(consolidatedtable[k][0][0]+'\t'+'+'.join(refpos)+'\t'+consolidatedtable[k][0][3]+'\t'+'+'.join(targetpos)+'\t'+'\t'.join(consolidatedtable[k][0][6:]))
	return(outputtable)

def write_seq_from_bed_txt(bed_txt,output_handle):
	for l in bed_txt:
		l=str(l).split()
		#write header
		d=output_handle.write(f">{l[7]}|{l[6]}|{l[2]}_{l[3]}\n")
		tg_pos=l[3].split('+')
		tg_seq=''
		#write seq, merge multiple regions if needed
		for i in tg_pos:
			start=int(i.split('-')[0])
			end=int(i.split('-')[1])
			if start<end:
				tg_seq=tg_seq+str(ref_recs[l[2]].seq[(start-1):end])
			else:
				tg_seq=tg_seq+str(ref_recs[l[2]].seq[(end-1):start].reverse_complement())
		d=output_handle.write(tg_seq+'\n')
					
def seq2seq_ortho_extraction(seed_file,targets_file):
	S='makeblastdb -in '+targets_file+' -out '+sp+'.temp -dbtype nucl >/dev/null'
	os.system(S)
	S='blastn -task dc-megablast -query '+seed_file+' -db '+sp+'.temp -outfmt 6 -evalue 1e-20 | sort -k2,2 -k11,11n> '+sp+'.temp.blast'
	os.system(S)
	hits=open(sp+'.temp.blast').readlines()
	if len(hits)==0:return('')
	new_order = [0, 6, 7, 1, 8, 9, 11, 10]
	reordered_hits = []
	for line in hits:
		parts = line.strip().split('\t')
		subject = parts[1]
		reordered = [parts[i] for i in new_order]
		# Add new columns from dictionaries
		reordered.append(id2sp.get(subject, 'NA'))
		reordered.append(id2fam.get(subject, 'NA'))
		reordered_hits.append('\t'.join(reordered) + '\n')
	hits_bed_txt=aln_scaffolder(reordered_hits,split_block=False)
	return(hits_bed_txt)

def dedup_pt_ids(bed_id):
	GROUP_COL = 3
	EVAL_COL = 7
	BITSCORE = 6
	data_lines = [seq_loc[i] for i in bed_id]
	data_lines = [l.split() for l in data_lines]
	df = pd.DataFrame(data_lines)
	df.rename(columns={GROUP_COL: 'GroupID', EVAL_COL: 'evalue', BITSCORE: 'bitscore'}, inplace=True)
	df_max_per_group = (
		df.sort_values(by=['GroupID', 'evalue'], ascending=[True, True])
		.drop_duplicates(subset=['GroupID'], keep='first')
	)
	filtered_bed = df_max_per_group.iloc[:, -1].tolist()
	return [i.strip() for i in filtered_bed]
			
def filter_blast_results(bed_id, family, top_n=100):
	GROUP_COL = 3
	EVAL_COL = 7
	BITSCORE = 6
	FAMILY = 9
	START = 1
	END = 2
	data_lines = [seq_loc[i] for i in bed_id]
	data_lines = [l.split() for l in data_lines]
	df = pd.DataFrame(data_lines)
	df.rename(columns={START:'start', END: 'end', GROUP_COL: 'GroupID', EVAL_COL: 'evalue', BITSCORE: 'bitscore', FAMILY: 'family'}, inplace=True)
	df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
	df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
	df["start"] = pd.to_numeric(df["start"], errors="coerce")
	df["end"] = pd.to_numeric(df["end"], errors="coerce")
	df_sorted = df.sort_values(
		by=["evalue", "bitscore"],
		ascending=[True, False]
	)
	df_filtered = df_sorted.head(top_n)
	total_len=df["end"].max()-df["start"].min()
	df_filtered['perc'] = (df_filtered['end']-df_filtered['start'])/float(total_len)
	number_of_long_hit = (df_filtered['perc'] > 0.6).sum()
	df_completness = number_of_long_hit/len(df_filtered)
	df_filtered = pd.concat([
		df_sorted.head(top_n),
		df_sorted[df_sorted["family"] == family]
	]).drop_duplicates()
	df_final = df[df['GroupID'].isin(df_filtered['GroupID'])]
	filtered_bed = df_final.iloc[:, -1].tolist()
	return (df_completness,number_of_long_hit,[i.strip() for i in filtered_bed])

def find_intersect_bed(bedA,bedB):
	a_hits = bedA.intersect(bedB, wa=True, f=0.7)
	b_hits = bedA.intersect(bedB, wa=True, F=0.7)
	combined = a_hits.cat(b_hits, postmerge=False).sort()
	unique_lines = list(dict.fromkeys(map(str, combined)))
	intersect_unique_bed = pybedtools.BedTool('\n'.join(unique_lines), from_string=True)
	return(intersect_unique_bed)

def build_support_signal(df, bin_size=50):
    region_start = df.start.min()
    region_end = df.end.max()
    bins = np.arange(region_start, region_end, bin_size)
    support = np.zeros(len(bins), dtype=int)
    for i, b in enumerate(bins):
        b_end = b + bin_size
        support[i] = np.sum(
            (df.start < b_end) & (df.end > b)
        )
    return bins, support

def fit_hmm(support):
    """
    Fit a 2-state Poisson HMM
    """
    X = support.reshape(-1, 1)
    model = PoissonHMM(
        n_components=2,
        n_iter=200,
        tol=1e-4,
        verbose=False
    )
    model.fit(X)
    # Use lambda_ for PoissonHMM instead of means_
    rates = model.lambdas_.flatten()
    block_state = np.argmax(rates)   # higher rate = BLOCK
    bridge_state = 1 - block_state
    return model, block_state, bridge_state

def decode_states(model, support):
    X = support.reshape(-1, 1)
    states = model.predict(X)
    return states

def extract_blocks(bins, states, block_state, bin_size,
                   min_block_len=80):
    blocks = []
    in_block = False
    for i, s in enumerate(states):
        if s == block_state and not in_block:
            start = bins[i]
            in_block = True
        elif s != block_state and in_block:
            end = bins[i]
            if end - start >= min_block_len:
                blocks.append((start, end))
            in_block = False
    # Handle tail
    if in_block:
        end = bins[-1] + bin_size
        if end - start >= min_block_len:
            blocks.append((start, end))
    return blocks


def hmm_split_intervals(bed_id,bin_size=50,min_block_len=80):
	data_lines = [seq_loc[i] for i in bed_id]
	data_lines = [l.split() for l in data_lines]
	df = pd.DataFrame(data_lines)
	df.rename(columns={0:'chrom', 1:'start', 2: 'end'}, inplace=True)
	df["start"] = pd.to_numeric(df["start"], errors="coerce")
	df["end"] = pd.to_numeric(df["end"], errors="coerce")
	bins, support = build_support_signal(df, bin_size)
	# Edge case: too few bins
	if len(bins) < 10:
		return [(df.start.min(), df.end.max())]
	model, block_state, bridge_state = fit_hmm(support)
	states = decode_states(model, support)
	blocks = extract_blocks(
		bins, states, block_state,
		bin_size, min_block_len
	)
	# If HMM fails to split, return original
	if len(blocks) == 0:
		return [(df.start.min(), df.end.max())]
	return blocks


def find_short_candidates(df, max_len=80):
    lengths = df.end - df.start
    return df.index[lengths < max_len].tolist()

def proportional_support(candidate, evidence_df):
    c_start, c_end = candidate.start, candidate.end
    left_total = 0
    right_total = 0
    for _, j in evidence_df.iterrows():
        # must overlap candidate to count
        if j.end <= c_start or j.start >= c_end:
            continue
        # left side contribution
        left = max(0, min(j.end, c_start) - j.start)
        # right side contribution
        right = max(0, j.end - max(j.start, c_end))
        left_total += left
        right_total += right
    return left_total, right_total

def choose_merge_direction(candidate, left, right, evidence_df):
    left_prop, right_prop = proportional_support(candidate, evidence_df)
    if left_prop > right_prop:
        return "left"
    if right_prop > left_prop:
        return "right"
    # tie-breaker: merge to larger neighbor
    left_len = left.end - left.start
    right_len = right.end - right.start
    return "left" if left_len >= right_len else "right"

def merge_short_intervals(primary_bed, evidence_ids, max_len=80):
	df = pd.DataFrame([[j.chrom,j.start,j.end] for j in primary_bed])
	df.rename(columns={0:'chrom', 1:'start', 2: 'end'}, inplace=True)
	df["start"] = pd.to_numeric(df["start"], errors="coerce")
	df["end"] = pd.to_numeric(df["end"], errors="coerce")
	df = df.copy().reset_index(drop=True)
	data_lines = [seq_loc[i] for i in evidence_ids]
	data_lines = [l.split() for l in data_lines]
	evidence_df = pd.DataFrame(data_lines)
	evidence_df.rename(columns={0:'chrom', 1:'start', 2: 'end'}, inplace=True)
	evidence_df["start"] = pd.to_numeric(df["start"], errors="coerce")
	evidence_df["end"] = pd.to_numeric(df["end"], errors="coerce")
	to_drop = set()
	short_idxs = find_short_candidates(df, max_len)
	for j in short_idxs:
		if j in to_drop:continue
		candidate = df.loc[j]
		if j == 0:
			# merge right
			df.at[j + 1, "start"] = candidate.start
			to_drop.add(j)
			continue
		if j == len(df) - 1:
			# merge left
			df.at[j - 1, "end"] = candidate.end
			to_drop.add(j)
			continue
		left = df.loc[j - 1]
		right = df.loc[j + 1]
		direction = choose_merge_direction(
			candidate, left, right, evidence_df
		)
		if direction == "left":
			df.at[j - 1, "end"] = candidate.end
		else:
			df.at[j + 1, "start"] = candidate.start
		to_drop.add(j)
	df = df.drop(index=list(to_drop)).reset_index(drop=True)
	return pybedtools.BedTool.from_dataframe(df)

############
#Classification related function
def ncbi_ref_root(tree,reference_phylo):
	all_family=[node.name for node in tree]
	all_family = [i.split('|')[0] for i in all_family if not i.startswith('query')]
	all_family=list(set(all_family))
	query_tip = [node.name for node in tree if node.name.startswith('query')]
	reference_family=[node.name for node in reference_phylo]
	overlapping_fam=[i for i in all_family if i in reference_family]
	if len(overlapping_fam)<2:
		tree.set_outgroup(tree.get_midpoint_outgroup())
		return(tree)
	reference_phylo.prune(overlapping_fam)
	best_root=[]
	best_root_fam_num=1000
	for child in reference_phylo.children:
		if len(child.get_leaves())<best_root_fam_num:
			best_root=[node.name for node in child]
			best_root_fam_num = len(child.get_leaves())
	#if all members of best_root form a clade, root with all of them
	best_root_all_tips = [node.name for node in tree if node.name.split('|')[0] in best_root]
	if len(best_root_all_tips) == 1:
		if tree.check_monophyly(values=best_root_all_tips+query_tip, target_attr="name")[0]:
			#query nested within this family, set them together as outgroup
			tree.set_outgroup(tree.get_common_ancestor(best_root_all_tips+query_tip))
		else:	
			tree.set_outgroup(tree&best_root_all_tips[0])
		return(tree)
	elif tree.check_monophyly(values=best_root_all_tips, target_attr="name")[0]:
		if tree.check_monophyly(values=best_root_all_tips+query_tip, target_attr="name")[0]:
			#query nested within this family, set them together as outgroup
			tree.set_outgroup(tree.get_common_ancestor(best_root_all_tips+query_tip))
		else:
			tree.set_outgroup(tree.get_common_ancestor(best_root_all_tips))
		return(tree)
	#else, root with one of the best_root family
	else:
		mono=0
		for one_family in best_root:
			best_root_fam_tips = [node.name for node in tree if node.name.split('|')[0] == one_family]
			if len(best_root_fam_tips)==1:
				tree.set_outgroup(tree&best_root_fam_tips[0])
				mono=1
				return(tree)
			else:
				#this is a node
				if tree.check_monophyly(values=best_root_fam_tips, target_attr="name")[0]:
					tree.set_outgroup(tree.get_common_ancestor(best_root_fam_tips))
					mono=1
					return(tree)
				elif tree.check_monophyly(values=best_root_fam_tips+query_tip, target_attr="name")[0]:
					#query nested within this family, set them together as outgroup
					tree.set_outgroup(tree.get_common_ancestor(best_root_fam_tips+query_tip))
					mono=1
					return(tree)
		if mono ==0:
			tree.set_outgroup(tree.get_midpoint_outgroup())
			return(tree)

def max_clade_support(query,tip_collection,tree):
	q_branch=tree&query
	best_support=0
	for nd in tree:
		if nd.name in tip_collection:
			nd.add_features(type="ingroup")
	for nd in tree.get_monophyletic(values=["ingroup"], target_attr="type"):
		if any(leaf.name.startswith('query|') for leaf in nd.get_leaves()):
			q_clade = nd
			if nd.support >best_support:
				best_support = nd.support
			break
	for child in q_clade.traverse('preorder'):
		if any(leaf.name.startswith('query|') for leaf in child.get_leaves()):
			allchildren_tips=[leaf.name for leaf in child]
			if child.support >best_support and (not child.is_leaf()) and (not all(j.startswith(("query", fam)) for j in allchildren_tips)):
				best_support = child.support
	return(best_support)

# Function to find the largest continuous block in a list around a specific element based on its name
#takes in a list of true or false value indicating where the element is a member of the group of interest
def find_max_block_around_element(lst, element_start, element_end, max_diff_count):
	diff_count = 0
	first_ingroup = element_start
	last_ingroup = element_end
	ingroup_id=[]
	#extend on the left side
	ingroup_id.append(element_start)
	for i in range(max(element_start-1,0),-1,-1):
		if lst[i]:  # When same as the target taxon rank
			#current_start = max(0, current_start - 1)  # Start searching from one element before
			diff_count=0
			first_ingroup=i
			ingroup_id.append(i)
		else:
			if diff_count>= max_diff_count:break
			diff_count += 1
	#extend on the right side
	ingroup_id.append(element_end)
	for i in range(min(element_end+1,len(lst)),len(lst)):
		if lst[i]:  # When same as the target taxon rank
			#current_end = min(len(lst), current_end+1)  # Start searching from one element before
			diff_count=0
			last_ingroup=i
			ingroup_id.append(i)
		else:
			if diff_count>= max_diff_count:break
			diff_count += 1
	return(first_ingroup, last_ingroup, list(set(ingroup_id)))

def GetSister(tree,target_sp,q_fam):
	q_branch=tree&target_sp
	#make sure the query branch is not the root
	if q_branch.get_ancestors()[0].is_root():
		ancestor=tree.get_midpoint_outgroup()
		tree.set_outgroup(ancestor)
	if not q_branch.get_ancestors()[0].is_root():
		for leaf in tree:
			leaf.add_features(family=leaf.name.split('|')[0])
		q_branch.add_features(family=q_fam)
		q_fam_branch=''
		for nd in tree.get_monophyletic(values=[q_fam], target_attr="family"):
			if target_sp in [leaf.name for leaf in nd]:
				q_fam_branch=nd
		receiver=[nd.name for nd in q_fam_branch]
		donor=[leaf.name for leaf in q_fam_branch.get_sisters()[0]]
		#get donor at family level
		donor_family=list(set([j.split('|')[0] for j in donor]))
		donor_genera=[j.split('|')[1] for j in donor]
		donor_genera=list(set([j.split('_')[0] for j in donor_genera]))
		bs=tree.get_common_ancestor(donor+[target_sp]).support
		return(receiver,donor_family,donor_genera,donor,str(bs))
	else:
		return('NA','NA','NA','NA','NA')


#take a tree and return true or false for VGT based on tip orders, this is a generous way to classify VGT and better accommodates topology error in a phylo tree
def VGTFromTipOrder(tree,target_sp,q_family):
	tips=[node.name for node in tree]
	q_index = tips.index(target_sp)
	ingroup_binary = [1 if j.split('|')[0] in ingroup else 0 for j in tips]
	ingroup_binary[q_index] = 1
	qfam_binary = [1 if j.split('|')[0] in [q_family] else 0 for j in tips]
	qfam_binary[q_index] = 1
	# find the query family block and the larger ingroup block
	qfam_start,qfam_end,qfam_index=find_max_block_around_element(qfam_binary, q_index, q_index, max_diff_count=1)
	sis_ingroup_start,sis_ingroup_end,sis_ingroup_index=find_max_block_around_element(ingroup_binary, q_index, q_index, max_diff_count=1)
	qfam_genera =[tips[j].split('|')[1] for j in qfam_index]
	qfam_genera = set([j.split('_')[0] for j in qfam_genera])
	sis_ingroup_tips=[tips[j] for j in sis_ingroup_index]
	sis_ingroup_fam=set([j.split('|')[0] for j in sis_ingroup_tips]+[q_family])
	#if this qfam cluster contains more than five genera of the query family, or the ingroup cluster contains more than three families, this is most likely a VGT
	if len(sis_ingroup_fam)>2 or len(qfam_genera)>4:
		return(1)
	#check if this cluster is nested within a larger ingroup cluster
	else:
		output=0
		if qfam_start!=0:
			if ingroup_binary[qfam_start-1]==1:output=1
		if qfam_end!=(len(tips)-1):
			if ingroup_binary[qfam_end+1]==1:output=1
		return(output)	

def correct_nd_support(tree):
	#FastTree does not report branch support for super short branches or when it is close to 0, ete3 store those value as node.name instead of node.support, needs to be corrected
	for nd in tree.traverse():
		if not nd.is_leaf():
			if nd.name:nd.support = float(nd.name)  # real support
			else:nd.support = 0.0 # or None
			nd.name = ""
	return(tree)

def coverage_fragmentation_metric(aln):
	coverage={}
	fragmentation={}
	for rec in aln:
		fragmentation[rec.id]=rec.id.count('+')+1
		coverage[rec.id]=1-(float(rec.seq.count('-'))/len(rec.seq))
	return(coverage, fragmentation)

   
####################################
#III. MTPT mode
###################################
if args.m =='mtpt':
	#blast
	if args.e:evalue=args.e
	else:evalue=1e-40
	print(str(datetime.datetime.now())+'\tPerforming BLAST to identify candidate MTPT with e-value threshold of '+str(evalue))
	#make the query sequence as the database otherwise many potential hits from the Viridiplantae get ignored by blastn if the other way around
	os.system(f"makeblastdb -in {query} -out {sp}.mtpt -dbtype nucl >/dev/null")
	S='blastn -task dc-megablast -query '+sp+'.pt_db.fas -db '+sp+'.mtpt -outfmt 6 -evalue '+str(evalue)+' >'+sp+'.mtpt.blast'
	os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
	#add taxon information to blast hit
	out=open(sp+'.mtpt.blast.taxon','w')
	id2sp={}
	id2fam={}
	taxon_ref=open(script_path+'/database/pt_Viridiplantae_taxonomy.tsv').readlines()
	for l in taxon_ref[1:]:
		id2sp[l.split()[0]]=l.split('\t')[7].strip()
		id2fam[l.split()[0]]=l.split('\t')[5]
	x=open(sp+'.mtpt.blast').readlines()
	for l in x:
		try:
			ncbi_id=l.split()[0]
			ncbi_id=l.split('.')[0]
			d=out.write(l.strip()+'\t'+id2sp[ncbi_id]+'\t'+id2fam[ncbi_id]+'\n')
		except KeyError:
			#try to get taxon info from custom reference species
			if '|' in l.split()[0]:d=out.write(l.strip()+'\t'+l.split('|')[1]+'\t'+l.split('|')[0]+'\n')
			else:d=out.write(l.strip()+'\t'+l.split()[0]+'\tNA\n')
	out.close()
	#sort blast results and give each row an unique id
	S="awk -v OFS='\\t' '{if ($9 <= $10) print $2, $9, $10, $1, $7, $8, $12, $11, $13, $14; else print $2, $10, $9, $1, $8, $7, $12, $11, $13, $14}' "+sp+".mtpt.blast.taxon| sort -k1,1 -k2,2n -k4,4n -k8,8g | awk -v FS='\\t' -v OFS='\\t' '{print $0, NR}' > "+sp+".mtpt.bed"
	os.system(S)
	os.system(f"rm {sp}.mtpt.blast.taxon")
	#Merge MTPT blocks
	hits=open(sp+'.mtpt.bed').readlines()
	loci=pybedtools.BedTool(''.join(hits), from_string=True).merge(c=11,o='collapse')
	out=open(sp+'.mtpt.merged.bed','w')
	d=out.write(''.join([str(i) for i in loci]))
	out.close()
	q_recs=SeqIO.index(query,'fasta')
	ref_recs=SeqIO.index(sp+'.pt_db.fas', 'fasta')
	seq_loc={}
	for l in hits:seq_loc[l.split()[-1]]=l
	order=1
	retained_order=[]
	sum_out=open(sp+'.mtpt.sum.tsv','w')
	sum_out.write('Locus_ID\tQuery\tStart\tEnd\tClassification\tSupport\tSister_family\tSister_genus\tSister_species\n')
	output_dir = sp+'_HGTscanner_supporting_files'
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	order = 1
	for l in loci:
		ids=l.fields[3]
		ids=ids.split(',')
		#remove duplicate loci from the inverted repeat region, retain only one per species
		ids=dedup_pt_ids(ids)
		#If too few plastid hits for this loci (unlikely because plastid synteny is well conserved across order), skip. In practice, I found a lot of these were ancestral mt transfers.
		if len(ids)<11:
			order=order+1
			continue
		out=open(output_dir+'/'+sp+'.mtpt.'+str(order)+'.fas','w')
		#write target
		#if too many blast hits, only retain the top 200 hits to save computational time. I also notice the IR region was included twice for a single species. Thus only one best hit will be selected from each species.
		hit_num = 200
		if args.hit:hit_num=args.hit
		if len(ids)>hit_num:
			_, _, ids = filter_blast_results(ids,fam,hit_num)
		for i in ids:
			line=seq_loc[i]
			id=line.split()[3]
			start=min(int(line.split()[4]),int(line.split()[5]))
			end=max(int(line.split()[4]),int(line.split()[5]))
			try:d=out.write('>'+id2fam[id.split('.')[0]]+'|'+id2sp[id.split('.')[0]]+'|'+id+'_'+str(start)+'_'+str(end)+'\n')
			except KeyError:d=out.write('>'+id+'_'+str(start)+'_'+str(end)+'\n')
			d=out.write(str(ref_recs[id].seq[start-1:end])+'\n')
		#write query
		start=min([int(seq_loc[i].split()[1]) for i in ids])-1
		end=max([int(seq_loc[i].split()[2]) for i in ids])
		q_rec_out = q_recs[l.chrom][start:end]
		q_rec_out.id = "query|" + q_rec_out.id
		q_rec_out.description = ""
		d = SeqIO.write(q_rec_out, out, 'fasta')
		#d=SeqIO.write(q_recs[l.split()[0]][(int(l.split()[1])-1):int(l.split()[2])],out,'fasta')
		out.close()
		sum_out.write(f"{order}\t{l.chrom}\t{start}\t{end}\tTBD\tTBD\tTBD\tTBD\tTBD\n")
		retained_order.append(order)
		order=order+1
	sum_out.close()
	print(str(datetime.datetime.now())+'\tFound '+ str(len(retained_order))+' potential MTPT sequences for '+sp+' from the original '+str(order-1)+' loci showing plastid seq homology')
	#alignment and phylogenetic reconstruction
	if args.notree:
		print(str(datetime.datetime.now())+'\tStart alignment for '+str(len(retained_order))+' loci with more than 10 taxa. May take a while...')
	else:
		print(str(datetime.datetime.now())+'\tStart alignment and phylogeny for '+str(len(retained_order))+' loci with more than 10 taxa. May take a while...')
	for i in retained_order:
		print(f"{str(datetime.datetime.now())}\tAnalyzing Loci #{i}", end='\r')
		S="mafft --quiet --adjustdirection "+output_dir+"/"+sp+".mtpt."+str(i)+".fas | sed 's/_R_//g' > "+output_dir+"/"+sp+".mtpt."+str(i)+".aln.fas"
		os.system(S)
		if args.notree:
			continue
		else:
			S=f"FastTree -quiet -nt {output_dir}/{sp}.mtpt.{i}.aln.fas >{output_dir}/{sp}.mtpt.{i}.aln.fas.treefile"
			#S=f"iqtree2 -T 1 -B 1000 -s {output_dir}/{sp}.mtpt.{i}.aln.fas"
			subprocess.run(S, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if args.notree:
		print(str(datetime.datetime.now())+'\tAlignment complete')
		print(str(datetime.datetime.now())+'\tCleaning intermediate files')
		os.system(f"rm {sp}.mtpt.n*")
		os.system(f"rm {sp}.pt_db.fas")
		sys.exit(str(datetime.datetime.now())+'\tYou have chosen to complete phylogeny separately. Exit...')
	print(str(datetime.datetime.now())+'\tAlignment and FastTree phylogeny complete')
	print(str(datetime.datetime.now())+'\tCleaning intermediate files')
	os.system(f"rm {sp}.mtpt.n*")
	os.system(f"rm {sp}.pt_db.fas")
	#######################################
	#classification based on the fasttree phylogeny
	current_time = datetime.datetime.now()
	print(f"{current_time}\tReading and processing phylogenies from {len(retained_order)} loci in {sp}.mtpt.sum.tsv")
	sum_txt=open(sp+'.mtpt.sum.tsv').readlines()
	new_sum_txt = []
	new_sum_txt.append(sum_txt[0])
	for line in sum_txt[1:]:
		id=line.split()[0]
		q='query|'+line.split()[1]
		try:
			t=Tree(output_dir+"/"+sp+'.mtpt.'+id+'.aln.fas.treefile',format=1)
			t=correct_nd_support(t)
			#check branch length of query first, if too long, certainly an ancestral mt transfer and not a young mtpt
			q_branch=t&q
			branch_lengths = [node.dist for node in t.traverse() if not node.is_root()]
			z_score = abs(q_branch.dist - stats.mean(branch_lengths)) / stats.stdev(branch_lengths)
			#print(id, z_score, q_branch.dist)
			if q_branch.dist>0.15 and z_score >10:
				raw_support = q_branch.get_ancestors()[0].support
				sister_nd = q_branch.get_sisters()[0]
				sister_tips=[leaf.name for leaf in sister_nd]
				sister_family = list(set([i.split('|')[0] for i in sister_tips]))
				sister_genera=[]
				for i in sister_tips:
					try:
						genus=i.split('|')[1]
						genus=genus.split('_')[0]
						sister_genera.append(genus)
					except IndexError:
						sister_genera.append(i.split('_')[0])
				new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (high confidence)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				continue
			#now these look like true mtpt
			#Rooting based on the ncbi common taxonomy
			ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
			t=ncbi_ref_root(t,ncbi_tree)
			#print(t.write())
			#Check sister clade of the quqery
			q_branch=t&q
			sister_nd = q_branch.get_sisters()[0]
			sister_tips=[leaf.name for leaf in sister_nd]
			sister_family = list(set([i.split('|')[0] for i in sister_tips]))
			try:sister_family.remove('NA')
			except ValueError:pass
			sister_genera=[]
			for i in sister_tips:
				try:
					genus=i.split('|')[1]
					genus=genus.split('_')[0]
					sister_genera.append(genus)
				except IndexError:
					sister_genera.append(i.split('_')[0])
			sister_genera = list(set([i.split('_')[0] for i in sister_genera]))
			raw_support = q_branch.get_ancestors()[0].support
			#Check the neighbors of (query+sister)
			alltips=[j.name for j in t]
			qs_index = [j for j, name in enumerate(alltips) if name in [q] + sister_tips]
			ingroup_binary = [1 if j.split('|')[0] in ingroup else 0 for j in alltips]
			ingroup_start,ingroup_end,ingroup_index=find_max_block_around_element(ingroup_binary, min(qs_index), max(qs_index), max_diff_count=1)
			nested_in_ingroup = 0
			if ingroup_start<min(qs_index) and ingroup_end>max(qs_index):nested_in_ingroup = 1
			if all(i in ingroup for i in sister_family):
				#native MTPT: sister clade is ingroup
				if len(sister_family)==1:
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tnative MTPT\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
			elif len([i for i in sister_family if i in ingroup])==0:
				#sister clade is outgroup
				if nested_in_ingroup == 1:
					#but nested within ingroup: inconclusive
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tinconclusive\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:
					#do not nest within ingroups
					if len(sister_family)==1:
						#check the neighbor of the q+sister
						donor_fam_binary = [1 if j.split('|')[0] in sister_family else 0 for j in alltips]
						donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
						nested_in_donor_fam = 0
						if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
						if nested_in_donor_fam:
							#alien MTPT: a single donor family is involved and the q+sister nest within other species of the same donor family
							if donor_fam_end==len(alltips):donor_fam_tips = alltips[donor_fam_start:]
							else:donor_fam_tips = alltips[donor_fam_start:donor_fam_end+1]
							fam_support = max_clade_support(q,donor_fam_tips,t)
							new_sum_txt.append('\t'.join(line.split()[0:4])+f"\thigh confidence alien MTPT\t{fam_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
						else:
							#does not nest with other species of the same family: putative mt transfer
							new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
					elif len(sister_family)<6:
						#putative alien MTPT: no more than five donor family is involved
						new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
					else:
						#ancestral mt transfer: too many donor family involved
						new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (high confidence)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
			else:
				#inconclusive: sister clade is ingroup+outgroup
				if len(sister_family)>5:
					#too many sister family likely ancestral mt transfer
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tinconclusive\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
		except ete3.parser.newick.NewickError:
			new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tTree not found\tNA\tNA\tNA\tNA\n")
	#update the summary file
	sum_out = open(sp+'.mtpt.sum.tsv','w')
	d=sum_out.write(''.join(new_sum_txt))
	#print(''.join(new_sum_txt))
	print(str(datetime.datetime.now())+'\tCompleted evaluation of MTPT source. See summary file in '+sp+'.mtpt.sum.tsv')
	
##################
elif args.m =='mtpt_eval':
	#Evaluate the source of the region and update the summary file
	current_time = datetime.datetime.now()
	sum_txt = open(f"{args.wd}/{sp}.mtpt.sum.tsv").readlines()
	print(f"{current_time}\tReading and processing phylogenies from {len(sum_txt)-1} loci in {args.wd}/{sp}.mtpt.sum.tsv")
	output_dir = f"{args.wd}/{sp}_HGTscanner_supporting_files"
	new_sum_txt = []
	new_sum_txt.append(sum_txt[0])
	for line in sum_txt[1:]:
		id=line.split()[0]
		q='query|'+line.split()[1]
		try:
			t=Tree(output_dir+"/"+sp+'.mtpt.'+id+'.aln.fas.treefile',format=1)
			t=correct_nd_support(t)
			#check branch length of query first, if too long, certainly an ancestral mt transfer and not a young mtpt
			q_branch=t&q
			branch_lengths = [node.dist for node in t.traverse() if not node.is_root()]
			z_score = abs(q_branch.dist - stats.mean(branch_lengths)) / stats.stdev(branch_lengths)
			#print(id, z_score, q_branch.dist)
			if q_branch.dist>0.15 and z_score >10:
				raw_support = q_branch.get_ancestors()[0].support
				sister_nd = q_branch.get_sisters()[0]
				sister_tips=[leaf.name for leaf in sister_nd]
				sister_family = list(set([i.split('|')[0] for i in sister_tips]))
				sister_genera=[]
				for i in sister_tips:
					try:
						genus=i.split('|')[1]
						genus=genus.split('_')[0]
						sister_genera.append(genus)
					except IndexError:
						sister_genera.append(i.split('_')[0])
				new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (high confidence)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				continue
			#now these look like true mtpt
			#Rooting based on the ncbi common taxonomy
			ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
			t=ncbi_ref_root(t,ncbi_tree)
			#print(t.write())
			#Check sister clade of the quqery
			q_branch=t&q
			sister_nd = q_branch.get_sisters()[0]
			sister_tips=[leaf.name for leaf in sister_nd]
			sister_family = list(set([i.split('|')[0] for i in sister_tips]))
			try:sister_family.remove('NA')
			except ValueError:pass
			sister_genera=[]
			for i in sister_tips:
				try:
					genus=i.split('|')[1]
					genus=genus.split('_')[0]
					sister_genera.append(genus)
				except IndexError:
					sister_genera.append(i.split('_')[0])
			sister_genera = list(set([i.split('_')[0] for i in sister_genera]))
			raw_support = q_branch.get_ancestors()[0].support
			#Check the neighbors of (query+sister)
			alltips=[j.name for j in t]
			qs_index = [j for j, name in enumerate(alltips) if name in [q] + sister_tips]
			ingroup_binary = [1 if j.split('|')[0] in ingroup else 0 for j in alltips]
			ingroup_start,ingroup_end,ingroup_index=find_max_block_around_element(ingroup_binary, min(qs_index), max(qs_index), max_diff_count=1)
			nested_in_ingroup = 0
			if ingroup_start<min(qs_index) and ingroup_end>max(qs_index):nested_in_ingroup = 1
			if all(i in ingroup for i in sister_family):
				#native MTPT: sister clade is ingroup
				if len(sister_family)==1:
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tnative MTPT\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
			elif len([i for i in sister_family if i in ingroup])==0:
				#sister clade is outgroup
				if nested_in_ingroup == 1:
					#but nested within ingroup: inconclusive
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tinconclusive\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:
					#do not nest within ingroups
					if len(sister_family)==1:
						#check the neighbor of the q+sister
						donor_fam_binary = [1 if j.split('|')[0] in sister_family else 0 for j in alltips]
						donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
						nested_in_donor_fam = 0
						if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
						if nested_in_donor_fam:
							#alien MTPT: a single donor family is involved and the q+sister nest within other species of the same donor family
							if donor_fam_end==len(alltips):donor_fam_tips = alltips[donor_fam_start:]
							else:donor_fam_tips = alltips[donor_fam_start:donor_fam_end+1]
							fam_support = max_clade_support(q,donor_fam_tips,t)
							new_sum_txt.append('\t'.join(line.split()[0:4])+f"\thigh confidence alien MTPT\t{fam_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
						else:
							#does not nest with other species of the same family: putative mt transfer
							new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
					elif len(sister_family)<6:
						#putative alien MTPT: no more than five donor family is involved
						new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
					else:
						#ancestral mt transfer: too many donor family involved
						new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (high confidence)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
			else:
				#inconclusive: sister clade is ingroup+outgroup
				if len(sister_family)>5:
					#too many sister family likely ancestral mt transfer
					new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tancestral mt transfer (putative)\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
				else:new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tinconclusive\t{raw_support}\t{','.join(sister_family)}\t{','.join(sister_genera)}\t{','.join(sister_tips)}\n")
		except ete3.parser.newick.NewickError:
			new_sum_txt.append('\t'.join(line.split()[0:4])+f"\tTree not found\tNA\tNA\tNA\tNA\n")
	#update the summary file
	sum_out = open(f"{sp}.mtpt.sum.tsv",'w')
	d=sum_out.write(''.join(new_sum_txt))
	sum_out.close()
	print(str(datetime.datetime.now())+f"\tCompleted evaluation of MTPT source based on supporting information from {args.wd}/{sp}_HGTscanner_supporting_files. See summary file in {sp}.mtpt.sum.tsv")

#####################################
#IV. MT mode for HGT detection in mito
elif args.m == 'mt':
	if args.b:
		mask_bed=args.b
		print(str(datetime.datetime.now())+'\tMasking query sequence '+query+' using the bed file: '+args.b)
		mask_fasta_with_bed(query, mask_bed, sp+'.mask.fas')
		#make blast db
		S='makeblastdb -in '+sp+'.mask.fas -out '+sp+'.mt -dbtype nucl >/dev/null'
		os.system(S)
	else:
		S='makeblastdb -in '+query+' -out '+sp+'.mt -dbtype nucl >/dev/null'
		os.system(S)
	###################
	#BLAST homology and scaffolding
	##################
	evalue=1e-20
	if args.e:evalue=args.e
	print(str(datetime.datetime.now())+'\tPerforming BLAST to identify candidate HGT with evalue threshold of '+str(evalue))
	S='blastn -task dc-megablast -query '+sp+'.mt_db.fas -db '+sp+'.mt -outfmt 6 -evalue '+str(evalue)+' >'+sp+'.mt.blast'
	#os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
	#Add taxonomic information to each hit at species and family level to help identify syntenic blocks
	x=open(sp+'.mt.blast').readlines()
	out=open(sp+'.mt.blast.taxon','w')
	id2sp={}
	id2fam={}
	taxon_ref=open(script_path+'/database/mt_Viridiplantae_taxonomy.tsv').readlines()
	for l in taxon_ref:
		id2sp[l.split('\t')[0]]=l.split('\t')[-1].strip()
		id2fam[l.split('\t')[0]]=l.split('\t')[6]
	for l in x:
		try:
			d=out.write(l.strip()+'\t'+id2sp[l.split()[0]]+'\t'+id2fam[l.split()[0]]+'\n')
		except KeyError:
			#try to get taxon info from custom reference species
			if '|' in l.split()[0]:d=out.write(l.strip()+'\t'+l.split('|')[1]+'\t'+l.split('|')[0]+'\n')
			else:d=out.write(l.strip()+'\t'+l.split()[0]+'\tNA\n')
	out.close()
	#sort blast results, give each row a unique id, only hits <20k are printed
	#S="cat "+sp+".taxon.blast | sort -k1,1 -k7,7n -k11,11n | awk -v OFS='\\t' '{if ($8-$7 <= 20000) print $1, $7, $8, $2, $9, $10, $13, $14, NR}' > "+sp+".taxon_sorted.bed"
	S="awk -v OFS='\\t' '{if ($9 <= $10) print $2, $9, $10, $1, $7, $8, $12, $11, $13, $14; else print $2, $10, $9, $1, $8, $7, $12, $11, $13, $14}' "+sp+".mt.blast.taxon| sort -k1,1 -k2,2n -k4,4n -k8,8g | awk -v FS='\\t' -v OFS='\\t' '{print $0, NR}' > "+sp+".mt.bed"
	os.system(S)
	os.system('rm '+sp+'.mt.blast.taxon')
	#define potential HGT blocks
	#classify these hits based on source, here, alignments from the ingroups will be filtered to merge synteny blocks
	#This is necessary because they create long synteny blocks that may be broken by HGT in the middle
	x=open(sp+".mt.bed").readlines()
	otherfam=[]
	samefam=[]
	for l in x:
		if not l.split()[9] in ingroup:
			otherfam.append(l)
		else:
			samefam.append(l)
	#merge BLAST hits, but require at least 50 bp overlap
	otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=11,o='collapse',d=-50)
	samefam_bed=pybedtools.BedTool(''.join(samefam), from_string=True)
	out=open(sp+'.mt.merged.bed','w')
	d=out.write(str(otherfam_merged))
	out.close()
	print(str(datetime.datetime.now())+'\tFound '+str(len(otherfam_merged))+' homologous genetic blocks for further examination')
	#extract sequences for each block
	q_recs=SeqIO.index(query,'fasta')
	ref_recs=SeqIO.index(sp+'.mt_db.fas', 'fasta')
	order=1
	num=1
	mapout=open(sp+'.alnmap.bed','w')
	hits=open(sp+'.mt.bed').readlines()
	seq_loc={}
	for l in hits:seq_loc[l.split()[-1]]=l
	if not os.path.isdir(sp+'_HGTscanner_supporting_files'):os.mkdir(sp+'_HGTscanner_supporting_files')
	for hit in otherfam_merged:
		#gather overlapping alignment from both other families and close relatives
		ids=hit.fields[3]
		ids=ids.split(',')
		#If too many hits (>300), select the best ones
		hit_num = 300
		if args.hit:hit_num=args.hit
		completeness, num_long_hit, new_ids = filter_blast_results(ids,fam,hit_num)
		###Evaluate whether blocks needs to be further divided up
		if completeness < 0.15 and num_long_hit < 6 and len(new_ids)>100 and hit.end-hit.start> 600:
			#This block needs to be split up
			to_spilt=1
			#use HMM to identify high confidence homology ranges
			major_fragments = hmm_split_intervals(new_ids)
			#print(order, major_fragments)
			if major_fragments[-1][1]>hit.end:major_fragments[-1]=(major_fragments[-1][0],hit.end)
			#if there is only 1 block and both ends close to the raw edges
			if len(major_fragments)==1:
				if major_fragments[0][0]-hit.start<100 and hit.end-major_fragments[0][1]<100:to_spilt=0
			if to_spilt:
				#find some small blocks
				new_individual_bed_txt = '\n'.join([f"{hit.fields[0]}\t{j[0]}\t{j[1]}" for j in major_fragments])
				new_individual_bed = pybedtools.BedTool(new_individual_bed_txt, from_string=True)  # internal blocks
				hit_bed = pybedtools.BedTool([hit])
				remaining_bed = hit_bed.subtract(new_individual_bed)
				all_new_bed = remaining_bed.cat(new_individual_bed, postmerge=False).sort()
				#merge small intervals <100 bp
				all_new_bed = merge_short_intervals(all_new_bed,new_ids,100)
				#rerank best hits within each individual block
				raw_beds_txt=[seq_loc[j] for j in ids]
				raw_beds=pybedtools.BedTool(''.join(raw_beds_txt), from_string=True)
				for block_i in all_new_bed:
					out=open(sp+'_HGTscanner_supporting_files'+'/'+sp+'.hgt.'+str(order)+'.fas','w')
					subhit=pybedtools.BedTool([block_i])
					#filter for region that either covers >70% of the query seq or have >70% of its own length within the target loci
					new_raw_beds = find_intersect_bed(raw_beds,subhit)
					new_raw_beds_ids = str(new_raw_beds).split('\n')
					new_raw_beds_ids = [j.split()[-1] for j in new_raw_beds_ids[:-1]]
					#filter for top hits
					_, _, new_raw_beds_ids = filter_blast_results(new_raw_beds_ids,fam,hit_num)
					#rescaffold
					if len(new_raw_beds_ids)>1:
						new_raw_beds_txt = [seq_loc[j] for j in new_raw_beds_ids]
						new_raw_beds=aln_scaffolder(new_raw_beds_txt,split_block=False)
					else:
						j=str(new_raw_beds).split()
						new_raw_beds=[f"{j[0]}\t{j[1]}-{j[2]}\t{j[3]}\t{j[4]}-{j[5]}\t{j[6]}\t{j[7]}\t{j[8]}\t{j[9]}\t{j[10]}"]
					#write other family
					#min_start=100000000000
					#max_end=0
					for l in new_raw_beds:
						l=str(l).split()
						#query_pos=l[1].split('-')
						#if int(query_pos[0])<min_start:min_start=int(query_pos[0])
						#if int(query_pos[-1])>max_end:max_end=int(query_pos[-1])
						#write sequence
						d=out.write(f">{l[7]}|{l[6]}|{l[2]}_{l[3]}\n")
						tg_pos=l[3].split('+')
						tg_seq=''
						for i in tg_pos:
							start=int(i.split('-')[0])
							end=int(i.split('-')[1])
							if start<end:
								tg_seq=tg_seq+str(ref_recs[l[2]].seq[(start-1):end])
							else:
								tg_seq=tg_seq+str(ref_recs[l[2]].seq[(end-1):start].reverse_complement())
						d=out.write(tg_seq+'\n')
					#write query
					q_pos=str(block_i).split()
					q_rec_out = q_recs[l[0]][(int(q_pos[1])-1):int(q_pos[2])]
					q_rec_out.id = "query|" + q_rec_out.id
					q_rec_out.description = ""
					d = SeqIO.write(q_rec_out, out, 'fasta')
					#write close relative
					samefam_hit = find_intersect_bed(samefam_bed,subhit)
					if not samefam_hit.count()==0:
						d=SeqIO.write(q_recs[l[0]][(int(q_pos[1])-1):int(q_pos[2])],sp+'.tempseed.fas','fasta')
						out2=open(sp+'.tempTarget.fas','w')
						for ll in samefam_hit:
							d=SeqIO.write(ref_recs[ll.fields[3]],out2,'fasta')
						out2.close()
						samefam_bed_txt=seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas')
						if not samefam_bed_txt=='':write_seq_from_bed_txt(samefam_bed_txt,out)
					d=mapout.write(hit.chrom+'\t'+q_pos[1]+'\t'+q_pos[2]+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
					order=order+1
					out.close()
			else:
				#HMM method was unable to divide it, output the whole chunk
				raw_beds_txt=[seq_loc[j] for j in new_ids]
				seqout_beds=aln_scaffolder(raw_beds_txt,split_block=False)
				out=open(sp+'_HGTscanner_supporting_files'+'/'+sp+'.hgt.'+str(order)+'.fas','w')
				#write other family
				write_seq_from_bed_txt(seqout_beds,out)
				#write query
				q_rec_out = q_recs[hit.chrom][(hit.start-1):hit.end]
				q_rec_out.id = "query|" + q_rec_out.id
				q_rec_out.description = ""
				d = SeqIO.write(q_rec_out, out, 'fasta')
				#write close relatives
				subhit=pybedtools.BedTool(str(hit),from_string=True)
				samefam_hit = find_intersect_bed(samefam_bed,subhit)
				if not samefam_hit.count()==0:
					#write sequences from close relatives if they are not empty
					d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],sp+'.tempseed.fas','fasta')
					out2=open(sp+'.tempTarget.fas','w')
					for l in samefam_hit:d=SeqIO.write(ref_recs[l.fields[3]],out2,'fasta')
					out2.close()
					samefam_bed_txt=seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas')
					if not samefam_bed_txt=='':write_seq_from_bed_txt(samefam_bed_txt,out)
				d=mapout.write(hit.chrom+'\t'+str(hit.start)+'\t'+str(hit.end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
				order=order+1
				out.close()
		else:
			#This block can be directly output
			raw_beds_txt=[seq_loc[j] for j in new_ids]
			seqout_beds=aln_scaffolder(raw_beds_txt,split_block=False)
			out=open(sp+'_HGTscanner_supporting_files'+'/'+sp+'.hgt.'+str(order)+'.fas','w')
			#write other family
			write_seq_from_bed_txt(seqout_beds,out)
			#write query
			q_rec_out = q_recs[hit.chrom][(hit.start-1):hit.end]
			q_rec_out.id = "query|" + q_rec_out.id
			q_rec_out.description = ""
			d = SeqIO.write(q_rec_out, out, 'fasta')
			#write close relatives
			subhit=pybedtools.BedTool(str(hit),from_string=True)
			samefam_hit = find_intersect_bed(samefam_bed,subhit)
			if not samefam_hit.count()==0:
				#write sequences from close relatives if they are not empty
				d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],sp+'.tempseed.fas','fasta')
				out2=open(sp+'.tempTarget.fas','w')
				for l in samefam_hit:d=SeqIO.write(ref_recs[l.fields[3]],out2,'fasta')
				out2.close()
				samefam_bed_txt=seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas')
				if not samefam_bed_txt=='':write_seq_from_bed_txt(samefam_bed_txt,out)
			d=mapout.write(hit.chrom+'\t'+str(hit.start)+'\t'+str(hit.end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
			order=order+1
			out.close()
		current_time = datetime.datetime.now()
		print(f"{current_time}\tExtracting sequences from homologous genetic block #{num}", end='\r')
		num=num+1
	mapout.close()		
	print(f"{current_time}\tA total of #{order-1} aligned sequences from #{num-1} merged homologous genetic blocks were extracted.", end='\r')
	#organize files
	os.system('rm '+sp+'.temp*')
	os.system('rm '+sp+'.mt.n*')
	os.system('rm '+sp+'.mt_db.fas')
	#os.system('mv '+sp+'.hgt.*.fas '+sp+'_HGTscanner_supporting_files')
	#############
	#alignment and phylogenetic reconstruction
	#############
	if args.notree:
		print(f"{current_time}\tStart alignment for #{order-1} loci")
	else:
		print(f"{current_time}\tStart alignment and FastTree phylogeny for #{order-1} loci")
	for i in range(1,order):
		print(f"{current_time}\tWorking on loci #{i}...", end='\r')
		#alignment
		fasta_in = f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.fas"
		temp1 = f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.temp1.fas"  # sequences starting with 'query|'
		temp2 = f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.temp2.fas"  # all other sequences
		aln_out = f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas"
		with open(temp1, 'w') as t1, open(temp2, 'w') as t2:
			for rec in SeqIO.parse(fasta_in, "fasta"):
				if rec.id.startswith("query|"):d=SeqIO.write(rec, t1, "fasta")
				else:d=SeqIO.write(rec, t2, "fasta")
		os.system(f"mafft --quiet --adjustdirection --6merpair --keeplength --addfragments {temp2} {temp1} | sed 's/_R_//g' > {aln_out}")
		#os.system(f"mafft --quiet --adjustdirection --6merpair --addfragments {temp2} {temp1} | sed 's/_R_//g' > {aln_out}")
		os.remove(temp1)
		os.remove(temp2)
		records = list(SeqIO.parse(aln_out, "fasta"))
		records = records[1:] + [records[0]]
		SeqIO.write(records, aln_out, "fasta")
		if args.notree:
			continue
		else:
			S=f"FastTree -quiet -nt {sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas >{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas.treefile"
			subprocess.run(S, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	if args.notree:
		print(str(datetime.datetime.now())+'\tCompleted alignment. See sequence file in '+sp+'_HGTscanner_supporting_files')
		sys.exit(str(datetime.datetime.now())+'\tYou have chosen to complete phylogeny separately. Exit...')
	else:
		print(str(datetime.datetime.now())+'\tCompleted alignment and FastTree phylogeny. See sequence and tree file in '+sp+'_HGTscanner_supporting_files')
	#######################
	#Classify based on blast, synteny, and phylogeny
	#########################
	print(f"{current_time}\tStart evaluating hypothesis of HGT...")
	#classify mt HGT
	lines=open(sp+'.alnmap.bed').readlines()
	out=open(f"{sp}.hgt.sum.tsv","w")
	d=out.write('Query\tStart\tEnd\tAlignment\tClassification\tRecepient\tDonor_Family\tDonor_genera\tDonor_species\tMethod\tBS\tIngroup_seq_coverage\tIngroup_seq_fragment\tDonor_seq_coverage\tDonor_seq_fragment\n')
	for l in lines:
		i = l.split()[-1]
		i = i.split('.')[2]
		recs=SeqIO.parse(f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas","fasta")
		coverage_metric, fragmentation_metric = coverage_fragmentation_metric(recs)
		allsp=[j for j in coverage_metric.keys()]
		target_tip=[j for j in allsp if j.startswith('query')]
		target_tip=target_tip[0]
		allsp.remove(target_tip)
		families=[j.split('|')[0] for j in allsp]
		families=list(set(families))
		try:families.remove(fam)
		except ValueError:pass
		if not families==['NA']:
			try:families.remove('NA')
			except ValueError:pass
		ingroup_sp = [j for j in allsp if j.split('|')[0] in ingroup]
		ingroup_coverage = [coverage_metric[j] for j in ingroup_sp]
		ingroup_fragmentation = [fragmentation_metric[j] for j in ingroup_sp]
		non_ingroup_sp = [j for j in allsp if not j.split('|')[0] in ingroup]
		genera=[j.split('|')[1] for j in allsp if not j.startswith(fam)]
		genera=list(set([j.split('_')[0] for j in genera]))
		same_fam_sp_id=[j for j in allsp if j.startswith(fam)]
		#same_fam_sp = list(set([j.split('|')[1] for j in same_fam_sp_id]))
		#BLAST case 1: only one or two families are in the blast result
		if len(families)==0:
			#VGT: only the query family
			d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','NA','BLAST: homology only found in ingroup','NA','NA','NA','NA','NA'])+'\n')
		elif len(families)==1:
			if families[0] in ingroup:
				#VGT: the only other family is an ingroup
				d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','NA','BLAST: homology only found in ingroup','NA','NA','NA','NA','NA'])+'\n')
			else:
				#HGT: the only other family is not ingroup
				receiver=';'.join([target_tip]+same_fam_sp_id)
				donor_fam=families[0]
				donor_gen=';'.join(genera)
				donor_sp=';'.join(non_ingroup_sp)
				donor_coverage = [coverage_metric[j] for j in non_ingroup_sp]
				donor_fragmentation = [fragmentation_metric[j] for j in non_ingroup_sp]
				d=out.write(f"{l.strip()}\tHigh confidence HGT\t{receiver}\t{donor_fam}\t{donor_gen}\t{donor_sp}\tBLAST: exclusive homology in two families\tNA\t{median(ingroup_coverage)}\t{median(ingroup_fragmentation)}\t{median(donor_coverage)}\t{median(donor_fragmentation)}\n")
		#>=2 other families
		else:
			#BLAST case 2: no other ingroup is in the tree, but there are at least two other families
			if len(ingroup_sp)==0:
				#check tree for donor
				try:
					t=Tree(f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas.treefile",format=1)
					t=correct_nd_support(t)
					#reroot
					ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
					t=ncbi_ref_root(t,ncbi_tree)
					receiver,donor_fam,donor_gen,donor_tips,bs=GetSister(t,target_tip,fam)
					donor_coverage = [coverage_metric[j] for j in donor_tips]
					donor_fragmentation = [fragmentation_metric[j] for j in donor_tips]
					if not receiver=='NA':
						#check if the query is well nested within the donor family
						if len(donor_fam)==1:
							alltips=[node.name for node in t]
							qs_index = [j for j, name in enumerate(alltips) if name in receiver]
							donor_fam_binary = [1 if j.split('|')[0] in donor_fam else 0 for j in alltips]
							donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
							nested_in_donor_fam = 0
							if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
							if nested_in_donor_fam:
								d=out.write(l.strip()+'\t'+'\t'.join(['High confidence HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'BLAST+Phylogeny: no ingroup shows homology and nested in donor family',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
							else:
								d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'BLAST: no other ingroup shows homology',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
						else:
							d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'BLAST: no other ingroup shows homology',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
					else:
						#the query was placed on root by mid-point rooting
						d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join([target_tip]+same_fam_sp_id),';'.join(families),';'.join(genera),';'.join(non_ingroup_sp),'BLAST: check rooting problem',bs,'NA','NA','NA','NA'])+'\n')
				except ete3.parser.newick.NewickError:
					receiver=target_tip
					donor_fam=''
					donor_gen=''
					bs=''
					if len(allsp)<4:
						#too few tips to run a tree
						donor_fam=families
						donor_gen=genera
						bs='NA'
					d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',receiver,';'.join(donor_fam),';'.join(donor_gen),';'.join(non_ingroup_sp),'BLAST: check rooting problem',bs,'NA','NA','NA','NA'])+'\n')
			#Phylogeny to assess classification	
			else:
				try:
					#examine if this can be a VGT without rerooting the tree
					t=Tree(f"{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas.treefile",format=1)
					t=correct_nd_support(t)
					#VGT based on tip order: the query if nested within an ingroup cluster
					if VGTFromTipOrder(t,target_tip,fam):
						d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','NA','Phylogeny','NA','NA','NA','NA','NA'])+'\n')
					else:
						#root tree
						ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
						t=ncbi_ref_root(t,ncbi_tree)
						#use tree topology to make decisions
						receiver,donor_fam,donor_gen,donor_tips,bs=GetSister(t,target_tip,fam)
						donor_coverage = [coverage_metric[j] for j in donor_tips]
						donor_fragmentation = [fragmentation_metric[j] for j in donor_tips]
						if VGTFromTipOrder(t,target_tip,fam):
							d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','NA','Phylogeny','NA','NA','NA','NA','NA'])+'\n')
							continue
						if len(donor_fam)==1:
							if donor_fam[0] in ingroup:
								#VGT: query sister to an ingroup
								d=out.write(l.strip()+'\t'+'\t'.join(['VGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'Phylogeny',bs,'NA','NA','NA','NA'])+'\n')
							else:
								#check if the query is well nested within the donor family
								alltips=[node.name for node in t]
								qs_index = [j for j, name in enumerate(alltips) if name in receiver + donor_tips]
								donor_fam_binary = [1 if j.split('|')[0] in donor_fam else 0 for j in alltips]
								donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
								nested_in_donor_fam = 0
								if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
								if nested_in_donor_fam:
									#High confidence HGT: a single donor family is involved and the q+sister nest within other species of the same donor family
									if donor_fam_end==len(alltips):donor_fam_tips = alltips[donor_fam_start:]
									else:donor_fam_tips = alltips[donor_fam_start:donor_fam_end+1]
									fam_support = max_clade_support(target_tip,donor_fam_tips,t)
									d=out.write(l.strip()+'\t'+'\t'.join(['High confidence HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'Phylogeny: nested in donor family',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
								else:
									#Putative HGT: one family that's not close relative, but phylogeny is not completely well nested within
									d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'Phylogeny',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
						else:
							#multiple donor families
							d=out.write(l.strip()+'\t'+'\t'.join(['inconclusive',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),';'.join(donor_tips),'Phylogeny',bs,str(median(ingroup_coverage)),str(median(ingroup_fragmentation)),str(median(donor_coverage)),str(median(donor_fragmentation))])+'\n')
				except ete3.parser.newick.NewickError:
					receiver=target_tip
					d=out.write(l.strip()+'\t'+'\t'.join(['tree_not_found',receiver,'NA','NA','NA','Phylogeny','NA','NA','NA','NA','NA'])+'\n')
	out.close()
	print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See summary file in '+sp+'.hgt.sum.tsv')

########################
#Standalone mt_eval
elif args.m =='mt_eval':
	#Evaluate the source of the region and update the summary file
	current_time = datetime.datetime.now()
	print(f"{current_time}\tReading and processing phylogenies from {len(sum_txt)-1} loci in {args.wd}/{sp}.mtpt.sum.tsv")
	lines=open(f"{args.wd}/{sp}.alnmap.bed").readlines()
	print(f"{current_time}\tReading and processing phylogenies from {len(lines)} loci in {args.wd}/{sp}.alnmap.bed")
	out=open(f"{sp}.hgt.sum.tsv","w")
	d=out.write('Query\tStart\tEnd\tAlignment\tClassification\tRecepient\tDonor_Family\tDonor_genera\tMethod\tBS\n')
	for l in lines:
		i = l.split()[-1]
		i = i.split('.')[2]
		recs=open(f"{args.wd}/{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.fas").readlines()
		allsp=[j[1:].strip() for j in recs if j.startswith('>')]
		target_tip=[j for j in allsp if j.startswith('query')]
		target_tip=target_tip[0]
		allsp.remove(target_tip)
		families=[j.split('|')[0] for j in allsp]
		families=list(set(families))
		try:families.remove(fam)
		except ValueError:pass
		if not families==['NA']:
			try:families.remove('NA')
			except ValueError:pass
		ingroup_sp = [j for j in allsp if j.split('|')[0] in ingroup]
		genera=[j.split('|')[1] for j in allsp if not j.startswith(fam)]
		genera=list(set([j.split('_')[0] for j in genera]))
		same_fam_sp_id=[j for j in allsp if j.startswith(fam)]
		#same_fam_sp = list(set([j.split('|')[1] for j in same_fam_sp_id]))
		#BLAST case 1: only one or two families are in the blast result
		if len(families)==0:
			#VGT: only the query family
			d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','BLAST: homology only found in ingroup','NA'])+'\n')
		elif len(families)==1:
			if families[0] in ingroup:
				#VGT: the only other family is an ingroup
				d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','BLAST: homology only found in ingroup','NA'])+'\n')
			else:
				#HGT: the only other family is not ingroup
				receiver=';'.join([target_tip]+same_fam_sp_id)
				donor_fam=families[0]
				donor_gen=';'.join(genera)
				d=out.write(f"{l.strip()}\tHigh confidence HGT\t{receiver}\t{donor_fam}\t{donor_gen}\tBLAST: exclusive homology in two families\tNA\n")
		#>=2 other families
		else:
			#BLAST case 2: no other ingroup is in the tree, but there are at least two other families
			if len(ingroup_sp)==0:
				#check tree for donor
				try:
					t=Tree(f"{args.wd}/{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas.treefile",format=1)
					t=correct_nd_support(t)
					#reroot
					ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
					t=ncbi_ref_root(t,ncbi_tree)
					receiver,donor_fam,donor_gen,donor_tips,bs=GetSister(t,target_tip,fam)
					if not receiver=='NA':
						#check if the query is well nested within the donor family
						if len(donor_fam)==1:
							alltips=[node.name for node in t]
							qs_index = [j for j, name in enumerate(alltips) if name in receiver]
							donor_fam_binary = [1 if j.split('|')[0] in donor_fam else 0 for j in alltips]
							donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
							nested_in_donor_fam = 0
							if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
							if nested_in_donor_fam:
								d=out.write(l.strip()+'\t'+'\t'.join(['High confidence HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'BLAST+Phylogeny: no ingroup shows homology and nested in donor family',bs])+'\n')
							else:
								d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'BLAST: no other ingroup shows homology',bs])+'\n')
						else:
							d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'BLAST: no other ingroup shows homology',bs])+'\n')
					else:
						#the query was placed on root by mid-point rooting
						d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join([target_tip]+same_fam_sp_id),';'.join(families),';'.join(genera),'BLAST: check rooting problem',bs])+'\n')
				except ete3.parser.newick.NewickError:
					receiver=target_tip
					donor_fam=''
					donor_gen=''
					bs=''
					if len(allsp)<4:
						#too few tips to run a tree
						donor_fam=families
						donor_gen=genera
						bs='NA'
					d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',receiver,';'.join(donor_fam),';'.join(donor_gen),'BLAST: check rooting problem',bs])+'\n')
			#Phylogeny to assess classification	
			else:
				try:
					#examine if this can be a VGT without rerooting the tree
					t=Tree(f"{args.wd}/{sp}_HGTscanner_supporting_files/{sp}.hgt.{i}.aln.fas.treefile",format=1)
					t=correct_nd_support(t)
					#VGT based on tip order: the query if nested within an ingroup cluster
					if VGTFromTipOrder(t,target_tip,fam):
						d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','Phylogeny','NA'])+'\n')
					else:
						#root tree
						ncbi_tree=Tree(script_path+'/database/ncbi_common_taxonomy.phy',format=1)
						t=ncbi_ref_root(t,ncbi_tree)
						#use tree topology to make decisions
						receiver,donor_fam,donor_gen,donor_tips,bs=GetSister(t,target_tip,fam)
						if VGTFromTipOrder(t,target_tip,fam):
							d=out.write(l.strip()+'\t'+'\t'.join(['VGT','NA','NA','NA','Phylogeny','NA'])+'\n')
							continue
						if len(donor_fam)==1:
							if donor_fam[0] in ingroup:
								#VGT: query sister to an ingroup
								d=out.write(l.strip()+'\t'+'\t'.join(['VGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'Phylogeny',bs])+'\n')
							else:
								#check if the query is well nested within the donor family
								alltips=[node.name for node in t]
								qs_index = [j for j, name in enumerate(alltips) if name in receiver + donor_tips]
								donor_fam_binary = [1 if j.split('|')[0] in donor_fam else 0 for j in alltips]
								donor_fam_start,donor_fam_end,donor_fam_index=find_max_block_around_element(donor_fam_binary, min(qs_index), max(qs_index), max_diff_count=1)
								nested_in_donor_fam = 0
								if donor_fam_start<min(qs_index) or donor_fam_end>max(qs_index):nested_in_donor_fam = 1
								if nested_in_donor_fam:
									#High confidence HGT: a single donor family is involved and the q+sister nest within other species of the same donor family
									if donor_fam_end==len(alltips):donor_fam_tips = alltips[donor_fam_start:]
									else:donor_fam_tips = alltips[donor_fam_start:donor_fam_end+1]
									fam_support = max_clade_support(target_tip,donor_fam_tips,t)
									d=out.write(l.strip()+'\t'+'\t'.join(['High confidence HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'Phylogeny: nested in donor family',bs])+'\n')
								else:
									#Putative HGT: one family that's not close relative, but phylogeny is not completely well nested within
									d=out.write(l.strip()+'\t'+'\t'.join(['Putative HGT',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'Phylogeny',bs])+'\n')
						else:
							#multiple donor families
							d=out.write(l.strip()+'\t'+'\t'.join(['inconclusive',';'.join(receiver),';'.join(donor_fam),';'.join(donor_gen),'Phylogeny',bs])+'\n')
				except ete3.parser.newick.NewickError:
					receiver=target_tip
					donor_fam=''
					donor_gen=''
					bs=''
					d=out.write(l.strip()+'\t'+'\t'.join(['tree_not_found',receiver,';'.join(donor_fam),';'.join(donor_gen),'Phylogeny',bs])+'\n')
	out.close()
	print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See summary file in '+sp+'.hgt.sum.tsv')
