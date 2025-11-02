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
A Python tool for detecting horizontal gene transfer\n')

try:
    # Attempt to import all required (non-standard) modules
    from Bio import SeqIO
    from ete3 import Tree
    import ete3
    import pybedtools
    from numpy import median
    from statistics import mode
    import pandas as pd
    
    # Import standard modules after the checks are passed
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


parser = argparse.ArgumentParser(description='HGTscanner is a tool to identify HGT blocks in organellar genomes.')
parser.add_argument('-q', metavar='query', help='Fasta file of the query genome', required=True)
parser.add_argument('-o', metavar='output_prefix', help='Output prefix', required=True)
parser.add_argument('-f', metavar='family', help='Family of the query for HGT classification', required=True)
parser.add_argument('-mtpt', action='store_true', help='Invoking the MTPT mode')
parser.add_argument('-pt_fix_id', metavar='id_file', help='A file of user-selected GenBank accession numbers for MTPT detection.')
parser.add_argument('-pt_add_seq', metavar='fatsa', help='A fasta file containing plastid references for MTPT detection.')
parser.add_argument('-pt_add_id', metavar='id_file', help='A file user-selected GenBank accession numbers for MTPT detection.')
parser.add_argument('-hit', metavar='number of hits', help='Number of best blast hits to be included.')
parser.add_argument('-mt_add_seq', metavar='fasta', help='A fasta file containing mitochondrial references for mt HGT detection.')
parser.add_argument('-b', metavar='bed_file', help='A bed file for regions to be masked')
parser.add_argument('-e', metavar='evalue', help='BLAST evalue threshold')

####################################
#I. pass argument values, check required arguments
################################

args = parser.parse_args()
script_path = os.path.abspath(sys.argv[0])
script_path = os.path.dirname(script_path)

#set family for the focal clade
if args.f:
	fam = args.f

def id2seq(ids,output_file):
	recs=SeqIO.parse(script_path+'/database/Viridiplantae_pt_aug2025.genome.fas','fasta')
	out=open(output_file,'a')
	for rec in recs:
		id=rec.id
		if id.split('.')[0] in ids:d=SeqIO.write(rec,out,'fasta')
		out.close()

if args.mtpt:
	#mtpt mode
	#assemble the plastid dataset for MTPT
	try:
		sp=args.o
		query=args.q
		print(str(datetime.datetime.now())+'\tDetecting MTPT in '+sp+' using the query sequence '+query)
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
				S='cat '+pt_reference+' '+script_path+'/database/Viridiplantae_pt_aug2025.representative.fas >'+sp+'.pt_db.fas'
				os.system(S)
				add_seq=1
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
python HGTScanner.py -mtpt -q <query sequence> -o <output prefix> -f <family> [optional] -e <e value> -pt_add_seq <fasta> -pt_add_id <list of id file> -pt_fix_id <list of id file>')
	except IOError as e:print(e.errno)
else:
	#default mode for mtHGT
	try:
		query=args.q
		sp=args.o
		print(str(datetime.datetime.now())+'\tDetecting HGT in '+sp+' using the query sequence '+query)
		if args.mt_add_seq:
			#add custom references to database
			mt_reference = args.mt_add_seq
			print(str(datetime.datetime.now())+'\tAdded custom mitochondrial reference '+mt_reference+' to the NCBI Viridiplantae mitochondrial database')
			S='cat '+mt_reference+' '+script_path+'/database/Viridiplantae_mt_aug2025.genome.fas >'+sp+'.mt_db.fas'
			os.system(S)
		else:
			S='cp '+script_path+'/database/Viridiplantae_mt_aug2025.genome.fas '+sp+'.mt_db.fas'
			os.system(S)
			print(str(datetime.datetime.now())+'\tNo custom mitochondrial reference provided. Will use the built-in Viridiplantae mitochondrial database')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -q <query sequence> -o <output prefix> -f <family> [optional] -mt_add_seq <reference fasta> -e <e value> -b <bed file for masking>')
	except IOError as e:print(e.errno)

	
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
			outputtable.append(consolidatedtable[k][0][0]+'\t'+';'.join(refpos)+'\t'+consolidatedtable[k][0][3]+'\t'+';'.join(targetpos)+'\t'+'\t'.join(consolidatedtable[k][0][6:]))
	return(outputtable)

def write_seq_from_bed_txt(bed_txt,output_handle):
	for l in bed_txt:
		l=str(l).split()
		#write header
		d=output_handle.write(f">{l[7]}|{l[6]}|{l[2]}_{l[3]}\n")
		tg_pos=l[3].split(';')
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
			
def filter_blast_results(bed_id, top_n=100,pt_screen=False):
	GROUP_COL = 3
	EVAL_COL = 7
	BITSCORE = 6
	data_lines = [seq_loc[i] for i in bed_id]
	data_lines = [l.split() for l in data_lines]
	df = pd.DataFrame(data_lines)
	df.rename(columns={GROUP_COL: 'GroupID', EVAL_COL: 'evalue', BITSCORE: 'bitscore'}, inplace=True)
	if pt_screen:
		df_max_per_group = (
			df.sort_values(by=['GroupID', 'evalue'], ascending=[True, True])
			.drop_duplicates(subset=['GroupID'], keep='first')
		)
	else:df_max_per_group=df
	df_max_per_group["evalue"] = pd.to_numeric(df_max_per_group["evalue"], errors="coerce")
	df_max_per_group["bitscore"] = pd.to_numeric(df_max_per_group["bitscore"], errors="coerce")
	df_final = df_max_per_group.sort_values(
		by=['evalue', 'bitscore'],
		ascending=[True, False]
	).head(top_n)
	filtered_bed = df_final.iloc[:, -1].tolist()
	return [i.strip() for i in filtered_bed]

def find_intersect_bed(bedA,bedB):
	a_hits = bedA.intersect(bedB, wa=True, f=0.7)
	b_hits = bedA.intersect(bedB, wa=True, F=0.7)
	combined = a_hits.cat(b_hits, postmerge=False).sort()
	unique_lines = list(dict.fromkeys(map(str, combined)))
	intersect_unique_bed = pybedtools.BedTool('\n'.join(unique_lines), from_string=True)
	return(intersect_unique_bed)
	
####################################
#III. MTPT mode
###################################
if args.mtpt:
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
	#sort blast results and give each row an uniq id
	S="awk -v OFS='\\t' '{if ($9 <= $10) print $2, $9, $10, $1, $7, $8, $12, $11, $13, $14; else print $2, $10, $9, $1, $8, $7, $12, $11, $13, $14}' "+sp+".mtpt.blast.taxon| sort -k1,1 -k2,2n -k4,4n -k8,8g | awk -v FS='\\t' -v OFS='\\t' '{print $0, NR}' > "+sp+".mtpt.bed"
	os.system(S)
	os.system(f"rm {sp}.mtpt.blast.taxon")
	#Merge MTPT blocks
	hits=open(sp+'.mtpt.bed').readlines()
	loci=pybedtools.BedTool(''.join(hits), from_string=True).merge(c=11,o='collapse')
	out=open(sp+'.mtpt.merged.bed','w')
	d=out.write('\n'.join([str(i) for i in loci]))
	out.close()
	q_recs=SeqIO.index(query,'fasta')
	ref_recs=SeqIO.index(sp+'.pt_db.fas', 'fasta')
	seq_loc={}
	for l in hits:seq_loc[l.split()[-1]]=l
	order=1
	retained_order=[]
	output_dir = sp+'_HGTscanner_supporting_files'
	if not os.path.isdir(output_dir):os.mkdir(output_dir)
	for l in loci:
		ids=l.fields[3]
		ids=ids.split(',')
		#If too few plastid hits for this loci (unlikely because plastid synteny is well conserved across order), skip. In practice, I found a lot of these were ancestral mt transfers.
		if len(ids)<3:
			order=order+1
			continue
		out=open(output_dir+'/'+sp+'.mtpt.'+str(order)+'.fas','w')
		#write query
		q_rec_out = q_recs[l.chrom][(l.start-1):l.end]
		q_rec_out.id = "query|" + q_rec_out.id
		q_rec_out.description = ""
		d = SeqIO.write(q_rec_out, out, 'fasta')
		#d=SeqIO.write(q_recs[l.split()[0]][(int(l.split()[1])-1):int(l.split()[2])],out,'fasta')
		#write target
		#if too many blast hits, only retain the top 200 hits to save computational time. I also notice the IR region was included twice for a single species. Thus only one best hit will be selected from each species.
		hit_num = 200
		if args.hit:hit_num=args.hit
		if len(ids)>hit_num:
			ids = filter_blast_results(ids,hit_num, pt_screen=True)
		for i in ids:
			line=seq_loc[i]
			id=line.split()[3]
			start=min(int(line.split()[4]),int(line.split()[5]))
			end=max(int(line.split()[4]),int(line.split()[5]))
			try:d=out.write('>'+id2fam[id.split('.')[0]]+'|'+id2sp[id.split('.')[0]]+'|'+id+'_'+str(start)+'_'+str(end)+'\n')
			except KeyError:d=out.write('>'+id+'_'+str(start)+'_'+str(end)+'\n')
			d=out.write(str(ref_recs[id].seq[start-1:end])+'\n')
		out.close()
		retained_order.append(order)
		order=order+1
	print(str(datetime.datetime.now())+'\tExatracted '+ str(order-1)+' potential MTPT sequences for '+sp)
	#alignment and phylogenetic reconstruction
	print(str(datetime.datetime.now())+'\tStart alignment for '+str(len(retained_order))+' loci with more than 3 taxa. May take a while...')
	for i in retained_order:
		print(f"{str(datetime.datetime.now())}\tAnalyzing Loci #{i}", end='\r')
		#S="nohup mafft --quiet --adjustdirection "+output_dir+"/"+sp+".temp."+str(i)+".fas | sed 's/_R_//g' > "+output_dir+"/"+sp+".mtpt."+str(i)+".aln.fas > /dev/null 2>&1 &"
		S="mafft --quiet --adjustdirection "+output_dir+"/"+sp+".mtpt."+str(i)+".fas | sed 's/_R_//g' > "+output_dir+"/"+sp+".mtpt."+str(i)+".aln.fas"
		#print(S)
		os.system(S)
		#b=SeqIO.index(output_dir+"/"+sp+".mtpt."+str(i)+".aln.fas",'fasta')
		#q=loci[i-1].split()[0]
		#if len(b[q].seq)<10000:
		#	S="nohup iqtree -B 1000 -T 2 --quiet -redo -s "+output_dir+"/"+sp+".mtpt."+str(i)+".aln.fas >/dev/null 2>&1"
		#	os.system(S)
		#else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')
	print(str(datetime.datetime.now())+'\tAlignment complete')
	print(str(datetime.datetime.now())+'\tCleaning intermediate files')
	os.system(f"rm {sp}.mtpt.n*")
	os.system('rm '+output_dir+"/"+sp+'*.bionj')
	os.system('rm '+output_dir+"/"+sp+'*.gz')
	os.system('rm '+output_dir+"/"+sp+'*.log')
	os.system('rm '+output_dir+"/"+sp+'*.iqtree')
	os.system('rm '+output_dir+"/"+sp+'*.mldist')
	os.system('rm '+output_dir+"/"+sp+'*.phy')
	os.system('rm '+output_dir+"/"+sp+'*.contree')
	os.system('rm '+output_dir+"/"+sp+'*.nex')
	##############
	#Evaluate the source of the region and output summary file
	out=open(sp+'.mtpt.sum.tsv','w')
	out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
	for i in retained_order:
		q=loci[i-1].split()[0]
		outgroup=[]
		try:
			t=Tree(output_dir+"/"+sp+'.mtpt.'+str(i)+'.aln.fas.treefile')
			#check branch length of query first, if too long, certainly an ancestral mt transfer and not a young mtpt
			ancestor=t.get_midpoint_outgroup()
			t.set_outgroup(ancestor)
			q_branch=t&q
			if q_branch.dist>1:
				out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\tNA\tNA\n')
			if q_branch.get_ancestors()[0].is_root():
				pass#handle root
			sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
			out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
		except ete3.parser.newick.NewickError:
			sisters=open(output_dir+"/"+sp+'.mtpt.'+str(i)+".fas").readlines()
			sisters=[i[1:].strip() for i in sisters if (i.startswith('>')) and (not i[1:].strip()==q)]
			out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+'ALL HITS: '+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
	out.close()
	#organizing files
	os.system('rm '+sp+'.pt_db.fas')
	print(str(datetime.datetime.now())+'\tCompleted evaluation of MTPT source. See summary file in '+sp+'.hgt.sum.tsv')

#####################################
#IV. Default mode for HGT detection
else:
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
	#BLAST
	evalue=1e-20
	if args.e:evalue=args.e
	print(str(datetime.datetime.now())+'\tPerforming BLAST to identify candidate HGT with evalue threshold of '+str(evalue))
	S='blastn -task dc-megablast -query '+sp+'.mt_db.fas -db '+sp+'.mt -outfmt 6 -evalue '+str(evalue)+' >'+sp+'.mt.blast'
	os.system(S)
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
		if not l.split()[9] in ['Orobanchaceae','Lamiaceae','Scrophulariaceae','Oleaceae','Phrymaceae','Acanthaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Bignoniaceae','Anacardiaceae']:
			otherfam.append(l)
		else:
			samefam.append(l)
	#merge BLAST hits, but require at least 20 bp overlap
	otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=11,o='collapse',d=-20)
	samefam_bed=pybedtools.BedTool(''.join(samefam), from_string=True)
	out=open(sp+'.merged.bed','w')
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
	for hit in otherfam_merged:
		#gather overlapping alignment from both other families and close relatives
		ids=hit.fields[3]
		ids=ids.split(',')
		#If too many hits (>300), select the best ones
		hit_num = 300
		if args.hit:hit_num=args.hit
		if len(ids)>hit_num:
			ids = filter_blast_results(ids,hit_num,pt_screen=False)
		raw_beds_txt=[seq_loc[j] for j in ids]
		seqout_beds=aln_scaffolder(raw_beds_txt,split_block=True)
		if seqout_beds[0]=='X':
			#this hit needs to be further divided into smaller chunks
			raw_beds=pybedtools.BedTool(''.join(raw_beds_txt), from_string=True)
			#print(''.join(raw_beds_txt)) #debug purposes
			for i in seqout_beds[1]:
				out=open(sp+'.hgt.'+str(order)+'.fas','w')
				subhit=pybedtools.BedTool(i, from_string=True)
				#filter for region that either covers >70% of the query seq or have >70% of its own length within the target loci
				new_raw_beds = find_intersect_bed(raw_beds,subhit)
				#rescaffold
				if new_raw_beds.count()>1:
					new_raw_beds=aln_scaffolder(str(new_raw_beds).split('\n')[:-2],split_block=False)
				else:
					j=str(new_raw_beds).split()
					new_raw_beds=[f"{j[0]}\t{j[1]}-{j[2]}\t{j[3]}\t{j[4]}-{j[5]}\t{j[6]}\t{j[7]}\t{j[8]}\t{j[9]}\t{j[10]}"]
				#write other family
				min_start=1000000
				max_end=0
				for l in new_raw_beds:
					l=str(l).split()
					query_pos=l[1].split('-')
					if int(query_pos[0])<min_start:min_start=int(query_pos[0])
					if int(query_pos[-1])>max_end:max_end=int(query_pos[-1])
					#write sequence
					d=out.write(f">{l[7]}|{l[6]}|{l[2]}_{l[3]}\n")
					tg_pos=l[3].split(';')
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
				q_rec_out = q_recs[l[0]][(min_start-1):max_end]
				q_rec_out.id = "query|" + q_rec_out.id
				q_rec_out.description = ""
				d = SeqIO.write(q_rec_out, out, 'fasta')
				#write close relative
				samefam_hit = find_intersect_bed(samefam_bed,subhit)
				if not samefam_hit.count()==0:
					d=SeqIO.write(q_recs[l[0]][(min_start-1):max_end],sp+'.tempseed.fas','fasta')
					out2=open(sp+'.tempTarget.fas','w')
					for ll in samefam_hit:
						d=SeqIO.write(ref_recs[ll.fields[3]],out2,'fasta')
					out2.close()
					samefam_bed_txt=seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas')
					if not samefam_bed_txt=='':write_seq_from_bed_txt(samefam_bed_txt,out)
				d=mapout.write(hit.chrom+'\t'+str(min_start)+'\t'+str(max_end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
				order=order+1
		else:
			out=open(sp+'.hgt.'+str(order)+'.fas','w')
			#write other family
			write_seq_from_bed_txt(seqout_beds,out)
			#write query
			q_rec_out = q_recs[hit.chrom][(hit.start-1):hit.end]
			q_rec_out.id = "query|" + q_rec_out.id
			q_rec_out.description = ""
			d = SeqIO.write(q_rec_out, out, 'fasta')
			#d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],out,'fasta')
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
	print(f"{current_time}\tA total of #{order} aligned sequences from #{num-1} merged homologous genetic blocks were extracted.", end='\r')
	#organize files
	os.system('rm '+sp+'.temp*')
	os.system('rm '+sp+'.mt.n*')
	os.system('rm '+sp+'.mt_db.fas')
	if not os.path.isdir(sp+'_HGTscanner_supporting_files'):os.mkdir(sp+'_HGTscanner_supporting_files')
	os.system('mv '+sp+'.hgt.*.fas '+sp+'_HGTscanner_supporting_files')
	print(f"{current_time}\tStart alignment for #{order} loci")
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
		os.system(f"mafft --quiet --adjustdirection --6merpair --addfragments {temp2} {temp1} | sed 's/_R_//g' > {aln_out}")
		os.remove(temp1)
		os.remove(temp2)
	print(str(datetime.datetime.now())+'\tCompleted alignment. See sequence file in '+sp+'_HGTscanner_supporting_files')
	#############
	#alignment and phylogenetic reconstruction
	print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')

