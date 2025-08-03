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

from Bio import SeqIO
import os,argparse
import datetime
from ete3 import Tree
import ete3
import pybedtools
from numpy import median
from statistics import mode
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description='HGTscanner is a tool to identify HGT blocks in organellar genomes.')
parser.add_argument('-q', metavar='query', help='fasta file of the target genome', required=True)
parser.add_argument('-mtpt', metavar='mode', help='invoke the MTPT mode')
parser.add_argument('-pt_ref', metavar='reference', help='fasta file containing custom plastid references of both close relatives and potential HGT donor. This will be combined with the NCBI Viridiplantae plastid database.')
parser.add_argument('-mt_ref', metavar='reference', help='fasta file containing custom mitochondrial references of both close relatives and potential HGT donor. This will be combined with the NCBI Viridiplantae mito database.')
parser.add_argument('-o', metavar='output', help='output prefix', required=True)
parser.add_argument('-f', metavar='family', help='family of the query for HGT classification', required=True)
parser.add_argument('-b', metavar='file', help='bed file for regions to be masked')
parser.add_argument('-e', metavar='evalue', help='BLAST evalue threshold')

####################################
#I. pass argument values, check required arguments
################################

args = parser.parse_args()
#optional
if args.f:
	fam = args.f

if args.mtpt:
	#mtpt mode
	try:
		sp=args.o
		query=args.q
		print(str(datetime.datetime.now())+'\tDetecting MTPT in '+sp+' using the query sequence '+query)
		if args.pt_ref:
			#add custom plastid references to database
			pt_reference = args.pt_ref
			print(str(datetime.datetime.now())+'\tAdd custom plastid reference '+pt_reference+' to the representative NCBI Viridiplantae plastid database')
			S='cat '+pt_reference+' Viridiplantae_pt.fasta >'+sp+'.pt_db.fas'
			os.system(S)
		else:
			S='cp Viridiplantae_pt.fasta '+sp+'.pt_db.fas'
			os.system(S)
			print(str(datetime.datetime.now())+'\tNo custom plastid reference provided. Will only identify MTPT locations without make inference for HGT')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -mtpt -q <query sequence> -o <output prefix> [optional] -pt_ref <reference fasta>')
	except IOError as e:print(e.errno)
else:
	#default mode
	try:
		query=args.q
		sp=args.o
		print(str(datetime.datetime.now())+'\tDetecting HGT in '+sp+' using the query sequence '+query)
		if args.mt_ref:
			#add custom references to database
			mt_reference = args.mt_ref
			print(str(datetime.datetime.now())+'\tAdded custom mitochondrial reference '+mt_reference+' to the NCBI Viridiplantae mitochondrial database')
			S='cat '+mt_reference+' Viridiplantae_mt.fasta >'+sp+'.mt_db.fas'
			os.system(S)
		else:
			S='cp Viridiplantae_mt.fasta '+sp+'.mt_db.fas'
			os.system(S)
			print(str(datetime.datetime.now())+'\tNo custom mitochondrial reference provided. Will use the built-in Viridiplantae mitochondrial database')
	except TypeError:
		print('############################################################\n\
#ERROR:Insufficient arguments!\n\
Usage:\n\
python HGTScanner.py -q <query sequence> -o <output prefix> [optional] -mt_ref <reference fasta> -b <bed file for masking> -f <query species family>')
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

def aln_scaffolder(bedtxt):
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
	if mode_consolidated_block_count<2:
		outputtable=[]
		for k in consolidatedtable.keys():
			refpos=[l[1]+'-'+l[2] for l in consolidatedtable[k]]
			targetpos=[l[4]+'-'+l[5] for l in consolidatedtable[k]]
			outputtable.append(consolidatedtable[k][0][0]+'\t'+';'.join(refpos)+'\t'+consolidatedtable[k][0][3]+'\t'+';'.join(targetpos)+'\t'+'\t'.join(consolidatedtable[k][0][6:]))
	else:
		#break this scaffold into smaller chunks
		for k in consolidated_block_count.keys():
			outputtable=['X']
			if consolidated_block_count[k]==mode_consolidated_block_count:
				outputtable.append(['\t'.join(i[:3]) for i in consolidatedtable[k]])
				break
	return(outputtable)

def seq2seq_ortho_extraction(seed_file,targets_file,output_handle):
	S='makeblastdb -in '+targets_file+' -out '+sp+'.temp -dbtype nucl >/dev/null'
	os.system(S)
	S='blastn -task dc-megablast -query '+seed_file+' -db '+sp+'.temp -outfmt 6 -evalue 1e-20 | sort -k2,2 -k11,11n> '+sp+'.temp.blast'
	os.system(S)
	hits=open(sp+'.temp.blast').readlines()
	#select only the best hit per target file
	cur_hit=''
	for l in hits:
		ncbiID=l.split()[1]
		if ncbiID!=cur_hit:
			try:
				d=out.write('>'+family[ncbiID]+'|'+ncbiID+'|'+species[ncbiID]+'\n')
			except KeyError:
				if ncbiID.startswith('Oro'):d=out.write('>Orobanchaceae|'+ncbiID+'\n')
				else:d=out.write('>NA|'+ncbiID+'\n')
			if int(l.split()[8])<int(l.split()[9]):
				d=out.write(str(ref_recs[l.split()[1]].seq[(int(l.split()[8])-1):int(l.split()[9])])+'\n')
			else:
				d=out.write(str(ref_recs[l.split()[1]].seq[(int(l.split()[9])-1):int(l.split()[8])])+'\n')
			cur_hit=ncbiID
		else:
			pass
			

####################################
#III. MTPT mode
###################################
if args.mtpt:
	#blast
	if args.e:evalue=args.e
	else:evalue=1e-40
	S='makeblastdb -in '+query+' -out '+sp+' -dbtype nucl >/dev/null'
	os.system(S)
	S='blastn -task dc-megablast -query '+sp+'.pt_db.fas -db '+sp+' -outfmt 6 -evalue '+evalue+' >'+sp+'.mtpt.blast'
	os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
	#sort blast results and give each row an uniq id
	S="awk -v OFS='\\t' '{if ($9 <= $10) print $2, $9, $10, $1, $7, $8, $3, S11, $12; else print $2, $10, $9, $1, $7, $8, $3, S11, $12}' "+sp+".mtpt.blast| sort -k1,1 -k2,2n -k4,4n | awk 'BEGIN{FS=OFS='\\t'} {print $0, NR}' > "+sp+".mtpt.bed"
	os.system(S)
	#define potential HGT blocks
	#otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=9,o='collapse')

	S="bedtools merge -i "+sp+".mtpt.bed -c 9 -o collapse >"+sp+'.temp.bed'
	os.system(S)
	#extract sequences for each block
	loci=open(sp+'.temp.bed').readlines()
	hits=open(sp+'.mtpt.bed').readlines()
	q_recs=SeqIO.index(query,'fasta')
	ref_recs=SeqIO.index(pt_reference, 'fasta')
	seq_loc={}
	for l in hits:
		seq_loc[l.split()[-1]]=l
	order=1
	for l in loci:
		out=open(sp+'.temp.'+str(order)+'.fas','w')
		d=SeqIO.write(q_recs[l.split()[0]][(int(l.split()[1])-1):int(l.split()[2])],out,'fasta')
		out.close()
		out=open(sp+'.temp.'+str(order)+'.fas','a')
		bed_ids=l.split()[3].strip()
		bed_ids=bed_ids.split(',')
		for i in bed_ids:
			line=seq_loc[i]
			id=line.split()[3]
			start=int(line.split()[4])-1
			end=int(line.split()[5])
			d=out.write('>'+id+'_'+str(start)+'_'+str(end)+'\n')
			d=out.write(str(ref_recs[id].seq[start:end])+'\n')
		out.close()
		order=order+1
	print(str(datetime.datetime.now())+'\tExatracted '+ str(order-1)+' potential MTPT sequences for '+sp)
	#alignment and phylogenetic reconstruction
	print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')
	for i in range(1,order):
		S="mafft --quiet --adjustdirection "+sp+".temp."+str(i)+".fas | sed 's/_R_//g' > "+sp+".gt."+str(i)+".aln.fas"
		os.system(S)
		b=SeqIO.index(sp+".gt."+str(i)+".aln.fas",'fasta')
		q=loci[i-1].split()[0]
		if len(b[q].seq)<10000:
			S="nohup iqtree -B 1000 -T 4 --quiet -redo -s "+sp+".gt."+str(i)+".aln.fas >/dev/null 2>&1"
			os.system(S)
			print(str(datetime.datetime.now())+'\tLoci #'+str(i))
		else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')
	print(str(datetime.datetime.now())+'\tCleaning intermediate files')
	os.system('rm '+sp+'*.bionj')
	os.system('rm '+sp+'*.gz')
	os.system('rm '+sp+'*.log')
	os.system('rm '+sp+'*.iqtree')
	os.system('rm '+sp+'*.mldist')
	os.system('rm '+sp+'*.phy')
	os.system('rm '+sp+'*.contree')
	os.system('rm '+sp+'*.nex')
	##############
	#Evaluate the source of the region and output summary file
	out=open(sp+'.mtpt.sum.tsv','w')
	out.write('ID\tTarget_scaffold\tStart\tEnd\tPhylo_source\tBlast_hit_ID\n')
	for i in range(1,order):
		q=loci[i-1].split()[0]
		outgroup=[]
		try:
			t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
			ancestor=t.get_midpoint_outgroup()
			t.set_outgroup(ancestor)
			q_branch=t&q
			if not q_branch.get_ancestors()[0].is_root():
				sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
			else:
				t=Tree(sp+'.gt.'+str(i)+'.aln.fas.treefile')
				outgroup=[leaf.name for leaf in t if leaf.name.startswith('Sorghum')]
				if len(outgroup)>0:
					t.set_outgroup(outgroup[0])
					q_branch=t&q
					sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
				else:
					outgroup=[leaf.name for leaf in t if leaf.name.startswith('Rehm')]
					if outgroup:
						t.set_outgroup(t&outgroup[0])
						q_branch=t&q
						sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
					else:
						#no sorghum no rehmannia
						sisters=[leaf.name for leaf in q_branch.get_sisters()[0]]
			out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
		except ete3.parser.newick.NewickError:
			sisters=open(sp+'.temp.'+str(i)+".fas").readlines()
			sisters=[i[1:].strip() for i in sisters if (i.startswith('>')) and (not i[1:].strip()==q)]
			out.write(str(i)+'\t'+loci[i-1].split()[0]+'\t'+loci[i-1].split()[1]+'\t'+loci[i-1].split()[2]+'\t'+'ALL HITS: '+','.join(sisters)+'\t'+loci[i-1].split()[3]+'\n')
	out.close()
	#organizing files
	os.system('rm '+sp+'.temp.*.fas')
	os.system('rm '+sp+'.temp.bed')
	os.system('rm '+sp+'.n*')
	if not os.path.isdir('HGTscanner_supporting_files'):os.mkdir('HGTscanner_supporting_files')
	os.system('mv '+sp+'.*.aln.fas* HGTscanner_supporting_files')
	print(str(datetime.datetime.now())+'\tCompleted evaluation of MTPT source. See summary file in '+sp+'.hgt.sum.tsv')

#####################################
#IV. Default mode for HGT detection

else:
	if args.b:
		mask_bed=args.b
		print(str(datetime.datetime.now())+'\tMasking query sequence '+query+' using the bed file: '+args.b)
		mask_fasta_with_bed(query, mask_bed, sp+'.mask.fas')
	#BLAST
	if args.ref:
		#add custom references to database
		print(str(datetime.datetime.now())+'\tAdded custom reference '+reference+' to the NCBI Viridiplantae mitochondrial database')
		S='cat '+reference+' Viridiplantae_mt.fasta >'+sp+'.mt_db.fas'
		os.system(S)
	else:
		S='cp Viridiplantae_mt.fasta '+sp+'.mt_db.fas'
		os.system(S)
	S='makeblastdb -in '+sp+'.mt_db.fas -out '+sp+'.mt -dbtype nucl >/dev/null'
	os.system(S)
	if args.b:S='blastn -task dc-megablast -query '+sp+'.mask.fas -db '+sp+'.mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
	else:S='blastn -task dc-megablast -query '+query+' -db '+sp+'.mt -outfmt 6 -evalue 1e-20 >'+sp+'.raw.blast'
	os.system(S)
	print(str(datetime.datetime.now())+'\tBLAST completed for '+sp)
	#Add taxonomic information to each hit at species and family level to help identify syntenic blocks
	x=open(sp+'.raw.blast').readlines()
	out=open(sp+'.taxon.blast','w')
	species={}
	family={}
	y=open('Viridiplantae_mt.taxonomy').readlines()
	for l in y:
		species[l.split('\t')[0]]=l.split('\t')[1]
		family[l.split('\t')[0]]=l.split('\t')[2].strip()

	for l in x:
		try:
			d=out.write(l.strip()+'\t'+species[l.split()[1]]+'\t'+family[l.split()[1]]+'\n')
		except KeyError:
			if l.split()[1].startswith('Oro'):d=out.write(l.strip()+'\t'+l.split()[1]+'\tOrobanchaceae\n')
			else:d=out.write(l.strip()+'\t'+l.split()[1]+'\tNA\n')
	out.close()
	#sort blast results, give each row a unique id, only hits <20k are printed
	S="cat "+sp+".taxon.blast | sort -k1,1 -k7,7n -k11,11n | awk -v OFS='\\t' '{if ($8-$7 <= 20000) print $1, $7, $8, $2, $9, $10, $13, $14, NR}' > "+sp+".taxon_sorted.bed"
	os.system(S)
	#define potential HGT blocks
	#classify these hits based on source, here, alignments from Orobanchaceae, Lamiaceae, Oleaceae, and Scrophulariaceae will be regarded as close relatives
	#This is necessary because they create long synteny blocks that may be broken by HGT in the middle
	x=open(sp+".taxon_sorted.bed").readlines()
	otherfam=[]
	samefam=[]
	for l in x:
		if not l.split()[7] in ['Orobanchaceae','Lamiaceae','Scrophulariaceae','Oleaceae','Phrymaceae','Acanthaceae','Verbenaceae','Plantaginaceae','Gesneriaceae','Bignoniaceae','Anacardiaceae']:
			otherfam.append(l)
		else:
			samefam.append(l)

	otherfam_merged=pybedtools.BedTool(''.join(otherfam), from_string=True).merge(c=9,o='collapse')
	samefam_bed=pybedtools.BedTool(''.join(samefam), from_string=True)
	out=open(sp+'.merged.bed','w')
	d=out.write(str(otherfam_merged))
	print(str(datetime.datetime.now())+'\tFound '+str(len(otherfam_merged))+' homologous genetic blocks for further examination')
	#extract sequences for each block
	q_recs=SeqIO.index(query,'fasta')
	ref_recs=SeqIO.index(sp+'.mt_db.fas', 'fasta')
	order=1
	num=1
	mapout=open(sp+'.alnmap.bed','w')
	for hit in otherfam_merged:
		#gather overlapping alignment from both other families and close relatives
		ids=hit.fields[3]
		raw_beds_txt=id2bed(ids.split(','),x)
		seqout_beds=aln_scaffolder(raw_beds_txt)
		if seqout_beds[0]=='X':
			#this hit needs to be further divided into smaller chunks
			raw_beds=pybedtools.BedTool(''.join(raw_beds_txt), from_string=True)
			for i in seqout_beds[1]:
				out=open(sp+'.hgt.'+str(order)+'.fas','w')
				subhit=pybedtools.BedTool(i, from_string=True)
				new_raw_beds=raw_beds.intersect(subhit,wa=True,f=0.4)
				#write other family
				min_start=1000000
				max_end=0
				for l in new_raw_beds:
					l=str(l).split()
					d=out.write('>'+l[7]+'|'+l[3]+'|'+l[4]+'-'+l[5]+'|'+l[6]+'\n')
					start=int(l[4])
					end=int(l[5])
					if int(l[1])<min_start:min_start=int(l[1])
					if int(l[2])>max_end:max_end=int(l[2])
					if start<end:
						d=out.write(str(ref_recs[l[3]].seq[(start-1):end])+'\n')
					else:
						d=out.write(str(ref_recs[l[3]].seq[(end-1):start].reverse_complement())+'\n')
				#write query
				d=SeqIO.write(q_recs[l[0]][(min_start-1):max_end],out,'fasta')
				#write close relative
				samefam_hit=samefam_bed.intersect(subhit,f=0.4)
				if str(samefam_hit)!='':
					d=SeqIO.write(q_recs[l[0]][(min_start-1):max_end],sp+'.tempseed.fas','fasta')
					out2=open(sp+'.tempTarget.fas','w')
					for ll in samefam_hit:
						d=SeqIO.write(ref_recs[ll.fields[3]],out2,'fasta')
					out2.close()
					seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas',out)
				d=mapout.write(hit.chrom+'\t'+str(min_start)+'\t'+str(max_end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
				order=order+1
		else:
			out=open(sp+'.hgt.'+str(order)+'.fas','w')
			#write other family
			for l in seqout_beds:
				d=out.write('>'+'|'.join([l.split()[5]]+l.split()[2:5])+'\n')
				secs=l.split()[3]
				sequence=''
				for sec in secs.split(';'):
					start=int(sec.split('-')[0])
					end=int(sec.split('-')[1])
					if start<end:
						sequence=sequence+str(ref_recs[l.split()[2]].seq[(start-1):end])	
					else:
						sequence=sequence+str(ref_recs[l.split()[2]].seq[(end-1):start].reverse_complement())
					d=out.write(sequence+'\n')	
			#write query
			d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],out,'fasta')
			#write close relatives
			hit_bed=pybedtools.BedTool(str(hit),from_string=True)
			samefam_hit=samefam_bed.intersect(hit_bed)
			d=SeqIO.write(q_recs[hit.chrom][(hit.start-1):hit.end],sp+'.tempseed.fas','fasta')
			out2=open(sp+'.tempTarget.fas','w')
			for l in samefam_hit:
				d=SeqIO.write(ref_recs[l.fields[3]],out2,'fasta')
			out2.close()
			seq2seq_ortho_extraction(sp+'.tempseed.fas',sp+'.tempTarget.fas',out)
			d=mapout.write(hit.chrom+'\t'+str(hit.start)+'\t'+str(hit.end)+'\t'+sp+'.hgt.'+str(order)+'.fas\n')
			order=order+1
			out.close()
		current_time = datetime.datetime.now()
		print(f"{current_time}\tExtracting sequences from homologous genetic block #{num}", end='\r')
		num=num+1	
	mapout.close()		
	print(f"{current_time}\tA total of #{order} aligned sequences from #{num} merged homologous genetic blocks were extracted.", end='\r')
	#organize files
	os.system('rm '+sp+'.raw.blast')
	os.system('rm '+sp+'.temp*')
	os.system('rm '+sp+'.mt.n*')
	os.system('rm '+sp+'.mt_db.fas')
	if not os.path.isdir(sp+'_HGTscanner_supporting_files'):os.mkdir(sp+'_HGTscanner_supporting_files')
	os.system('mv '+sp+'.hgt.*.fas '+sp+'_HGTscanner_supporting_files')
	print(str(datetime.datetime.now())+'\tCompleted evaluation of HGT source. See sequence file in '+sp+'_HGTscanner_supporting_files')
	#############
	#alignment and phylogenetic reconstruction
	print(str(datetime.datetime.now())+'\tStart alignment and phylogenetic reconstruction with mafft and iqtree for '+str(order-1)+' regions. May take a while...')

#for i in range(1,order):
#	current_time = datetime.datetime.now()
#	print(f"{current_time}\t Sequence alignment and IQTREE for alignment #{i}", end='\r')
#	S="timeout 20m mafft --genafpair --maxiterate 1000 --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
#	os.system(S)
	#recs=list(SeqIO.parse(sp+".hgt."+str(i)+".fas",'fasta'))
	#max_len=max([len(rec.seq) for rec in recs])
	#if max_len<1500:
	#	S="mafft --localpair --maxiterate 1000 --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
	#	S="mafft --adjustdirection --6merpair --addfragments othersequences referencesequence > output"
	#	os.system(S)
	#	S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".hgt."+str(i)+".aln.fas >/dev/null 2>&1"
	#	os.system(S)
	#elif max_len<5000:
	#	S="mafft --quiet --adjustdirection "+sp+".hgt."+str(i)+".fas | sed 's/_R_//g' > "+sp+".hgt."+str(i)+".aln.fas"
	#	os.system(S)
	#	S="nohup iqtree -B 1000 -T 4 --quiet -m GTR+F -redo -s "+sp+".hgt."+str(i)+".aln.fas >/dev/null 2>&1"
	#	os.system(S)
	#else:print(str(datetime.datetime.now())+'\tLoci #'+str(i)+' is longer than 10kb. Skip tree building. Check manually.')

#os.system('rm '+sp+'*.bionj')
#os.system('rm '+sp+'*.gz')
#os.system('rm '+sp+'*.log')
#os.system('rm '+sp+'*.iqtree')
#os.system('rm '+sp+'*.mldist')
#os.system('rm '+sp+'*.phy')
#os.system('rm '+sp+'*.contree')
#os.system('rm '+sp+'*.nex')

