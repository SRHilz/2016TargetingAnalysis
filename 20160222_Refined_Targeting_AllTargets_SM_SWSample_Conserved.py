#2016.02.01
#S. Hilz
#Grimson Lab

#Description: Takes lists of seeds (foreground and background), and first uses targetscan to predict their targets, then performs a ks-test on the RPKMs of their targets
#Input: <TargetScan_Express_File> File build by BuildTargetScanFile_BySeed_forGeneID.py
#       <NNN_is_seed_de.txt> is_seed_de.txt file for every comparison made
#       <NNN_RPKM.txt> file containing RPKM info for every gene for every comparison made
#Output:<NNN_analysis.txt> file for every comparison made with the following columns: seed, is_de, in_ts, target_num, target_num_RPKM, background_target_num, ks test p-value, highest_cpm, log2fc
#       Column info: seed = the seed, is_de = whether or not the seed was differentially expressed for that condition, in_ts = if it was in the TS database, can_test = if there were at least 5 predicted targets, target_num = how many targets for either the experimental or background set, ks test p-value = p-value from the ks test performed on targets, highest_cpm = average cpm value for condition with the highest exp of that seed, log2fc = log2 fc for that seed 
#Usage: To run, call python <name_of_this_program> <TargetScan_Express_File> <NNN_is_seed_de.txt> <NNN_RPKM.txt> <Comparison> <PCT Cutoff> <Context Cutoff>


#Other: Column Identifiers in TS Express File: Refseq; 1: Symbol; 2: Seed; 3: species filter; 4: #con; 5: #con8mer; 6: #con7merm8; 7: #con7merA1; 8:#noncon; 9:#noncon8mer; 10: #noncon7merm8; 11: #noncon7merA1
    #12: miR family member; #13: cumulative context score; 14:#cumulative PCT
    #output, same format as input, but only outputs lines matching species ID
    #SCORE: more negative implies more expected repression
    #PCT: probability of significant conservation

#NOTE: The big differences between this program and all others are that it assumes a background distribution that is for all miRNA targets, and not just targets of not changing miRNAs. This distribution is made by taking a random sample of 10% of all targets of miRNAs
#USERINPUT (non-command line - will make command line later)
RNASeqCPMCutoff = 1
KSTestListSize = 20
miRNACPMCutoff = 250
BackgroundSetSampleSize = .05

#Updates: The only update since 2015.05.25 was that I changed the gene identifying from the gene symbol (column 2 of TargetScan file) to the gene ID (column 1); I also made it so you can enter the TS PCT and Conservation scores at the command line
import sys
import os
import numpy as np
import scipy
import scipy.stats
import random
from scipy.stats import ks_2samp
from scipy.stats import mannwhitneyu
from collections import defaultdict

def get_targets(Seed, TSDic, ConsCutoff, ContextCutoff):#ConsCutoff=cutoff for conservation score, must be greater than value given; ContextCutoff is for context score, must be less than value given
    TargetList=[]
    N8merList = []
    N7merm8List = []
    N7mer1aList = []
    ContextList = []
    PCTList = []
    if Seed in TSDic:
        for y in TSDic[Seed]:
            ConsPass=0
            ContextPass=0
            Entry=y.split('\t')
            if ConsCutoff!='NA':
                if Entry[14]!='NULL':
                    if float(Entry[14])>float(ConsCutoff):#conservation cutoff
                        ConsPass=1
            
            else: ConsPass=1
            if ContextCutoff!='NA':
                if Entry[13]!='NULL':
                    if float(Entry[13])<float(ContextCutoff):#context score cutoff
                        ContextPass=1
            else: ContextPass=1
            if ConsPass==1 and ContextPass==1:
                TargetList.append(Entry[0])
                N8merList.append(str(int(Entry[5])+ int(Entry[9])))#appends sums of conserved and non-conserved
                N7merm8List.append(str(int(Entry[6])+ int(Entry[10])))
                N7mer1aList.append(str(int(Entry[7])+ int(Entry[11])))
                PCTList.append(str(Entry[14]))
                ContextList.append(str(Entry[13]))
                
    else: print "Seed not in Dic",Seed
    output = [TargetList, N8merList, N7merm8List, N7mer1aList, ContextList, PCTList]#bundles output into list
    return output

def ReadTargetScanFile_BySeed(SummaryCountsExpressFile):#reads in file built by BuildTargetScanFile_BySeed_forGeneID.py into dic
    to_return = {}#this dic will be indexed first by seed, and then genesymbol, and contain an array of all entries per gene (the best transcript by context score is chosen for those with redundant transcript entries)
    infile = open(SummaryCountsExpressFile,'r')
    while 1:
        line=infile.readline()#goes through every line of summary_counts.txt
        if not line: break
        line=line.rstrip()
        if '>' in line:
            seed=line[1:]
            to_return[seed]=[]
        else: to_return[seed].append(line)
    return to_return       

##Part 0: Import List of Conserved Seeds
infile = open('/local/workdir/sh745/Databases/miR_Family_Info_v6.2.txt','r')
conserved = []
while 1:
    line=infile.readline()
    if not line: break
    line=line.rstrip()
    sline = line.split('\t')
    seed = sline[1]
    con_category = str(sline[5])
    if con_category == '2':
        if seed not in conserved:
            conserved.append(seed)
print "Only looking at these conserved seeds: ",len(conserved), conserved
infile.close()

##Part 1: Import DE Seeds
try:
    infile = open(sys.argv[2],'r')
except IndexError:
    print "No input file provided"
    sys.exit()

seed_info = {}

while 1:
    line=infile.readline()
    if not line: break
    line=line.rstrip()
    sline = line.split('\t')
    if float(sline[1])>=miRNACPMCutoff and sline[0] in conserved:#will only consider seeds with an expression over threshold
        seed_info[sline[0]]=sline[1:]

infile.close()

condition = sys.argv[4]   

##Part 2: Import RNASeq
try:
    infile = open(sys.argv[3],'r')
except IndexError:
    print "No input file provided"
    sys.exit()

rna_seq = {}    

infile.readline()#remove header
all_rna = []
while 1:
    line=infile.readline()
    if not line: break
    line=line.rstrip()
    sline=line.split('\t')
    if float(sline[1])>=int(RNASeqCPMCutoff) or float(sline[2])>=int(RNASeqCPMCutoff):#this is the filtering step
        rna_seq[sline[0]]=sline[1:]#indexed by gene name; will thus only contain RNASeq data expressed in at least one sample over the cutoff threshold
        all_rna.append(float(sline[3]))

##Part 3: Import TS Dictionary from Express File and Identify Best Transcript for Each Gene
ts_dic_redundant = ReadTargetScanFile_BySeed(sys.argv[1])

sites_by_transcripts = {}#dic indexed first by gene symbol, then by transcript ID, giving counts for number of sites for each transcript

for seed in ts_dic_redundant:
    for line in ts_dic_redundant[seed]:
        sline = line.split('\t')
        gene_symbol = sline[1]
        transcript_id = sline[0]
        total_sites = int(sline[4]) + int(sline[8])
        if gene_symbol not in sites_by_transcripts:
            sites_by_transcripts[gene_symbol] = {}
        if transcript_id not in sites_by_transcripts[gene_symbol]:
            sites_by_transcripts[gene_symbol][transcript_id] = 0
        sites_by_transcripts[gene_symbol][transcript_id] += total_sites
        
optimal_transcripts = []#list of the best transcript to use per gene, determined by firstly if this transcript is expressed, and secondly if this transcript has the most target sites

for gene in sites_by_transcripts:
    expressed_transcripts = []
    for transcript in sites_by_transcripts[gene]:
        if transcript in rna_seq:
            expressed_transcripts.append(transcript)
    best_transcript = ['empty',0]#name and # sites for "optimal" transcript; starts out empty, we hope to find something better than this
    for transcript in expressed_transcripts:
        if sites_by_transcripts[gene][transcript] > best_transcript[1]:#if we have an actual context score entry
            best_transcript = [transcript,sites_by_transcripts[gene][transcript]]
    if best_transcript[0] != 'empty':
        optimal_transcripts.append(best_transcript[0])#in the end, will have only one transcript per gene

ts_dic = {}#will now only contain non-redundant entries, and only those in optimal transcripts. Means won't contain non-expressed RNASeq

for seed in ts_dic_redundant:
    for line in ts_dic_redundant[seed]:
        sline = line.split('\t')
        transcript_id = sline[0]
        if transcript_id in optimal_transcripts:
            if seed not in ts_dic:
                ts_dic[seed] = []
            ts_dic[seed].append(line)

##Part 4: Build target lists for each seed (EDIT IF NEED TO CHANGE HOW TARGETS ARE FILTERED)

for seed in seed_info:
    if seed in ts_dic:
        seed_info[seed].append('TRUE')
        target_info = get_targets(seed, ts_dic, sys.argv[5],sys.argv[6])#will return a list of 6 lists. Index of lists is important as index links info among lists
        i = 0
        while i<len(target_info):
            seed_info[seed].append(target_info[i])
            i+=1  
    else: 
        seed_info[seed].append('FALSE')
        seed_info[seed].append([])
    
##Part 5: Divide target_dic into experimental and background (EDIT IF NEED TO CHANGE WHAT IS EXPERIMENTAL AND BACKGROUND)
exp_background_targets = []#will for this type of analysis contain all targets of miRNAs for that context score, weighted by number of genes with each type of seed
site_dic = {}
for seed in seed_info:
    candidate_background_set = []
    for target in seed_info[seed][4]:
        if target in rna_seq:
            candidate_background_set.append(target)
            if target not in site_dic:
                site_dic[target] = 0#will have all genes that are in RNAseq data as indexes
            site_dic[target] += 1
    if len(candidate_background_set) >= KSTestListSize:
        percent = int(round(len(candidate_background_set)*BackgroundSetSampleSize))
        sample = random.sample(candidate_background_set,percent)
        for target in sample:
            exp_background_targets.append(target)#this is where the "weighting" happens - we don't make the final background list non-redundant. Instead, we let it contain multiple instances of the same genes - the more time a gene appears, the more sites it has.

outfile = open('sanity_check.txt','w')
d = defaultdict(int)
for target in exp_background_targets:
    d[target] += 1
for target in d:
    outfile.write(target+'\t'+str(d[target])+'\n')
outfile.close()

print "All expressed targets of miRNAs, ",len(list(set(exp_background_targets)))
                
##Part 6: Sort RNASeq data into Experimental and Background targets for each seed, perform ks test, and output
output_context_tag = str(sys.argv[6])
output_context_tag = output_context_tag.replace('-.0','p')#will get it either way, if the user uses a zero or not
output_context_tag = output_context_tag.replace('-.','p')
output_CPM_tag = str(RNASeqCPMCutoff)
main_outfile = open('20160221_'+condition+'_analysis_Context_'+output_context_tag+'_RNACPM'+output_CPM_tag+'.txt','w')  

for seed in seed_info:
    experimental_rpkms, background_rpkms = [],[]
    unique_background_targets = exp_background_targets[:]#creates a copy of background targets; we will subtract from this list all things that are also targets of the seed we are looking at right now
    seed_output = []#this will collect all of the potential targeting output for that particular seed; will only be output if a ks test can be performed in the end
    i = 0#counter to go through each target for that seed
    while i<len(seed_info[seed][4]):#we go through each target in order (important because info for each target is stored in order
        if seed_info[seed][4][i] in rna_seq:#if target i is in our expressed rna_seq dic built in Part 5
            seed_output.append(seed_info[seed][4][i]+'\tExperimental\t'+'\t'.join(rna_seq[seed_info[seed][4][i]])+'\t'+str(site_dic[seed_info[seed][4][i]])+'\t'+seed_info[seed][5][i]+'\t'+seed_info[seed][6][i]+'\t'+seed_info[seed][7][i]+'\t'+seed_info[seed][8][i]+'\t'+seed_info[seed][9][i]+'\n')
            experimental_rpkms.append(float(rna_seq[seed_info[seed][4][i]][2]))
        if seed_info[seed][4][i] in  unique_background_targets:
            unique_background_targets.remove(seed_info[seed][4][i])
        i+=1
    for gene in unique_background_targets:
        seed_output.append(gene+'\tBackground\t'+'\t'.join(rna_seq[gene])+'\t'+str(site_dic[gene])+'\tNA\tNA\tNA\tNA\tNA\n')
        background_rpkms.append(float(rna_seq[gene][2]))
    if len(experimental_rpkms)>=20 and len(background_rpkms)>=20:#this decides finally if we can perform our test
        test_stat = float(mannwhitneyu(experimental_rpkms, background_rpkms)[1])
        if sys.argv[7]=='-v':
            outfile = open('20160215_'+condition+'_'+seed+'_analysis_Context_'+output_context_tag+'_RNACPM'+output_CPM_tag+'.txt','w')#seed specific output; needed to look at the distributions in R
            outfile.write('gene\tlist\ts1_CPM\ts2_CPM\tlog2(s2/s1)fc\tsitecount\t8mers\t7merm8s\t7mer1as\tcontext++\tpct\n')
            for entry in seed_output:
                outfile.write(entry)
            outfile.close()
    else: test_stat = 'NA'
    main_outfile.write(seed+'\t'+seed_info[seed][2]+'\t'+seed_info[seed][3]+'\t'+str(len(seed_info[seed][4]))+'\t'+str(len(experimental_rpkms))+'\t'+str(len(background_rpkms))+'\t'+str(test_stat)+'\t'+str(np.mean(experimental_rpkms))+'\t'+str(np.mean(background_rpkms))+'\t'+str(seed_info[seed][0])+'\t'+str(seed_info[seed][1])+'\n')

main_outfile.close()    




    
    
