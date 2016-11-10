#2015.03.11
#S. Hilz
#Grimson Lab

#Description: Module for taking TargetScan's Summary_Counts.txt file and turning it into an express file indexed by seed containing redundant entries for each seed
#Input: <File_Name> Summary_Counts.txt file
#Output: <To_Return> An express file containing non-redundant entries for each gene symbol, indexed by seed.
#NOTES: This script generates redundant TargetScan Express files. In my targeting analysis,
#        I then first get a list of all expressed transcripts, then determine which of the transcripts in TargetScan are expressed. Finally, if multiple transcripts for the
#        same gene are expressed, I pick the one with the highest context score.
def BuildTargetScanFile_BySeed(SummaryCountsFile,Species_filter,OutputFileName):#Builds a file indexed by seed and then genesymbol, for a particular species
    ToReturn={}#this dic will be indexed first by seed, and then gene symbol, and contain an array of all entries per gene (the best transcript by context score is chosen for those with redundant transcript entries)
    DataIn=open(SummaryCountsFile,'r')
    #input
    #0: Refseq; 1: Symbol; 2: Seed; 3: species filter; 4: #con; 5: #con8mer; 6: #con7merm8; 7: #con7merA1; 8:#noncon; 9:#noncon8mer; 10: #noncon7merm8; 11: #noncon7merA1
    #12: miR family member; #13: cumulative context score; 14:#cumulative PCT
    #output, same format as input, but only outputs lines matching species ID
    #SCORE: more negative implies more expected repression
    #PCT: probability of significant conservation
    tempL=DataIn.readline()#removes header in SummaryCounts.txt
    while 1:
        tempL=DataIn.readline()#goes through every line of summary_counts.txt
        if not tempL: break
        tempL=tempL.rstrip()
        StempL=tempL.split('\t')
        if len(StempL)!=15: print StempL
        Species=StempL[3]
        Seed=StempL[2]
        if str(Species)==str(Species_filter):#appends only if of the right species
            if Seed not in ToReturn:#if we already have an entry in ToReturn for a given Seed
                ToReturn[Seed] = []
            ToReturn[Seed].append(tempL)
    DataIn.close()
    outfile = open(OutputFileName,'w')
    for key in ToReturn:
        outfile.write('>'+key+'\n')
        for line in ToReturn[key]:
            outfile.write(line+'\n')
    print 'Entries for ',len(ToReturn.keys()),' seeds output.'

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
