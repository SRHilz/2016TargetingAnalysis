#2016.02.01
#S. Hilz
#Grimson Lab

#Description: Takes edgeR RNASeq CPM output file, group info (same as edgeR), and comparison information (same as is_seed_de)and creates individual de files for each condition to use in targetting analysis

#Libraries
import numpy as np
import math

#Variables (User-defined)
cpm_file = '20160201_RNA_CPM_Values.csv'
de_comparison_array = ['D1vLZ','D3vLZ','D3vPach','D7vPach']#must correspond with de_file_array!
cpm_dic={}#this one is complicated - its basically a dic that links each comparison with CPM values, one set for each of the two conditions being compared
cpm_dic['D1vLZ']=[[1,2],[7,8]]
cpm_dic['D3vLZ']=[[3,4],[7,8]]
cpm_dic['D3vPach']=[[3,4],[9,10]]
cpm_dic['D7vPach']=[[5,6],[9,10]]

#Part 1: Read in cpm file
infile = open(cpm_file,'r')
remove_header = infile.readline()
cpm_info_dic = {}
while 1:
    line = infile.readline()
    if not line: break
    line = line.rstrip()
    sline = line.split(',')
    cpm_info_dic[sline[0].replace('"','')] = [float(x) for x in sline[1:]]#indexed by gene (or transcript) name
infile.close()

#Part 2: Output data
i = 0
while i <len(de_comparison_array):#for each comparison
    outfile = open(de_comparison_array[i]+'_CPM_FC.txt','w')#create output
    outfile.write('Gene\tS1\tS2\tS2/S1_log2FC\n')
    for transcript in cpm_info_dic:#now output the average CPM of S1, the average CPM of S2, and the log2FCS2/S1 for each transcript
        values = []#to average
        for index in cpm_dic[de_comparison_array[i]][0]:#grabs each index associated with S1 in the comparison
            #print cpm_info_dic[seed][index-1]
            values.append(cpm_info_dic[transcript][index-1])
        avg_cpm_condition1 = np.mean(values)#S1 average
        values = []#to average
        for index in cpm_dic[de_comparison_array[i]][1]:#grabs each index associated with S2 in the comparison
            values.append(cpm_info_dic[transcript][index-1])
        avg_cpm_condition2 = np.mean(values)
        log2fc = math.log((avg_cpm_condition2+1)/(avg_cpm_condition1+1),2)#calculates log2FC, with psuedocounts added
        outfile.write(transcript+'\t'+str(avg_cpm_condition1)+'\t'+str(avg_cpm_condition2)+'\t'+str(log2fc)+'\n')
    outfile.close()
    i+=1 
    
        
    
    
