#2015.04.15
#S. Hilz
#Grimson Lab

#Description: Takes a CPM file plus an infinite number of DE files and creates both a master de file and individual de files for each condition to use in targetting analysis

#Libraries
import numpy as np

#Variables (User-defined)
cpm_file = '20160201_SmallRNA_CPM_Values.csv'
de_file_array = ['20160216_D1vD3_edgeR_results_250.csv','20160216_D3vD7_edgeR_results_250.csv','20160216_D7vLZ_edgeR_results_250.csv','20160216_LZvPach_edgeR_results_250.csv','20160216_D1vD7_edgeR_results_250.csv','20160216_D1vPach_edgeR_results_250.csv', '20160216_D7vPach_edgeR_results_250.csv', '20160217_D1vLZ_edgeR_results_250.csv']
de_comparison_array = ['D1vD3','D3vD7','D7vLZ','LZvPach','D1vD7','D1vPach','D7vPach','D1vLZ']#[',,'D1vPach',,'D7vD8','D7vPach','D8vLZ',]#must correspond with de_file_array!
cpm_dic={}#this one is complicated - its basically a dic that links each comparison with CPM values, one set for each of the two conditions being compared
cpm_dic['D1vD3']=[[1,2],[3,4]]
cpm_dic['D3vD7']=[[3,4],[5,6]]
cpm_dic['D7vLZ']=[[5,6],[7,8]]
cpm_dic['LZvPach']=[[7,8],[9,10]]
cpm_dic['D1vD7']=[[1,2],[5,6]]
cpm_dic['D1vPach']=[[1,2],[9,10]]
cpm_dic['D7vPach']=[[5,6],[9,10]]
cpm_dic['D1vLZ']=[[1,2],[7,8]]


#Part 1: Read in each de_file
comparison_info_dic = {}
i = 0
while i < len(de_file_array):
    comparison_info_dic[de_comparison_array[i]] = {}
    infile = open(de_file_array[i],'r')
    remove_header = infile.readline()
    while 1:
        line = infile.readline()
        if not line: break
        line = line.rstrip()
        sline = line.split(',')
        if float(sline[4])>=.05:
            comparison_info_dic[de_comparison_array[i]][sline[0].replace('"','')] = [str(sline[1]),'NotDE']
        else:
            comparison_info_dic[de_comparison_array[i]][sline[0].replace('"','')] = [str(sline[1]),'DE']
    infile.close()
    i+=1

#Part 2: Read in cpm file
infile = open(cpm_file,'r')
remove_header = infile.readline()
cpm_info_dic = {}
while 1:
    line = infile.readline()
    if not line: break
    line = line.rstrip()
    sline = line.split(',')
    cpm_info_dic[sline[0].replace('"','')] = [float(x) for x in sline[1:]]
infile.close()

#Part 3: Output data
master_dic = {}#will hold each seed and whether or not it was DE, with comparisons indexed in order of the de_comparision_array
i = 0
while i <len(de_comparison_array):
    outfile = open(de_comparison_array[i]+'_is_seed_de.txt','w')
    for seed in comparison_info_dic[de_comparison_array[i]]:
        #print seed
        values = []
        for index in cpm_dic[de_comparison_array[i]][0]:
            #print cpm_info_dic[seed][index-1]
            values.append(cpm_info_dic[seed][index-1])
        avg_cpm_condition1 = np.mean(values)
        values = []
        for index in cpm_dic[de_comparison_array[i]][1]:
            values.append(cpm_info_dic[seed][index-1])
        avg_cpm_condition2 = np.mean(values)
        if avg_cpm_condition1>avg_cpm_condition2:
            highest_cpm = avg_cpm_condition1
        else:
            highest_cpm = avg_cpm_condition2
        outfile.write(seed+'\t'+str(highest_cpm)+'\t'+comparison_info_dic[de_comparison_array[i]][seed][0]+'\t'+comparison_info_dic[de_comparison_array[i]][seed][1]+'\n')
        if seed in master_dic:
            master_dic[seed].append(comparison_info_dic[de_comparison_array[i]][seed][1])
        else:
            master_dic[seed]=[comparison_info_dic[de_comparison_array[i]][seed][1]]
    outfile.close()
    i+=1

master_outfile = open('is_seed_de.txt','w')
header = []
i=0
while i < len(de_comparison_array):
    header.append(de_comparison_array[i])
    i+=1
    
master_outfile.write('\t'+'\t'.join(header)+'\n')

for seed in master_dic:
    master_outfile.write(seed+'\t'+'\t'.join(master_dic[seed])+'\n')

master_outfile.close()    
    
        
    
    
