#!/usr/bin/env python3

import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import sys

#counting the kmers
def kmers_counted(seq, k):
    counted = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer not in counted:
            counted[kmer] = 0
        counted[kmer] += 1
    return counted

#find the false DNA sequence
def sequence_DNA(seq):
    false_DNA = {}
    for i in seq:
        if i not in 'AGCT':
            if i in false_DNA: false_DNA[ i ] += 1
            else: false_DNA[ i ] = 1
    if false_DNA != {}: return false_DNA

#export the table(kmers, observed, and possible) into .csv format for each sequence name
def create_file(file_name, detail_data, seq_name):
    data = pd.DataFrame([[i[0], i[1], i[2]] for i in detail_data], columns=['k-mers', 'observed', 'possible'])
    #print('data-frame',detail_data)
    data.to_csv('output/data/'+ file_name + '_' + seq_name+'.csv', index=False)
    #return data

#show a linguistic complexity graph for all sequence name of DNA
def generate_complex_graph(file_name, seq_name_tab, linguistic_complexity_tab):       
    plt.plot(seq_name_tab, linguistic_complexity_tab, 'ro')
    plt.xlabel('Sequence Name')
    plt.ylabel('Linguistic Complexity')
    #plt.show()
    graph_name = file_name + '_complexity.png'
    plt.savefig('output/image/'+graph_name)
    #return data

#compare the possible and observed kmers of all sequence name
def create_kmers_graph(file_name, seq_name_tab, possible_total_tab, observed_total_tab):
    fig, ax = plt.subplots()
    index = np.arange(len(seq_name_tab))
    bar_width = 0.3
    opacity = 0.55
    rects1 = plt.bar(index, possible_total_tab, bar_width,
                     alpha=opacity,
                     color='y',
                     label='Possible Kmers')
    rects2 = plt.bar(index + bar_width, observed_total_tab, bar_width,
                     alpha=opacity,
                     color='r',
                     label='Observed Kmers')
    plt.xlabel('Sequence Names')
    plt.ylabel('Kmers')
    plt.title('Kmers per Sequence Name')
    plt.xticks(index + bar_width, seq_name_tab)
    plt.legend()
    plt.tight_layout()
    #plt.show()
    graph_name = file_name + '_kmers.png'
    plt.savefig('output/image/'+graph_name)

#show false DNA from each sequence name 
def graph_false_dna(file_name, data, seq_name):
    plt.clf()
    plt.bar(range(len(data)), list(data.values()), align='center')
    plt.xticks(range(len(data)), list(data.keys()))
    graph_name = 'false_DNA_'+ file_name + '_' + seq_name + '.png'
    plt.savefig('output/image/'+graph_name)

#main function
if __name__ == "__main__":
    #get an input variable as a file name
    file_name=sys.argv[1]
    if file_name in glob.glob('*.fasta'):
        f = open(file_name,'r')
        seq = f.readlines()
        file_name = file_name.replace('.', '_')
        seq_name_tab = []
        observed_total_tab = []
        possible_total_tab = []
        linguistic_complexity_tab = []
        for line_num, line in enumerate(seq[0:len(seq)]):
            if len(line) > 1 :
                if '>' in line :
                    line = line.replace(">", "")
                    seq_name = line.rstrip()
                    seq_name_tab.append(seq_name)
                else:
                    seq = line.rstrip()
                    #check each sequence name and format
                    seq_format = sequence_DNA(seq) 
                    if seq_format != None :						
                        wrong_seq_graph = graph_false_dna(file_name, seq_format, seq_name)
                        print('you have wrong DNA sequence in ', seq_name, ' detail:', seq_format)
                    k_list = []
                    possible_tab = []
                    observed_tab = []
                    for k in range(1,len(seq)+1):
                        #check all possible kmers
                        if 4**k < len(seq):
                            possible = 4**k 
                        else:
                            possible = len(seq) - k + 1                         
                        #get the kmers
                        counted = kmers_counted(seq, k)
                        #print(counted)
                        #get the observed kmers  
                        observed = len(counted)
                        k_list.append(k)
                        possible_tab.append(possible)                          
                        observed_tab.append(observed)                        
                    #get the total of all possible kmers                        
                    possible_total = sum(possible_tab)  
                    possible_tab.append(possible_total)
                    possible_total_tab.append(possible_total)
                    #get the total of all observed kmers
                    observed_total = sum(observed_tab)
                    observed_tab.append(observed_total)
                    observed_total_tab.append(observed_total)
                    k_list.append('Total'); 
                    #merge kmers-data per each sequence name
                    detail_data = list(zip(k_list, observed_tab, possible_tab))
                    #print('list',detail_data)
                    #create file from each sequence name                      
                    detail = create_file(file_name, detail_data, seq_name)
                    #get all linguistic complexity
                    linguistic_complexity = observed_total/possible_total
                    linguistic_complexity_tab.append(linguistic_complexity)
        #print(linguistic_complexity_tab) 
        summary_linguistic_complexity = generate_complex_graph(file_name, seq_name_tab, linguistic_complexity_tab)                 
        summary_kmers = create_kmers_graph(file_name, seq_name_tab, possible_total_tab, observed_total_tab)           
        print('succeed to create files and graphs')         
    else:
        print('file shall in fasta format')