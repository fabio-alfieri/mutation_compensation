'''
@author: fabio alfieri
'''

import pandas as pd
from Bio import SeqIO
import re
import os

# download the whole proteome fasta file from https://www.uniprot.org/uniprotkb?facets=reviewed:true&query=(taxonomy_id:9606)
proteins = SeqIO.parse(os.getcwd() + "humanProteome.fasta", "fasta")

# produce the mutation file with mutations that provoke amino acid changes
mutations = pd.read_csv(os.getcwd() + "PANCANCER_mutations.csv")
# the mutation file should have at least this two columns:
#     - SWISSPROT
#     - HGVSp_Short

mutatedProteins = 0 
nonMutatedProteins = 0
n = 0

os.system("mkdir " + os.getcwd() + "/sequences")
os.system("mkdir " + os.getcwd() + "/tangoResults")


# print one sequence each loop and mutate it
for seq_record in proteins:
    seq_name = seq_record.id
    prot_id = seq_name.split('|')[2]
    mutations_filt = mutations[(mutations.SWISSPROT == prot_id)] # select protein-specific mutations (based on UniProtID)
    mut = mutations_filt.HGVSp_Short # mutation position and aminoacid change
    
    if len(mutations_filt) != 0:
        file = open(os.getcwd() + "/sequences/protein_" + prot_id + ".txt", "a")
        file.write(str(prot_id) + "_WT N N 7 298 0.1 ")
        file.write(str(seq_record.seq))
        file.write("\n")
        
        for mut_record in mut:
            n = n + 1
            
            if mut_record != '\\N': 
                mut_r = mut_record.split('.')[1]
                mutPosition = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", mut_r)
                mutWT = mut_r[0]
                mutCancer = mut_r[-1]
                s = seq_record.seq
                
                if int(tuple(mutPosition)[0]) <= len(seq_record.seq):
                    
                    if mutWT == seq_record.seq[int(tuple(mutPosition)[0])-1]:
                        sm = s.tomutable()
                        sm[int(tuple(mutPosition)[0])-1] = mutCancer # mutate the WT sequence
                        mutatedProteins = mutatedProteins + 1
                        file.write(str(prot_id) + "_" + mut_r + " N N 7 298 0.1 ")
                        for base in sm:
                            file.write(base)
                        file.write("\n")
                        
                    else:
                        nonMutatedProteins = nonMutatedProteins +1 
                        
                else:
                    nonMutatedProteins = nonMutatedProteins + 1
                    
        file.close()
        os.chdir(os.getcwd() "/tangoResults/")
        os.system("cd " + os.getcwd() + "/tangoResults/")
        cmd = "TANGO/tango --inputfile=" + "sequences/protein_" + prot_id + ".txt"
        os.system(cmd)
        chmod = "chmod 777 *_aggregation*"
        os.system(chmod)
        cmd_mv = "mv *_aggregation* tangoResults/mutProt_aggrScore_" + prot_id + ".txt"
        os.system(cmd_mv)
        
        
        