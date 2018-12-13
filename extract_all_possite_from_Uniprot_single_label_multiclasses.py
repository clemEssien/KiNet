import numpy as np
import re
import pandas as pd

f1="pos_seq_file.fasta";
f2="pos_seq_file.fasta";
f3="pos_seq_file.fasta";

def read_fasta(fasta_file):
    fp = open(fasta_file)
    lines = fp.readlines()
    
    fasta_dict = {} #record seq for one id
    positive_dict={} #record positive positions for one id
    idlist=[] #record id list sorted
    gene_id = ""
    for line in lines:
        if line[0] == '>':
            if gene_id != "":
                fasta_dict[gene_id] = seq
                idlist.append(gene_id)
            seq = ""
            gene_id = line.strip('\n') #  line.split('|')[1] all in > need to be id
        else:
            seq += line.strip('\n')
    
    fasta_dict[gene_id] = seq #last seq need to be record
    idlist.append(gene_id)
    
    for gene_id in fasta_dict:
       seq=fasta_dict[gene_id];
       posnum=0; #record positive position num
       for pos in range(len(seq)):
           mid_aa=seq[pos];
           if(mid_aa=='#'):
              if(posnum==0):
                 positive_dict[gene_id]=[pos-1]
              else:
                 positive_dict[gene_id]+=[pos-1-posnum]
              posnum+=1;
        
       fasta_dict[gene_id]=fasta_dict[gene_id].replace('#','') #delete all #
    return fasta_dict,positive_dict,idlist

def get_frag(fasta_dict,positive_dict,idlist,nb_windows, empty_aa,focus):
    seq_list_2d = []
    
    for id in idlist: #for sort
        seq = fasta_dict[id]
        if(id in positive_dict):
            positive_list=positive_dict[id]
        else:
            positive_list=[]
        for pos in range(len(seq)):
            mid_aa=seq[pos];
            #if mid_aa != "S":
            if not(mid_aa in focus):
               continue
            #print(id)
            #print(pos)
            #print(mid_aa)
            start = 0
            if pos-nb_windows>0:
                start = pos-nb_windows 
            left_seq = seq[start:pos]
            
            end = len(seq)
            if pos+nb_windows<end:
               end = pos+nb_windows+1
            right_seq = seq[pos+1:end]
            
            if len(left_seq) < nb_windows:
                if empty_aa is None:
                     continue
                nb_lack = nb_windows - len(left_seq)
                left_seq = ''.join([empty_aa for _count in range(nb_lack) ]) + left_seq
            
            if len(right_seq) < nb_windows:
              if empty_aa is None:
                continue
              nb_lack = nb_windows - len(right_seq)
              right_seq = right_seq + ''.join([empty_aa for _count in range(nb_lack) ])
            
            final_seq = left_seq + mid_aa + right_seq
            
            if(pos in positive_list):
                pid=id.split('|')[1]
                final_seq_list = [pid+","+str(pos)]+["0"] + [ AA for AA in final_seq]
                seq_list_2d.append(final_seq_list)
            
    df = pd.DataFrame(seq_list_2d)
    
    
    return df

fasta_dict,positive_dict,idlist = read_fasta(f1) #{id} seq  fasta_file="train_test_fasta" fasta_file="./cross-validation-protein/validation_proteins_0.fasta"
frag1 = get_frag(fasta_dict,positive_dict,idlist, 16, '-',["S"])
fasta_dict,positive_dict,idlist = read_fasta(f2) #{id} seq  fasta_file="train_test_fasta" fasta_file="./cross-validation-protein/validation_proteins_0.fasta"
frag2 = get_frag(fasta_dict,positive_dict,idlist, 16, '-',["T"])
fasta_dict,positive_dict,idlist = read_fasta(f3) #{id} seq  fasta_file="train_test_fasta" fasta_file="./cross-validation-protein/validation_proteins_0.fasta"
frag3 = get_frag(fasta_dict,positive_dict,idlist, 16, '-',["Y"])

frags_all=np.row_stack([frag1,frag2,frag3])  #1512 749

import os
folder="./all_kinase_seqs/"
#list_dirs = os.listdir(folder)
list_dirs = ['KIS_fragments.txt', 'Ack_fragments.txt', 'LISK_fragments.txt', 'DAPK_fragments.txt', 'DMPK_fragments.txt', 'EGFR_fragments.txt', 'WNK_fragments.txt', 'Fer_fragments.txt', 'Haspin_fragments.txt', 'Alpha_fragments.txt', 'STKR_fragments.txt', 'TSSK_fragments.txt', 'Ret_fragments.txt', 'PDK1_fragments.txt', 'PKC_fragments.txt', 'NKF4_fragments.txt', 'Csk_fragments.txt', 'MLCK_fragments.txt', 'CLK_fragments.txt', 'JAK_fragments.txt', 'STE11_fragments.txt', 'LRRK_fragments.txt', 'RAF_fragments.txt', 'Met_fragments.txt', 'CAMKL_fragments.txt', 'PDGFR_fragments.txt', 'IRAK_fragments.txt', 'ALK_fragments.txt', 'Tie_fragments.txt', 'RSK_fragments.txt', 'CAMK2_fragments.txt', 'CAMK1_fragments.txt', 'PIM_fragments.txt', 'GRK_fragments.txt', 'Bud32_fragments.txt', 'TTBK_fragments.txt', 'FGFR_fragments.txt', 'RAD53_fragments.txt', 'MAPKAPK_fragments.txt', 'PKN_fragments.txt', 'NAK_fragments.txt', 'PKA_fragments.txt', 'CK1_fragments.txt', 'TLK_fragments.txt', 'Abl_fragments.txt', 'SRPK_fragments.txt', 'NDR_fragments.txt', 'TTK_fragments.txt', 'NEK_fragments.txt', 'BAZ_fragments.txt', 'Akt_fragments.txt', 'CDK_fragments.txt', 'MLK_fragments.txt', 'CDC7_fragments.txt', 'Syk_fragments.txt', 'PEK_fragments.txt', 'InsR_fragments.txt', 'CAMKK_fragments.txt', 'Axl_fragments.txt', 'PDHK_fragments.txt', 'PLK_fragments.txt', 'DYRK_fragments.txt', 'CK2_fragments.txt', 'PIKK_fragments.txt', 'Src_fragments.txt', 'MAPK_fragments.txt', 'VRK_fragments.txt', 'Trk_fragments.txt', 'DDR_fragments.txt', 'RCK_fragments.txt', 'WEE_fragments.txt', 'PKD_fragments.txt', 'Eph_fragments.txt', 'STE20_fragments.txt', 'IKK_fragments.txt', 'STE7_fragments.txt', 'GSK_fragments.txt', 'FAK_fragments.txt', 'PKG_fragments.txt', 'Tec_fragments.txt', 'PHK_fragments.txt', 'SGK_fragments.txt', 'VEGFR_fragments.txt', 'Aur_fragments.txt', 'STE-Unique_fragments.txt', 'TAF1_fragments.txt']

kinase_clusters=list(range(len(list_dirs)))
kinase_list=list_dirs
num=0
maxnum=50000

for kf in list_dirs:
    fastaFileInput=folder+kf
    file = open(fastaFileInput,'r')
    lines = file.readlines()
    lnum=0;
    labels={}
    seqs={}
    for line in lines:
              id=line.strip().split("\t")[1] # pid and pos used as id
              if(id not in labels.keys()):
                  labels[id]=list_dirs.index(kf)+1 #start from 1
              else:
                  labels[id]=list_dirs.index(kf)+1 #start from 1
              
              seqs[id]='\t'.join([x for x in line.strip().split("\t")[2:]])
    
    all_seqs=[]
    tempc=0
    for (k,v) in labels.items():
       strj=""
       line=k+"\t"+str(int(v))+"\t"+seqs[k];
       all_seqs.append(line);
       if tempc >=maxnum:
           break;
       tempc+=1;
    
    kinase_clusters[num]=all_seqs
    num+=1;
    
    file.close()


all_data={} #3721->41123 merging kinase and all phosphorylation data

for frags in frags_all:
    all_data[frags[0]]="\t".join([x for x in frags])

for i in range(len(kinase_clusters)):   
         for j in range(len(kinase_clusters[i])):
              temp=kinase_clusters[i][j]
              all_data[temp.split("\t")[0]]=temp


###all_data is 41123

f=open("All_uniprot_pos_fragments_multiclass.txt",'w')
for each in all_data:
   #temp="\t".join([x for x in each])
   f.write(all_data[each]+"\n")

