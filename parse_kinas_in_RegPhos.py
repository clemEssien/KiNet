import numpy as np
import re
import os

fd1=open("C:/Users/ucexbw/Documents/Phos/RegPhos_kinase_human_cluster.txt",'r')
focus_type=fd1.readline().split("\t") #head
focus_num=focus_type.index('Family') ## can be 'Kinase', 'Group', 'Family', 'Subfamily' Family is CK1 CK2.... try this first

kinase_list={}
for l in fd1.readlines():
    kinase_list[l.split("\t")[1]]=l.split("\t")[focus_num]; # focus_type: kinase

print("number of kinases "+str(len(kinase_list))) #Family 141 Group 12

fd2=open("C:/Users/ucexbw/Documents/Phos/RegPhos_Phos_human.txt",'r')
fd2.readline() # remove header
kinase_information={}
for l in fd2.readlines():
    kinase=l.split("\t")[4] # kinase;  xxx_group could be family/subfamily/group; xxx(xxx) front is kinase
    if kinase == '': # no kinase information
        continue
    
    if '_group' in kinase:
        search=re.match('(.*)_group',kinase).group(1)
        if search not in kinase_list.values():
               print(search+" not in kinase_list family")
               continue
               #break
        
    elif '(' in kinase:
        if re.match('(.*)\(.*\)',kinase).group(1) not in kinase_list.keys():
             print(re.match('(.*)\(.*\)',kinase).group(1)+" 1 not in kinase_list kinase")
             continue
        
        search=kinase_list[re.match('(.*)\(.*\)',kinase).group(1)]
    else:
        if kinase not in kinase_list.keys():
            print(kinase+" 2 not in kinase_list kinase")
            continue
            #break
        
        search=kinase_list[kinase]
    
    pid=l.split("\t")[1]
    pos=l.split("\t")[2] #start from 1 in python should be start from 0 later will be -1 
    residue=l.strip().split("\t")[7]
    if search not in kinase_information.keys():
        kinase_information[search]=[]
        kinase_information[search].append(str(pid)+','+str(pos)+','+str(residue))
    else:
        kinase_information[search].append(str(pid)+','+str(pos)+','+str(residue))
    

len(kinase_information) #86
############extract seqs####################
#fd3=open("/home/dlg/wangdu/deeplearningforgenomic/DATA/case_raw/Uniprot-reviewed-all/uniprot-all-reviewed_2017_8_31.fasta",'r');
fd3=open("C:/Users/ucexbw/Documents/Phos/pos_seq_file.fasta",'r');
lines=fd3.readlines()
seqs={}
lnum=0
for line in lines:
    if line[0] == '>':
        if lnum !=0:
              seqs[pid]=seq
           
        pid=line.strip().split('|')[1]
        seq = ""
    else:
        seq +=line.strip()
        
    lnum+=1;
    
    #for last
    seqs[pid]=seq
    fd3.close()

##########################################
#######extract all positive fragments#####
windows=16
fout2=open("record_statistics.txt",'w')
pos_all=0;
for (k,v) in kinase_information.items():
    print("kinase family is "+k);
    fout=open("./all_kinase_seqs/"+str(k)+"_fragments.txt",'w')
    pos_num=0;
    for each in kinase_information[k]:
        pid=each.split(",")[0]
        pos=int(each.split(",")[1])-1
        type=each.split(",")[2]
        frags=[];
        if pid in seqs: # is uniprot data
           pos_num+=1;
           pos_all+=1;
           start=0;
           if pos-windows>0:
              start=pos-windows
           left_seq=seqs[pid][start:pos]
           end=len(seqs[pid])
           if pos+windows<end:
              end=pos+windows+1
           right_seq=seqs[pid][pos+1:end]
           if len(left_seq) < windows:
                nb_lack = windows - len(left_seq)
                left_seq=''.join(['-' for _count in range(nb_lack)]) + left_seq
           if len(right_seq) < windows:
                nb_lack = windows - len(right_seq)
                right_seq = right_seq + ''.join(['-' for _count in range(nb_lack) ])
           
           final_seq = left_seq+seqs[pid][pos]+right_seq
           if(seqs[pid][pos] == "S"): #s 1 T 2 Y 3 other 4
                 print_seq = "1\t"+pid+','+str(pos)+'\t'+'\t'.join([x for x in final_seq])
           elif(seqs[pid][pos] == "T"):
                 print_seq = "2\t"+pid+','+str(pos)+'\t'+'\t'.join([x for x in final_seq])
           elif(seqs[pid][pos] == "Y"):
                 print_seq = "3\t"+pid+','+str(pos)+'\t'+'\t'.join([x for x in final_seq])
           else: 
                 print_seq = "4\t"+pid+','+str(pos)+'\t'+'\t'.join([x for x in final_seq])
           fout.write(print_seq+"\n");
        else:
           continue
    
    #print("number of positive sites of "+k+" "+str(pos_num)+"\n");
    fout2.write(k+" "+str(pos_num))
    fout.close()

fout2.close()
print(pos_all) #5727 
 
######################################################################    