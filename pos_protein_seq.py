# fetch list of prot ids from GFF to be used to filter FASTA file for positive sequence
downloaded_gff_file="uniprot-filtered-organism_Human_.gff";

def pos_protein_seq():
  gff_file = open(downloaded_gff_file)
  positive_seq_ids = [];
  count = 0;
  for line in gff_file:
  	ignore_line_case = line.lower();
  	if("Modified residue".lower() in ignore_line_case  and ("Note=phosphoserine".lower() in ignore_line_case or "Note=phosphothreonine".lower() in ignore_line_case or "Note=phosphotyrosine".lower() in ignore_line_case)):
  	  positive_seq_ids.append(line.split()[0]);
  	  count = count +1;
  positive_seq_set_ids = set(positive_seq_ids);
  print("list count: "+str(len(positive_seq_ids)));
  print("set count: "+str(len(list(positive_seq_set_ids))));
  print(list(positive_seq_set_ids));
  positive_seq(list(positive_seq_set_ids))
  gff_file.close();


#filter GFF file for positive sequences using list of prot id(s)
def positive_seq(prot_id_list): 
  downloaded_fasta_file="file.fasta";
  positive_seq_fasta_file = open("pos_seq_file.fasta","w"); 
  pos_seq = [];

  with open(downloaded_fasta_file, 'r') as fasta_file:
    file_data=fasta_file.read();
    seq_list = file_data.split('>sp|');

    for i in range(len(prot_id_list)):
    	print("elem: "+str(i)+" - "+ str(prot_id_list[i]));
    	for j in range(len(seq_list)):
    	    if(prot_id_list[i] in seq_list[j]):
    	   	  positive_seq_fasta_file.write('>sp|'+seq_list[j]);
    	   	  print('>sp|'+seq_list[j]);
    	   	  break;  
  positive_seq_fasta_file.close();


pos_protein_seq()
  