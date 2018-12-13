"# KiNet" 
To train, run the command:

python train_models.py -input input-fasta-file -output-prefix [prefix of pre-trained model] -residue-types specify residue types

for example;

python train_models.py -input training_file.fasta -output-prefix example -residue-types S,T,Y


To predict on test dataset, run the command

python predict.py -input input-fasta-file -model-prefix [prefix of pre-trained model] -output [specify output file] 

for example;

python predict.py -input test_file.fasta -model-prefix example -output output 
