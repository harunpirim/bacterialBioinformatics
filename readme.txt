General instructions for using this repository:

    File uses:

        1) combined.ipynb: takes json file and fasta file downloaded from uniprot and compiles features to be used in training
        2) training_prep.ipynb: converts json file of features to a tsv file that can be trainied on (e.g. encodes qualitative features)
        3) ml_genome.ipynb: performs training, validation, prediction on prepared tsv file
        4) other_code folder: older files no longer necessary 
    
    Before using repository:
        1) select protein entries from Uniprot and download both the fasta and json file for the selection of protein entries.
        2) use iFeature libary to extract protein features such as DPC, CTriad, Moran, etc. (full list of required features found at top of combined.ipynb)
        3) adjust GO terms to be used as lables in both combined.ipynb and training_prep.ipynb
    
    During use of repository:
        1) adjust file location names to fit your machine
        2) be aware that ipynb files may need to be run in smaller chunks if protein selection is large due to memory 

How to use repository:
    1) combined.ipynb:
        -begin with this file.
        -after ensuring file location variables are set correctly, run this file.
        -file output of this file should be train_protein_information.json and train_proteins.fasta
    
    2) training_prep.ipynb:
        -notice that the file can prepare training and test sets, but test sets are commented out.
        -run this file on the files created from combined.ipynb to get a tsv file of data ready for training in ml_genome.ipynb

    3) ml_genome.ipynb:
        -running this file should produce predictions for each of the three models        