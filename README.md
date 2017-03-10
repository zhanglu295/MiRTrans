# MiRTrans
MiRTrans is a framework to incorporate trans-omics sequencing data to predict miRNA targets.

MicroRNA targets prediction by Trans-omics data

Input file(examples are in the ./input):

Four input files are required for detecting miRNA targets by MiRTrans，examples are in the ./input.


1. Degradome sequencing: The raw sequencing data should be processed by Cleaveland4(http://sites.psu.edu/axtell/software/cleaveland4/). Three columns: miRNA name, transcript name and p-value with '\t' as delimiter.


2. Sequence-based prediction: The union result combined from several sequence-based prediction software, such as PITA, TargetScan, psRNATarget etc. Two columns: miRNA name and transcript name with '\t' as delimiter.


3. Transcript expression: The mRNA expression matrix as (m+1)*(n+1), where m is the number of mRNA and n is the sample size, with one row header and one column with mRNA name.


4. MiRNA expresssion:The microRNA expression matrix as (m+1)*(n+1), where m is the number of microRNA and n is the sample size, with one row header and one column with microRNA name.


Run:
python3 MiRTrans.py config1.txt
The config1.txt is the configuration file with the path of inputs.
