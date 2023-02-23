##################################################################################################### dependencies

import os
import sys
import pysam
import shutil
import argparse
import subprocess
import pandas as pd
from os import path
from Bio import SeqIO
from random import sample

##################################################################################################### arguments

parser = argparse.ArgumentParser(prog='mapping summary', description='mapping summary')

parser.add_argument('--input_baits', "-inb",  help='input baits generated by gotcha')
parser.add_argument('--input_dbase', "-ind",  help='input sequence database from which baits have been generated')
parser.add_argument('--output', "-out",  help='output folder')
parser.add_argument('--reads_number', "-rn",  help='number of simulated reads')
parser.add_argument('--reads_length', "-rl",  help='length of simulated reads')

args=parser.parse_args()

##################################################################################################### folder system

if path.exists(args.output):
	shutil.rmtree(args.output + "/tmp")
	
os.makedirs(args.output + "/tmp")
os.chdir(args.output + "/tmp")

input_dbase =  "../../" + args.input_dbase
shutil.copy(input_dbase,"input_dbase.fa")

input_baits =  "../../" + args.input_baits
shutil.copy(input_baits,"input_baits.fa")

##################################################################################################### simulate reads from input sequence database

reads_number =  "reads=" + args.reads_number
reads_length =  "length=" + args.reads_length

with open('input_dbase.fa', 'r') as sampled, open('bbmap.log', 'w') as bbmap_log:
	subprocess.call(["bbmap.sh",  "-Xmx30g", "ref=input_dbase.fa"], stdout=bbmap_log, stderr=bbmap_log)
	subprocess.call(["randomreads.sh", "-Xmx10g", "in=input_dbase.fa", "out=simulated_reads.fq", reads_number, reads_length], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

##################################################################################################### reformat baits so that they have a univocal header

input_baits=open('input_baits.fa','r')
input_baits_ref=open('input_baits_ref.fa', 'w')

records = set()
of = open("output.fa", "w")
for record in SeqIO.parse(input_baits, "fasta"):
    ID = record.id
    num = 1
    while ID in records:
        ID = "{}#{}".format(record.id, num)
        num += 1
    records.add(ID)
    record.id = ID
    record.name = ID
    record.description = ID
    SeqIO.write(record, input_baits_ref, "fasta")
    
input_baits_ref.close()

##################################################################################################### map reads to input baits

with open('align.bam', 'w') as align_bam, open('bowtie.log', 'w') as bowtie_log:
	subprocess.call(["bowtie2-build", "input_baits_ref.fa", "input_index"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	subprocess.call(["bowtie2", "-x", "input_index", "-U", "simulated_reads.fq", "--very-sensitive-local" ], stdout=align_bam, stderr=bowtie_log)
	
##################################################################################################### stats part 1

samfile = pysam.AlignmentFile("align.bam", "rb")

with open('mapping.log', 'w') as mapping_log:
	
	i=0
	
	mapping_log.write("ref_seq\tbait_seq\taln_lenght\tmissmatches\tcigar\n")
	
	for seq in samfile :
	    if(seq.flag!=4): 
	        mapping_log.write(seq.get_reference_sequence()+"\t"+
	        seq.query_sequence[seq.query_alignment_start:seq.query_alignment_end]+"\t"+
	        str(len(seq.get_reference_sequence()))+"\t"+
	        str(seq.get_tag("XM:i"))+"\t"+
	        seq.cigarstring + "\n")
	    else:
	        mapping_log.write("NA"+"\t"+
	        seq.query_sequence[seq.query_alignment_start:seq.query_alignment_end]+"\t"+
	        "NA"+"\t"+
	        "NA"+"\t"+
	        "NA" + "\n")
	
##################################################################################################### stats part 2

df=pd.read_csv('mapping.log', sep='\t')
tot_reads=int(args.reads_number)

lesstthan05mismatches=len(df[df.missmatches > 6])
lesstthan05mismatchespercentage=round((lesstthan05mismatches/tot_reads)*100)
print("sequences with >05 mismatches with baits")
print("total =" , lesstthan05mismatches, "\t", "percentage =", str(lesstthan05mismatchespercentage) + "%\n")

lesstthan10mismatches=len(df[df.missmatches > 11])
lesstthan10mismatchespercentage=round((lesstthan10mismatches/tot_reads)*100)
print("sequences with >10 mismatches with baits")
print("total =" , lesstthan10mismatches, "\t", "percentage =", str(lesstthan10mismatchespercentage) + "%\n")