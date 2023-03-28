#! /usr/bin/env python

#####################################################################################
#                                                                                   #
#                    gotha catch'em all (envronimental DNA baits)                   #
#                                                                                   #
#####################################################################################

#notes Kevin;
        # Taxon filtering in combination with clustering sequences
# bug notes
        # python3 ../gotcha_git/gotcha/gotcha.3.0.1.py -bp ./bold-cli -t "Apis,Eucera,Ceratina" -m COI-5P -c 5 -o Bees_beta_3.1 -ft species -e -v -cl 0.99
        # size selection bug

# new
        # tree output with branch length.
        # Fixed gap problem with selecting baits (no bait should have a gap)
        # Fixed problem with node selection. There is a problem that sequence close to the root are placed on tips.
        # Fixed, for each tip sequence look back into the tree to the closest ancestral node.

import re
import os
import re
import sys
import glob
import shutil
import os.path
import datetime
import argparse
import operator
import subprocess
from os import path
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
import ete3 #new v2.2.2
from ete3 import Tree  #new v2.2.2
from statistics import mean #new v2.2.2
import statistics #new v2.3.2
from collections import Counter #new v2.2.2

##################################################################################################### parsing arguments and setting variables

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

def check_positive(value):
    try:
        value = int(value)
        if value <= 0:
            raise argparse.ArgumentTypeError("{} is not a positive integer".format(value))
    except ValueError:
        raise Exception("{} is not an integer".format(value))
    return value

parser = argparse.ArgumentParser(prog='gotcha', description='gotcha version 1.0.0 built Feb 2023\ndeveloped by Kevin Nota & Giobbe Forni\n',add_help=False)

parser.formatter_class = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=150)
formatter_class=argparse.RawTextHelpFormatter

man = parser.add_argument_group('MARKER')
man.add_argument('-m', '--marker', metavar='\b', choices=['COI-5P', 'COI-3P', 'rbcL', 'matK', 'ITS', 'coding' , 'noncoding'], help='can be COI/rbcL/matK/ITS or coding/noncoding using custom .fna/.aln')
man.add_argument('-c', '--code', metavar='\b', help='genetic code - e.g. 1 for plastid and nuclear, 5 for mitochondrial invertebrate')

bld = parser.add_argument_group('BOLD DOWNLOAD')
bld.add_argument('-t', '--taxa', metavar='\b', help='input taxa')
bld.add_argument('-g', '--geo', metavar='\b', default="", help='geographic location of samples')
bld.add_argument('-ft', '--filter_tax', metavar='\b', default="", choices=['species', 'genus', 'family'], help='can be species/genus/family - defeault is none')
bld.add_argument("-bp", "--boldcli_path", metavar='\b', default="bold-cli", help="path to bold-cli")

ctm = parser.add_argument_group('CUSTOM INPUTS')
ctm.add_argument('-cf', '--custom_fna', metavar='\b', default="", help='custom nucleotide sequences file (.fasta format)')
ctm.add_argument('-ca', '--custom_aln', metavar='\b', default="", help='custom nucleotide alignment file (.fasta format) - gotcha will jump to size selection')
ctm.add_argument('-cn', '--custom_nwk', metavar='\b', default="", help='custom newick file - tree can be multifurcating and can lack taxa')

par = parser.add_argument_group('BAITS SEARCH PARAMETERS')
par.add_argument('-fl', '--filter_len', metavar='\b', default=100, help='minimum amminoacids / nucleotides length - default is 100', type=check_positive)
par.add_argument('-sc', '--seq_collapse', metavar='\b', type=float, default=1, help='percent identity to collapse sequences - default is 1')
par.add_argument('-bc', '--bts_collapse', metavar='\b', type=float, default=1, help='percent identity to collapse baits - defeault is 1')
par.add_argument('-bl', "--baitlength", metavar='\b', type=int, default=80, help='lenght of the bait sequences - default is 80')
par.add_argument('-tl', "--tiling", metavar='\b', default=30, type=int, help='kmer tiling for ancestral node selection - default is 30') #changed the default from 10 to 30 (think most cases ~3x coverage is alright and for a 80mer, that would be ~30)
par.add_argument('-ds', "--distance", metavar='\b', type=float, default=0.09, help='maximal distance to the ancestral node - default is 0.09')
par.add_argument('-so', "--trim_seqoverlap", metavar='\b', default="", help='trimal param between 0 and 100 - if not specified will skip')
par.add_argument('-ro', "--trim_resoverlap", metavar='\b', default="", help='trimal param between 1 and 0 - if not specified will skip')
par.add_argument('-dt', "--dust_threshold", metavar='\b', default=1, help='score threshold for low complexity region masking - default is 1')
par.add_argument('-mf', "--max_farg_lenght", metavar='\b', default=3000, type=int, help='this value can restrain the selected fragment size by choosing the optimal fragment below this value, \
                                                                                                                                                                                \ndefault will take the highest value if number of sequences multiplied by fragment size') #new in version 2.3.2
par.add_argument('-fs', '--fast', action='store_true', help='performs the fast baits search w/out ancestral state reconstruction')
par.add_argument('-mc', '--manual_check', action='store_true', help='stops prior to tree inference for a manual check of the alignment')

otr = parser.add_argument_group('OTHER')
otr.add_argument("-h", "--help", action="help", help="show this help message and exit")
otr.add_argument("-path", default=".", help=argparse.SUPPRESS)
otr.add_argument('-o', '--out', metavar='\b', default="probes", help='basename of output file and folders - default is probes')
otr.add_argument('-v', '--verbose', action='store_false', help='keeps temporary folder and files')
otr.add_argument('-e', '--erase', action='store_true', help='erases and rewrites a pre existing output folder')
otr.add_argument('-th', '--threads', metavar='\b', type=check_positive, default=1, help='number of threads used for tree inference - default is 1')
otr.add_argument('--version', action='store_true', help="show program's version number")

args=parser.parse_args()

#args.out = args.out.replace("/","")

def check_if_all_is_installed():

    if args.version == 1:
        print('\ngotcha version 5.0.0 built March 2023\ndeveloped by Kevin Nota & Giobbe Forni\n')
        exit()

    if args.boldcli_path != "bold-cli" :
        args.boldcli_path = args.boldcli_path
    rc = subprocess.call(['which', args.boldcli_path],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0 and not args.custom_fna:
        print('bold-cli\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'baseml'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('baseml\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'transeq'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('transeq\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'cd-hit'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('cdhit\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'iqtree'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('iqtree\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'trimal'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('trimal\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'translatorx'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('translatorx\tmissing in path')
        exit()

    rc = subprocess.call(['which', 'muscle'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if rc != 0:
        print('muscle\tmissing in path')
        exit()

    ##################################################################################################### errors !!!

    if not args.custom_fna and not args.custom_aln and (args.taxa is None):
        print("\n WARNING! Either a custom alignment (--custom_fna) or a lineage to download from BOLD (--taxa) is required! \n")
        quit()

    if  args.custom_fna and args.custom_aln:
        print("\n WARNING! both custom sequences and alignment have been specified! \n")
        quit()

    if not args.custom_fna and (args.marker is None):
        print("\n WARNING! When not using custom alignment --marker flags is required! \n")
        quit()

    if  args.custom_nwk and (args.custom_fna is None):
        print("\n WARNING! When using a custom tree (--custom_nwk) either custom sequences (--custom_fna) or alignment (--custom_aln) has to be specified! \n")
        quit()

    if args.erase == False and path.exists(args.out):
        print("\n WARNING! An output folder with the same name already exists! Use --erase to overwrite \n")
        quit()
    elif args.erase == True and path.exists(args.out):
        shutil.rmtree(args.out)
        os.makedirs(args.out)
    else:
        os.makedirs(args.out)

    if (args.code != None):
        if int(args.code) not in range (1, 33):
            print("\n WARNING! Unknown geneitc code has been specified - please refer to NCBI! \n")
            quit()

    if not 0.01 <= args.distance <= 0.2:
        print("\n WARNING! A value between 0.01 and 0.2 should be specified for parameter distance! \n")
        quit()

    if not 1 >= args.seq_collapse >= 0.7:
        print("\n WARNING! A value between 1 and 0.7 should be specified for parameter seq_collapse! \n")
        quit()    
    
    if not 1 >= args.bts_collapse >= 0.7:
        print("\n WARNING! A value between 1 and 0.7 should be specified for parameter bts_collapse! \n")
        quit()  

    if args.trim_seqoverlap != "" and args.trim_resoverlap == "":
        print("\n WARNING! You have to specify both seqoverlap and resoverlap parameters for the trimal step \n")
        quit()

    if args.trim_seqoverlap == "" and args.trim_resoverlap != "":
        print("\n WARNING! You have to specify both seqoverlap and resoverlap parameters for the trimal step \n")
        quit()

    if args.trim_seqoverlap != "":
        if not 100 >= float(args.trim_seqoverlap) >= 1:
            print("\n WARNING! A value between 100 and 1 should be specified for parameter seqoverlap! \n")
            quit()

    if args.trim_resoverlap != "":
        if not 1 >= float(args.trim_resoverlap) >= 0:
            print("\n WARNING! A value between 1 and 0 should be specified for parameter resoverlap! \n")
            quit()

    os.makedirs(args.out + "/tmp")
    os.chdir(args.out + "/tmp")

    ##################################################################################################### define markers

    if ( args.marker == "coding" ):
        coding = True
    if ( args.marker == "non-coding" ):
        coding = False

    if ( args.marker == "COI-5P" ):
        marker_list=[ "COI-5P" ]
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    if ( args.marker == "COI-3P" ):
        marker_list=[ "COI-3P" ]
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    elif ( args.marker == "rbcL" ):
        marker_list=[ "rbcL", "rbcl" , "RBCL" , "Rbcl" , "rbcla" , "rbcl-a" ]
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    elif ( args.marker == "matK" ):
        marker_list=[ "matK" , "matk" , "MATK" ]
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    elif ( args.marker == "ITS" ):
        marker_list=[ "ITS" , "its" ]
        coding = False
        if ( args.code != None ):
            print("\n WARNING! a genetic code has been specified when using a non-coding marker! \n")
            quit()

    if ( args.custom_fna != "" ) and ( args.marker == "coding" ):
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    elif ( args.custom_fna != "" ) and ( args.marker == "noncoding" ):
        coding = False
        if ( args.code != None ):
            print("\n WARNING! a genetic code has been specified for a non-coding marker! \n")
            quit()

    if ( args.custom_aln != "" ) and ( args.marker == "coding" ):
        coding = True
        if ( args.code == None ):
            print("\n WARNING! a genetic code has to be specified when using a coding marker! \n")
            quit()
    elif ( args.custom_aln != "" ) and ( args.marker == "noncoding" ):
        coding = False
        if ( args.code != None ):
            print("\n WARNING! a genetic code has been specified for a non-coding marker! \n")
            quit()
    print("\nanalysis started on" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "\n")
    
    return coding, marker_list

def write_log_file_01():
    with open("log.txt", "a+") as log:

        log.write("# analysis started on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")

    # INPUTS

        log.write("\n\n\t INPUTS:" + "\n\n")

        if ( args.custom_fna != "" ):
            log.write("\t BOLD TAXONOMY \t\tN" + "\n")
            log.write("\t CUSTOM SEQUENCES \tY" + "\n")
            log.write("\t CUSTOM ALIGNMENT \tN" + "\n")
        elif ( args.custom_aln != "" ):
            log.write("\t BOLD TAXONOMY \t\tN" + "\n")
            log.write("\t CUSTOM SEQUENCES \tN" + "\n")
            log.write("\t CUSTOM ALIGNMENT \tY" + "\n")
        else:
            if ( args.geo == "" ):
                log.write("\t BOLD TAXONOMY \t\t" + str(args.taxa) + "\n")
                log.write("\t GEOGRAPHY \t\tN" + "\n")
                log.write("\t CUSTOM SEQUENCES \tN" + "\n")
                log.write("\t CUSTOM ALIGNMENT \tN" + "\n")
            else:
                log.write("\t BOLD TAXONOMY \t\t" + str(args.taxa) + "\n")
                log.write("\t GEOGRAPHY \t\t" + args.geo + "\n")
                log.write("\t CUSTOM SEQUENCES \tN" + "\n")
                log.write("\t CUSTOM ALIGNMENT \tN" + "\n")

        if ( args.custom_nwk == "" ):
            log.write("\t CUSTOM TREE \t\tN" + "\n")
        else:
            log.write("\t CUSTOM TREE \t\tY" + "\n")

        log.write("\t MARKER \t\t" + str(args.marker) + "\n")

        if (coding == True):
            log.write("\t GEN CODE \t\t" + str(args.code) + "\n")
        else:
            log.write("\t GEN CODE \t\t" + "-" + "\n")

        if (args.fast == False):
            log.write("\t FAST MODE \t\t" + "N" + "\n")
        else:
            log.write("\t FAST MODE \t\t" + "Y" + "\n")

    # PARAMETER

        log.write("\n\n\t PARAMETERS:" + "\n\n")
        if ( args.filter_tax == "" ):
            log.write("\t TAXONOMIC FILTER \tN" + "\n")
        else:
            log.write("\t TAXONOMIC FILTER \t" + str(args.filter_tax) + "\n")
        log.write("\t MIN MARKER LENGTH \t" + str(args.filter_len) + "\n")
        log.write("\t MAX MARKER IDENTITY \t" + str(args.seq_collapse) + "\n")
        log.write("\t BAITS LENGTH \t\t" + str(args.baitlength) + "\n")
        log.write("\t BAITS TILING \t\t" + str(args.tiling) + "\n")
        if ( args.trim_seqoverlap != "" ) and ( args.trim_resoverlap != "" ):
            log.write("\t SEQ OVERLAP \t\t" + str(args.trim_seqoverlap) + "\n")
            log.write("\t RES OVERLAP \t\t" + str(args.trim_resoverlap) + "\n")
        else:
            log.write("\t TRIMAL \t\tN" + "\n")
        log.write("\t DUST TRESHOLD \t\t" + str(args.dust_threshold) + "\n")

def bold_download():

    if args.boldcli_path != "bold-cli" :
        args.boldcli_path = "../../" + args.boldcli_path

    if ( args.geo == "" ):
        print("downloading" , args.marker , "sequences for taxa" , args.taxa)
        subprocess.run([args.boldcli_path , "-taxon" , args.taxa , "-marker" , args.marker , "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
        if os.path.getsize("tmp.bold") == 0 :
            print("\n WARNING! The BOLD search for marker" , args.marker , "and taxa" , args.taxa , "returned nothing! \n")
            quit()
    else:
        print("downloading" , args.marker , "sequences for taxa" , args.taxa, "from" , args.geo)
        subprocess.run([args.boldcli_path , "-taxon" , args.taxa , "-marker" , args.marker , "-geo" , args.geo, "-output" , "tmp.bold"] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)
        if os.path.getsize("tmp.bold") == 0 :
            print("\n WARNING! The BOLD search for marker" , args.marker , "and taxa" , args.taxa , "in" , args.geo , "returned nothing! \n")
            quit()

    num0 = sum(1 for line in open('tmp.bold', errors='ignore'))
    print("#        downloaded: ", num0)

    with open("log.txt", "a+") as log:
        log.write("\n\n\t FILTERING:" + "\n\n")
        log.write("\t " + str(num0) + "\t starting sequences" + "\n")
    log.close()

    with open("tmp.bold", 'r', errors='ignore') as file:
        if 'Fatal Error' in file.read():
            print("\n WARNING! The taxonomy specified appears to be not correct! \n")
            quit()

##################################################################################################### filter marker and taxonomy

    tmp_1_fasta=[]

    if (args.filter_tax != ""):
        print("filtering sequences with" , args.filter_tax , "identification")
    else:
        print("skipping taxonomic filtering")

    with open("tmp.bold", 'r', errors='ignore') as file:
        header = file.readline()
        for l in file :
            sl = l.split('\t')
            if (args.filter_tax == "family"):
                if sl[15]:
                    header = (">" + sl[0])
                    seq = sl[71].replace('-','')
                    if sl[69] in marker_list:
                        tmp_1_line=[header, seq]
                        tmp_1_fasta.append(tmp_1_line)
            elif (args.filter_tax == "genus"):
                if sl[19]:
                    header = (">" + sl[0])
                    seq = sl[71].replace('-','')
                    if sl[69] in marker_list:
                        tmp_1_line=[header, seq]
                        tmp_1_fasta.append(tmp_1_line)
            elif (args.filter_tax == "species"):
                if sl[21]:
                    header = (">" + sl[0])
                    seq = sl[71].replace('-','')
                    if sl[69] in marker_list:
                        tmp_1_line=[header, seq]
                        tmp_1_fasta.append(tmp_1_line)
            else:
                header = (">" + sl[0])
                seq = sl[71].replace('-','')
                if sl[69] in marker_list:
                    tmp_1_line=[header, seq]
                    tmp_1_fasta.append(tmp_1_line)

    with open('tmp1.fna', 'w') as tmp_1:
        for line in tmp_1_fasta:
            for element in line:
                tmp_1.write(str(element) + '\n')

    if os.path.getsize("tmp1.fna") == 0 :
        print("\n WARNING! No sequence passed the taxonomic filter! \n")

    num1 = len([1 for line in open("tmp1.fna") if line.startswith(">")])                                                    # count tmp1.fna
    print("#        kept: ", num1)
    
    with open("log.txt", "a+") as log:
        log.write("\t " + str(num1) + "\t after taxonomic filtering" + "\n")
    log.close()

def seq_quality_filters():
    if (args.custom_fna != ""):
        print("custom sequences detected\n")
        custom_fna =  "../../" + args.custom_fna
        shutil.copy(custom_fna,"tmp1.fna")

    if (args.custom_aln == ""):
        print("collapsing identical sequences with a percent identitiy of", args.seq_collapse)
        subprocess.run(["cd-hit","-i","tmp1.fna","-o","tmp2a.fna","-c", str(args.seq_collapse)] , stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

        num2 = len([1 for line in open("tmp2a.fna") if line.startswith(">")])                                                   # count tmp2.fna
        print("#        kept: ", num2)
        with open("log.txt", "a+") as log:
            log.write("\t " + str(num2) + "\t after identity collapse" + "\n")
        log.close()
    ##################################################################################################### remove sequences with non ATGC nucleotides

        print("removing sequences with non ATGC nucleotides")

        selected_seqs_ns = list()

        tmp2b_fna=open("tmp2b.fna",'w')

        for record in SeqIO.parse("tmp2a.fna", "fasta"):
    #               if record.seq.count('N') == 0:
            matches = ['B','D','E','F','H','I','J','K','L','M','N','O','P','Q','R','S','U','V','W','X','Y','Z','b','d','e','f','h','i','j','k','l','m','n','o','p','q','r','s','u','v','w','x','y','z']
            if not any(x in record.seq for x in matches):
                selected_seqs_ns.append(record)
        nseq=(len(selected_seqs_ns))
        SeqIO.write(selected_seqs_ns, tmp2b_fna , "fasta")
        tmp2b_fna.close()
        print("#        kept: ", nseq)

        with open("log.txt", "a+") as log:
            log.write("\t " + str(nseq) + "\t removing sequences with non ATGC nucleotides" + "\n")
        log.close()
    ##################################################################################################### length cutoff

        print("filtering sequences with a minimum length of", args.filter_len)

        selected_seqs_nt = list()
        tmp2_fna=open("tmp2.fna",'w')

        for record in SeqIO.parse("tmp2b.fna", "fasta"):
            if len(record.seq) >= int(args.filter_len):
                selected_seqs_nt.append(record)
        nseq=(len(selected_seqs_nt))
        SeqIO.write(selected_seqs_nt, tmp2_fna , "fasta")

        print("#        kept: ", nseq)
        tmp2_fna.close()

        with open("log.txt", "a+") as log:
            log.write("\t " + str(nseq) + "\t after length filter" + "\n")
        log.close()
    ##################################################################################################### too few sequences to proceed?

        nseq = len([1 for line in open("tmp2.fna") if line.startswith(">")])
        if nseq < 4:
            print("\n WARNING! Less than 4 sequences passed the filtering steps! \n" )
            quit()

    ##################################################################################################### if the marker is coding

    if (coding == True and args.custom_aln == ""):
        #print("aligning sequences")

        # adjust direction with mafft
        with open('tmp3.fna', 'w') as tmp3_fna, open('tmp2.fna', 'r') as tmp2_fna:
            if 'MAFFT_BINARIES' in os.environ:
                os.environ.pop('MAFFT_BINARIES')
            subprocess.call(["mafft" , "--thread", str(args.threads), "--adjustdirection" , "tmp2.fna"], stdout=tmp3_fna, stderr=subprocess.DEVNULL)

        # remove "N" and "-"
        with open("tmp4.fna", "w") as tmp4_fna:
            for record in SeqIO.parse("tmp3.fna", "fasta"):
                record.seq = record.seq.ungap("-")
                record.seq = record.seq.ungap("n")
                SeqIO.write(record, tmp4_fna, "fasta")

        # protein-guided nucleotide alignment with translatorx
        with open('tmp4.fna', 'r') as tmp4_fna, open('translatorx.log', 'w') as translatorx_log:
            subprocess.call(["translatorx" , "-i" , "tmp4.fna", "-t", "T", "-o", "tmp", "-c", args.code], stdout=translatorx_log, stderr=translatorx_log)
        os.rename("tmp.nt_ali.fasta", "tmp5.fna")

        # translate sequences to find stop codons
        with open('tmp5.fna', 'r') as tmp5_fna:
            subprocess.call(["transeq" , "-sequence" , "tmp5.fna", "-outseq", "translated4stop.fna", "-table", args.code], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # find sequences with stopcodons
        identifiers = set()
        with open('translated4stop.fna', 'r') as translated4stop_fna:
            for seq_record in SeqIO.parse(translated4stop_fna, "fasta"):
                if "*" in (seq_record.seq[0:(len(seq_record.seq)-1)]):
                    identifiers.add(seq_record.id[:-2])

        print("removing sequences with stop codons")

        # find sequences with stopcodons
        with open('tmp5.fna', 'r') as tmp5_fna, open('tmp6.fna', 'w') as tmp6_fna:
            records = SeqIO.parse(tmp5_fna, 'fasta')
            for record in records:
                if record.id not in identifiers:
                    SeqIO.write(record, tmp6_fna, 'fasta')

        num6 = len([1 for line in open("tmp6.fna") if line.startswith(">")])                                                    # count tmp6.fna
        print("#        kept: ", num6)

        with open("log.txt", "a+") as log:
            log.write("\t " + str(num6) + "\t the removal of sequences with stop codons" + "\n")
        log.close()

        os.rename("tmp6.fna", "tmp.aln")

    ##################################################################################################### if the marker is not coding

    if (coding == False and args.custom_aln == ""):

        print("aligning non-coding sequences")

        with open('tmp.aln', 'w') as tmp_aln, open('tmp2.fna', 'r') as tmp2_fna:
            if 'MAFFT_BINARIES' in os.environ:
                os.environ.pop('MAFFT_BINARIES')
            subprocess.call(["mafft" , "--thread", str(args.threads) , "--adjustdirection" , "tmp2.fna"], stdout=tmp_aln, stderr=subprocess.DEVNULL)

    ##################################################################################################### reformat after mafft adjustdirection

    os.rename('tmp.aln', 'tmptmp.aln')

    tmp_aln=open("tmp.aln",'w')

    for record in SeqIO.parse("tmptmp.aln", "fasta"):
        record.description = (record.id).replace('_R_', '')
        record.id=record.description
        SeqIO.write(record, tmp_aln , "fasta")

    tmp_aln.close()

    os.remove('tmptmp.aln')

    ##################################################################################################### can enter custom .aln

    if (args.custom_aln != ""):
        print("custom alignment detected")
        custom_aln =  "../../" + args.custom_aln
        shutil.copy(custom_aln,"tmp.aln")
        num = len([1 for line in open("tmp.aln") if line.startswith(">")])

def size_selection():

    input_sequences=list(SeqIO.parse("tmp.aln", "fasta")) # path to the fasta file
    fragment_selections= pd.DataFrame()
    #counting the number of gaps in the alignment file
    gaps = pd.DataFrame()

    print("\nselecting optimal fragment size based on", len(input_sequences),":")# "sequences with max length of",len(input_sequences[0].seq), "\n")

    #print("\nnumber of input sequences =", len(input_sequences))
    #print("max sequences lenght =", len(input_sequences[0].seq))

    #print("calculating missing data")
    #progress=tqdm(total=None, desc="Alignment scanner ")
    for p in range(len(input_sequences[0].seq)):
        gap_positions = 0
        for i in range(len(input_sequences)):
            if (input_sequences[i].seq[p] == "-"):
                gap_positions += 1
            data = {'position':[p], 'proportion gaps':[gap_positions/len(input_sequences)]}

        gaps=gaps.append(pd.DataFrame(data), ignore_index=True)
       #progress.update()
    #progress.close()

    #progress=tqdm(total=None, desc="Calculating distance = ")
    for treshold in range(0, 100):
        threshold_found=False
        i=0

        # find the first and last position where thresholds values are met
        # tressholds that are tested are 0-1 with 0.01 increments (101 loops)
        while not (threshold_found==True) :
            if (gaps.loc[i,"proportion gaps"]<=float(treshold/100)):
                threshold_found=True
                start_position=gaps.loc[i,"position"]
            i=i+1
            if(i==len(gaps)):
                threshold_found=True
                start_position=0

        threshold_found=False
        i=len(gaps)-1
        while not (threshold_found==True) :
            if (gaps.loc[i,"proportion gaps"]<=float(treshold/100)):
                threshold_found=True
                finish_position=gaps.loc[i,"position"]
            i=i-1
            if(i==0):
                threshold_found=True
                finish_position=0
        count=0

        lenght_no_gaps=[]
        # this counts the number of sequences withing the threshold with no gaps (-)
        for i in range(len(input_sequences)):
            if(len(re.findall("^-|-$|^N|N$", str(input_sequences[i].seq[start_position:finish_position])))==0):
                count+=1
                lenght_no_gaps.append(len(re.sub("-", "", str(input_sequences[i].seq[start_position:finish_position]))))
        if len(lenght_no_gaps) == 0 :
            continue
        data = {'threshold':[str(treshold)+"%"],
            'number of sequences':[count],
            'seq_length':max(Counter(lenght_no_gaps).most_common())[0], #statistics.mode(lenght_no_gaps), 
            'number * lenght':[count*len(input_sequences[i].seq[start_position:finish_position])],
            'start position':[start_position],
            'finish_position':[finish_position]}
        fragment_selections = fragment_selections.append(pd.DataFrame(data), ignore_index=True)
        #progress.update()
    fragment_selections.to_csv("size_selection_summary.df", sep="\t") # new version 2.3.2
    #progress.close()    

    # loop to find the most optimal fragment lenght
    find_optimal_fragment=False
    i=0

    fragment_selections = fragment_selections.loc[fragment_selections['seq_length'] <= args.max_farg_lenght] # new version 2.3.2 (gaps is a problem here)

    while not (find_optimal_fragment==True):
        if((fragment_selections.loc[i,"number * lenght"]==fragment_selections['number * lenght'].max())==True):
            #print(fragment_selections.loc[i])
            optimal_threshold=fragment_selections.loc[i, 'threshold']
            find_optimal_fragment=True
            start_position=fragment_selections.loc[i,'start position']
            finish_position=fragment_selections.loc[i,'finish_position']
            selected_size=fragment_selections.loc[i,'seq_length']
        seq_count=fragment_selections.loc[i, 'number of sequences']
        i=i+1


    while not (((finish_position-start_position)/3.0).is_integer()):
        finish_position-=1
        lenght=((finish_position-start_position)/3.0)

    #correcting size for the triplet cutoff (prevent sequences being lost due to gaps)
    count_eng_gap=0
    frag_count=0
    while not (seq_count<=(frag_count)):
        count_eng_gap=0
        frag_count=0
        for i in range(len(input_sequences)):
            if(len(re.findall("-$", (str(input_sequences[i].seq[start_position:finish_position]))))!=0):
                count_eng_gap+=1
            if(len(re.findall("^-|-$|^N|N$", (str(input_sequences[i].seq[start_position:finish_position]))))==0):
                frag_count+=1
        finish_position=finish_position-3

    #Creating new fasta file with original labels but new sequence
    new_fasta=list()
    frag_count=0
    for i in range(len(input_sequences)):
        if(len(re.findall("^-|-$|^N|N$", str(input_sequences[i].seq[start_position:finish_position])))==0):  #removes all sequences that start with w gap, and sequences that contain N's
            new_seq=SeqRecord(Seq(str(input_sequences[i].seq[start_position:finish_position])),
            id=input_sequences[i].id,
            description="")
            new_fasta.append(new_seq)
            frag_count+=1

    #Writing the new fasta to tmp folder
    SeqIO.write(new_fasta, "size_selected.fasta", "fasta") #path to the output file
    os.rename("tmp.aln", "no_size_selection.fasta")
    shutil.copy('size_selected.fasta', 'tmp.aln')
    #os.rename("size_selected.fasta", "tmp.aln")

    #print("optimal missing data threshold = "+str(optimal_threshold)) #removed v2.3.2
    print("#\tkept:\t\t"+str(frag_count))
    print("#\tomitted:\t "+str(len(input_sequences)-frag_count))
    print("#\tfragment size:\t "+str(selected_size)+"\n")
    #print("start position = "+str(start_position)) #removed v2.3.2
    #print("end position = "+str(finish_position)) #removed v2.3.2

    num7=str(frag_count)
    with open("log.txt", "a+") as log:
        log.write("\t " + str(num7) + "\t after size selection" + "\n")
    log.close()

def extra_trimal_fiter():
    print("Misaligned sequences removal step with seqoverlap", args.trim_seqoverlap, "and resoverlap", args.trim_resoverlap)
    with open('trimal.aln', 'w') as trimal_aln, open('tmp.aln', 'r') as tmp_aln, open('trimal.log', 'w') as trimal_log:
        subprocess.call(["trimal" , "-in" , "tmp.aln", "-out", "trimal.aln", "-noallgaps", "-seqoverlap",args.trim_seqoverlap, "-resoverlap", args.trim_resoverlap], stdout=trimal_log, stderr=trimal_log)
    os.rename('trimal.aln' , 'tmp.aln')

    if os.stat('tmp.aln').st_size == 0:
        print("\nThe misaligned sequence removal step eliminated all sequences - use different sequoverlap and resoverlap parameters!")
        exit()

    num8 = len([1 for line in open("tmp.aln") if line.startswith(">")])                                                     # count 8
    print("#        kept: ", num8)
        
    with open("log.txt", "a+") as log:
        log.write("\t " + str(num8) + "\t after the removal of misaligned sequences" + "\n")
    log.close()

def standard_trimal_filter():
    with open('trimal.aln', 'w') as trimal_aln, open('tmp.aln', 'r') as tmp_aln, open('trimal.log', 'w') as trimal_log:
        subprocess.call(["trimal" , "-in" , "tmp.aln", "-out", "trimal.aln", "-noallgaps"], stdout=trimal_log, stderr=trimal_log)
    os.rename('trimal.aln' , 'tmp.aln')

def tree_inverence():

    if (args.custom_nwk != ""):

        custom_tre = '../../' + args.custom_nwk
        shutil.copy(custom_tre, './custom.tre')
        custom_tre = './custom.tre'

        brlen = re.findall(":[0-9]*.[0-9]*", custom_tre)

        with open(custom_tre, 'r') as tre_file:

            if brlen != "":
                content = tre_file.read()
                tre_file = re.sub(":[0-9]*.[0-9]*", '', content )

            nwk_sp_list=tre_file.replace(',', ' ').replace('(', '').replace(')', '').replace(';', '').split()

        aln_sp_list = []
        for record in SeqIO.parse("tmp.aln", "fasta"):
            aln_sp_list.append(record.id)

        check = all(sp in aln_sp_list for sp in nwk_sp_list)

        if  check == False:
            print("custom nwk detected with more species in the tree than sequences - these tips will be pruned.")
            toprune_list = list(set(nwk_sp_list) - set(aln_sp_list))
            with open('toprune.lst', 'w') as f:
                for line in toprune_list:
                    f.write(line)
                    f.write('\n')

            tmp_rscript_pruning=[]
            tmp_rscript_pruning.append("library(phytools)")
            tmp_rscript_pruning.append("tree <- read.newick(\"custom.tre\")")
            tmp_rscript_pruning.append("spec <- read.table(\"toprune.lst\")")
            tmp_rscript_pruning.append("ptre<-drop.tip(tree,tree$tip.label[match(spec$V1, tree$tip.label)])")
            tmp_rscript_pruning.append("write.tree(ptre, file=\"custom.tre\")")

            with open('tmp_rscript_pruning.R', 'w') as tmp_rscript_pruning_file:
                for line in tmp_rscript_pruning:
                    tmp_rscript_pruning_file.write(line + '\n')

            with open('tmp_rscript_pruning.R_log', 'w') as tmp_rscript_pruning_log:
                subprocess.run(["Rscript" , "tmp_rscript_pruning.R"], stdout=tmp_rscript_pruning_log, stderr=tmp_rscript_pruning_log)

        elif len(nwk_sp_list) == len(aln_sp_list) and check == True:
            print("custom nwk detected with 1:1 species match to sequences - skipping topology inference.")
        elif len(nwk_sp_list) < len(aln_sp_list)  and check == True:
            print("custom nwk detected with fewer species than sequences - using constrained topology.")

##################################################################################################### if the custom .nwk is coding

        if (coding == True):

            tmp_partitions=[]
            filename = "tmp.aln"
            format = "fasta"
            tmp_aln = AlignIO.read("tmp.aln", "fasta")
            tmp_partitions.append("DNA, st = 1-%i\\3" % tmp_aln.get_alignment_length())
            tmp_partitions.append("DNA, nd = 2-%i\\3" % tmp_aln.get_alignment_length())
            tmp_partitions.append("DNA, rd = 3-%i\\3" % tmp_aln.get_alignment_length())

            with open('tmp.partitions', 'w') as tmp_partitions_file:
                for line in tmp_partitions:
                    tmp_partitions_file.write(line + '\n')

            subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.threads) , "-spp" , "tmp.partitions" , "-g" , "custom.tre", "-fast"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            def_tre_file = args.out + ".nwk"
            os.rename('tmp.partitions.treefile' , def_tre_file)
            def_aln_file = args.out + ".aln"
            os.rename('tmp.aln' , def_aln_file)

##################################################################################################### if the custom .nwk is not coding

        if (coding == False):

            with open('tmp.nwk', 'w') as tmp_nwk :
                subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.threads) , "-g" , custom_tre, "-fast"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            def_tre_file = args.out + ".nwk"
            os.rename('tmp.aln.treefile' , def_tre_file)

            def_aln_file = args.out + ".aln"
            os.rename('tmp.aln' , def_aln_file)

##################################################################################################### inferring a coding marker tree

    elif (args.custom_nwk == ""):

        print("\ninferring tree\n")

        if (coding == True):

            tmp_partitions=[]
            filename = "tmp.aln"
            format = "fasta"
            tmp_aln = AlignIO.read("tmp.aln", "fasta")
            tmp_partitions.append("DNA, st = 1-%i\\3" % tmp_aln.get_alignment_length())
            tmp_partitions.append("DNA, nd = 2-%i\\3" % tmp_aln.get_alignment_length())
            tmp_partitions.append("DNA, rd = 3-%i\\3" % tmp_aln.get_alignment_length())

            with open('tmp.partitions', 'w') as tmp_partitions_file:
                for line in tmp_partitions:
                    tmp_partitions_file.write(line + '\n')

            with open('tmp.nwk', 'w') as tmp_nwk :
                subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.threads) , "-spp" , "tmp.partitions", "-fast"], stdout=tmp_nwk, stderr=subprocess.DEVNULL)

            def_tre_file = args.out + ".nwk"
            os.rename('tmp.partitions.treefile' , def_tre_file)

            def_aln_file = args.out + ".aln"
            os.rename('tmp.aln' , def_aln_file)

##################################################################################################### inferring a not coding marker tree

        if (coding == False):

            with open('tmp.nwk', 'w') as tmp_nwk :
                subprocess.run(["iqtree" , "-s" , "tmp.aln" , "-nt" , str(args.threads), "-fast"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            def_tre_file = args.out + ".nwk"
            os.rename('tmp.aln.treefile' , def_tre_file)

            def_aln_file = args.out + ".aln"
            os.rename('tmp.aln' , def_aln_file)

##################################################################################################### ancestral sequences inference

    if (coding == True):
        print("inferring ancestral sequences using genetic code" , args.code , "\n")
    if (coding == False):
        print("inferring ancestral sequences\n")

    tmp_ctl=[]

    aln_line="seqfile = " + def_aln_file
    tmp_ctl.append(aln_line)
    tmp_ctl.append("outfile = tmp_baseml.out")
    tre_line="treefile = " + def_tre_file
    tmp_ctl.append(tre_line)
    tmp_ctl.append("noisy = 3")
    tmp_ctl.append("verbose = 1")
    tmp_ctl.append("runmode = 0")
    tmp_ctl.append("model = 7")
    tmp_ctl.append("Mgene = 0")
    tmp_ctl.append("clock = 0")
    tmp_ctl.append("fix_kappa = 0")
    tmp_ctl.append("kappa = 2.5")
    tmp_ctl.append("fix_alpha = 1")
    tmp_ctl.append("alpha = 0.")
    tmp_ctl.append("Malpha = 0")
    tmp_ctl.append("ncatG = 5")
    tmp_ctl.append("fix_rho = 1")
    tmp_ctl.append("rho = 0.")
    tmp_ctl.append("nparK = 0")
    tmp_ctl.append("nhomo = 0")
    tmp_ctl.append("getSE = 0")
    tmp_ctl.append("RateAncestor = 1")
    tmp_ctl.append("Small_Diff = 1e-6")
    tmp_ctl.append("cleandata = 0")
    if coding == True :
        icode = int(args.code) - 1
    #       icode_line="icode = " + str(icode)
        icode_line="icode = 11"
        tmp_ctl.append(icode_line)
    tmp_ctl.append("fix_blength = 2")
    tmp_ctl.append("method = 0")

    with open('tmp.ctl', 'w') as tmp_ctl_file:
        for line in tmp_ctl:
            tmp_ctl_file.write(line + '\n')

    with open('tmp.nwk', 'w') as tmp_nwk, open('tmp_baseml.log', 'w') as tmp_baseml_log:
        subprocess.run(["baseml" , "tmp.ctl"], stdout=tmp_baseml_log, stderr=subprocess.DEVNULL)

        tmp_baseml_log = './tmp_baseml.log'
        baseml_error = re.findall("stop codon TAA", tmp_baseml_log)
        if(len(baseml_error) >= 1):
            print("\n A VERY RARE EVENT HAPPENED! baseml crashed because it inferred an ancestral stop codon ... quite funny right? \n")
            quit()

##################################################################################################### extract baseml trees

    with open("rst", 'r') as rst, open("brlen.phylo.txt", 'w') as brlen_phylo, open("nodes.phylo.txt", 'w') as nodes_phylo:
        brlen = [10]
        nodes = [17]
        for position, line in enumerate(rst):
            if position in brlen:
                brlen_phylo.write(line.replace(" ", ""))
            if position in nodes:
                nodes_phylo.write(line.replace(" ", ""))

    return def_tre_file, def_aln_file

def bait_finder():
    print("finding baits of", args.baitlength, "nt with max distance of", args.distance,"% and tiling",args.tiling)

    input_seq = list(SeqIO.parse(args.out + ".aln", "fasta"))

    max_n_miss = int(args.distance*args.baitlength)
    anc_seq=[]
    fasta_dict = {}
    for seq in input_seq:
        fasta_dict[str(seq.name)]=str(seq.seq)

    file = open("rst")
    patern=r" {3,}"
    node_seq_dict = {}
    c=0
    for row in file :
        if "List of extant and reconstructed sequences" in row:
            c=c+1
        if(c==2):
            if "node #" in row :

                node_name=re.split(patern, row)[0]
                node_seq=re.sub("\n", "", re.sub(" ", "", re.split(patern, row)[1]))
                node_seq_dict[str(node_name)] = str(node_seq)
                anc_seq.append(">"+node_name+"\n"+node_seq+"\n")
    node_seq_dict.update(fasta_dict)

    fasta_file = open("Anc_seqs.fasta", "w")
    fasta_file.writelines(anc_seq)
    fasta_file.close()
    
    t=Tree("nodes.phylo.txt",1)
    node_content = t.get_cached_content(store_attr='name')

    pos=0
    i=1
    seq_len = len(seq)
    use_node="NO"
    bait_fasta=[]
    summary_stats=["Pos\tMean\tmin\tmax\tuniq\tgc\tmin_gc\tmax_gc\n"]
    progress=tqdm(total=None, desc="Scanning tile ")
    while (pos+args.baitlength) < seq_len:
        tips_included = []
        miss_list_sum_bait = []
        baitmer_list = []
        tip_no_node = []
        used_nodes_it=[]
        c_baits=0
        for node in t.traverse():
            if not node.is_leaf():
                node_name ="node #"+str(node.name)
                level="NODE"
            else:
                node_name =re.sub("^\d*_", "", str(node.name))
                node_name_tip =str(node.name)
                level="TIP"

            anc_str= node_seq_dict[node_name][pos:(args.baitlength+pos)]

            clade=node_content[node]
            miss_list=[]
            gap_counter_all=[]
            for seq in list(clade):
                #print("sequence = " + seq)
                gap_counter=0
                tip_seq=fasta_dict[re.sub("^\d*_", "", seq)][pos:(args.baitlength+pos)]
                while (tip_seq.count("-") != 0) and (len(re.sub("-", "",tip_seq))!=80):
                    gap_counter = gap_counter+(80-len(re.sub("-", "",tip_seq)))
                    tip_seq=fasta_dict[re.sub("^\d*_", "", seq)][pos:(args.baitlength+pos+gap_counter)]
                    if (args.baitlength+pos+gap_counter) >= seq_len:
                        break
                miss=0
                gap_counter_all.append(gap_counter)

                for nt in range(len(tip_seq)):
                    anc_str = node_seq_dict[node_name][pos:(args.baitlength+pos+gap_counter)]
                    if tip_seq[nt].upper()!=anc_str[nt].upper() and tip_seq[nt] != "-" :
                        miss=miss+1
                miss_list.append(miss)

            if(max(miss_list)<=max_n_miss) and len(Counter(gap_counter_all).keys())==1 :
                if level == "NODE":
                    clade_tips=node.get_leaf_names()
                elif level == "TIP":
                    clade_tips=[node_name_tip]
                for tip in clade_tips:
                    if tip in tips_included :
                        continue
                    else:
                        tips_included.append(str(tip))
                    if level == "NODE":
                        use_node="YES"
                    else:
                        tip_no_node.append(str(tip_seq))

            if(use_node=="YES"):
                if gap_counter != 0:
                    Anc_gap_correct=str()
                    anc_str=node_seq_dict[node_name][pos:(args.baitlength+pos+gap_counter)]
                    for nuc in range(len(tip_seq)):
                        if(tip_seq[nuc]!="-"):
                            Anc_gap_correct=Anc_gap_correct+anc_str[nuc]
                    anc_str=Anc_gap_correct

                bait_fasta.append(">"+node_name+ " "+ str(pos)+":"+str(pos+args.baitlength) +"\n" +anc_str+"\n")
                baitmer_list.append(anc_str)
                use_node="NO"
                miss_list_sum_bait=miss_list_sum_bait+miss_list
                c_baits=c_baits+1
                used_nodes_it.append(node_name)

        if len(tip_no_node)!=0:
            for seq in tip_no_node:
                miss_list=[]
                miss_dic={}
                for node in used_nodes_it:
                    miss=0
                    if len(seq)!=args.baitlength:
                        anc_str=node_seq_dict[node][pos:pos+len(seq)]
                        for nt in range(len(seq)):
                            if seq[nt].upper()!=anc_str[nt].upper() and seq[nt] != "-":
                                miss=miss+1
                    else:
                        anc_str=node_seq_dict[node][pos:pos+args.baitlength]
                        for nt in range(len(seq)):
                            if seq[nt].upper()!=anc_str[nt].upper():
                                 miss=miss+1

                    miss_list.append(miss)
                if len(miss_list)==0:
                    print("Error in node fiding, tip sequences only has gaps, check allingment")    
                    exit()
                if min(miss_list)<= max_n_miss :
                    continue
                #loop if the nodes already present are okay to catch the tips.

                miss_list=[]
                for key in node_seq_dict.keys():
                    if not "node" in key  :
                        continue
                    miss=0
                    if len(seq)!=args.baitlength:
                        anc_str=node_seq_dict[key][pos:(pos+len(seq))]
                        for nt in range(len(seq)):
                            if seq[nt].upper()!=anc_str[nt].upper() and seq[nt] != "-":
                                miss=miss+1
                    else:
                        anc_str=node_seq_dict[key][pos:pos+args.baitlength]
                        for nt in range(len(seq)):
                            if seq[nt].upper()!=anc_str[nt].upper():
                                
                                miss=miss+1	                        
                    miss_dic.update({key:miss})

            if(len(miss_dic)!=0):
                needed_key = [key for key in miss_dic.keys() if miss_dic[key] == min(miss_dic.values())]
                bait_fasta.append(">"+str(needed_key[0])+ " "+ str(pos)+":"+str(pos+args.baitlength) +"\n" +node_seq_dict[needed_key[0]][pos:(args.baitlength+pos)]+"\n")
                baitmer_list.append(anc_str)
                c_baits=c_baits+1
                used_nodes_it.append(needed_key)

        bait_gc=[]
        bait_gc_uniq=[]
        for bait in baitmer_list :
            if not bait in bait_gc_uniq:
                bait_gc.append((bait.upper().count("G")+bait.upper().count("C"))/args.baitlength)
                bait_gc_uniq.append(bait)


        summary_stats.append(str(pos)+":"+str(args.baitlength+pos)+"\t"+\
              str(mean(miss_list_sum_bait))+"\t"+ \
              str(min(miss_list_sum_bait))+"\t"+ \
              str(max(miss_list_sum_bait))+"\t"+\
              str(len(Counter(baitmer_list).keys()))+"\t"+\
              str(mean(bait_gc))+"\t"+\
              str(min(bait_gc))+"\t"+\
              str(max(bait_gc))+"\n")
        pos=pos+int(args.tiling)
        progress.update()
    progress.close()    
    fasta_file = open("bait_tmp.fasta", "w")
    fasta_file.writelines(bait_fasta)
    fasta_file.close()
    
    summary_stats_file = open("node_selection_stats.tsv", "w")
    summary_stats_file.writelines(summary_stats)
    summary_stats_file.close()

def bait_gc_complexity_flagger():
   ## do GC filteing
    baits=[]
    for record in (SeqIO.parse("bait_tmp.fasta", "fasta")):
        bait=record.seq.upper()
        gc=((bait.count("G")+bait.count("C"))/args.baitlength)
        if(gc<0.30):
            flag="GC="+str(gc)+"; "+"GC_flag=low;"
        elif(gc>0.70):
            flag="GC="+str(gc)+"; "+"GC_flag=high;"
        else:
            flag="GC="+str(gc)+"; "+"GC_flag=good;"

        baits.append(SeqRecord(Seq(str(bait)), id=record.id, description=record.description+" "+flag))
    SeqIO.write(baits, "tmp_f.fasta", "fasta")

    subprocess.run(["dustmasker", "-in", "tmp_f.fasta", "-out", "tmp.dust.fasta" ,"-outfmt", "fasta", "-window", str(args.baitlength), "-level", "15"], stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

    flagged_baits=[]
    file=open("tmp.dust.fasta")
    for record in (SeqIO.parse(file, "fasta")):
        seq=record.seq
        if not str(seq).islower() and not str(seq).isupper():
            count=0
            for nucl in seq:
                if(nucl.islower()):
                    count=count+1
            if(count>=int(args.dust_threshold)):
                record.description=record.description+"; dust_flag=failed; low_complex_lenght="+str(count)+";"
            elif(count<int(args.dust_threshold)):
                record.description=record.description+"; dust_flag=pass; low_complex_lenght="+str(count)+";"
        else:
            record.description=record.description+"; dust_flag=pass;"
        flagged_baits.append(record)

    SeqIO.write(flagged_baits, "bait.fasta", "fasta")

def summary_info_r_script():
    tmp_rscript=[]
    tmp_rscript.append("library(ggplot2)")
    tmp_rscript.append("library(cowplot)")
    tmp_rscript.append("library(ggtree)")
    tmp_rscript.append("library(treeio)")
    tmp_rscript.append("library(Biostrings)")
    tmp_rscript.append("library(data.table)")
    tmp_rscript.append("library(stringr)")
    tmp_rscript.append("")
    tmp_rscript.append("fasta <- readDNAStringSet(\""+ "bait_tmp.fasta" +"\")") # read the alignment file
    tmp_rscript.append("node_sum<- data.table()")
    tmp_rscript.append("node_sum<- node_sum[,.(label=gsub(\" [0-9]{1,}:[0-9]{1,}.*\", \"\", fasta@ranges@NAMES),")
    tmp_rscript.append("                       position=str_split_fixed(gsub(\"node #\", \"node_#\", fasta@ranges@NAMES), \" \", 3)[,2])]")
    tmp_rscript.append("")
    tmp_rscript.append("node_sum.dt<- node_sum[,.N,by=(label)]")
    tmp_rscript.append("node_sum.dt[, label:=(gsub(\"^[a-z]{1,} #\", \"\", label))]")
    tmp_rscript.append("")
    tmp_rscript.append("tree <- read.tree(\"" + "nodes.phylo.txt" + "\")") # read the tree with node lavels nodes.phylo.txt
    tmp_rscript.append("tree.df <-  as_tibble(tree)")
    tmp_rscript.append("tree.df$label <- gsub(\"[0-9]{1,}_\", \"\", tree.df$label)")
    tmp_rscript.append("")
    tmp_rscript.append("brlen.df <- as_tibble(read.tree(\"brlen.phylo.txt\"))")
    tmp_rscript.append("tree.df$branch.length <- brlen.df$branch.length")
    tmp_rscript.append("")
    tmp_rscript.append("for (node in node_sum.dt$label){")
    tmp_rscript.append("tree.df$N[tree.df_original$node==node] <- node_sum.dt$N[node_sum.dt$label==node]")
    tmp_rscript.append("}")
    tmp_rscript.append("")
    tmp_rscript.append("if(sum(list.files()==\"tmp.bold\")==1){")
    tmp_rscript.append("")
    tmp_rscript.append("bold<- read.delim(\"tmp.bold\")")
    tmp_rscript.append("for (id in tree.df$label){")
    tmp_rscript.append("  tree.df$label[tree.df$label==id] <- paste(id, \" species = \", bold$species_name[bold$processid==id][1])")
    tmp_rscript.append("}")
    tmp_rscript.append("}")
    tmp_rscript.append("")
    tmp_rscript.append("write.tree(as.phylo(tree.df), \"node_tree_branch_lenght.nwk\")")
    tmp_rscript.append("tree_plot <- as.treedata(tree.df)")
    tmp_rscript.append("e <- ggtree(tree_plot, layout=\"fan\")+")
    tmp_rscript.append("  geom_tippoint(aes(size=N), col=\"red\")+")
    tmp_rscript.append("  geom_tiplab(align=T, size=(4/log(length(tree_plot@phylo$tip.label))+1))+")
    tmp_rscript.append("  geom_nodepoint(aes(size=N),col=\"black\")+ggtitle(\"Tree with node counts\")")
    tmp_rscript.append("  #geom_nodelab(nudge_x = 0.1, nudge_y = 0.1, colour=\"darkgreen\")")
    tmp_rscript.append("")
    tmp_rscript.append("n.df <- read.delim(\""+ "node_selection_stats.tsv" + "\")") # read summary table from node selection .py
    tmp_rscript.append("    ")
    tmp_rscript.append("plot_theme <- theme(panel.background = element_rect(fill = \"white\", colour = \"black\"),")
    tmp_rscript.append("                    strip.background = element_blank(),")
    tmp_rscript.append("                    panel.grid.minor = element_blank(),")
    tmp_rscript.append("                    panel.grid.major = element_blank(),")
    tmp_rscript.append("                    strip.text = element_text(face = \"bold\", size = 11, angle = 0, hjust = 0.5),")
    tmp_rscript.append("                    strip.background.x = element_blank(),")
    tmp_rscript.append("                    text = element_text(size = 12),")
    tmp_rscript.append("                    axis.text.x = element_text(size = , angle = 0),")
    tmp_rscript.append("                    axis.title = element_text(size=12),")
    tmp_rscript.append("                    axis.text.y = element_text(size=12),")
    tmp_rscript.append("                    plot.margin = margin(0.5, 0.5, 0.5, 0.5, \"cm\"),")
    tmp_rscript.append("                    legend.position = \"None\",")
    tmp_rscript.append("                    legend.key = element_blank(),")
    tmp_rscript.append("                    legend.text = element_text(size=12),")
    tmp_rscript.append("                    plot.title = element_text(hjust = 0.5))")
    tmp_rscript.append("")
    tmp_rscript.append("n.df$Pos <- factor(n.df$Pos, levels = n.df$Pos)")
    tmp_rscript.append("")
    tmp_rscript.append("a <- ggplot()+plot_theme+ggtitle(\"Range and mean number missmatches between baits and tip\")+")
    tmp_rscript.append("  theme(axis.text.x = element_blank())+")
    tmp_rscript.append("        #plot.margin = margin(t=0.5, r=0.5, b=-1, l=0.5, \"cm\"))+")
    tmp_rscript.append("  geom_errorbar(data=n.df, aes(Pos, ymin=min, ymax=max), alpha=0.3)+")
    tmp_rscript.append("  geom_point(data=n.df, aes(Pos, Mean), size=3)+")
    tmp_rscript.append("  scale_y_continuous(limits = c(0,7))+")
    tmp_rscript.append("  xlab(\"\")+ylab(\"Mean distance bait (in nt.)\")")
    tmp_rscript.append("  ")
    tmp_rscript.append("b <- ggplot()+plot_theme+ggtitle(\"Number of uniq baits\")+")
    tmp_rscript.append("  theme(axis.text.x = element_blank())+")
    tmp_rscript.append("        #plot.margin = margin(t=-1, r=0.5, b=-1, l=0.5, \"cm\"))+")
    tmp_rscript.append("  geom_point(data=n.df, aes(Pos, uniq), size=3)+")
    tmp_rscript.append("  xlab(\"\")+")
    tmp_rscript.append("  ylab(\"Number of uniq baits\")")
    tmp_rscript.append("")
    tmp_rscript.append("c <- ggplot()+plot_theme+ggtitle(\"Range and mean GC-content baits\")+")
    tmp_rscript.append("  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))+")
    tmp_rscript.append("        #plot.margin = margin(t=-1, r=0.5, b=0.5, l=0.5, \"cm\"))+")
    tmp_rscript.append("  geom_errorbar(data=n.df, aes(Pos, ymin=min_gc, ymax=max_gc), alpha=0.3)+")
    tmp_rscript.append("  geom_point(data=n.df, aes(Pos, gc), size=3)+")
    tmp_rscript.append("  scale_y_continuous(limits = c(0,1))+")
    tmp_rscript.append("  geom_hline(yintercept = c(0.75, 0.30), col=\"red\", linetype=\"dashed\")+")
    tmp_rscript.append("  ylab(\"Mean GC\")+xlab(\"Bait position in the selected fragment\")")
    tmp_rscript.append("")
    tmp_rscript.append("pdf(\"summary_stats.pdf\", width = 8.27, height = 11.69)")
    tmp_rscript.append("plot_grid(a,b,c, ncol = 1, align = \"vh\", axis = \"rltb\")")
    tmp_rscript.append("plot_grid(e)")
    tmp_rscript.append("dev.off()")
    tmp_rscript.append("svg(\"node_tree.svg\",  height = 20, width = 20)")
    tmp_rscript.append("plot(e)")
    tmp_rscript.append("dev.off()")


    with open('tmp.R', 'w') as tmp_rscript_file:
        for line in tmp_rscript:
            tmp_rscript_file.write(line + '\n')
    with open('tmp.R', 'w') as tmp_rscript_file:
        for line in tmp_rscript:
            tmp_rscript_file.write(line + '\n')

    with open('tmp.R_log', 'w') as tmp_rscript_log:
        subprocess.run(["Rscript" , "tmp.R"], stdout=tmp_rscript_log, stderr=subprocess.DEVNULL)

def fast_mode():
    print("\n\n###finding baits with the fast approach###")

    rc = subprocess.call(['which', 'dustmasker'],stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    baits=[]
    baitsize=int(args.baitlength)

    run="TRUE"
    if(run=="TRUE"):
        for record in (SeqIO.parse("tmp.aln", "fasta")):
            seq=record.seq.ungap("-")
            for t in range(0, len(seq), int(args.tiling)):
                if(len(seq[t:(t+baitsize)])!=baitsize):
                    break
                    print(done)
                bait=seq[t:(t+baitsize)].upper()
                gc=((bait.count("G")+bait.count("C"))/baitsize)
                if(gc<0.30):
                    flag="GC="+str(gc)+"; "+"GC_flag=low;"
                elif(gc>0.70):
                    flag="GC="+str(gc)+"; "+"GC_flag=high;"
                else:
                    flag="GC="+str(gc)+"; "+"GC_flag=good;"

                baits.append(SeqRecord(Seq(str(bait)), id=record.id+"_"+str(t)+str(t+baitsize), description=flag))

        SeqIO.write(baits, "tmp_f.fasta", "fasta")

    subprocess.run(["dustmasker", "-in", "tmp_f.fasta", "-out", "tmp.dust.fasta" ,"-outfmt", "fasta", "-window", str(baitsize), "-level", "15"], stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

    flagged_baits=[]
    file=open("tmp.dust.fasta")
    for record in (SeqIO.parse(file, "fasta")):
        seq=record.seq
        if not str(seq).islower() and not str(seq).isupper():
            count=0
            for nucl in seq:
                if(nucl.islower()):
                    count=count+1
            if(count>=int(args.dust_threshold)):
                record.description=record.description+"; dust_flag=failed; low_complex_lenght="+str(count)+";"
            elif(count<int(args.dust_threshold)):
                record.description=record.description+"; dust_flag=pass; low_complex_lenght="+str(count)+";"
        else:
            record.description=record.description+"; dust_flag=pass;"
        flagged_baits.append(record)
    SeqIO.write(flagged_baits, "bait.fasta", "fasta")

def bait_gc_complexity_filter():
    filtered_baits = []

    for record in (SeqIO.parse("bait.fasta", "fasta")):
        flags = record.description.split(" ")
        if "GC_flag=good;" and "dust_flag=pass;" in flags:
            filtered_baits.append(record)

    SeqIO.write(filtered_baits, "filtered_baits.fasta", "fasta")

    subprocess.run(["cd-hit","-i","filtered_baits.fasta","-o","filtere_uniq_baits.fasta","-c", "1"], stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

def final_clustering():

    subprocess.run(["cd-hit","-i","filtere_uniq_baits.fasta","-o","final_clusterd_baits.fasta","-c", str(args.bts_collapse)], stdout=subprocess.DEVNULL , stderr=subprocess.DEVNULL)

def write_log_file_02():

    with open("log.txt", "a+") as log:
        log.write("\n\n\t BAITS:" + "\n\n")

        baits = len([1 for line in open("bait.fasta") if line.startswith(">")])
        log.write("\t " + str(baits) + "\t baits"+ "\n")
        baits_filterd = len([1 for line in open("filtered_baits.fasta") if line.startswith(">")])
        log.write("\t " + str(baits_filterd) + "\t baits passing gc+complexity filter"+ "\n")
        baits_colapse = len([1 for line in open("filtere_uniq_baits.fasta") if line.startswith(">")])
        log.write("\t " + str(baits_colapse) + "\t unique baits"+ "\n\n") 
        baits_colapse = len([1 for line in open("final_clusterd_baits.fasta") if line.startswith(">")])
        log.write("\t " + str(baits_colapse) + "\t unique clustered baits"+ "\n\n") 
        log.write("# analysis finished on " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
        print("\n\nGotcha run successfully\n\na totall of " + str(baits_colapse)+ " baits are needed to capture the selected taxa\n")

def clean_up():
    #def_prb_file = args.out + "raw_.fasta"
    os.rename('bait.fasta' , args.out + "_raw_baits.fasta")
    os.rename('filtered_baits.fasta' , args.out + "_filtered_baits.fasta")
    os.rename('filtere_uniq_baits.fasta' , args.out + "_filtered_uniq_baits.fasta")
    os.rename('final_clusterd_baits.fasta' , args.out + "_filtered_uniq_clustered_baits.fasta")

    #shutil.copy(def_prb_file, '..')
    shutil.copy("log.txt", '..')
    shutil.copy(args.out + "_raw_baits.fasta", '..')
    shutil.copy(args.out + "_filtered_baits.fasta", '..')
    shutil.copy(args.out + "_filtered_uniq_baits.fasta", '..')


    if args.fast != True :
        shutil.copy('summary_stats.pdf', '..')
        #shutil.copy(def_tre_file, '..')
        shutil.copy("node_tree_branch_lenght.nwk", '..')
        os.rename("node_tree_branch_lenght.nwk", args.out + "_node_tree_branch_lenght.nwk")
        shutil.copy(def_aln_file, '..')
        shutil.copy("node_tree.svg", '..')

    os.chdir('..')

    if args.verbose == False :
        pass
    else :
        shutil.rmtree("tmp/", ignore_errors=True) # this is not sovling the problem with folder not being removed
        #files = glob.glob("tmp/*")
        #for f in files:
        #       os.remove(f)
        #os.rmdir("tmp")

    print("\nanalysis finished on" , datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

if __name__ == '__main__':
    print("Run Gotcha Main")
    coding, marker_list = check_if_all_is_installed()

    #####################################################################
    #                         starts the log file                       #
    #####################################################################
    
    write_log_file_01()

    #####################################################################


    #####################################################################
    #  downloads sequences from bold if no costume alingent or fasta is given
    #####################################################################
    
    if (args.custom_fna == "") and (args.custom_aln == ""):
         bold_download()

    #####################################################################
    #                 set of quality filters fo any input               #
    #####################################################################
  
    seq_quality_filters()

    #####################################################################
    #            size selection of the reference sequences              #
    #####################################################################
   
    size_selection()
    
    #####################################################################


    #####################################################################
    # an aditonal optional timal filter, otherwise standard filer is applied
    #####################################################################
   
    if args.trim_seqoverlap != "" and args.trim_seqoverlap != "" :
        extra_trimal_fiter()
    else:
        standard_trimal_filter()
    

    #####################################################################
    #           stop the gotch for a manual quality check               #
    #####################################################################
  
    if (args.manual_check == True):
        shutil.copy('tmp.aln', '..')
        print("\nstopping prior to tree inference for a manual check of the alignment\n")
        quit()


    #####################################################################
    #               Phylogentic backbone bait finding                   #
    #####################################################################
    
    if (args.fast == False):
        def_tre_file, def_aln_file = tree_inverence()
        bait_finder()
        bait_gc_complexity_flagger()
        summary_info_r_script()


    #####################################################################
    #                           fast mode                               #
    #####################################################################

    if (args.fast == True):
        fast_mode()

    #####################################################################
    
    
    #####################################################################
                        # final filtering and clean up
    #####################################################################
    bait_gc_complexity_filter()
    final_clustering()
    write_log_file_02()
    clean_up()
    #####################################################################