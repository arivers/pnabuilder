#/env/bin/python
# Scripts to test design the optimal PNA sequence from subsequence
# Adam R. Rivers, DOE Joint Genome Institute/ Lawrence Berkeley National Laboratory
# September 10, 2016

import argparse
import primer3
import os
import subprocess
import re
import csv
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
from itertools import product


#Functions

def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   d = IUPAC.IUPACData.ambiguous_dna_values
   return [ list(map("".join, product(*map(d.get, seq)))) ]


def tmok(probe,hi,lo):
    """check the tm of the primer returns true if within range false otherwise"""
    probeseq = Seq(str(probe), IUPAC.ambiguous_dna) # convert text to seqobjact
    probelist = extend_ambiguous_dna(probeseq) #  create a list of all combinations
    tm = [] 
    for i in probelist[0]: # for each ambiguous primer:
        tm.append(primer3.calcTm(str(i))) # calculate the tm
    avetm =  float(sum(tm))/float(len(tm)) # average the tm's
    if (avetm < float(hi)) and (avetm > float(lo)): #check they are within the target temp
        return True
    else:
        return False
  
# TODO: the Tm's for regular primers and PNA's are different. 
# Update the function to access the PNA calculator on the PNABio site instead of using primer3-py      
def pnabiotemp(probe):
	"""Connect to PNABIO and use their calculator for Tm"""
	pass
	       

def subprobes(probe,probeminlen,hi,lo):
    """Scans through all subsequences grater than the minimum length in a sequence appending subsequences meeting tm criteria to a list"""
    subp = []
    length = len(probe) # check probe length
    for i in range(int(probeminlen),length): #for each set of substrings:
        for j in range(0,(length-i)): # and for each substring in a set:
            ss = probe[j:j+i+1] #save the sub-sequnce to a list
            if tmok(ss,hi,lo): # test if the average temp of the probe is in temp range
                subp.append(ss)
    return subp # return list of passing subsequences

def find_amped(usearch, path, primem, primerfile, out):
    """call Usearch to find sequences that amplifiy with primers"""
    params = [usearch, "-search_pcr", path, "-db", primerfile, "-strand", "both", "-maxdiffs", primem, "-minamp", "300", "-maxamp", "1000", "-ampout", out]
    p0 = subprocess.Popen(params, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p0.communicate()

def split_euks(fin, tmp):
    """write fasta files containing Eukaryotic and prokaryotic reads, return a list of [total seqs, euk seqs]"""
    with open(fin, 'r') as f:
        with open(os.path.join(tmp ,"ssu.euks.fasta"), 'w') as euk_file:
            with open(os.path.join(tmp ,"ssu.proks.fasta"), 'w') as prok_file:
                sequences = SeqIO.parse(f,'fasta')
                euks = [] # list of euk sequence ids
                proks = [] # list of prok sequence ids
                i = 0 # count of all sequences
                e = 0 # count of euks
                for record in sequences:
                    i += 1
                    #print(record.description)
                    if re.search( "Eukaryota",record.description): # if euk is in the description:
                        e += 1  # increment euk counter 
                        euks.append(record) #  add to euk list
                    else:
                        proks.append(record) # add to prok list
            	SeqIO.write(euks, euk_file, "fasta")
            	SeqIO.write(proks, prok_file, "fasta")
    return [i,e]
        
def find_probehits(usearch, inpath, outpath, probem, probes ):
    """Use Usearch to find sequences that amplifiy with primers"""
    params = [usearch, "-search_oligodb", inpath, "-db", probes, "-strand", "both", "-maxdiffs", probem, "-userout", outpath, "-userfields", "target+query+mism"]
    p1 = subprocess.Popen(params, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p1.communicate()

def parse_probehits(report):
    probecnt = Counter()
    with open(report, 'r') as hits:
        for line in hits:
            ln = line.split()
            probe = ln[0]
            ssu = ln[1]
            mism = ln[2]
            probecnt[probe] += 1
    return probecnt.most_common()

def main():
    parser = argparse.ArgumentParser(description='A script to design an optimal PNA sequence to block 18S amplification')
    parser.add_argument('--sequence', help="a long target probe sequence for testing", default='ACTTTCGTTCTTGATYRA')
    parser.add_argument('--ssu', help="SSU sequences of interest", default='/global/dna/shared/databases/silva/123/Exports/SILVA_123_SSURef_Nr99_tax_silva_trunc.fasta.gz')
    parser.add_argument('--usearch', help="usearh path", default='/global/dna/projectdirs/MEP/tools/bin/usearch64')
    parser.add_argument('--primers', help="path to primers in fasta fromat", default="primers.fa")
    parser.add_argument('--tmhi', help="highest tm",default='80')
    parser.add_argument('--tmlo', help="lowest tm",default='40')
    parser.add_argument('--primem', help="The gratest allowed primer mismatch", default = "2")
    parser.add_argument('--probem', help="The gratest allowed primer mismatch", default = "2")
    parser.add_argument('--probeminlen', help="probe minimum length", default = "14")
    parser.add_argument('--tmp', help="root directory to write temp files in", default = ".")
    parser.add_argument('--outeuks', help="the euks output report file", type=argparse.FileType('w'), default = 'euks.out.txt')
    parser.add_argument('--outproks', help="the proks output report file", type=argparse.FileType('w'), default = 'proks.out.txt')
    args = parser.parse_args()
    
    # Identify all subprobes with correct mt
    print("Identifing all subprobes with an accptable Tm.") 
    subprobelist = subprobes(probe=args.sequence,probeminlen=args.probeminlen,hi=args.tmhi,lo=args.tmlo)
    print("{} Subsequences were identified with an acceptable Tm.".format(len(subprobelist)))
    # Write probes to file
    print("Saving probes to file")
    with open(os.path.join(args.tmp , "probes.fa"), 'w') as probefile:
        n = 0
        for i in subprobelist:
            n += 1
            probefile.write(">" + str(n) + "\n")
            probefile.write(i + "\n")
    probefile.close()
    
    # split euk and prok sequences into two files, recording total sequences and euks copied
    print("Splitting the Eukaryotic and Prokaryotic sequences") 
    counts = split_euks(fin=args.ssu, tmp=args.tmp)
    print("Total sequences read: {}".format(counts[0]))
    print("Eukaryotic sequences: {}".format(counts[1]))
    
    # to ePCR to find Euks who have sequences amplified by the v4-v5 primers
    print("Running ePCR on eukaryotic sequences")
    find_amped(usearch= args.usearch, path = os.path.join(args.tmp ,"ssu.euks.fasta"), primem = args.primem, primerfile = args.primers, out =  os.path.join(args.tmp , "euks.amplicons.fasta"))
    
    # to ePCR to find Proks who have sequences amplified by the v4-v5 primers
    print("Running ePCR on prokaryotic sequences")
    find_amped(usearch= args.usearch, path = os.path.join(args.tmp ,"ssu.proks.fasta"), primem = args.primem, primerfile = args.primers, out =  os.path.join(args.tmp , "proks.amplicons.fasta"))
   
    # identify euks amplicons that the probes bind to
    print("Testing which probes bind to eukaryotic amplicons")
    find_probehits(usearch = args.usearch, inpath = os.path.join(args.tmp ,"euks.amplicons.fasta"), outpath = os.path.join(args.tmp , "euks.probehits.txt"), probem= args.probem, probes=os.path.join(args.tmp, "probes.fa"))
    
    # identify proks amplicons that the probes bind to
    print("Testing which probes bind to prokaryotic amplicons")
    find_probehits(usearch = args.usearch, inpath = os.path.join(args.tmp ,"proks.amplicons.fasta"), outpath = os.path.join(args.tmp , "proks.probehits.txt"), probem= args.probem, probes=os.path.join(args.tmp, "probes.fa"))
    
    # Parse the euk probe report
    if os.path.getsize(os.path.join(args.tmp , "euks.probehits.txt")) > 0:   
        print("Counting the probes that bound to eukayotic SSU")
        euklist = parse_probehits(report= os.path.join(args.tmp , "euks.probehits.txt"))
    else:
        print("The probe did not match any eukaryotic sequences")
    
    # Parse the prok probe report 
    if os.path.getsize(os.path.join(args.tmp , "proks.probehits.txt")) > 0:
        print("Counting the probes that bound to prokaryotic SSU")
        proklist = parse_probehits(report= os.path.join(args.tmp , "proks.probehits.txt"))
    else:
        print("The probe did not match any prokaryotic sequences")
    
    #print euk probe hits
    print("Creating euk output report")
    csv_out=csv.writer(args.outeuks)
    csv_out.writerow(['probe','number of euks removed'])
    for row in euklist:
        csv_out.writerow(row)
    try:
        print("Creating prok output report")
        csv_out=csv.writer(args.outproks)
        csv_out.writerow(['probe','number of proks removed'])
        for row in proklist:
            csv_out.writerow(row)
    except:
    	pass
    	
if __name__ == '__main__':
    main()