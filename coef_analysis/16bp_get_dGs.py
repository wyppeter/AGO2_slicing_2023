# Generate bp-ing energies for segments of miRNA
# Bartel Lab
# Peter Y Wang 2023

import sys
import subprocess
from numpy import mean, median

WC = {"G":"C", "C":"G",
      "A":"U", "U":"A"}
def revcomp(s):
    """ Reverse complement RNA """
    return "".join([WC[nt] for nt in s[::-1]])
def numList2Str(l):
    """ List of floats to list of strs for export """
    return [str(n) for n in l]
def numList2StrRound(l, dp):
    """ List of floats to list of strs for export, WITH ROUNDING """
    return [str(round(n, dp)) for n in l]

#### Thermodynamic param + purine content parsing for miRNA domains
divNames = [
    "full", #0
    "seed", #1
    "cent", #2
    "supp", #3
    "tail"  #4
    ]
def miRDiv(seq):
    """ Extract specific miR seq domains/regions """
    return {
        divNames[0]:seq[1:],    #2  - end
        divNames[1]:seq[1:8],   #2  - 8
        divNames[2]:seq[8:12],  #9  - 12
        divNames[3]:seq[12:16], #13 - 16
        divNames[4]:seq[16:]    #17 - end
        }
def evalInput(seq):
    """ Get params for RNAeval """
    constr  = "("*len(seq) + "&" + ")"*len(seq)
    combSeq =         seq  + "&" + revcomp(seq)
    return combSeq + "\n" + constr
def evalCall(callInput, tempParam = []):
    """ Wrapper to call RNAeval """
    # Outputs the dG energy value
    foldOutput = subprocess.run(
            "RNAeval",
            input = callInput,
            capture_output = True, text = True
            ).stdout
    return float(foldOutput.split("\n")[1].split()[-1].strip("()"))

def DomParse(seq):
    """ Call getDomParams on each region defined in miRDiv """
    domOut = []
    div = miRDiv(seq)
    for dom in divNames:
        regseq = div[dom]
        dG = evalCall(evalInput(regseq)) - 4.09  # correct for initiation penalty
        if regseq == revcomp(regseq):
            # Need to correct for symm penalty given
            dG = dG - 0.43
        domOut.append(dG)
    return domOut  # order preserved at Python v3.6+


#### MAIN ####
# Load miRNA sequences
with open("miR-seqs.csv", "r") as infile:
    miRNAs = {}
    infile.readline()
    for l in infile:
        info = l.strip().split(",")
        miRNAs[info[0]] = info[1]

# Prep outfile and write header
outfile = open("16bp_analysis-dGvals.tsv", "w")
header = ["miR","seq"] + \
    ["%s_dG" % dID for dID in divNames]
outfile.write("\t".join(header) + "\n")

# Parse and obtain thermodynamics data
for miR in miRNAs.keys():
    miRseq = miRNAs[miR]
    print("Processing %s..." % miR)

    dGout = numList2StrRound(DomParse(miRseq), 2)

    # Ready for output
    outfile.write("\t".join([miR, miRseq] + dGout) + "\n")

outfile.close()
