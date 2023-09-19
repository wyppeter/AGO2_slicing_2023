# Calculate the dG of structure at seed-match site of targets
# Bartel Lab
# Peter Y Wang 2023

import subprocess
import re

def RNAfold(rnaseq):
    """ Wrapper to call RNAfold, with dotbracket output """
    # Returns db and dG
    foldOutput = subprocess.run(
            ["RNAfold", "--noPS"],
            input = rnaseq,
            capture_output = True, text = True
            ).stdout
    out = foldOutput.split("\n")[1].split()
    return out[0], float(out[-1].strip("()"))

# Use 8mer as seed
t1 = re.compile("[AUGC][augc]")  # Target site is capitalized
seedSize = 8
pseudoC = 1e-20  # pseudocount
def fracFoldSeed(seq, db):
    """ How much of base-pairing is within the seed? Assuming no internal bp within seed itself """
    t1pos = t1.search(seq).start()
    dbseed = db[t1pos-seedSize+1:t1pos+1]
    foldseed = len(dbseed) - dbseed.count(".") + pseudoC
    foldfull = (len(seq) - db.count("."))/2 + pseudoC
    return foldseed/foldfull

# Main
with open("targ-seqs.csv", "r") as inF, open("targ-seqs+seedStruct.csv", "w") as outF:
    header = inF.readline()

    outF.write(",".join([
            "rxnID",
            "targID",
            "targ.seq",
            "db",
            "targ.dG",
            "seed.dG"
            ]) + "\n")

    for l in inF:
        items = l.rstrip("\n").split(",")
        rxnID  = items[0]
        targID = items[1]
        seq    = items[2]
        targdb, targdG = RNAfold(seq)
        seedfrac = fracFoldSeed(seq, targdb)

        outF.write(",".join([
            rxnID,
            targID,
            seq,
            targdb,
            "%.8f"%targdG,
            "%.8f"%(seedfrac*targdG)
            ]) + "\n")
