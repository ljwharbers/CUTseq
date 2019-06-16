"""

Get the cut site locations per enzyme

"""

from Bio import SeqIO
from Bio.Restriction import *
import re
import time

#define variables
enzyme = {"BfaI", "CviAII", "FatI", "NlaIII", "MseI", "AflII", "AseI", "BsrGI", "HindIII", "NcoI", "NdeI", "NsiI",
          "PciI", "PstI", "SacI", "SpeI", "SphI", "EcoRI"}

#specify reference location
reference = "E:\SciLife Lab\References\Homo_sapiens.GRCh37.75.dna.fa"

#create restrictionbatch class of the enzymes
rb = RestrictionBatch(enzyme)

start = time.time()
#loop through fasta sequences (chromosomes)
for record in SeqIO.parse(reference, "fasta"):
    sites = rb.search(record.seq)                           #search for cut sites of all the enzymes in the batch
    for x in sites:
        with open("E:\SciLife Lab\CUTseq\Cut_site_analysis\DistanceDistribution/cutsiteLocations_" + str(x) + ".bed",
                  "a") as output:
            for value in sites[x]:
                output.write("chr" + re.sub(" .*", "", record.description) + "\t" + str(value) + "\t" +
                             str(value + 1) + "\n")
    print("Chromosome " + re.sub(" .*", "", record.description) + ": Done")
end = time.time()
print("Time elapsed:", end - start)