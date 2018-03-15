# -*- encoding: utf-8 -*-
#Generating MultiFASTA using geneId instead of seqname
#Using a .gtf or .gff files, create a new file in which each line is a gene 

from Bio import SeqIO
import sys

############################# PARAMETERS #############################
num_args = 4

if len(sys.argv) != num_args:
    print("\nERROR calling the script. Number of arguments must be %d\nSyntax: python %s [progenitor.GFF] [progenitor.FASTA] [output.FASTA]\n" %(num_args-1, sys.argv[0]))
    sys.exit(1)

progenitor = sys.argv[1]
fastaScaffolds = sys.argv[2]
multiFasta = sys.argv[3]
#progenitor = "PleosPC9_1_GeneModels_Filteredmodels2.gff"
#fastaScaffolds = "PleosPC9_1_Assembly_scaffolds.fasta"
#multiFasta = "multiFastaGFF.fasta"

if (progenitor[-3:].upper()  != "GFF" and progenitor[-3:].upper()  != "GTF") or (fastaScaffolds[-5:].upper() != "FASTA"):
    print("\nERROR calling the script. Provide .GFF and .FASTA files\nSyntax: python %s [progenitor.GFF] [progenitor.FASTA] [output.FASTA]\n" %( sys.argv[0] ) )
    sys.exit(2)

############################# VARIABLES #############################

#Initiailzing values
currentScaffoldNumber = "" 
currentScaffoldName = ""
currentScaffoldIndex = ""
currentSeqId = ""
newSeqId = ""

geneStart = 99999999999999999999999999999
geneEnd = 0

#Do NOT change names like 'scaffold' 'geneIdGTF'. Change instead number if needed (or change as well names below)
headers = {"scaffold":0, "start":3, "end":4, "geneIdGTF":10, "geneIdGFF":11}
separator = "|"
fileWidth = 50
scaffoldName =  ""

############################# FUNCTIONS #############################

def updateGeneSection(newStart,  geneStart,  newEnd,  geneEnd):
    #All parameters (arguments) in the Python language are passed by reference
    if newStart < geneStart:
        geneStart = newStart
    if newEnd > geneEnd:
        geneEnd = newEnd
    return geneStart,  geneEnd
    
def writegene( SeqId, ScaffoldName, currentScaffoldIndex,  geneStart, geneEnd,  multiFasta,  fastaScaffolds,  separator):
    gene = ">" + SeqId + separator + ScaffoldName + separator + str( geneStart ) + separator + str( geneEnd )
    multiFasta.write( gene + "\n" )
    #Saving gene seqcuence
    gene = fastaScaffolds[  currentScaffoldIndex + geneStart : currentScaffoldIndex + geneEnd+1 ]
    #Print all data of current gene
    while len(gene) > fileWidth: 
        multiFasta.write( gene[ 0:fileWidth ] + "\n" )
        gene = gene[ fileWidth: ]
    #In case <geneCompleteSeqcuence> still has some data to print
    if len( gene ) > 0:
        multiFasta.write( gene + "\n" )
    return 
    
############################# PROGRAM #############################

if (progenitor[-3:].upper()  == "GFF"):
    scaffoldName = ">chr"
    geneIdSearch = 1
elif (progenitor[-3:].upper()  == "GTF"):
    scaffoldName = ">scaffold_"
    geneIdSearch = 2
progenitor = open(progenitor, "r")
fastaScaffolds = open(fastaScaffolds, "r")
multiFasta = open(multiFasta,"w")

fastaScaffoldsRed = fastaScaffolds.read().replace("\n", "")

for line in progenitor:
    divided = line.split()
    newScaffoldName = divided[ headers["scaffold"] ]
    newScaffoldNumber = newScaffoldName[ newScaffoldName.index('_')+1: ]
    newStart = int( divided[ headers["start"] ]  )
    newEnd =  int( divided[ headers["end"] ]  )
    #Some lines has no geneId. And some has : after it.
    if geneIdSearch == 1:
        if len(divided) == 14:
            newSeqId = divided[ headers["geneIdGFF"] ]
            newSeqId = newSeqId[:-1]
        elif len(divided) == 12:
            newSeqId = divided[ headers["geneIdGFF"] ]
    else:
        newSeqId = divided[ headers["geneIdGTF"] ]
    #Check if gene has changed
    if currentSeqId != newSeqId:
        if geneEnd > 0:
            writegene(currentSeqId, currentScaffoldName,  currentScaffoldIndex,  geneStart, geneEnd,  multiFasta,  fastaScaffoldsRed,  separator)
            geneStart = 99999999999999999999999999999
            geneEnd = 0
        #Changing current gene which we will work now
        currentSeqId = newSeqId
        #Check if scaffold has changed (if seq hasn't changed either scaffold)
        if currentScaffoldNumber != newScaffoldNumber:
            currentScaffoldNumber = newScaffoldNumber
            currentScaffoldName = newScaffoldName
            currentScaffoldIndex = fastaScaffoldsRed.find( scaffoldName + str(currentScaffoldNumber) ,  0) + len( scaffoldName + str(currentScaffoldNumber) ) -1
        #Updating gene start and then    
        geneStart,  geneEnd = updateGeneSection(newStart,  geneStart,  newEnd,  geneEnd)
    else:
        #Updating gene start and then    
        geneStart,  geneEnd = updateGeneSection(newStart,  geneStart,  newEnd,  geneEnd)
#Must write last gene (file ends without changing <currentSeqId>)
writegene( currentSeqId, currentScaffoldName, currentScaffoldIndex,  geneStart, geneEnd,  multiFasta,  fastaScaffoldsRed,  separator)
#Closing files
progenitor.close()
fastaScaffolds.close()
multiFasta.close()

############################# END #############################
