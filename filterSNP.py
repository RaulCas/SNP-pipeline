# -*- encoding: utf-8 -*-

import sys


############################# PARAMETERS #############################
num_args = 4

if len(sys.argv) != num_args:
    print("\nERROR calling the script. Number of arguments must be %d\nSyntax: python %s sortedVCFXML.vcf\n" %(num_args-1, sys.argv[0]))
    sys.exit(1)
modes = ["_Strict", "_Half"]
GATKoutputfile = sys.argv[1]
filterMode = int(sys.argv[2])
window = int(sys.argv[3])
output = open(GATKoutputfile[:GATKoutputfile.rfind(".")] + "_Cleared_" + str(window) + modes[int(filterMode) - 1] + GATKoutputfile[GATKoutputfile.rfind("."):], 'w')
#GATKoutputfile = "output.cvs"
#filterMode = 1


############################# PROGRAM #############################

snpFile = open(GATKoutputfile, "r")
oldLine = ""
currentLine = ""
nextLine = ""
i = 0
oldId = 0
currentId = 0
nextId = 0
oldPos = 0
currentPos = 0
nextPos = 0
for line in snpFile:
    oldLine = currentLine
    currentLine = nextLine
    nextLine = line
    flag = 0
    if i == 2:
        currentLineAux = currentLine.split()
        currentId = currentLineAux[2][:currentLineAux[2].find("-")]
        nextLineAux = nextLine.split()
        nextId = nextLineAux[2][:nextLineAux[2].find("-")]
        if currentId == nextId:
            currentPos = int(currentLineAux[1])
            nextPos = int(nextLineAux [1])
            if nextPos - currentPos > window:
                output.write(currentLine)
    elif i > 2:
        oldLineAux = oldLine.split()
        oldId = oldLineAux[2][:oldLineAux[2].find("-")]
        currentLineAux = currentLine.split()
        currentId = currentLineAux[2][:currentLineAux[2].find("-")]
        nextLineAux = nextLine.split()
        nextId = nextLineAux[2][:nextLineAux[2].find("-")]
        if oldId == currentId:
            oldPos = int(oldLineAux[1])
            currentPos = int(currentLineAux[1])
            if currentPos - oldPos < window:
                flag = flag + 1
        if currentId == nextId:
            currentPos = int(currentLineAux[1])
            nextPos = int(nextLineAux [1])
            if nextPos - currentPos < window:
                flag = flag + 1
        if flag < filterMode:
            output.write(currentLine)
    else:
        output.write(line)
    i = i + 1
if currentId == nextId:
    if nextPos - currentPos < window:
        output.write(nextLine)
snpFile.close()
