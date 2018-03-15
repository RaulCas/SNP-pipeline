# -*- encoding: utf-8 -*-

import sys


############################# PARAMETERS #############################
num_args = 2

if len(sys.argv) != num_args:
    print("\nERROR calling the script. Number of arguments must be %d\nSyntax: python %s sortedVCFXML.vcf\n" %(num_args-1, sys.argv[0]))
    sys.exit(1)

sortedVCFXMLFile = sys.argv[1]

#sortedVCFXMLFile = "resultsRef_PC9_VCFXML.vcf"


############################# PROGRAM #############################

sortedVCF = open(sortedVCFXMLFile, "r")
oldPos = -1
oldId = -1
lastLine = ""
for line in sortedVCF:
    campos = line.split()
    if line[:1] != "#":
        newPos = int(campos[1])
        newId = int(campos[2][:campos[2].find('-')])
        if oldPos == newPos and oldId == newId:
                lastLine = (line[:line.rfind("\t")])
        else:
            if lastLine != "":
                print (lastLine)
            lastLine = (line[:line.rfind("\t")])
            oldPos = newPos
            oldIde = newId
    else:
        print (line[:-1])
sortedVCF.close()
