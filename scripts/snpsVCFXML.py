# -*- encoding: utf-8 -*-

from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import sys
#num_args = 3
#
#if len(sys.argv) != num_args:
#    print("\nERROR calling the script. Number of arguments must be %d\nSyntax: python %s [BLAST1.txt] [BLAST2.txt]\n" %(num_args-1, sys.argv[0]))
#    sys.exit(1)

blastPC15 = sys.argv[1]
blastPC9 = sys.argv[2]
#blastPC15 = "PC15*PC9.xml"
#blastPC9 = "PC9*PC15.xml"
alelleLimit = 80 #At least 80% of concordance
semiAlelleLimit = 50 #At least 50% of concordance

############################# FUNCTIONS #############################
def readBlast(blastFile):
    #Parse and List all records in order to access in an easy way to data
    #Note: not parsing data is like reading plain-text. Not listing records will only allow to access once to each record
    blast_records = list(NCBIXML.parse(blastFile))
    genes = {}
    for record in blast_records:
        query_title = record.query.split('|')
        query_id = query_title[0]
        query_scaffold = query_title[1]
        query_start = query_title[2]
        query_end = query_title[3]
#        print ("\n" + str(query_id) + " == " + str(query_scaffold) + " == " + str(query_start) + " == " + str(query_end))
        for alignment in record.alignments:
            num_alignment = 0
            alignments = {}
            sbjct_title = alignment.title.split()
            sbjct_title  = sbjct_title[1].split('|')
            sbjct_id = sbjct_title[0]
            sbjct_scaffold = sbjct_title[1]
            sbjct_start = sbjct_title[2]
            sbjct_end = sbjct_title[3]
#            print (str(sbjct_id) + " == " + str(sbjct_scaffold) + " == " + str(sbjct_start) + " == " + str(sbjct_end))
            for hsp in alignment.hsps:
                 alignments[num_alignment] = [hsp.query_start,  hsp.sbjct_start,  hsp.frame,  hsp.query,  hsp.match,  hsp.sbjct,  hsp.identities,  hsp.align_length,  hsp.query_end,  hsp.sbjct_end,  hsp.score]
                 num_alignment =  num_alignment + 1
#                 print ("\t" + str(hsp.query_start) + " -- " + str(hsp.sbjct_start) + " -- " + str(hsp.frame) + " -- " +str(hsp.identities))
            gen = [query_id, query_scaffold, query_start, query_end, sbjct_id, sbjct_scaffold, sbjct_start, sbjct_end, alignments,  alignment.length]
            genes[query_id] = gen
    return genes

def getAlelles(PC15, PC9,  alelleLimit,  semiAlelleLimit):
    alelos = []
    semialelos = []
    wrong = []
    girados = []
    noreciprocos = []
    noconcordantes = []

    alelos_cont = 0
    semialelos_cont = 0
    wrong_cont = 0 
    girados_cont = 0
    no_reciprocos_cont = 0
    no_concordantes_cont = 0
    reciprocos_cont = 0
    alellesAndSemiAlelles ={}

    for gen in PC15:
        #First we check if all of them exists
        if gen and PC15[gen][4] and PC15[gen][4] in PC9:
            #Are PC15 query == PC9 subject AND PC9 query == PC15 subject?
            PC15Query = gen
            PC15Subject = PC15[gen][4]
            PC9Query = PC15Subject
            PC9Subject = PC9[PC9Query][4]
            #If for both the best matching is each other
            if PC15Query == PC9Subject and PC15Subject == PC9Query:
#                alellesAndSemiAlelles[PC15Query]=PC15[PC15Query]
                #Chech PC15's and PC9's matches
                alingnments = PC15[gen][8]
                matchPC15 = str(alingnments[0][2])
                matchPC15 = matchPC15[matchPC15.find(",")+1:-1]
                matchPC15 = int(matchPC15)
                identitiesPC15 = 0
                errMatchPC15 = 0
                errMatchPC9 = 0
                for hsp in alingnments:
                    match = str(alingnments[hsp][2])
                    match = match[match.find(",")+1:-1]
                    match = int(match)
                    if match != matchPC15:
                        errMatchPC15 = -1
                    identitiesPC15 = identitiesPC15 + alingnments[hsp][6]
                alingnments = PC9[PC15Subject][8]
                matchPC9 = str(alingnments[0][2])
                matchPC9 = matchPC9[matchPC9.find(",")+1:-1]
                matchPC9 = int(matchPC9)
                identitiesPC9 = 0
                for hsp in alingnments:
                    match = str(alingnments[hsp][2])
                    match = match[match.find(",")+1:-1]
                    match = int(match)
                    if match != matchPC9:
                        errMatchPC9 = -1
                    identitiesPC9 = identitiesPC9 + alingnments[hsp][6]
                #Both matches are OK (maybe not the same)
                if errMatchPC15 != -1 and errMatchPC9 != -1:
                    PC15Length = PC15[gen][9]
                    PC9Length = PC9[PC15Subject][9]
                    #IdentitiesPC15/PC15 length AND IdentitiesPC15/PC9 length (and viceversa)
                    rel = min( (identitiesPC15 * 100 / PC15Length),  (identitiesPC15 * 100 / PC9Length), (identitiesPC9 * 100 / PC15Length),  (identitiesPC9 * 100 / PC9Length))
                    if rel>=alelleLimit:
                        alellesAndSemiAlelles[PC15Query]=PC15[PC15Query]
                        alelos.append(PC15Query)
                        alelos_cont = alelos_cont + 1
                    elif rel>=semiAlelleLimit:
                        alellesAndSemiAlelles[PC15Query]=PC15[PC15Query]
                        semialelos.append(PC15Query)
                        semialelos_cont = semialelos_cont + 1
                    else:
                        wrong.append(PC15Query)
                        wrong_cont = wrong_cont + 1
                else:
                    girados_cont = girados_cont + 1
                    girados.append(PC15Query)
                reciprocos_cont = reciprocos_cont + 1
            else:
                no_reciprocos_cont = no_reciprocos_cont + 1
                noreciprocos.append(PC15Query)
        else:
            no_concordantes_cont = no_concordantes_cont + 1
            noconcordantes.append(PC15Query)
            
#    print ("Reciprocos:\t" + str(reciprocos_cont))
#    print ("\tAlelos:\t\t" +  str(alelos_cont))
#    print ("\tSemialelos:\t" +  str(semialelos_cont))
#    print ("\tAlelos y semialelos:\t" +  str(alelos_cont+semialelos_cont))
#    print ("\tWrong:\t\t" +  str(wrong_cont))
#    print ("\tGirados:\t" +  str(girados_cont))
#    print ("No reciprocos:\t" + str(no_reciprocos_cont))
#    print ("No concordantes:" + str(no_concordantes_cont))
#    print ("Genes PC15:\t" + str(len(PC15)))
#    print ("Genes PC9:\t" + str(len(PC9)))
#    print ("")

    return alellesAndSemiAlelles

def getSnps(genes):
#alignments[num_alignment] = [hsp.query_start,  hsp.sbjct_start,  hsp.frame,  hsp.query,  hsp.match,  hsp.sbjct,  hsp.identities,  hsp.align_length]
#gen = [query_id, query_scaffold, query_start, query_end, sbjct_id, sbjct_scaffold, sbjct_start, sbjct_end, alignments,  alignment.length]
    print ("##fileformat=VCFv4.2")
    print ("##reference=PC15v2 & PC9v1")
    print ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    chrom = "1"
    quality = "."
    filter = "PASS"
    info = "."
    counterGen = 0
    for gen in genes:
        subject = genes[gen][1]
        alignment = genes[gen][8]
        geneQueryStart = int(genes[gen][2])
        geneSubjectStart = int(genes[gen][6])
        scaffold = genes[gen][1]
        chrom = str(scaffold)
        #Get only number and remove initial zeros
        chrom = str(int(chrom[chrom.find("_")+1:]))
        counterGen = 0
#         ##########################################################
#        if int(gen) == 165005:
#            print ("--->\n")
#            print (str(subject) + "\n")
#            for align in alignment:
#                queryStart = alignment[align][0]
#                subjectStart = alignment[align][1]
#                strand =alignment[align][2]
#                querySeq = alignment[align][3]
#                snpSeqcuence = alignment[align][4]
#                subjectSeq = alignment[align][5]
#                print (str(queryStart) + "\t" + str(subjectStart) + "\t" + str(strand) + "\t" + str(querySeq) + "\t" + str(snpSeqcuence) + "\t" + str(subjectSeq))
#            print ("\n")
#            print ("<---\n")
#        ##########################################################
#        if subject = 1080027:
        if 1:
            for align in alignment:
                queryStart = int(alignment[align][0])
                subjectStart = int(alignment[align][1])
                strand =alignment[align][2]
                querySeq = alignment[align][3]
                snpSeqcuence = alignment[align][4]
                subjectSeq = alignment[align][5]
                queryEnd = int(alignment[align][8])
                subjectEnd = int(alignment[align][9])
                score = str(alignment[align][10])
    #            print (str(strand))
                #Depending on strand, index must be added or subtracted
                queryStrand = 1
                subjectStrand = 1
                if queryEnd<queryStart:
                    queryStrand = -1
                    my_dna_query = Seq(querySeq, generic_dna)
                    querySeq = my_dna_query.complement()
                    subjectStrand = -1                
                    my_dna_subject = Seq(subjectSeq, generic_dna)
                    subjectSeq = my_dna_subject.complement()
                #1) Check if snp is in first position
                if snpSeqcuence[:2] == " |": 
                    if querySeq[0] != "-" and querySeq[0] != "n" and subjectSeq[0] != "-" and subjectSeq[0] != "n":
                        counterGen = counterGen + 1
    #                    print str(gen) + " " + str(subject)  + " " + str(queryStart)  + " " + str(subjectStart)  + " " + str(geneQueryStart + queryStart - 1)  + " " + str(geneSubjectStart + subjectStart - 1)  + " " + str(querySeq[0]) + " " + str(subjectSeq[0])
                        print  (chrom + "\t" + str(geneQueryStart + queryStart - 1) + "\t" +  str(gen) + "-" + str(counterGen) + "\t" + str(querySeq[0]) + "\t" + str(subjectSeq[0]) + "\t" + quality + "\t" + filter + "\t" + info + "\t" + score)
                #2) Check rest of snp in seqcuence
                snpIndex = snpSeqcuence.find("| |")
                snpPos = snpIndex + 1
                queryPos = queryStart + snpPos - querySeq.count('-', 0, snpIndex)  
                subjectPos = subjectStart + snpPos - subjectSeq.count('-', 0, snpIndex)  
                if queryStrand == -1:
                    queryPos = queryStart - snpPos + querySeq.count('-', 0, snpIndex)  
                if subjectStrand == -1:
                    subjectPos = subjectStart - snpPos + subjectSeq.count('-', 0, snpIndex)  
                while snpIndex >= 0:
                    #Indicate position, CAUTION: may be '-', in this case they MUST NOT be taken in count
                    #CAUTION: snpIndex says <| |> index in string, which is actually it's begining. Index of snp is snpIndex + 2 (+1 cause python goes from 0 to (length - 1) and +1 cause <|> is first character searched for index)

                    if querySeq[snpPos] != "-" and querySeq[snpPos] != "n" and subjectSeq[snpPos] != "-" and subjectSeq[snpPos] != "n":
                        counterGen = counterGen + 1
    #                    print str(gen) + " " + str(subject)  + " " + str(queryPos)  + " " + str(subjectPos)  + " " + str(geneQueryStart + queryPos - 1)  + " " + str(geneSubjectStart + subjectPos - 1)   + " " + str(querySeq[snpPos]) + " " + str(subjectSeq[snpPos])
                        print (chrom + "\t" + str(geneQueryStart + queryPos - 1) + "\t" +  str(gen) + "-" + str(counterGen) + "\t" + str(querySeq[snpPos]) + "\t" + str(subjectSeq[snpPos]) + "\t" + quality + "\t" + filter + "\t" + info + "\t" + score)
                    snpIndex = snpSeqcuence.find("| |",  snpIndex + 1)
                    
                    snpPos = snpIndex + 1
                    queryPos = queryStart + snpPos - querySeq.count('-', 0, snpIndex)  
                    subjectPos = subjectStart + snpPos - subjectSeq.count('-', 0, snpIndex)  
                    queryRel = snpPos - querySeq.count('-', 0, snpIndex)  
                    subjectRel = snpPos - subjectSeq.count('-', 0, snpIndex)  
                    if queryStrand == -1:
                        queryPos = queryStart - snpPos + querySeq.count('-', 0, snpIndex)  
                        queryRel = snpPos + querySeq.count('-', 0, snpIndex)  
                    if subjectStrand == -1:
                        subjectPos = subjectStart - snpPos + subjectSeq.count('-', 0, snpIndex)  
                        subjectRel = snpPos + subjectSeq.count('-', 0, snpIndex)  

                snpPos = len(snpSeqcuence) - 1
                queryPos = queryStart + snpPos - querySeq.count('-', 0, snpIndex)  
                subjectPos = subjectStart + snpPos - subjectSeq.count('-', 0, snpIndex)  
                if queryStrand == -1:
                    queryPos = queryStart - snpPos + querySeq.count('-', 0, snpIndex)  
                if subjectStrand == -1:
                    subjectPos = subjectStart - snpPos + subjectSeq.count('-', 0, snpIndex)  
                #3) Check if snp is in last position
                if snpSeqcuence[-2:] == "| ":
                    if querySeq[snpPos] != "-" and querySeq[snpPos] != "n" and subjectSeq[snpPos] != "-" and subjectSeq[snpPos] != "n":
                        counterGen = counterGen + 1
    #                    print str(gen) + " " + str(subject)  + " " + str(queryPos)  + " " + str(subjectPos)  + " " + str(geneQueryStart + len(snpSeqcuence) - 1 - 1)  + " " + str(geneSubjectStart + len(snpSeqcuence) - 1 - 1)  + " " + str(querySeq[snpPos]) + " " + str(subjectSeq[snpPos])
                        print  (chrom + "\t" + str(geneQueryStart + len(snpSeqcuence) - 1 - 1) + "\t" +  str(gen) + "-" + str(counterGen) + "\t" + str(querySeq[snpPos]) + "\t" + str(subjectSeq[snpPos]) + "\t" + quality + "\t" + filter + "\t" + info + "\t" + score)
############################# PROGRAM #############################

if __name__ == '__main__':
    
    blastPC15 = open(blastPC15, "r")
    blastPC9 = open(blastPC9, "r")
    
    genesPC15 = readBlast(blastPC15)
    genesPC9 = readBlast(blastPC9)

    filteredGenes = getAlelles(genesPC15,  genesPC9, alelleLimit,  semiAlelleLimit)
    
    snps = getSnps(filteredGenes)
    
#    for gen in filteredGenes:
#        print ("\n" + str(filteredGenes[gen][0]) + " == " + str(filteredGenes[gen][1]) + " == " + str(filteredGenes[gen][2]) + " == " + str(filteredGenes[gen][3]))
#        print (str(filteredGenes[gen][4]) + " == " + str(filteredGenes[gen][5]) + " == " + str(filteredGenes[gen][6]) + " == " + str(filteredGenes[gen][7]))
#        alingnments = filteredGenes[gen][8]
#        for hsp in alingnments:
#            print ("\t" + str(alingnments[hsp][0]) + " -- " + str(alingnments[hsp][1]) + " -- " + str(alingnments[hsp][2]) + " -- " +str(alingnments[hsp][6]) + " -- " +str(alingnments[hsp][7]))
#            print (str(len(alingnments[hsp][3])) + " = " + str(alingnments[hsp][3]))
#            print (str(len(alingnments[hsp][4])) + " = " + str(alingnments[hsp][4]))
#            print (str(len(alingnments[hsp][5])) + " = " + str(alingnments[hsp][5]))
    #Closing files
    blastPC15.close()
    blastPC9.close()
####ESTABAS REVISANDO LOS MATCHES, PARA CALCULAR SI ESTAN DADOS LA VUELTA O NO
