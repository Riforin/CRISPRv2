
# coding: utf-8

# In[10]:

#!/usr/bin/env python
#title           : yCRISPRv3.py
#description     : defines several methods to generate oligos for CRISPR
#author          : Leo d'Espaux <leodespaux@gmail.com>, help from Xingkai Li
#date            : 20150423
#version         : 3.0
#notes           : includes options for 
#notes           : delGene replGene cutOligos
#notes           : at some point there will be one interface
#==============================================================================

from Bio.Seq import * 
from Bio.SeqUtils import MeltingTemp
import unittest

"""
Created on Thu Sep 11 16:01:49 2014
@author: Leo leodespaux@gmail.com
"""

def cutOligos(GeneName, cutsite, DNA):
    
# This function asks user for a target cut site, and generates oligos to clone that cut site into 
# a CRISPR plasmid, and to screen for positive clones containing that cut site. It may also generate
# universal homology regions to integrate or otherwise alter that target region. 
#
# The current cut site architecture is:  >>>>> YtRNAp-HDV ribozyme- >20nt< -gRNA <<<<<

    GeneName=input("Name, using quotes: ")
    cutsite=input("20-mer cut sequence, using quotes: ")
    DNA=input("Locus sequence +/- a few kb, using quotes: ")

    if DNA.find(cutsite)==-1:                   # If cutiste sequence found in ANTISENSE
        DNA=Seq(DNA).reverse_complement()       # then reverse DNA, and turn it into a string
    
    index=DNA.find(cutsite)+16                  # index gives the start position of the string, e.g., 0. 
                                                # we add 16 since index+0=start of 20-mer, so index+16=cut site, 
                                                # 3 nt before last of 20mer
        
    Lup=DNA[index-520:index-490]                # This primer binds 500bp upstream of cut site

    
    cutSequence=Seq("cgggtggcgaatgggacttt")+cutsite+Seq("gttttagagctagaaatagc")
    seqprimer=Seq("gacttt")+cutsite
    

    print("cut" + GeneName + "  " + cutSequence)
    print("Lcolony" + GeneName + "  " + seqprimer)
    print("Lup" + GeneName + "  " + Lup)

    
    
# Parameters are here for testing purposes; will be removed when complete. 
def delGene(geneName, cutsite):

# This function asks user for a chromosomal locus, a region to be deleted, a suitable CRIPSR cutsite 
# and outputs oligos for cloning of a pL308 Cas9-gRNA vector, and ones for generating a donor DNA
# to delete the unwanted chromosomal region. Primers Lup+Rdown produce a 1kb band if deletion was
# successful. 
# part of yCRISPRv3 by leodespaux@gmail.com

    #GeneName=input("Name, using quotes: ")
    #cutsite=input("20-mer cut sequence, using quotes: ").upper()
    locus = genomicData[geneName][0]    
    deletion = genomicData[geneName][1]

    deletion = Seq(deletion)
     
    if deletion.find(cutsite)==-1:
        if deletion.reverse_complement().find(cutsite)==-1:
            print ("WARNING: Guide 20-mer sequence not found in deletion region.")
        
    locus=Seq(locus)
    
    index=locus.find(deletion)                  
   
    # index gives the start position within locus of the string deletion.
    # now we delete the deletion region to redefine a newlocus:
    
    newlocus=locus[0:index]+locus[index+len(deletion):]

    # note that since index starts at 0, a value of n points to, in the newlocus,
    # the first nt after the deletion. So we define the newlocus as above. Note too
    # that a string of len=40 ends at an index of 39--so we pick up at index+len-1. 
    
    Lup=newlocus[index-500:index-470]
    Rdown=newlocus[index+469:index+499].reverse_complement()

    Rtemp1 = newlocus[:index].reverse_complement()
    Rtemp2 = newlocus[index:].reverse_complement()

    rPrimer, rLength = getPrimer(Rtemp1)
    lPrimer, lLength = getPrimer(newlocus[index:])

    Rup = getOverhang(Rtemp2, rLength) + rPrimer
    Ldown = getOverhang(newlocus[:index], lLength) + lPrimer

    cutSequence=Seq("cgggtggcgaatgggacttt")+cutsite+Seq("gttttagagctagaaatagc")
    seqprimer=Seq("gacttt")+cutsite
    
    print("cut" + GeneName + "  " + cutSequence)
    print("seq" + GeneName + "  " + seqprimer)
    print("Lup" + GeneName + "del" + " " + Lup)
    print("Rup" + GeneName + "del" + " " + Rup)
    print("Ldown" + GeneName + "del" + " " + Ldown)
    print("Rdown" + GeneName + "del" + " " + Rdown)

    return Ldown, Rup
    
# TODO: Test insertGene
def insertGene(locusName, cutsite, insertName):

# This function asks user for a target cut site, a target locus, and a DNA sequence to be inserted.
# The function outputs primers to allow one to PCR the insert with homology to the cut site. The final
# integration happens 17nt from the cut site in both direction, such that a 34nt region is deleted.
# The current cut site architecture is:  >>>>> YtRNAp-HDV ribozyme- >20nt< -gRNA <<<<<
    """
    locusName=input("Locus name, using quotes: ")
    cutsite=input("20-mer cut sequence, using quotes: ").upper()
    insertName=input("Name of gene to be inserted, using quotes: ")
    """
    locus = genomicData[locusName][0]
    insert = genomicData[locusName][1]

    if locus.find(cutsite)==-1:                   # If cutsite sequence found in ANTISENSE
        locus=Seq(locus).reverse_complement()     # then reverse DNA, and turn it into a string
    
    index=locus.find(cutsite)                     # index gives the start position of cutsite within locus

    Rtemp1 = locus[:index].reverse_complement()
    Rtemp2 = insert.reverse_complement()

    rPrim, rLength = getPrimer(Rtemp1)
    lPrim, lLength = getPrimer(locus[index+34:])

    Rup = getOverhang(Rtemp2, rLength) + rPrim        # This is the sequence to the left of the insertion
    Ldown = getOverhang(insert, lLength)  + lPrim                         # This is the sequence to the right of the insertion 

    Lhom=locus[index-60:index]
    Rhom=Seq(locus[index+34:index+94]).reverse_complement()
    
    print("L_" + insertName + "_:" + locusName + " " + Lprimer)
    print("R_" + insertName + "_:" + locusName + " " + Rprimer)
    print("Lext:" + locusName + " " + Lhom)
    print("Rext:" + locusName + " " + Rhom)

    return Rup, Ldown

def getPrimer(seq):
    
    Tmax = 57
    maxLen = 60
    overhang = 0.8

    mp = 0
    length = 0
    primer = Seq("")

    while mp <= Tmax and length <= maxLen / (1 + overhang):
        primer = primer + seq[length]
        mp = MeltingTemp.Tm_NN(primer)
        length += 1

    lastThree = primer[-3:]
    while 'A' not in lastThree and 'T' not in lastThree:
        primer = primer[:-1]
        lastThree = primer[-3:]

    return primer, length

def getOverhang(seq, pLength):

    maxLen = 60
    overhang = 0.8

    length = 1
    oh = Seq("")
    while length <= overhang * pLength:
        oh = seq[-length] + oh
        length += 1
    return oh

# Fills the database. 
def getSequences(geneName):
    from intermine.webservice import Service

    template = service.get_template('Gene_GenomicDNA')

    rows = template.rows(
        E = {"op": "LOOKUP", "value": geneName, "extra_value": "S. cerevisiae"}
    )

    count = 0
    for row in rows:
        geneSeq = Seq(row["sequence.residues"])
        locusSeq = Seq(row["chromosome.residues.locus"])

        index = locusSeq.find(geneSeq)
        locusSeq = locusSeq[index-1000:locusSeq]

        # Reduce locusSize so it is only +/- 1 kbp of geneSeq
        break

    return geneSeq, locusSeq


class test_yCRISPRv4(unittest.TestCase):
    import unittest
    # TODO: Find a test case for Insertion
    def test_delGene(self):
        geneName = "LRO1"
        cutSite = "ATTATTTCATCGCAGGGTAT"
        locusSequence = "ATTACATTAGCATGCTTACATGGATGAATGGAAACAGATAACACTGAATTCTTGTAAAAAGGTAGCTTTTCAATGGTGGCTGTCTTTGTTCTATAATCTGCTGATATATCTTCGAACATCTGCTCAGGGCTTAGTGGTGAACCATTAGAATTAAACCCCACTATATACATTTTAGGGACCCTATAAGATGTCGAGTACGCAATATAAAGGTCGTAATACCTTTCTTGCGCCATATCTTTGGCTAGGCCACCTTTAGCATTAAATTCTTCTGTATCGTCATTTTCATCCTCCTCTTTTATTTCCATGTCTTGTATTAGCTCATCAATATCATCGATAGATGAGTCCTTGGTCCCCCCCGCAGGCGTACTTTGCACATGTTCAGTTTCAGATCCTATGTATTCCAGAACATCATCTTCATCTCCATCTTCTGCAAAACCTTTCATGATTACATCTGGCCCTTCAACTTCGACACATTGCTCGGCACGTTTATCACAGGGAACTTTTCTAATTATCAGAAACTGCTTATTCTTCGGTAAAAAATCTCTGTAACTAATATCTGACGACTCTTCGTTCCACTTCCAGGTGGGAAACATATGACATAAATAATCACCTGCTTGTACAAACTCCTCAGGAGTTATTTGACCTGTGGTTAAAAAGGTAGATTTGTGCGTTATGGGGGTAAGATATTCTCTCCAACTACTTAGTGTAGATCTAATCATGATAATAAAATTGTATTTACTCCTTTGTACTTCTTTGTTCCTAACTTCTAGCTAGCTTGATATACGTCAAGTGTACGTAATCATTCAAGTATCTACTTTCCTTTAAATAGCCCTTCGCTCTGCTTCCCTTATAATAATAATGTCCGTGTACAAAATTGAATGCCCAAGAAGGTGAAATTTGAAAAATAAGAAATCTACTAAATACCGATACGAAGAAGCGTATAGTAACAGCCATTACAAAAGGTTCTCTACCAACGAATTCGGCGACAATCGAGTAAAAAATGGGCACACTGTTTCGAAGAAATGTCCAGAACCAAAAGAGTGATTCTGATGAAAACAATAAAGGGGGTTCTGTTCATAACAAGCGAGAGAGCAGAAACCACATTCATCATCAACAGGGATTAGGCCATAAGAGAAGAAGGGGTATTAGTGGCAGTGCAAAAAGAAATGAGCGTGGCAAAGATTTCGACAGGAAAAGAGACGGGAACGGTAGAAAACGTTGGAGAGATTCCAGAAGACTGATTTTCATTCTTGGTGCATTCTTAGGTGTACTTTTGCCGTTTAGCTTTGGCGCTTATCATGTTCATAATAGCGATAGCGACTTGTTTGACAACTTTGTAAATTTTGATTCACTTAAAGTGTATTTGGATGATTGGAAAGATGTTCTCCCACAAGGTATAAGTTCGTTTATTGATGATATTCAGGCTGGTAACTACTCCACATCTTCTTTAGATGATCTCAGTGAAAATTTTGCCGTTGGTAAACAACTCTTACGTGATTATAATATCGAGGCCAAACATCCTGTTGTAATGGTTCCTGGTGTCATTTCTACGGGAATTGAAAGCTGGGGAGTTATTGGAGACGATGAGTGCGATAGTTCTGCGCATTTTCGTAAACGGCTGTGGGGAAGTTTTTACATGCTGAGAACAATGGTTATGGATAAAGTTTGTTGGTTGAAACATGTAATGTTAGATCCTGAAACAGGTCTGGACCCACCGAACTTTACGCTACGTGCAGCACAGGGCTTCGAATCAACTGATTATTTCATCGCAGGGTATTGGATTTGGAACAAAGTTTTCCAAAATCTGGGAGTAATTGGCTATGAACCCAATAAAATGACGAGTGCTGCGTATGATTGGAGGCTTGCATATTTAGATCTAGAAAGACGCGATAGGTACTTTACGAAGCTAAAGGAACAAATCGAACTGTTTCATCAATTGAGTGGTGAAAAAGTTTGTTTAATTGGACATTCTATGGGTTCTCAGATTATCTTTTACTTTATGAAATGGGTCGAGGCTGAAGGCCCTCTTTACGGTAATGGTGGTCGTGGCTGGGTTAACGAACACATAGATTCATTCATTAATGCAGCAGGGACGCTTCTGGGCGCTCCAAAGGCAGTTCCAGCTCTAATTAGTGGTGAAATGAAAGATACCATTCAATTAAATACGTTAGCCATGTATGGTTTGGAAAAGTTCTTCTCAAGAATTGAGAGAGTAAAAATGTTACAAACGTGGGGTGGTATACCATCAATGCTACCAAAGGGAGAAGAGGTCATTTGGGGGGATATGAAGTCATCTTCAGAGGATGCATTGAATAACAACACTGACACATACGGCAATTTCATTCGATTTGAAAGGAATACGAGCGATGCTTTCAACAAAAATTTGACAATGAAAGACGCCATTAACATGACATTATCGATATCACCTGAATGGCTCCAAAGAAGAGTACATGAGCAGTACTCGTTCGGCTATTCCAAGAATGAAGAAGAGTTAAGAAAAAATGAGCTACACCACAAGCACTGGTCGAATCCAATGGAAGTACCACTTCCAGAAGCTCCCCACATGAAAATCTATTGTATATACGGGGTGAACAACCCAACTGAAAGGGCATATGTATATAAGGAAGAGGATGACTCCTCTGCTCTGAATTTGACCATCGACTACGAAAGCAAGCAACCTGTATTCCTCACCGAGGGGGACGGAACCGTTCCGCTCGTGGCGCATTCAATGTGTCACAAATGGGCCCAGGGTGCTTCACCGTACAACCCTGCCGGAATTAACGTTACTATTGTGGAAATGAAACACCAGCCAGATCGATTTGATATACGTGGTGGAGCAAAAAGCGCCGAACACGTAGACATCCTCGGCAGCGCGGAGTTGAACGATTACATCTTGAAAATTGCAAGCGGTAATGGCGATCTCGTCGAGCCACGCCAATTGTCTAATTTGAGCCAGTGGGTTTCTCAGATGCCCTTCCCAATGTAAATGACCGACATTGACTCACTATCCATCCGTGTATTATTTCAAAGAGCGAAAAGAAGGCGCGTCGCGTCGACGCGCCTTTTTAGGCTAGAAAATAAACAGAAAACAAAAACAAAAACAAAAAAAGGCGAAAAAACAAACGAAAAAACAAACGACAGTAGATAGAGGAGAAGGTTTTTGACAGGTTTGTGTAATTGGTTATTTGCTTTTATAGATTATATATACACAAACTTAAAACACTGTTTTCTAAGTAGGCTGCGGAGAGGGCCAATTGCAGCAGCGACGTACGCATAGAGGAAAAACAACCGTTCATGTCCATTATGAAGCAGAGGCTACCACTGGGGGAGTTTTCCAGCTCGAAGATCAACAAACTGGCAATTGCCAATATTGCAGACGCTAGCGAGCCTAGAAACCATGGAGAAAACAACGTCGGGACTGTGTGTCTTCCCTCTATCAAGAGTTTGATGGTTAGTCCTGAAGTATACGAAAACACGAAGAGTCTTCCTGTTCCCCTAATGAGGAGCAGTGGCGGAGGTATGGCTTGTGCCAGTAAGTCCTCATGCCAGGATGGTATTAGTACGAAAACAACGTCTAGAGATTACTCTGAGCTTTCCAAGAAGTTGCAAATTCGTTTACAGTTCGCCTATTATAAATACAAGACTAAGCAAACAGACAAAAATTTCACAGACTTGAAATCAAAGCATAGCATTACGAGGCCTTCCAAAGTCGCGACTCACAGTAAATCAGAACCTTTAACTAGGAGAAGAAAATTAGTACTATCTCAGGGTCATTATAAAACCCCTGCTAGGTCTAAGATCAAAACTCCATCTTCGATTTGCTCTCACGATAATACATCTTCCTTTACATCTTTCCGCGGTGTCAGCGAAAGCTCGAGCACTACCGCAGATATGAACGTCGCAGATACCACCACACCAATACGCAACAACATAAACACAAAGCATTCAAACAGTCATAATCGCACATTGTATCAGAGACAAGAAAC"
        delSeq = "ATGGGCACACTGTTTCGAAGAAATGTCCAGAACCAAAAGAGTGATTCTGATGAAAACAATAAAGGGGGTTCTGTTCATAACAAGCGAGAGAGCAGAAACCACATTCATCATCAACAGGGATTAGGCCATAAGAGAAGAAGGGGTATTAGTGGCAGTGCAAAAAGAAATGAGCGTGGCAAAGATTTCGACAGGAAAAGAGACGGGAACGGTAGAAAACGTTGGAGAGATTCCAGAAGACTGATTTTCATTCTTGGTGCATTCTTAGGTGTACTTTTGCCGTTTAGCTTTGGCGCTTATCATGTTCATAATAGCGATAGCGACTTGTTTGACAACTTTGTAAATTTTGATTCACTTAAAGTGTATTTGGATGATTGGAAAGATGTTCTCCCACAAGGTATAAGTTCGTTTATTGATGATATTCAGGCTGGTAACTACTCCACATCTTCTTTAGATGATCTCAGTGAAAATTTTGCCGTTGGTAAACAACTCTTACGTGATTATAATATCGAGGCCAAACATCCTGTTGTAATGGTTCCTGGTGTCATTTCTACGGGAATTGAAAGCTGGGGAGTTATTGGAGACGATGAGTGCGATAGTTCTGCGCATTTTCGTAAACGGCTGTGGGGAAGTTTTTACATGCTGAGAACAATGGTTATGGATAAAGTTTGTTGGTTGAAACATGTAATGTTAGATCCTGAAACAGGTCTGGACCCACCGAACTTTACGCTACGTGCAGCACAGGGCTTCGAATCAACTGATTATTTCATCGCAGGGTATTGGATTTGGAACAAAGTTTTCCAAAATCTGGGAGTAATTGGCTATGAACCCAATAAAATGACGAGTGCTGCGTATGATTGGAGGCTTGCATATTTAGATCTAGAAAGACGCGATAGGTACTTTACGAAGCTAAAGGAACAAATCGAACTGTTTCATCAATTGAGTGGTGAAAAAGTTTGTTTAATTGGACATTCTATGGGTTCTCAGATTATCTTTTACTTTATGAAATGGGTCGAGGCTGAAGGCCCTCTTTACGGTAATGGTGGTCGTGGCTGGGTTAACGAACACATAGATTCATTCATTAATGCAGCAGGGACGCTTCTGGGCGCTCCAAAGGCAGTTCCAGCTCTAATTAGTGGTGAAATGAAAGATACCATTCAATTAAATACGTTAGCCATGTATGGTTTGGAAAAGTTCTTCTCAAGAATTGAGAGAGTAAAAATGTTACAAACGTGGGGTGGTATACCATCAATGCTACCAAAGGGAGAAGAGGTCATTTGGGGGGATATGAAGTCATCTTCAGAGGATGCATTGAATAACAACACTGACACATACGGCAATTTCATTCGATTTGAAAGGAATACGAGCGATGCTTTCAACAAAAATTTGACAATGAAAGACGCCATTAACATGACATTATCGATATCACCTGAATGGCTCCAAAGAAGAGTACATGAGCAGTACTCGTTCGGCTATTCCAAGAATGAAGAAGAGTTAAGAAAAAATGAGCTACACCACAAGCACTGGTCGAATCCAATGGAAGTACCACTTCCAGAAGCTCCCCACATGAAAATCTATTGTATATACGGGGTGAACAACCCAACTGAAAGGGCATATGTATATAAGGAAGAGGATGACTCCTCTGCTCTGAATTTGACCATCGACTACGAAAGCAAGCAACCTGTATTCCTCACCGAGGGGGACGGAACCGTTCCGCTCGTGGCGCATTCAATGTGTCACAAATGGGCCCAGGGTGCTTCACCGTACAACCCTGCCGGAATTAACGTTACTATTGTGGAAATGAAACACCAGCCAGATCGATTTGATATACGTGGTGGAGCAAAAAGCGCCGAACACGTAGACATCCTCGGCAGCGCGGAGTTGAACGATTACATCTTGAAAATTGCAAGCGGTAATGGCGATCTCGTCGAGCCACGCCAATTGTCTAATTTGAGCCAGTGGGTTTCTCAGATGCCCTTCCCAATGTAA"

        Ldown = "TTCGGCGACAATCGAGTAAAAAATGACCGACATTGACTCACTATCCATCC"
        Rup = "TAGTGAGTCAATGTCGGTCATTTTTTACTCGATTGTCGCCGAATTCGT"

        LdownOut, RupOut = delGene(geneName, cutSite, locusSequence, delSeq)
        self.assertEqual(str(LdownOut), Ldown)
        self.assertEqual(str(RupOut), Rup)

    @unittest.skip("Test not yet implemented - no test case.")
    def test_insertGene(self):
        #TODO: Fill this out
        pass


if __name__ == '__main__':
    unittest.main()
