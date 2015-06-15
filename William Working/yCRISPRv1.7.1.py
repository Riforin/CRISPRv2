#!/usr/bin/env python
#title           : yCRISPRv1.7 
#description     : defines several methods to generate oligos for CRISPR
#author          : Leo d'Espaux <leodespaux@gmail.com>
#date            : 9 June 2015
#==============================================================================



# import libraries we're using
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from intermine.webservice import Service


# define global variables
HomologyLength=1000 



def askUser():
    print("Hi baby, what are you in the mood for today? Using quotes enter")
    print("\"1\" to insert to an empty locus") 
    print("\"2\" to edit an existing gene")
    print(" ")

    Action= input("Your answer: ")
    if Action == "1":
        print(" ")
        print("Do you want to construct a cassette?")
        print("\"1\" to construct using standard promoter and terminator names, and a custom sequence")
        print("\"2\" to construct using standard promoter and terminator names, and a custom sequence")
        print(" ")

        typeEdit=input("Your answer: ")
        if typeEdit=="1":
            print("Coming")
        inserttoEmpty(typeEdit)
    elif Action == "2":
        editExisting()
    
                # promt and ask to pick between these two options
                # 1 == "insert a gene into a new locus"
                # 2 == "edit exiting locus
                # and go to either of the next two functios
                     
                    
                    
def buildCassette():
    print("Coming soon.")        

     
    
    
def inserttoEmpty(Sequence):
    print("Which of these characterized loci is going to get lucky: ")
    print(" ")
    print("208a, the first we used, defined as 100% ")    
    print("1021b, 105% ")
    print(" ")
    print("You know the drill, all answers in quotes.")
    
    filledLocus=input("Your answer: ")
    
    # find locus info to primers
    # promt a table  
    
    
    
    
def editExisting():
    print(" ")
    print("OK sweetie. Which locus do you want to edit, e.g., \"OAF1\": ")
    print("Right now I only have chromosome 1 for testing.")

    
    GeneName= input("Your answer: ")
    
    OrigGeneRecord=fetchGene(GeneName)
    #note that this returns a seqrecord
    
    print(" ")
    print("I'm fetching upstream and downstream homology reghhions for you...")
    print(" ")

    
    # let's load the appropriate chromosome file. The record of the gene we looked up
    # contains in the "features" the systematic name, wherein the second letter
    # corresponds to chromosome number, e.g., 1=A etc
    if OrigGeneRecord.features[1]=="A":
        origChromosomeRec=SeqIO.read("Scer01.fasta", "fasta")
    
    # let's explicitely name the sequences from the seq record
    OrigGeneSeq=OrigGeneRecord.seq
    OrigChromosomeSeq=origChromosomeRec.seq
    
    # flip the sequence to orient with respect to the old gene
    if OrigChromosomeSeq.find(OrigGeneSeq)==-1:
        OrigChromosomeSeq=OrigChromosomeSeq.reverse_complement()

        print ("OK. I'm done. You might want to know that your gene is in the Crick strand.")
        print ("Don't worry about who Crick is, though. I flipped your chromosome so you don't get confused.")
        
    # Note that fetchHomology takes seq as inputs for gene and locus
    UpHomSeq=fetchHomology(OrigGeneSeq,OrigChromosomeSeq,"upstream",HomologyLength)
    DownHomSeq=fetchHomology(OrigGeneSeq,OrigChromosomeSeq,"downstream",HomologyLength)

    
    UpHomRec = SeqRecord(UpHomSeq, id="UpHom")
    DownHomRec = SeqRecord(DownHomSeq, id="DownHom")
    
    print("OK I found them.")

    print(" ")
    print("UpstreamHomology reghhion is")
    print(UpHomRec.seq)
    
    print(" ")
    print("DownstreamHomology reghhion is")
    print(DownHomRec.seq)
    
    print(" ")
    print("What do you want to do to this gene?")
    print("\"1\" delete it cleanly")
    print("\"2\" replace it with an ORF/pre-made fragment")
    print("\"3\" replace it with a standard cassette I will build for you")
    print("\"4\" replace a specified region near your target gene")
    
    Action=input("Your answer: ")
                        
    #note that in all the below, we want to have fragments be records
    if Action=="1":
        fragments=[UpHomRec,DownHomRec]
    
    if Action=="2":
        print(" ")
        NewGeneName=input("What's the name of the gene you're inserting?")
        NewGeneSeq=Seq(input("What's the sequence of your new gene?"))
        InsertRec = SeqRecord(NewGeneSeq, id=NewGeneName)
        fragments=[UpHomRec, InsertRec, DownHomRec]

    if Action=="3":
        print(" ")
        print("I got you baby. I don't want you to strain your eyes looking up stuff.")
        print("Actually I haven't yet gotten to that so no insert.")

        InsertRec=buildCassette()
        fragments=[UpHomRec, DownHomRec] #we need to finish buildcassette to add InsertRec here

    if Action=="4":
        print("Actually I haven't yet gotten to that so no insert.")
        fragments=[UpHomRec, DownHomRec] #we need to finish buildcassette to add InsertRec here

        
        
    # note that fragments are seqrecords

    stitchTogether(fragments)
    
    print("")
    
    
"""
Stitches together fragments and returns primers to do so. 
@param fragments SeqRecord List of [up insert down] or [up down] 
"""    
def stitchTogether(fragments):
    #this function takes seq records and prints primers
    #at some point maybe a .gb file should be written
    print(" ")
    print("Hmm.. that's a bold choice. But OK, I'll design it for you.")
    print("These are what you want stitched together, in order, I think:")

    print(" ")
    for i in fragments:
        print(i.id)

    print(" ")
    print("Now we're cooking with gas! I have some primers for you:")
    
    up = fragments[0]
    insert = fragments[1]
    # If deletion is occuring, reassign down as insert
    if len(fragments) == 2:
        down = insert
    else:
        down = fragments[2]


    Rtemp1 = up.reverse_complement()
    Rtemp2 = insert.reverse_complement()

    rPrim, rLength = getPrimer(Rtemp1)
    lPrim, lLength = getPrimer(down)

    Rup = getOverhang(Rtemp2, rLength) + rPrim
    if len(fragments) == 2:
        Ldown = getOverhang(up, lLength) + lPrim
    else:
        Ldown = getOverhang(insert, lLength) + lPrim
    
 
    Lup = up[:30]
    Rdown = down[-30:].reverse_complement()

    # Lup and Rdown are actually SeqRecords and must be converted to seq and then str 
    # to print nicely. 

    print("Lup" + " " + str(Lup.seq))
    print("Rup" + " " + Rup)
    print("Ldown" + " " + Ldown)
    print("Rdown" + " " + str(Rdown.seq))
    
    print(" ")
    
"""
Helper method for stitchTogether
Gets the primer from a given sequence based on melting temperature and max length. 
@param seq A sequence (or string) that represents the 
       sequence from which to derive the primer. The primer
       is generated from left to right
"""

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


"""
Helper method for stitchTogether
Gets the overhang from a given sequence based on non-overhang length. 
@param seq Sequence from which to derive overhang. overhang is 
       generated from right to left
@param pLength Length of the non-overhang part of the primer. 
"""
def getOverhang(seq, pLength):

    maxLen = 60
    overhang = 0.8

    length = 1
    oh = Seq("")
    while length <= overhang * pLength:
        oh = seq[-length] + oh
        length += 1
    return oh


    #------------------------------ FETCH FUNCTIONS -------------------------------------

    
    
    
def fetchGene(GeneName):
    
    service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")
    template = service.get_template('Gene_GenomicDNA')

    rows = template.rows(
        E = {"op": "LOOKUP", "value": GeneName, "extra_value": "S. cerevisiae"}
    )
    
    # this service seems to return multiple similar genes but we want the first one only, so count
    # and it returns information about the gene you want
    count=0
    for row in rows:
        
        count=count+1
        if count==1:
            descr= row["description"]
            GeneSeq=Seq(row["sequence.residues"])
            GeneSysName=row["secondaryIdentifier"]
            print(" ")
            print("I think you want...... "+row["secondaryIdentifier"])
            print(row["description"])
            print(" ")
            print(row["sequence.residues"])
            print(" ")
            print("Good choice! I have a feeling you're going to get lucky with this one.")
            print(" ")
            print("Give me a second to put some of my ducks in a circle...")
       

            
    #let's create a record for the oldGene
    GeneRecord = SeqRecord(GeneSeq, id=GeneSysName)
    
    #now let's add some more information to make it useful
    GeneRecord.name=GeneName
    GeneRecord.features=GeneSysName

    return GeneRecord
    
  
    

    
def fetchHomology(OrigGeneSeq,OrigChromosomeSeq,direction,HomologyLength):
    # This function takes seq inputs for gene and locus, and also returns a seq

    StartIndex=OrigChromosomeSeq.find(OrigGeneSeq)
    EndIndex=StartIndex+len(OrigGeneSeq)
    
    if direction=="upstream":
        HomologySeq=OrigChromosomeSeq[StartIndex-HomologyLength:StartIndex]
    if direction=="downstream":
        HomologySeq=OrigChromosomeSeq[EndIndex:EndIndex+HomologyLength]

    return HomologySeq
    # note that these are returned as an array of .seq 

    
    
    #allRecords=[SeqIO.read("Yali0A_contig.fasta", "fasta"), SeqIO.read("Yali0B_contig.fasta", "fasta"),
      #  SeqIO.read("Yali0C_contig.fasta", "fasta"), SeqIO.read("Yali0D_contig.fasta", "fasta"),
    #    SeqIO.read("Yali0E_contig.fasta", "fasta"), SeqIO.read("Yali0F_contig.fasta", "fasta")]
        
       # for records in allRecords:
            # record.find(input)
            # return chromosome, id, direction. if direction is different,
            # return relative to the original ID, i.e., flip the chromosome view
            
     
    
    
    
    
    
    
    #------------------------------------ CONSTRUCTING STUFF --------------------------------------

    

        
def makeCassette(construction):
    print("I can help you with that. Don't worry, I won't take your job. You are smart and beautiful.")
    print("Do you want to type the name of the promoter? e.g. \"TEF1p\"? If not, type \"custom\"")
    return SeqRecord(Seq("ATCG",id="test"))

  
    
    



    
askUser()
