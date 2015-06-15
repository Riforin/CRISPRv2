from Bio.Seq import * 
from Bio.SeqUtils import MeltingTemp
import unittest
"""
Stitches together the three fragments generated from 
insertion/deletion/replacement. Fragments is a list of 
[up, insert, down] 
"""
def stitchTogether(fragments):
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
	Rdown = down[-30:0].reverse_complement()

	print("Lup" + " " + Lup)
	print("Rup" + " " + Rup)
	print("Ldown" + " " + Ldown)
	print("Rdown" + " " + Rdown)

	return Lup, Ldown, Rup, Rdown

"""
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
