import subprocess
import sys

length = sys.argv[1]
quantity = int(sys.argv[2])

n = 0
while n < quantity:

	subprocess.call(["python3 rand_seq.py "+length+"|RNAsubopt -e 6 -s --noLP|barriers --minh=5 > out.txt"], shell=True)

	f = open("out.txt", "r")
	liste = []
	for x in f:
		liste.append(x.split(" "))

	#print(liste)
	seq = liste[0][5][:-1]

	if len(liste) < 3:
		continue
	struc1 = liste[1][4]
	struc2 = liste[2][4]
	# save results
	with open('base.txt', 'a+') as f:
		f.writelines(seq + " "+ struc1+"\n")
		f.writelines(seq + " "+ struc2+"\n")

	n += 1



	#saveStep = 0
	#liste = []
	#for x in f:
#		if x.find("# Distribution") != -1:
#			break
#		if saveStep == 1:
#			liste.append(x.split(" "))
		#if saveSeq == 1:
		#	seq = x[x.find(' ')+1:]
		#	seq = seq[:-2]
		#	n = len(seq)
		#	saveSeq = 0

		#if x.find("# >NoName") != -1:
		#	saveSeq = 1

#		if x.find("# Transcription Step") != -1:
#			saveStep = 1

		# clean list
#	for i in range(0,len(liste)):
#		liste[i] = list(filter(None, liste[i]))
#		for j in range(0,len(liste[i])):
#			if liste[i][j].find(']') != -1:
#				liste[i][j] = liste[i][j][:-1]
#			if liste[i][j].find('\n') != -1:
#				liste[i][j] = liste[i][j][:-1]
