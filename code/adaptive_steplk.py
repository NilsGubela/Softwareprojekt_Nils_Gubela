import RNA
import random

# read in base file
f = open("base.txt", "r")
meta = []
for x in f:
	if x.find('\n') != -1:
		meta.append(x[:-1])
	else:
		meta.append(x.split(" "))

meta[4] = meta[4].split(" ")
print(meta)
# save pausing sides
pause = []
for i in range(1,len(meta[4])):
	pause.append(int(meta[4][i][:-2]))


# read in log file
f = open("NoName.log", "r")

# build list from log file
saveSeq = 0
saveStep = 0
liste = []
for x in f:
	if saveStep == 1:
		liste.append(x.split(" "))
	if saveSeq == 1:
		seq = x[x.find(' ')+1:]
		seq = seq[:-2]
		n = len(seq)
		saveSeq = 0

	if x.find("# >NoName") != -1:
		saveSeq = 1

	if x.find("# Tanscription Step") != -1:
		saveStep = 1


# clean list
for i in range(0,len(liste)):
	liste[i] = list(filter(None, liste[i]))
	for j in range(0,len(liste[i])):
		if liste[i][j].find(']') != -1:
			liste[i][j] = liste[i][j][:-1]
		if liste[i][j].find('\n') != -1:
			liste[i][j] = liste[i][j][:-1]


liste_end = liste[len(liste)-2:len(liste)]
liste.pop()
liste.pop()
liste.pop()
liste.pop()


# find structure rate at end of transcritption
for i in range(0,len(liste)):
	if int(liste[i][0]) != n:
		continue
	
	if liste[i][2] == meta[2]:
		new_ratio = float(liste[i][6])


# select new pausing side for next step
new_pause = 0
while new_pause == 0 or new_pause in pause:
	new_pause = random.randint(1,n)

# evaluate result from last step
if new_ratio <= float(meta[3]):
	# no benefit observed from introducing pausing side
	# remove pausing side
	meta[4].pop()
else:
	# new highscore, change ratio in base.txt
	meta[3] = new_ratio

if new_ratio >= float(meta[0]):
	# threshold is reached end program
	print("Threshold is reached, find pausing sides in base.txt")
	quit()

# add new pausing side
meta[4].append(str(new_pause)+"=1")
meta[4] = ' '.join(meta[4])
# write to base.txt
with open('base.txt', 'w') as f:
	f.writelines(str(line) + "\n" for line in meta)





