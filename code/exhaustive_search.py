import RNA
import random
import subprocess
import math
import numpy as np
import time
start_time = time.time()
#import matplotlib.pyplot as plt

# Simulation parameter
_RT = 0.61632077549999997

# function that transforms DrTransformer log into list
def dr_list(filename):
	# read result from log file
	f = open(filename, "r")

	# build list from log file
	#saveSeq = 0
	saveStep = 0
	saveEnd = 0
	liste = []
	liste_end = []
	for x in f:
		if x.find("# Distribution") != -1:
			saveEnd = 1
			continue
		if x.find("# Transcription Step") != -1:
			saveStep = 1
			continue
		if x.find("#") != -1:
			continue
		if saveStep == 1:
			liste.append(x.split(" "))
		if saveEnd == 1:
			liste_end.append(x.split(" "))


		# clean list
	for i in range(0,len(liste)):
		liste[i] = list(filter(None, liste[i]))
		for j in range(0,len(liste[i])):
			if liste[i][j].find(']') != -1:
				liste[i][j] = liste[i][j][:-1]
			if liste[i][j].find('\n') != -1:
				liste[i][j] = liste[i][j][:-1]

			# clean list_end
	for i in range(0,len(liste_end)):
		liste_end[i] = list(filter(None, liste_end[i]))
		for j in range(0,len(liste_end[i])):
			if liste_end[i][j].find(']') != -1:
				liste_end[i][j] = liste_end[i][j][:-1]
			if liste_end[i][j].find('\n') != -1:
				liste_end[i][j] = liste_end[i][j][:-1]

	return(liste, liste_end)




# length of pause is set to 20 seconds
#pause_len = 100

# obtain sequence and structure from base.txt
# read in base file
f = open("base.txt", "r")
meta = []
for x in f:
	if x.find('\n') != -1:
		meta.append(x[:-1].split(" "))
	else:
		meta.append(x.split(" "))

# perform exhaustive search for all sequence structure pairs in base.txt
for line in meta:
	if line == " ":
		break
	start_time = time.time()
	seq = line[0]
	struc = line[1]
	print("Structure "+struc+ " processing now")
	n = len(seq)
	# call DrTransformer without pausing site, controll
	subprocess.call(["echo "+seq+" | DrTransformer --logfile"], shell=True)
	# read in list from DrTransformer logfile
	liste, liste_end = dr_list("NoName.log")
	new_ratio = 0
	# find structure rate at end of transcritption
	for i in range(0,len(liste_end)):
		if liste_end[i][2] == struc:
			new_ratio = float(liste_end[i][6])
	# set highscore for all other calculations
	highscore = new_ratio
	old_highscore = highscore
	# search for all steps with no structure alternative in the beginning
	no_alt = 0
	while no_alt < n:
		if int(liste[no_alt][1]) == 1:
			no_alt += 1
		else:
			break
	
	if no_alt == n:
		print("There are no alternative structures possible")
		continue

	step = 0
	#highscore = 0
	best_pause = 0
	possible_length = [5,10,20,50,100]
	step_set = list(range(no_alt, n-5))
	closest_dist = n # find closest structure if pause is not available
	closest_call = 0
	closest_pause = 0
	closest_ratio = 0
	closest_struc = ""
	closest_length = 0

	# lists to store the positions for pauses and the length
	pause_list = []
	length_list = []
	# adaptive walk for one seq struc pair
	# stop after 2n tries for one pausing side -> exp to visit all sides twice
	while step < len(step_set) and len(pause_list) < 3:
	# alternative:
	#while step < n:

		# select pausing side for step
		new_pause = step_set[step]

		# check if equilibrium is reached
		id = new_pause - 1
		step_Z_list = []
		step_occupancy_list = []
		while(id < n):
			if int(liste[id][0]) == new_pause:
				step_Z_list.append(math.e**(-float(liste[id][6])/_RT))
				step_occupancy_list.append(float(liste[id][6]))
			if int(liste[id][0]) == new_pause + 1:
				break
			id +=1
		
		step_Z_list = np.array(step_Z_list)
		step_occupancy_list = np.array(step_occupancy_list)

		Z = sum(step_Z_list)
		myp8 = np.array(step_Z_list/Z)
		if(np.allclose(step_occupancy_list, myp8)):
			step += 1
			continue

		
		for i in possible_length:
			pause_len = i

			# build call for pausing sites
			pause_call = ""
			if len(pause_list)>=1:
				for l in range(len(pause_list)):
					pause_call = pause_call + " "+str(pause_list[l])+"="+str(length_list[l])
			pause_call = pause_call + " "+str(new_pause)+"="+str(pause_len)
			# call DrTransformer with pausing site
			subprocess.call(["echo "+seq+" | DrTransformer --logfile --pause-sites "+pause_call], shell=True)

			# read in list from DrTransformer logfile
			liste, liste_end = dr_list("NoName.log")

			# set new ratio to zero in case DrTransformer does not find strucutre
			new_ratio = 0
			# find structure rate at end of transcritption
			for i in range(0,len(liste_end)):
				if liste_end[i][2] == struc:
					new_ratio = float(liste_end[i][6])
					break 	
				elif RNA.bp_distance(struc, liste_end[i][2]) < closest_dist:
					closest_struc = liste_end[i][2]
					closest_dist = RNA.bp_distance(struc, liste_end[i][2])
					closest_ratio = float(liste_end[i][6])
					closest_pause = new_pause
					closest_length = pause_len

			# evaluate step
			if new_ratio > highscore:
				pause_list.append(new_pause)
				length_list.append(pause_len)
				highscore = new_ratio
				step_set.remove(new_pause)
				step = 0
				break


			# check if equilibrium is reached after this pause site length
			#id = new_pause - 1
			step_Z_list = []
			step_occupancy_list = []
			for id in range(0, len(liste_end)):
				if int(liste_end[id][0]) == new_pause:
					step_Z_list.append(math.e**(-float(liste_end[id][6])/_RT))
					step_occupancy_list.append(float(liste_end[id][6]))
				if int(liste_end[id][0]) == new_pause + 1:
					break
				
		
			step_Z_list = np.array(step_Z_list)
			step_occupancy_list = np.array(step_occupancy_list)

			Z = sum(step_Z_list)
			myp8 = np.array(step_Z_list/Z)
			if(np.allclose(step_occupancy_list, myp8)):
				break

		
		step += 1

	# build call for pausing sites for printing
	pause_call = ""
	if len(pause_list)>=1:
		for l in range(len(pause_list)):
			pause_call = pause_call + " "+str(pause_list[l])+"="+str(length_list[l])
			closest_call="NA,NA,NA,NA"
	else:
		pause_call = "no improvement found with pausing"
		closest_call =str(closest_dist)+","+closest_struc+","+str(closest_pause)+"="+str(closest_length)+","+str(closest_ratio)

	# save results
	#with open('result.txt', 'a+') as f:
	#	f.writelines(seq + " "+ struc+ " "+str(highscore)+ " "+ pause_call + "\n")
	with open('result.csv', 'a+') as f:
		# sequence/structure/old occupancy/new occupancy/pauses/closest distance/closest structure/closest pause/closest occupancy
		f.writelines(seq + ","+ struc+ ","+str(old_highscore)+","+str(highscore)+ ","+ pause_call +"," + closest_call + "\n")
	print(str(highscore)+ ","+ pause_call + "\n")

	print("--- %s seconds ---" % (time.time() - start_time))

	# save steps for later processing
	#with open('log/'+struc+'pause.txt', 'a+') as f:
	#	for i in range(0, len(pause_list)):
	#		f.writelines(str(pause_list[i]) + " " + str(pause_res[i])+"\n")
