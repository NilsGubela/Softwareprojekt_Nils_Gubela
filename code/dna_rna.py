import sys

seq = sys.argv[1]
n = len(seq)

out = ''

for i in range(0,n):
	if seq[i] == "A":
		out = out + "U"
	elif seq[i] == "T":
		out = out + "A"
	elif seq[i] == "G":
		out = out + "C"
	else:
		out = out + "G"

print(out)