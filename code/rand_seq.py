import random
import sys

n = int(sys.argv[1])

out = ''

while len(out) < n:
	r = random.randint(1,4)

	if r == 1:
		out = out + "A"
	elif r == 2:
		out = out + "C"
	elif r == 3:
		out = out + "G"
	else:
		out = out + "U"

print(out)