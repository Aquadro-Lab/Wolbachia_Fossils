import pandas as pd
import sys


with open(sys.argv[1]) as file:
	for line in file:
		if ("Wolbachia" in line):
			print(line.rstrip())
