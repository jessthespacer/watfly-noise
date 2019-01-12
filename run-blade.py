import os, csv
import numpy as np

# FOR USE WITH NAFNOISE

# --- Rotor characteristics ---
numBlades = 5

directory = "Inputs"	# Set input directory

# --- Program run flags ---
runCalcs = True
postProcess = True
writeToCSV = True

# --- Program constants ---
nFreq = 34		# 34 frequencies considered

# Run NAFNoise for all input files
if runCalcs:
	print("BEGINNING AERACOUSTIC CALCULATIONS\n")
	i = 0
	for file in os.listdir(os.fsencode(directory)):
		filename = os.fsdecode(file)
		if filename.endswith(".ipt"):
			i += 1
			print("SECTION " + str(i))		# Print section number
			os.system("nafnoise Inputs\\" + filename)	# Run NAFNoise
			print("---\n")

noise = []

# Compile output data
if postProcess:
	print("COMPILING OUTPUT DATA\n")
	for file in os.listdir(os.fsencode(directory)):
		filename = os.fsdecode(file)
		if filename.endswith(".out"):
			with open(directory + "\\" + filename) as f:
				data = f.readlines()[13:]	# Drop first 13 lines

				# Store acoustic data for one blade
				noise.append([
					list(		# Convert map object to list
						map(	# Convert each element to float
							float, x.split()	# Remove spaces
						)
					) for x in data])

	noise = np.array(noise)		# NumPy makes my life easier
	
	# Find total blade noise
	totNoise = noise[:, :, -1]			# Get totals column only
	bladeNoise = 10 * np.log10(np.sum(10**(totNoise / 10), 0))	# Add up SPLs
	# Make matrix of SPL vs. frequency
	bladeNoise = [list(noise[0, :, 0]), list(bladeNoise)]

if writeToCSV:
	print("WRITING TO CSV\n")
	head = ["Frequency", "Pressure side TBL [dB]", "Suction side TBL [dB]", \
		"Separation side TBL [dB]", "Laminar [dB]", "Bluntness [dB]", \
		"Turbulent inflow [dB]", "Total [dB]"]

	with open(directory + "\\" + "noiseBreakdown.csv", "w+", \
		newline = "\n", encoding = "utf-8") as f:
		i = 0
		writer = csv.writer(f, dialect = "excel")
		for section in noise:
			i += 1
			writer.writerow(["SECTION " + str(i)])
			writer.writerow(head)
			writer.writerows(section)
			writer.writerow("\n")
		
	with open(directory + "\\" + "bladeNoise.csv", "w+", \
		newline = "\n", encoding = "utf-8") as f:
		writer = csv.writer(f, dialect = "excel")
		writer.writerow(["Frequency [Hz]", "SPL [dB]"])
		# Transpose lists
		for i in range(len(bladeNoise[0])):
			writer.writerow([x[i] for x in bladeNoise])
