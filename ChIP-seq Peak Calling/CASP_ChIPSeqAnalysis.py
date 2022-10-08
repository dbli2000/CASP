#Import helpful packages and set params
import re, csv
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (12,6)

# Parse input files
def csv_parse(file):
	with open(file) as csv_file:
	    csv_reader = csv.reader(csv_file, delimiter=',')
	    position_count = 0
	    coverage_count = 0
	    counts = []
	    for row in csv_reader:
	    	if position_count == 0:
	    		position_count = 1
	    	else:
	    		#Filter out the LacI Artifact (high coverage due to high copy number)
	    		if (position_count <= 368000) and (position_count >= 366000):
	    			counts.append(float(0))
	    		else:
		    		counts.append(float(row[1]))
		    	coverage_count += float(row[1])
		    	position_count += 1 
	print(f'Processed {position_count} lines of '+file+f' with an average coverage of {coverage_count/position_count} reads.')
	# Note that line count = size of genome and read counts overcounts reads due to sum
	return counts

IP = csv_parse("CASP-Sigma_IPs.csv")
IN = csv_parse("CASP-Sigma_Inputs.csv")
noHAIP = csv_parse("noHA_IPs.csv")
noHAIN = csv_parse("noHA_Inputs.csv")

# Median normalize files
IPn = IP/np.median(IP)
INn = IN/np.median(IN)
noHAIPn = noHAIP/np.median(noHAIP)
noHAINn = noHAIN/np.median(noHAIN)

#Find hits
hit = IPn.copy()
hit[IPn<=5] = 0
hit[noHAIPn>=3] = 0

#Aggregate peaks and output max coverage positions
hits = []
hit_coverage = []
for x in range(len(hit)):
	if hit[x] != 0:
		if len(hits) == 0:
			hits.append(x)
			hit_coverage.append(hit[x])
		elif (hits[-1]<=x-500 and hit[x-500]==0):
			hits.append(x)
			hit_coverage.append(hit[x])
		elif (hit[x]>hit[hits[-1]]):
			hits[-1] = x
			hit_coverage[-1]=hit[x]
print("Initial hits = ", hits)
print("Initial hit coverage = ", hit_coverage)
# This finds hits at [63578, 221462, 457786, 564759, 697052, 1721008, 1733453, 1799887,
# 1848116, 2286250, 2466283, 2482621, 2707068, 2812163, 2914011, 2978872, 3867466, 3880807, 
# 4175577, 4245812, 4279904, 4399421]

# To further boost confidence, peaks are visual filtered for a triangular peak
# as this would be expected for a true ChIP-seq peak. This reduces hits down to
hits = [63578, 221462, 457786, 1721008, 1733453, 1799887, 1848116, 2482621, 2707068, 2812163, 2914011, 2978872, 3880807]

# Build a simple filter for coloring genome track figure
color_track = [0]*(len(IPn))
for region in hits:
	for x in range(1000):
		color_track[region-500+x] = 1
thresh_track = IPn.copy()
thresh_track[thresh_track<=5] = 0
thresh_track[thresh_track>=5] = 1
final_filter = np.array([int(a and b) for a, b in zip(thresh_track, color_track)])
IPb = np.ma.masked_where(final_filter == 1, IPn) #Below Threshold of 5 or not at Hit
IPa = np.ma.masked_where(final_filter == 0, IPn) #Above Threshold of 5 and at Hit

# Plot Median Normalized Read Coverage across Genome (for SI)
f, (ax11, ax12, ax13, ax14) = plt.subplots(4, 1, sharex='col')
# Subplot 1
p11, = ax11.plot(IPb, label = "CASP-σ IP", color="#BCBDBF")
ax11.plot(IPa, color="#3292B9")
ax11.legend(handles=[p11], loc = "upper right")
# Subplot 2
p12, = ax12.plot(INn, label = "CASP-σ Input", color="#BCBDBF")
ax12.legend(handles=[p12], loc = "upper right")
# Subplot 3
p13, = ax13.plot(noHAIPn, label = "Control IP", color="#BCBDBF")
ax13.legend(handles=[p13], loc = "upper right")
# Subplot 4
p14, = ax14.plot(noHAINn, label = "Control Input", color="#BCBDBF")
ax14.legend(handles=[p14], loc = "upper right")
plt.ylabel("Coverage")
plt.xlabel("Position (in bp)")
f.savefig("GenomeWideChIP-seqTracks_SuppFig.svg", format='svg' , dpi=1200)

# Look at max coverage of each of the hits
for hit_pos in hits:
	print(f"Coverage of {IPn[hit_pos]} with {IP[hit_pos]} reads at {hit_pos}")
# The top 4 hits are at indices 4, 5, 6, and 11

# Plot 4 Main Peaks Figure (for main figure)
f2, axs = plt.subplots(2, 4, sharex='col', sharey='row')
# Subplot 0,0 LHR 
axs[0, 0].plot(INn[hits[4]-500:hits[4]+500], label = "Input", color="#BCBDBF")
axs[0, 0].plot(IPn[hits[4]-500:hits[4]+500], label = "CASP-σ IPs", color="#3292B9")
axs[0,0].set_ylabel("Read Coverage")
axs[0,0].set_ylim(0, 20)
# Subplot 1,0 LHR
axs[1, 0].plot(noHAINn[hits[4]-500:hits[4]+500], label = "Control Input", color="#BCBDBF")
axs[1, 0].plot(noHAIPn[hits[4]-500:hits[4]+500], label = "Control IP", color="#3292B9")
axs[1,0].set_ylabel("Read Coverage")
axs[1,0].set_ylim(0, 20)
axs[1,0].set_xlabel(f"{hits[4]+1}")

# Subplot 0,1 RpmI
p5, = axs[0, 1].plot(INn[hits[5]-500:hits[5]+500], label = "Input", color="#BCBDBF")
p6, = axs[0, 1].plot(IPn[hits[5]-500:hits[5]+500], label = "CASP-σ IPs", color="#3292B9")
# Subplot 1,1 RpmI
p7, = axs[1, 1].plot(noHAINn[hits[5]-500:hits[5]+500], label = "Input", color="#BCBDBF")
p8, = axs[1, 1].plot(noHAIPn[hits[5]-500:hits[5]+500], label = "Control IP", color="#3292B9")
axs[1,1].set_xlabel(f"{hits[5]+1}")
#axs[1,1].legend(handles=[p7, p8], loc = "upper right")

# Subplot 0,2 SelD
p9, = axs[0, 2].plot(INn[hits[6]-500:hits[6]+500], label = "Input", color="#BCBDBF")
p10, = axs[0, 2].plot(IPn[hits[6]-500:hits[6]+500], label = "CASP-σ IPs", color="#3292B9")
# Subplot 1,2 SelD
p11, = axs[1, 2].plot(noHAINn[hits[6]-500:hits[6]+500], label = "Input", color="#BCBDBF")
p12, = axs[1, 2].plot(noHAIPn[hits[6]-500:hits[6]+500], label = "Control IP", color="#3292B9")
axs[1,2].set_xlabel(f"{hits[6]+1}")

# Subplot 0,3 LysA
p13, = axs[0, 3].plot(INn[hits[11]-500:hits[11]+500], label = "Input", color="#BCBDBF")
p14, = axs[0, 3].plot(IPn[hits[11]-500:hits[11]+500], label = "CASP-σ IPs", color="#3292B9")
axs[0,3].legend(handles=[p14, p13], loc = "upper right")
# Subplot 1,3 LysA
p15, = axs[1, 3].plot(noHAINn[hits[11]-500:hits[11]+500], label = "Input", color="#BCBDBF")
p16, = axs[1, 3].plot(noHAIPn[hits[11]-500:hits[11]+500], label = "Control IP", color="#3292B9")
axs[1,3].legend(handles=[p16, p15], loc = "upper right")
axs[1,3].set_xlabel(f"{hits[11]+1}")
f2.savefig("ChIP-seqPeaks_4tracks_MainFig.svg", format='svg' , dpi=1200)
plt.show()