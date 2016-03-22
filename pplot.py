from numpy import *
from matplotlib.pyplot import *

def load():
	result = []
	for i in range(2):
		result.append([])
	with open('graph', 'r') as f:
		for line in f:
			val = map(float, line.split())
			for i in range(2):
				result[i].append(val[i])
	return result

# radar green, solid grid lines
rc('grid', color='#316931', linewidth=1, linestyle='-')
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)

# force square figure and square axes looks better for polar, IMO
width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size, size))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='#ffffff')
another_bg = '#d5de9c'

line = load()
rad = line[0]
R = line[1]
n = len(line[0])
#print(n)
#R = array(R) - 20

def label(r):
	return str(r) + "dB"

ax.plot(rad, R, color='#ee8d18', lw=3)
ax.set_yticks(range(0, 80, 20))
labels = map(str, range(-40, 40, 20))
labels[3] = labels[3] + 'dB'
ax.set_yticklabels(labels)
ax.set_rmax(80)
grid(True)
show()


