#from numpy import *
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

line = load()
plot(line[0], line[1])
show()

