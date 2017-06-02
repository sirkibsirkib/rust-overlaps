class Sol:
	def __init__(self, idA, idB, O, OHA, OHB, OLA, OLB, K):
		self.idA = idA
		self.idB = idB
		self.O = O
		self.OHA = OHA
		self.OHB = OHB
		self.OLA = OLA
		self.OLB = OLB
		self.K = K
import sys

def load(filename):
	lines = [line.rstrip('\n') for line in open(filename)][1:]
	return {convert(line) for line in lines}

def convert(line):
	Sol(*tuple(line.split('\t')))

a = load(sys.argv[1])
b = load(sys.argv[2])

a_only = a.difference(b)
ab = a.intersection(b)
ab.remove(None)
b_only = b.difference(a)

print(len(a_only), len(ab), len(b_only))
print(ab)