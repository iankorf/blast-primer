import random

def random_seq(size, alph='ACGT'):
	seq = []
	for j in range(size):
		seq.append(random.choice(alph))
	return ''.join(seq)
