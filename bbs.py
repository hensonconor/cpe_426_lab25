import random
import sympy

def generateN():
	while True:
		random.seed(96)
		p_seed = random.randint(1, 1024)
		q_seed = random.randint(1, 1024)
		p = sympy.nextprime(p_seed)
		q = sympy.nextprime(q_seed)
		N = p * q
		if (N < 65536):
			return N
	
	
def bbs(width=16):
	random.seed(24)
	x0 = random.randint(0, 65535)
	N = generateN()
	bits = ""
	bits += str(x0%2)
	x = x0
	for i in range(width-1):
		xi = x * x
		bits += str((xi% N)%2)
		x = xi
	return bits
	
for _ in range(7):
	print(bbs(), end='')
