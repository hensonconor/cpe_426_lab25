import random
import sympy

def generateN():
	while True:
		random.seed()
		p_seed = random.randint(1, 1024)
		q_seed = random.randint(1, 1024)
		p = sympy.nextprime(p_seed)
		q = sympy.nextprime(q_seed)
		N = p * q
		if (N < 65536):
			print(N)
			return N
	
	
def bbs(width=16):
	random.seed()
	x0 = random.randint(0, 65535)
	N = generateN()
	bits = ""
	bits += str(x0%2)
	print(bits)
	x = x0
	for i in range(width):
		xi = x * x
		bits += str((xi% N)%2)
		x = xi
	return bits
	