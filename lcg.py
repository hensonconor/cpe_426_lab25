import random

def lcg(n=2):
	random.seed()
	a = random.randint(1, 255) 
	b = random.randint(1, 255)
	x0 = random.randint(1, 255)
	M = random.randint(1, 255)
	nums = ""
	for i in range(n):
		xi = (a * x0 + b) % M
		num = bin(xi).replace("0b", "")
		l = 8 - len(num)
		print(l)
		pad = ""
		while l > 0:
			pad += "0"
			l -= 1
		new = pad + num
		nums += new
		x0 = xi
	return nums