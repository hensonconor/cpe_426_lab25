#https://gist.github.com/StuartGordonReid
import scipy.special as spc
import math

def monobit_freq(input_str):
	s_n = 0
	for bit in input_str:
		if bit == '0':
			s_n -= 1
		else:
			s_n += 1

	s_obs = abs(s_n) / math.sqrt(len(input_str))
	p_val = spc.erfc(math.fabs(s_obs) / math.sqrt(2))

	return p_val

def runs(input_str):
	ones_count = 0
	n = len(input_str)
	for bit in input_str:
		if bit == '1':
			ones_count += 1
	p = float(ones_count / float(n))
	v_obs = 1
	tau = 2 / math.sqrt(n)
	if abs(p - 0.5) > tau:
		return 0.0
	for i in range(1, n):
		if input_str[i] != input_str[i-1]:
			v_obs += 1
	num = abs(v_obs - 2.0 * n * p * (1.0 - p))
	den = 2.0 * math.sqrt(2.0 * n) * p * (1.0 - p)
	p_val = spc.erfc(float(num / den))
	return p_val

def non_overlapping_pattern(input_str, pattern="000000001", num_blocks=8):
	n = len(input_str)
	pattern_size = len(pattern)
	block_size = math.floor(n / num_blocks)
	pattern_counts = numpy.zeros(num_blocks)
	for i in range(num_blocks):
		block_start = i * block_size
		block_end = block_start + block_size
		block_data = input_str[block_start:block_end]

		j = 0
		while j < block_size:
			sub_block = block_data[j:j + pattern_size]
			if sub_block == pattern:
				pattern_counts[i] += 1
				j += pattern_size
			else:
				j += 1
	mean = (block_size - pattern_size + 1) / pow(2, pattern_size)
	var = (block_size * ((1 / pow(2, pattern_size)) - (((2 * pattern_size) - 1) / 
		(pow(2, pattern_size * 2)))))
	chi_squared = 0
	for i in range(num_blocks):
		chi_squared += pow(pattern_counts[i] - mean, 2.0 ) / var
	p_val = spc.gammaincc(num_blocks / 2, chi_squared / 2)
	return p_val

def overlapping_pattern(input_str, pattern_size=9, block_size=1032):
	n - len(input_str)
	pattern = ""
	for i in range(pattern_size):
		pattern += "1"
	num_blocks = math.floor(n / block_size)
	lambda_val = float(block_size - pattern_size + 1) / pow(2, pattern_size)
	eta = lambda_val / 2.0

	piks = [self.get_prob(i, eta) for i in range(5)]
	diff = float(numpy.array(piks).sum())
	piks.append(1.0 - diff)

	pattern_counts = numpy.zeros(6)
	for i in range(num_blocks):
		block_start = i * block_size
		block_end = block_start + block_size
		block_data = input_str[block_start:block_end]

		pattern_count = 0
		j = 0
		while j < block_size:
			sub_block = block_data[j:j + pattern_size]
			if sub_block == pattern:
				pattern_count += 1
			j += 1
		if pattern_count <= 4:
			pattern_counts[pattern_count] += 1
		else:
			pattern_counts[5] += 1

	chi_squared = 0.0
	for i in range(len(pattern_counts)):
		chi_squared += pow(pattern_counts[i] - num_blocks * piks[i], 2.0) / (num_blocks * piks[i])
	return spc.gammaincc(5.0 / 2.0, chi_squared / 2.0)

def get_prob(u, x):
  out = 1.0 * numpy.exp(-x)
  if u != 0:
    out = 1.0 * x * numpy.exp(2 * -x) * (2 ** -u) * spc.hyp1f1(u + 1, 2, x)
  return out

def linear_complexity(self, bin_data, block_size=500):
  dof = 6
  piks = [0.01047, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833]

  t2 = (block_size / 3.0 + 2.0 / 9) / 2 ** block_size
  mean = 0.5 * block_size + (1.0 / 36) * (9 + (-1) ** (block_size + 1)) - t2

  num_blocks = int(len(bin_data) / block_size)
  if num_blocks > 1:
    block_end = block_size
    block_start = 0
    blocks = []
    for i in range(num_blocks):
      blocks.append(bin_data[block_start:block_end])
      block_start += block_size
      block_end += block_size

    complexities = []
    for block in blocks:
      complexities.append(self.berlekamp_massey_algorithm(block))

      t = ([-1.0 * (((-1) ** block_size) * (chunk - mean) + 2.0 / 9) for chunk in complexities])
      vg = numpy.histogram(t, bins=[-9999999999, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 9999999999])[0][::-1]
      im = ([((vg[ii] - num_blocks * piks[ii]) ** 2) / (num_blocks * piks[ii]) for ii in range(7)])

      chi_squared = 0.0
      for i in range(len(piks)):
        chi_squared += im[i]
      p_val = spc.gammaincc(dof / 2.0, chi_squared / 2.0)
      return p_val
  else:
    return -1.0

def berlekamp_massey_algorithm(self, block_data):
    n = len(block_data)
    c = numpy.zeros(n)
    b = numpy.zeros(n)
    c[0], b[0] = 1, 1
    l, m, i = 0, -1, 0
    int_data = [int(el) for el in block_data]
    while i < n:
      v = int_data[(i - l):i]
      v = v[::-1]
      cc = c[1:l + 1]
      d = (int_data[i] + numpy.dot(v, cc)) % 2
      if d == 1:
        temp = copy.copy(c)
        p = numpy.zeros(n)
        for j in range(0, l):
          if b[j] == 1:
            p[j + i - m] = 1
        c = (c + p) % 2
        if l <= 0.5 * i:
          l = i + 1 - l
          m = i
          b = temp
      i += 1
    return l

def approximate_entropy(bin_data, pattern_length=10):
  n = len(bin_data)
  # Add first m+1 bits to the end
  # NOTE: documentation says m-1 bits but that doesnt make sense, or work.
  bin_data += bin_data[:pattern_length + 1:]

  # Get max length one patterns for m, m-1, m-2
  max_pattern = ''
  for i in range(pattern_length + 2):
    max_pattern += '1'

  # Keep track of each pattern's frequency (how often it appears)
  vobs_one = numpy.zeros(int(max_pattern[0:pattern_length:], 2) + 1)
  vobs_two = numpy.zeros(int(max_pattern[0:pattern_length + 1:], 2) + 1)

  for i in range(n):
    # Work out what pattern is observed
    vobs_one[int(bin_data[i:i + pattern_length:], 2)] += 1
    vobs_two[int(bin_data[i:i + pattern_length + 1:], 2)] += 1

  # Calculate the test statistics and p values
  vobs = [vobs_one, vobs_two]
  sums = numpy.zeros(2)
  for i in range(2):
    for j in range(len(vobs[i])):
      if vobs[i][j] > 0:
        sums[i] += vobs[i][j] * math.log(vobs[i][j] / n)
  sums /= n
  ape = sums[0] - sums[1]
  chi_squared = 2.0 * n * (math.log(2) - ape)
  p_val = spc.gammaincc(pow(2, pattern_length-1), chi_squared/2.0)
  return p_val

def random_excursions(bin_data):
  # Turn all the binary digits into +1 or -1
  int_data = numpy.zeros(len(bin_data))
  for i in range(len(bin_data)):
    if bin_data[i] == '0':
      int_data[i] = -1.0
    else:
      int_data[i] = 1.0

  # Calculate the cumulative sum
  cumulative_sum = numpy.cumsum(int_data)
  # Append a 0 to the end and beginning of the sum
  cumulative_sum = numpy.append(cumulative_sum, [0])
  cumulative_sum = numpy.append([0], cumulative_sum)

  # These are the states we are going to look at
  x_values = numpy.array([-4, -3, -2, -1, 1, 2, 3, 4])

  # Identify all the locations where the cumulative sum revisits 0
  position = numpy.where(cumulative_sum == 0)[0]
  # For this identify all the cycles
  cycles = []
  for pos in range(len(position) - 1):
    # Add this cycle to the list of cycles
    cycles.append(cumulative_sum[position[pos]:position[pos + 1] + 1])
  num_cycles = len(cycles)

  state_count = []
  for cycle in cycles:
    # Determine the number of times each cycle visits each state
    state_count.append(([len(numpy.where(cycle == state)[0]) for state in x_values]))
  state_count = numpy.transpose(numpy.clip(state_count, 0, 5))

  su = []
  for cycle in range(6):
    su.append([(sct == cycle).sum() for sct in state_count])
  su = numpy.transpose(su)

  piks = ([([self.get_pik_value(uu, state) for uu in range(6)]) for state in x_values])
  inner_term = num_cycles * numpy.array(piks)
  chi = numpy.sum(1.0 * (numpy.array(su) - inner_term) ** 2 / inner_term, axis=1)
  p_values = ([spc.gammaincc(2.5, cs / 2.0) for cs in chi])
  return p_values

def get_pik_value(k, x):
  if k == 0:
    out = 1 - 1.0 / (2 * numpy.abs(x))
  elif k >= 5:
    out = (1.0 / (2 * numpy.abs(x))) * (1 - 1.0 / (2 * numpy.abs(x))) ** 4
  else:
    out = (1.0 / (4 * x * x)) * (1 - 1.0 / (2 * numpy.abs(x))) ** (k - 1)
  return out


print(runs("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"))
