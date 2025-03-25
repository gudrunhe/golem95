# vim: ts=3:sw=3

from gmpy import mpq

F90_WIDTH = 80

TENS_DEBUG = False

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
   31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
   73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
   127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
   179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
   233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
   283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
   353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
   419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
   467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
   547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
   607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
   661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
   739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
   811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
   877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
   947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
   1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
   1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
   1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
   1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
   1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373]

zero = mpq(0)
one  = mpq(1)

precision_golem = "precision_golem"
array_golem = "array"
matrice_golem = "matrice_s"
ff_golem = "form_factor_type"
ffp_golem = "form_factor_%dp"

ff_suffix = ""
# ff_suffix = "_b"

golem_minlegs = 1

PAT = {
		'globsolve': 'solve%d',
		'solve': 'solve%d_%d',
		'globtenseval': 'tenseval%d',
		'cglobtenseval': 'ctenseval%d',
		'tenseval': 'tenseval%d_%d',
		'tenseval': 'tenseval%d_%d',
		'ctenseval': 'ctenseval%d_%d',
		'reconstruct': 'reconstruct%d',
		'coefftype': 'coeff_type_%d',
		'printcoeffs': 'print_coeffs_%d',
		'evaluate': 'evaluate',
		'contract': 'contract%d_%d',
		'mucontract': 'contract_%s',
		'tenscontracta': 'contract_a_tensor_%d',
		'tenscontractb': 'contract_b_tensor_%d',
		'tenscontractc': 'contract_c_tensor_%d',
		'tenscontractd': 'contract_d_tensor_%d'
		}

def combinat(n, k):
	"""
		Calculates the binomial coefficient (n atop k).
	"""
	if k < 0 or k > n:
		return 0
	else:
		num = 1
		den = 1
		for i in range(1, k+1):
			num *= n-i+1
			den *= i
		return num/den

def generate_mapping(R, k):
	"""
		Generates a mapping from tensor components \hat{C}(a_1, ..., a_k)
		into a one dimensional array.

		PARAMETER

		R  -- rank
		k  -- number of non-zero components of q

		RETURN

		(lst, dic)

		lst -- list of (a_1, ..., a_k)
		dic -- mapping from (a_1, ..., a_k) -> int

		lst[dic[X]] = X if X in dic
	"""

	def rec_generator(k, R):
		if k == 0:
			yield []
		elif k <= R:
			for a_1 in range(1, R - (k - 1) + 1):
				if k > 1:
					for tail in rec_generator(k - 1, R - a_1):
						yield [a_1] + tail
				else:
					yield [a_1]
	
	lst = []
	dic = {}
	i = 0
	for indices in rec_generator(k, R):
		t = tuple(indices)
		lst.append(t)
		dic[t] = i
		i += 1

	assert i == combinat(R, k), \
			"len(%s) != %d, R=%d,k=%d" % (lst,combinat(R, k),R,k)
	return lst, dic

def generate_equations(R, k):
	"""
		Generates a set of equations for a given number of non-zero
		components and fixed maximum rank.
	
		PARAMETER

		R  -- rank
		k  -- number of non-zero components of q

		RETURN

		(LHS, RHS)

		LHS -- a matrix (i.e. list of lists) of coefficients
		RHS -- a list of values of q
	"""

	lst, dic = generate_mapping(R, k)
	l = len(lst)

	LHS = []
	RHS = []

	for num_eq in range(l):
		q = map(lambda i: PRIMES[i], lst[num_eq])
		coeffs = [
			reduce(lambda x,y: x*y,	map(lambda (b,e): b**e, zip(q, term)), 1)
			for term in lst]

		LHS.append(coeffs)
		RHS.append(q)

	return LHS, RHS, lst, dic

def matrix_inverse(M_int):
	"""
	Matrix inversion by Gauss elimination.

	PARAMETER

	M_int -- square matrix, represented as list of lists

	RETURN

	square matrix, represented as list of lists which is the inverse of M_int
	"""
	M = []
	R = []
	n = len(M_int)
	i = 0
	for row in M_int:
		M.append(map(mpq, row))
		R.append([zero]*i + [one] + [zero]*(n-i-1))
		i += 1

	# Forward substitution
	for i in range(n):
		# Find pivot element
		max_el = abs(M[i][i])
		max_row  = i
		for j in range(i+1, n):
			if abs(M[j][i]) > max_el:
				max_row = j
				max_el = abs(M[j][i])

		# Swap pivot row to the top
		if max_row != i:
			tmp = M[i]
			M[i] = M[max_row]
			M[max_row] = tmp

			tmp = R[i]
			R[i] = R[max_row]
			R[max_row] = tmp

		if M[i][i] == zero:
			raise ZeroDivisionError, "Cannot invert singular matrix."

		pivot = M[i][i]
		for k in range(i, n):
			M[i][k] /= pivot
		for k in range(0, n):
			R[i][k] /= pivot

		for j in range(i+1, n):
			f = M[j][i]
			if f == zero:
				continue

			for k in range(i, n):
				M[j][k] -= f * M[i][k]
			for k in range(n):
				R[j][k] -= f * R[i][k]
				
	# Backward substitution
	for i in range(n-1, -1, -1):
		for j in range(i):
			f = M[j][i]
			if f == zero:
				continue

			for k in range(n):
				R[j][k] -= f * R[i][k]

	return R

def fmt_mpq_f90(q):
	"""
	Formats a rational number according to our Fortran 90 conventions,
	as real(kind=ki).

	PARAMETER

	q -- a number of type gmpy.mpq

	RETURN

	a string representing the rational number in Fortran 90.
	"""
	num = q.numer()
	den = q.denom()

	if num == zero:
		return "0.0_ki"
	elif den == one:
		return "%s.0_ki" % num
	else:
		return "%s.0_ki/%s.0_ki" % (num, den)

def write_matrix_f90(f, indent, lname, rname, k, LHS, RHS):
	"""
	Writes a matrix and a vector to a Fortran90 file.
	The matrix and the vector form a linear system in the sense,
	that LHS[i] * x = N(q_i), where N is the numerator function,
	and the non-zero entries of q_i are in RHS[i].

	PARAMETER

	f      -- a file object open for writing
	indent -- number of blank characters to be added at the beginning
	          of the line.
	lname  -- name of the matrix LHS
	rname  -- name of the matrix RHS
	k      -- number of non-zero entries in q
	LHS    -- a matrix
	RHS    -- vector
	"""

	MI = matrix_inverse(LHS)

	n = len(RHS)

	f.write(" "*indent +
			"real(ki), dimension(%d,%d), parameter, private :: %s = &\n"
			% (n, n, lname))
	f.write(" "*indent + "& reshape((/&\n")
	f.write(" "*indent + "&")
	col = indent + 1
	first = True
	for row in MI:
		for element in row:
			if first:
				first = False
			else:
				f.write(",")
				col += 1
			s = fmt_mpq_f90(element)
			l = len(s)
			if col + l + 2 > F90_WIDTH:
				f.write(" &\n" + " "*indent + "&")
				col = indent + 1
			f.write(s)
			col += l
	f.write(" "*indent + "/),&\n")
	f.write(" "*indent + "& (/%d,%d/), order=(/2,1/))\n" % (n, n))

	f.write(" "*indent +
			"real(ki), dimension(%d,%d), parameter, private :: %s = &\n"
			% (n, k, rname))
	f.write(" "*indent + "& reshape((/&\n")
	f.write(" "*indent + "&")
	col = indent + 1
	first = True
	for row in RHS:
		for element in row:
			if first:
				first = False
			else:
				f.write(",")
				col += 1
			s = "%s.0_ki" % element
			l = len(s)
			if col + l + 2 > F90_WIDTH:
				f.write(" &\n" + " "*indent + "&")
				col = indent + 1
			f.write(s)
			col += l
	f.write(" "*indent + "/),&\n")
	f.write(" "*indent + "& (/%d,%d/), order=(/2,1/))\n" % (n, k))

def write_numeval_interface(f, indent, name, d):
	"""
	Writes an interface defining the expected signature of the
	numerator function.

	We use a real q and mu^2 in this implementation and allow the function
	to return a complex value.

	PARAMETER
	f      -- a file open for writing
	indent -- number of blank characters to be added at the beginning
	          of the line
	name   -- name of the function in the interface
	d      -- number of space-time dimensions

	OUTPUS

	For name=NAME and d=4, this function produces the following output in f:

	interface
	   function NAME(Q, mu2)
			use precision_golem, only: ki
			implicit none
			real(ki), dimension(0:3), intent(in) :: Q
			real(ki), intent(in) :: mu2
			complex(ki) :: NAME
		end function NAME
	interface
	"""
	f.write(" "*indent+"interface\n")
	nin = indent + 3
	f.write(" "*nin+"function %s(Q, mu2)\n" % name)
	nin+=3
	f.write(" "*nin+"use %s, only: ki\n" % precision_golem)
	f.write(" "*nin+"implicit none\n")
	f.write(" "*nin+"real(ki), dimension(0:%d), intent(in) :: Q\n" % (d-1))
	f.write(" "*nin+"real(ki), intent(in) :: mu2\n")
	f.write(" "*nin+"complex(ki) :: %s\n" % name)
	nin-=3
	f.write(" "*nin+"end function %s\n" % name)
	f.write(" "*indent+"end interface\n")

def write_subroutine_solve(f, indent, name, k, dim, lname, rname, type_name,
		recon_name, d, recon_name_extra=None, type_name_extra=None):
	"""
	Writes a subroutine to solve for a certain set of coefficients.

	PARAMETER

	f          -- file open for writing
	indent     -- indentation level (number of blanks prefixing the line)
	name       -- name of this subroutine
	k          -- number of non-zero components in q
	dim        -- number of coefficients to solve for
	lname      -- name of the LHS matrix (inverse of the original LHS)
	rname      -- name of the RHS (values for q)
	type_name  -- name of the type containing the coefficints
	recon_name -- name of the function used to subtract known bits
	d          -- number of space time dimensions
	recon_name_extra -- name of the function used to subtract
	              a second set of coefficients
	type_name_extra -- name of the type containing the extra coefficients
	"""

	if recon_name_extra is not None:
		xtra = ", coeffs2"
	else:
		xtra = ""

	for line in DOC.subroutine_solve(name,xtra,lname,rname,d-1,type_name):
		f.write(" "*indent + line + "\n")

	f.write(" "*indent
			+"subroutine     %s(numeval, indices, mu2, coeffs, idx%s)\n"
			% (name, xtra))
	nin = indent + 3
	f.write(" "*nin+"! generated by: write_subroutine_solve\n")
	f.write(" "*nin+"implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin+"integer, dimension(%d), intent(in) :: indices\n" % k)
	f.write(" "*nin+"real(ki), intent(in) :: mu2\n")
	f.write(" "*nin+"type(%s), intent(inout) :: coeffs\n" % type_name)
	f.write(" "*nin+"integer, intent(in) :: idx\n")
	if recon_name_extra is not None:
		f.write(" "*nin+"type(%s), intent(in), optional :: coeffs2\n"
				% type_name_extra)
	f.write(" "*nin+"complex(ki), dimension(%d) :: xnum\n" % dim)
	f.write(" "*nin+"real(ki), dimension(0:%d) :: Q\n" % (d-1))
	f.write(" "*nin+"integer :: i\n")
	if k < d:
		f.write(" "*nin+"Q(:)=0.0_ki\n")
	if recon_name_extra is not None:
		f.write(" "*nin+"if (present(coeffs2)) then\n")
		nin += 3
		f.write(" "*nin+"do i=1,%d\n" % dim)
		nin += 3
		for j in range(1, k+1):
			f.write(" "*nin+"Q(indices(%d)) = %s(i,%d)\n" % (j,rname,j))
		f.write(" "*nin+"xnum(i) = numeval(Q, mu2) &\n")
		f.write(" "*(nin+3) + "& - %s(Q, coeffs, %d) &\n" % (recon_name, k-1))
		f.write(" "*(nin+3) + "& - %s(Q, coeffs2, %d)\n"
				% (recon_name_extra, k))
		nin -= 3
		f.write(" "*nin+"end do\n")
		nin -= 3
		f.write(" "*nin+"else\n")
		nin += 3

	f.write(" "*nin+"do i=1,%d\n" % dim)
	nin += 3
	for j in range(1, k+1):
		f.write(" "*nin+"Q(indices(%d)) = %s(i,%d)\n" % (j,rname,j))
	f.write(" "*nin+"xnum(i) = numeval(Q, mu2) - %s(Q, coeffs, %d)\n"
			% (recon_name, k-1))
	nin -= 3
	f.write(" "*nin+"end do\n")
	if recon_name_extra is not None:
		nin -= 3
		f.write(" "*nin+"end if\n")
	f.write(" "*nin+"coeffs%%c%d(idx,:) = matmul(%s,xnum)\n" % (k, lname))
	f.write(" "*indent+"end subroutine %s\n" % name)

def write_function_recon(f, indent, name, k, dic, d, qtype = "real"):
	"""
	Writes a function which computes N(q) from the reconstructed
	coefficients for a given value of q.

	PARAMETER

	f             -- file open for writing
	indent        -- indentation level
	name          -- name of the function to be written
	k             -- number of non-zero components of q
	dic           -- dictionary mapping powers of components of q to
	                 indices in the coefficient arrays
	d             -- number of space-time dimensions
	qtype         -- type of q (real or complex) in this routine.
	"""
	dim = len(dic)
	for line in DOC.function_recon(name,k,d-1,qtype,dim):
		f.write(" "*indent + line + "\n")
	f.write(" "*indent+"pure function %s(Q, indices, coeffs)\n" % name)
	nin = indent + 3
	f.write(" "*nin+"! generated by: write_function_recon\n")
	f.write(" "*nin+"implicit none\n")
	f.write(" "*nin+"integer, dimension(%d), intent(in) :: indices\n" % k)
	f.write(" "*nin+"complex(ki), dimension(%d), intent(in) :: coeffs\n" % dim)
	f.write(" "*nin+"%s(ki), dimension(0:%d), intent(in) :: Q\n"
			% (qtype, (d-1)))
	f.write(" "*nin+"complex(ki) :: %s\n" % name)
	for i in range(k):
		f.write(" "*nin+"%s(ki) :: q%d\n" % (qtype, i))
	bin = [tuple(lst) for lst in binary(k)]
	lines = []
	regs = []
	calculate_polynomial(lines, regs, nin, bin[1:], dic)
	regs = set(regs)
	for reg in regs:
		f.write(" "*nin+"complex(ki) :: %s\n" % reg)
	f.write(" "*nin+"complex(ki) :: acc\n")

	for i in range(k):
		f.write(" "*nin+"q%d = Q(indices(%d))\n" % (i, i+1))

	for line in lines:
		f.write(line + "\n")

	f.write(" "*nin+"%s = acc\n" % name)
	f.write(" "*indent+"end  function %s\n" % name)

def binary(k):
	"""
	Yields the binary representations of the numbers up to 2^k:

	>>> [x for x in binary(3)]
	[[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]
	"""
	if k == 0:
		yield []
	else:
		for l in binary(k-1):
			yield l + [0]
			yield l + [1]

def poly_includes(monomial, exponents):
	"""
	Checks if a monomial contains at least the given set of powers
	of variables.

	EXAMPLE

	>>> m = [4, 0, 1]
	>>> poly_includes(m, [3, 0, 0])
	True
	>>> poly_includes(m, [3, 1, 1])
	False
	"""
	for p, q in zip(monomial, exponents):
		if p < q:
			return False
	return True

def poly_subtract(monomial, exponents):
	"""
	Divides out powers of variables from a monomial

	EXAMPLE

	>>> m = [4, 0, 1]
	>>> poly_subtract(m, [3, 0, 0])
	(1, 0, 1)
	>>> poly_subtract(m, [1, 0, 1])
	(3, 0, 0)
	"""
	result = []
	for p, q in zip(monomial, exponents):
		result.append(p - q)
	return tuple(result)

def format_monomial(exponents, coeff):
	"""
	Formats a monomial for printing in Fortran style.

	EXAMPLE

	>>> m = [4, 0, 1]
	>>> format_monomial(m, "5.0")
	'5.0*q0**4*q2'
	"""
	n = len(exponents)
	factors = [coeff]
	for i in range(n):
		if exponents[i] == 1:
			factors.append("q%d" % i)
		elif exponents[i] != 0:
			factors.append("q%d**%s" % (i, exponents[i]))
	return "*".join(factors)

def calculate_polynomial(lines, regs, indent, bin, dic, path=[]):
	"""
	Produces a program fragment to compute a polynomial using
	a multivariate Horner scheme

	PARAMETER

	lines  -- list of strings, used as an output buffer
	regs   -- list of temporary variables introduced by the algorithm
	indent -- indentation level
	bin    -- list of all possibilities selecting variables
	dic    -- dictionary mapping between indices of the coefficients
	          and powers of q_i
	path   -- location of the function call in the tree of recursive calls
	"""
	l = len(bin)
	count = [0 for i in range(l)]
	max_count = 0
	max_pos = 0
	for exponents in dic.iterkeys():
		for i in range(l):
			if poly_includes(exponents, bin[i]):
				count[i] += 1
				if count[i] >= max_count:
					max_count = count[i]
					max_pos = i
	pivot = bin[max_pos]

	if max_count <= 1:
		first = True
		col = indent + 6
		line = ""
		for exponents, coeff in dic.iteritems():
			term = format_monomial(exponents, "coeffs(%d)" % (coeff+1))
			lt = len(term)
			if first:
				line =  " "*indent+"acc = "+term
				col = indent + 6 + lt
				first = False
			elif col + lt + 3 <= F90_WIDTH:
				line = line + " + " + term
				col += lt + 3
			else:
				lines.append(line)
				line = (" "*indent) + "acc = acc + " + term
				col = indent + 12 + lt
		if not first:
			lines.append(line)
	else:
		A1 = {}
		A0 = {}

		for exponents, coeff in dic.iteritems():
			if poly_includes(exponents, pivot):
				A1[poly_subtract(exponents, pivot)] = coeff
			else:
				A0[exponents] = coeff

		if len(A1) > 0 and len(A0) > 0:
			calculate_polynomial(lines, regs, indent, bin, A1, [1]+path)
			reg = "reg%d" % len(path)
			regs.append(reg)
			term = format_monomial(pivot, "acc")
			lines.append(" "*indent + "%s = %s" % (reg, term))
			calculate_polynomial(lines, regs, indent, bin, A0, [0]+path)
			lines.append(" "*indent + "acc = acc + %s" % reg)
		elif len(A1) > 0:
			calculate_polynomial(lines, regs, indent, bin, A1, [1]+path)
			term = format_monomial(pivot, "acc")
			lines.append(" "*indent+"acc = %s" % term)
		elif len(A0) > 0:
			calculate_polynomial(lines, regs, indent, bin, A0, [0]+path)

def write_coeff_type(f, indent, name, R, block_info, d):
	"""
	Writes the definition of the derived type for a certain set of
	tensor coefficients to a file.

	PARAMETER

	f          -- file object open for writing
	indent     -- indentation level
	name       -- name of the derived type in Fortran
	R          -- rank of the tensor integral
	block_info -- information about the block structure of the system
	              and the corresponding tensor coefficients
	d          -- number of space time dimensions
	"""
	for line in DOC.coeff_type(name,R,min(R,d)):
		f.write(" "*indent + line + "\n")
	f.write(" "*indent+"type %s\n" % name)
	nin = indent + 3
	f.write(" "*nin + "complex(ki) :: c0\n")
	count = 1
	for k in range(1,min(R,d)+1):
		dim = combinat(d, k)
		lst, dic, lname, rname = block_info[(R,k)]
		dimk = len(lst)
		f.write(" "*nin + "complex(ki), dimension(%d,%d) :: c%d\n"
				% (dim, dimk, k))
		if TENS_DEBUG:
			lab = 0
			for indices in select(range(d), k):
				lab += 1
				for j in range(len(lst)):
					coeff = zip([[e] for e in indices],  [e for e in lst[j]])
					coeff = map(lambda (x, y): x*y, coeff)
					coeff = reduce(lambda x, y: x+y, coeff, [])

					perms = momenta_symm(coeff) - 1
					if perms == 0:
						perms = ""
					elif perms == 1:
						perms = " + 1 permutation"
					else:
						perms = " + %d permutations" % perms

					f.write(" "*nin+"! c%d(%d,%d) = %s%s\n"
							% (k, lab, j+1, coeff, perms))
		count += dim * dimk
	if TENS_DEBUG:
		if count == 1:
			entries = "entry"
		else:
			entries = "entries"
		f.write(" "*nin + "! This record has %d %s.\n" % (count, entries))
	f.write(" "*indent+"end type %s\n" % name)

def select(items, k):
	"""
	Iterator over all selections of k elements from a given list.

	PARAMETER

	items  --  list of elements to choose from (no repetitions)
	k      --  number of elements to select.
	"""
	n = len(items)
	# We use the fact that
	# (n choose k) = (1 choose 1)(n-1 choose k-1)+(1 choose 0)(n-1 choose k)
	if k == n:
		yield items[:]
	elif k == 0:
		yield []
	elif 0 < k and k < n:
		head = items[0:1]
		tail = items[1:]
		for result in select(tail, k-1):
			yield head + result

		for result in select(tail, k):
			yield result

def write_function_glob_recon(f, indent, name, type_name, proc_names, R, d):
	"""
	Writes a function that retrieves the value of N(Q) from
	reconstructed coefficients for a given real Q.
	"""
	for line in DOC.function_glob_recon(name,type_name,R,d-1):
		f.write(" "*indent + line + "\n")
	f.write(" "*indent + "pure function %s(Q, coeffs, max_k)\n" % name)
	nin = indent + 3
	f.write(" "*nin + "! generated by: write_function_glob_recon\n")
	f.write(" "*nin + "implicit none\n")
	f.write(" "*nin + "real(ki), dimension(0:%d), intent(in) :: Q\n" % (d-1))
	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % type_name)
	f.write(" "*nin + "integer, intent(in), optional :: max_k\n")
	f.write(" "*nin + "complex(ki) :: %s\n" % name)
	f.write(" "*nin + "integer :: maxk\n")
	f.write(" "*nin + "if (present(max_k)) then\n")
	f.write(" "*(nin+3) + "maxk = max_k\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "maxk = %d\n" % min(R, d))
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "%s = coeffs%%c0\n" % name)
	for k in range(1, min(R, d) + 1):
		f.write(" "*nin + "if (%d .le. maxk) then\n" % k)
		nin += 3

		lab = 0
		for indices in select(range(d), k):
			lab += 1
			f.write(" "*nin + "%s = %s + %s(Q, (/%s/), coeffs%%c%d(%d,:))\n"
					% (name, name, proc_names[(R,k)],
						",".join(map(str,indices)), k, lab))

		nin -= 3
		f.write(" "*nin + "end if\n")
	f.write(" "*indent + "end  function %s\n" % name)

def write_function_glob_recon_complex(f, indent, name,
		type_name, proc_names, R, d):
	"""
	Writes a function that retrieves the value of N(Q) from
	reconstructed coefficients for a given complex Q.
	"""
	for line in DOC.function_glob_recon_complex(name,type_name,R,d-1):
		f.write(" "*indent + line + "\n")
	f.write(" "*indent + "pure function %s(Q, coeffs)\n" % name)
	nin = indent + 3
	f.write(" "*nin + "! generated by: write_function_glob_recon_complex\n")
	f.write(" "*nin + "implicit none\n")
	f.write(" "*nin + "complex(ki), dimension(0:%d), intent(in) :: Q\n" % (d-1))
	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % type_name)
	f.write(" "*nin + "complex(ki) :: %s\n" % name)
	f.write(" "*nin + "%s = coeffs%%c0\n" % name)
	for k in range(1, min(R, d) + 1):
		lab = 0
		for indices in select(range(d), k):
			lab += 1
			f.write(" "*nin + "%s = %s + %s(Q, (/%s/), coeffs%%c%d(%d,:))\n"
					% (name, name, proc_names[(R,k)],
						",".join(map(str,indices)), k, lab))
	f.write(" "*indent + "end  function %s\n" % name)

def write_subroutine_glob_solve(f, indent, name, type_name, proc_names, R, d,
		extra_type_name=None):
	"""
	Writes a subroutine for the determination of the tensor coefficients 
	at a fixed value of mu^2 and for given maximum rank

	PARAMETER

	f            -- file open for writing
	indent       -- indentation level
	type_name    -- name of the type to be determined
	proc_names   -- lookup for the procedure names used for the determination
	                of the single parts of the problem.
	R            -- maximum rank
	d            -- number of space time dimensions
	extra_type_name -- coefficients to be subtracted

	"""
	if extra_type_name is not None:
		xtra = ", coeffs2"
	else:
		xtra = ""
	for line in DOC.subroutine_glob_solve(name,type_name,R,d-1,xtra):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent + "subroutine     %s(numeval, mu2, coeffs%s)\n"
			% (name, xtra))
	nin = indent + 3
	f.write(" "*nin + "! generated by: write_subroutine_glob_solve\n")
	f.write(" "*nin + "implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin + "real(ki), intent(in) :: mu2\n")
	f.write(" "*nin + "type(%s), intent(inout) :: coeffs\n" % type_name)
	if extra_type_name is not None:
		f.write(" "*nin + "type(%s), intent(in), optional :: coeffs2\n"
				% extra_type_name)
	if extra_type_name is not None:
		f.write(" "*nin + "if (present(coeffs2)) then\n")
		nin += 3
		f.write(" "*nin + "coeffs%c0 = numeval(null_vec, mu2) - coeffs2%c0\n")
		for k in range(1, min(d, R) + 1):
			lab = 0
			for indices in select(range(d), k):
				lab += 1
				f.write(" "*nin
						+ "call %s(numeval, (/%s/), mu2, coeffs, %d, coeffs2)\n" %
						(proc_names[(R, k)], ",".join(map(str, indices)), lab))
		nin -= 3
		f.write(" "*nin + "else\n")
		nin += 3
	f.write(" "*nin + "coeffs%%c0 = numeval((/%s/), mu2)\n"
			% ",".join(["0.0_ki"]*d))
	for k in range(1, min(d, R) + 1):
		lab = 0
		for indices in select(range(d), k):
			lab += 1
			f.write(" "*nin+"call %s(numeval, (/%s/), mu2, coeffs, %d)\n" %
					(proc_names[(R, k)], ",".join(map(str, indices)), lab))

	if extra_type_name is not None:
		nin -= 3
		f.write(" "*nin + "end if\n")
	f.write(" "*indent + "end subroutine %s\n" % name)

def write_print_coeffs(f, indent, name, type_name, block_info, R, d):
	"""
	Writes the code to print out a set of coefficients in human readable
	form.

	PARAMETER

	f           -- file open for writing
	indent      -- intentation level
	type_name   -- name of the type that contains the coefficients to print
	block_info  -- information about the block structure of the system
	R           -- maximum rank
	d           -- number of space time dimensions
	"""

	def format_power(pow):
		b, e = pow

		if e == 1:
			return str(b)
		else:
			return "%s^%d" % (b, e)

	for line in DOC.subroutine_print_coeffs(name,R,type_name):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent+"subroutine print_coeffs_%d(coeffs, unit)\n" % R)
	nin = indent + 3
	f.write(" "*nin+"! generated by: write_print_coeffs\n")
	f.write(" "*nin+"implicit none\n")
	f.write(" "*nin+"type(%s), intent(in) :: coeffs\n" % type_name)
	f.write(" "*nin+"integer, intent(in), optional :: unit\n")
	f.write(" "*nin+"integer :: ch\n")
	f.write(" "*nin+"if (present(unit)) then\n")
	f.write(" "*(nin+3)+"ch = unit\n")
	f.write(" "*nin+"else\n")
	f.write(" "*(nin+3)+"ch = 6\n")
	f.write(" "*nin+"end if\n")

	f.write(" "*nin
		+"write(ch,'(A4,G24.16,1x,G24.16,A1)') '   (', coeffs%c0, ')'\n")
	for k in range(1,min(R,d)+1):
		lst, dic, lname, rname = block_info[(R,k)]
		dim = len(lst)
		lab = 0
		for indices in select(range(d), k):
			lab += 1
			sindices = map(lambda i: "q(%d)" % i, indices)
			for i in range(dim):
				tail = "*".join([")"] + map(format_power, zip(sindices, lst[i])))
				l = len(tail)
				value = "coeffs%%c%d(%d,%d)" % (k, lab, i+1)
				f.write(" "*nin
					+ ("write(ch,'(A4,G24.16,1x,G24.16,A%d)')" % l)
					+ (" ' + (',  %s, '%s'\n" % (value, tail)))

	f.write(" "*indent+"end subroutine print_coeffs_%d\n" % R)

def write_subroutine_reconstruct(f, indent, name, R, d,
		solve_name,	type_name,
		solve_name_extra, type_name_extra):
	"""
	Writes a routine that determines all coefficients of an integral of
	maximum rank R, including the mu2 bits.

	PARAMETER

	f                 -- file open for writing
	indent            -- indentation level
	name              -- name of the subroutine
	R                 -- maximum rank
	d                 -- number of space time dimensions
	solve_name        -- function name used to solve for the mu2=0
	                     bits of the integral
	type_name         -- name of the type of the coefficients of the
	                     'normal' integral
	solve_name_extra  -- functions for solving for the mu2, mu2^2 and mu2^3 integrals
	type_name_extra   -- name of the type of the mu2, mu2^2 and mu2^3 coefficients

	"""
	# For lower rank numerators we use the fact that N(q, mu2) is at most quadratic
	# in mu2. After
	# subtracting N(q, 0) we have the form N(q, mu2) = a mu2^2 + b mu2.
	# N(q, +1) = a + b
	# N(q, -1) = a - b
	# 2 a = N(q, +1) + N(q, -1)
	# 2 b = N(q, +1) - N(q, -1)
	#
	# For higher rank numerators (rank-6 pentagons), we have the form
	# N~(q, mu2) = a mu2^3 + b mu2^2 + c mu^2 with N~(q, mu2)= N(q, mu2) - N(q, 0).
	# N~(q, +1) = a + b + c
	# N~(q, -1) = -a + b - c
	# N~(q, 2) = 8 a + 4 b + 2 c
	#   => a = -1/2 N~(q, +1) - 1/6 N~(q, -1) + 1/6 N~(q, 2)
	#      b = 1/2  N~(q, +1) + 1/2 N~(q, -1)
	#      c =      N~(q, +1) - 1/3 N~(q, -1) - 1/6 N~(q, 2)
	if type_name_extra is None:
		t = "complex(ki)"
	else:
		t = "type(%s)" % type_name_extra

	if R>=6:
		for line in DOC.subroutine_reconstruct_extended(name,type_name,t,R):
			f.write(" "*indent+line+"\n")
	else:
		for line in DOC.subroutine_reconstruct(name,type_name,t,R):
			f.write(" "*indent+line+"\n")
	if R>=6:
		f.write(" "*indent + "subroutine     %s(numeval, cm0, cm1, cm2, cm3)\n" % name)
	else:
		f.write(" "*indent + "subroutine     %s(numeval, cm0, cm1, cm2)\n" % name)
	nin = indent + 3
	f.write(" "*nin+"! generated by: write_subroutine_reconstruct\n")
	f.write(" "*nin + "implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin + "type(%s), intent(out) :: cm0\n" % type_name)

	f.write(" "*nin + "%s, intent(out), optional :: cm1\n" % t)
	f.write(" "*nin + "%s, intent(out), optional :: cm2\n" % t)
	if R>=6:
		f.write(" "*nin + "%s, intent(out), optional :: cm3\n" % t)
	f.write(" "*nin + "%s :: ca, cb\n" % t)
	if R>=6:
		f.write(" "*nin + "%s :: cc\n" % t)
	f.write(" "*nin + "call %s(numeval, 0.0_ki, cm0)\n" % solve_name)
	f.write(" "*nin + "if (present(cm1)) then\n")
	nin += 3
	f.write(" "*nin + "if (present(cm2)) then\n")
	nin += 3
	if R>=6:
		f.write(" "*nin + "if (present(cm3)) then\n")
		nin += 3
		if type_name_extra is None:
			f.write(" "*nin + "ca = numeval(null_vec, +1.0_ki) - cm0\n")
			f.write(" "*nin + "cb = numeval(null_vec, -1.0_ki) - cm0\n")
			f.write(" "*nin + "cc = numeval(null_vec, +2.0_ki) - cm0\n")
			f.write(" "*nin + "cm3= -0.5_ki*ca+(cc-cb)/6._ki\n")
			f.write(" "*nin + "cm2= 0.5_ki*(ca+cb)\n")
			f.write(" "*nin + "cm1= ca - cb/3._ki - cc/6.ki\n")
		else:
			f.write(" "*nin + "call %s(numeval, +1.0_ki, ca, cm0)\n"
					% solve_name_extra)
			f.write(" "*nin + "call %s(numeval, -1.0_ki, cb, cm0)\n"
					% solve_name_extra)
			f.write(" "*nin + "call %s(numeval, +2.0_ki, cc, cm0)\n"
					% solve_name_extra)
			f.write(" "*nin + "cm3%c0= -0.5_ki*ca%c0+(cc%c0-cb%c0)/6._ki\n")
			f.write(" "*nin + "cm2%c0= 0.5_ki*(ca%c0+cb%c0)\n")
			f.write(" "*nin + "cm1%c0= ca%c0 - cb%c0/3._ki - cc%c0/6._ki\n")

			for k in range(1,min(d,R-2)+1):
				f.write(" "*nin + "cm3%%c%d= -0.5_ki*ca%%c%d+(cc%%c%d-cb%%c%d)/6._ki\n"
						% (k,k,k,k))
				f.write(" "*nin + "cm2%%c%d= 0.5_ki*(ca%%c%d+cb%%c%d)\n"
						% (k,k,k))
				f.write(" "*nin + "cm1%%c%d= ca%%c%d - cb%%c%d/3._ki - cc%%c%d/6._ki\n"
						% (k,k,k,k))
		nin -= 3
		f.write(" "*nin + "else\n")
		nin += 3

	if type_name_extra is None:
		f.write(" "*nin + "ca = numeval(null_vec, +1.0_ki) - cm0%c0\n")
		f.write(" "*nin + "cb = numeval(null_vec, -1.0_ki) - cm0%c0\n")
		f.write(" "*nin + "cm1= 0.5_ki * (ca - cb)\n")
		f.write(" "*nin + "cm2= 0.5_ki * (ca + cb)\n")
	else:
		f.write(" "*nin + "call %s(numeval, +1.0_ki, ca, cm0)\n"
				% solve_name_extra)
		f.write(" "*nin + "call %s(numeval, -1.0_ki, cb, cm0)\n"
				% solve_name_extra)
		f.write(" "*nin + "cm1%c0= 0.5_ki * (ca%c0 - cb%c0)\n")
		f.write(" "*nin + "cm2%c0= 0.5_ki * (ca%c0 + cb%c0)\n")
		for k in range(1,min(d,R-2)+1):
			f.write(" "*nin + "cm1%%c%d = 0.5_ki * (ca%%c%d - cb%%c%d)\n"
					% (k,k,k))
			f.write(" "*nin + "cm2%%c%d = 0.5_ki * (ca%%c%d + cb%%c%d)\n"
					% (k,k,k))
	if R>=6:
		nin -= 3
		f.write(" "*nin + "end if\n")
	nin -= 3
	f.write(" "*nin + "else\n")
	nin += 3
	if type_name_extra is None:
		f.write(" "*nin + "cm1 = numeval(null_vec, +1.0_ki) - cm0%c0\n")
	else:
		f.write(" "*nin
				+ "call %s(numeval, +1.0_ki, cm1, cm0)\n" % solve_name_extra)
	nin -= 3
	f.write(" "*nin + "end if\n")
	nin -= 3
	f.write(" "*nin + "end if\n")
	f.write(" "*indent + "end subroutine %s\n" % name)

def write_subroutine_reconstruct_dummy(f, indent, name, R, d,
		solve_name,	type_name):
	"""
	Writes a routine that determines all coefficients of an integral of
	maximum rank R if the rank is smaller than 2.

	PARAMETER

	f                 -- file open for writing
	indent            -- indentation level
	name              -- name of the subroutine
	R                 -- maximum rank
	d                 -- number of space time dimensions
	solve_name        -- function name used to solve for the mu2=0
	                     bits of the integral
	type_name         -- name of the type of the coefficients of the
	                     'normal' integral
	"""
	for line in DOC.subroutine_reconstruct_dummy(name,type_name,R):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent + "subroutine     %s(numeval, cm0)\n" % name)
	nin = indent + 3
	f.write(" "*nin+"! generated by: write_subroutine_reconstruct_dummy\n")
	f.write(" "*nin + "implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin + "type(%s), intent(out) :: cm0\n" % type_name)
	f.write(" "*nin + "call %s(numeval, 0.0_ki, cm0)\n" % solve_name)
	f.write(" "*indent + "end subroutine %s\n" % name)


def write_module_solve(name, max_rank, d=4):
	"""
	Writes the module which is usually called 'tens_rec.f90'.

	PARAMETER

	name     -- the name of the module
	max_rank -- the maximum rank that should be implemented
	d        -- number of space time dimensions
	"""
	f = open("%s.f90" % name, "w")
	for line in DOC.module_solve():
		f.write(line + "\n")
	f.write("module %s\n" % name)
	f.write("use %s, only: ki\n" % precision_golem)
	f.write("implicit none\n")
	f.write("private :: ki\n")

	f.write("real(ki), dimension(0:%d), parameter, private :: null_vec = &\n"
			% (d-1))
	f.write(" & (/%s/)\n" % ",".join(["0.0_ki"]*d))
	block_info = {}

	for R in range(1,max_rank+1):
		for k in range(1,min(R,d)+1):
			LHS, RHS, lst, dic = generate_equations(R, k)
			lname = "mat%d_%d" % (R, k)
			rname = "q%d_%d" % (R, k)
			write_matrix_f90(f, 0, lname, rname, k, LHS, RHS)
			block_info[(R,k)] = (lst, dic, lname, rname)

	for R in range(1,max_rank+1):
		write_coeff_type(f, 0, PAT['coefftype'] % R, R, block_info, d)

	f.write("interface print_coeffs\n")
	for R in range(1,max_rank+1):
		f.write("   module procedure print_coeffs_%d\n" % R)
	f.write("end interface\n")

	f.write("contains\n")

	solve_names = {}
	recon_names = {}
	crecon_names = {}
	for R in range(1,max_rank+1):
		for k in range(1,min(R,d)+1):
			lst, dic, lname, rname = block_info[(R,k)]
			proc_name = PAT['solve'] % (R, k)
			solve_names[(R, k)] = proc_name
			if R <= max_rank - 2:
				extra_type = PAT['coefftype'] % (R + 2)
				extra_recon = PAT['globtenseval'] % (R + 2)
			else:
				extra_type = None
				extra_recon = None
			write_subroutine_solve(f, 0, proc_name, k, len(lst), lname, rname,
					PAT['coefftype'] % R, PAT['globtenseval'] % R, d,
					extra_recon, extra_type)
			proc_name = PAT['tenseval'] % (R, k)
			recon_names[(R,k)] = proc_name
			cproc_name = PAT['ctenseval'] % (R, k)
			crecon_names[(R,k)] = cproc_name
			write_function_recon(f, 0, proc_name, k, dic, d)
			write_function_recon(f, 0, cproc_name, k, dic, d, "complex")

		type_name = PAT['coefftype'] % R
		proc_name = PAT['globsolve'] % R
		if R <= max_rank - 2:
			extra_type_name = PAT['coefftype'] % (R+2)
		else:
			extra_type_name = None
		write_subroutine_glob_solve(f, 0, proc_name, type_name, solve_names,
				R, d, extra_type_name)
		proc_name = PAT['globtenseval'] % R
		cproc_name = PAT['cglobtenseval'] % R
		write_function_glob_recon(f, 0, proc_name, type_name, recon_names, R, d)
		write_function_glob_recon_complex(f, 0, cproc_name, type_name, crecon_names, R, d)
		write_print_coeffs(f, 0, PAT['printcoeffs'] % R, type_name, block_info,
				R, d)
		if R > 2:
			write_subroutine_reconstruct(f, 0, PAT['reconstruct'] % R, R, d,
				PAT['globsolve'] % R, PAT['coefftype'] % R,
				PAT['globsolve'] % (R-2), PAT['coefftype'] % (R-2))
		elif R == 2:
			write_subroutine_reconstruct(f, 0, PAT['reconstruct'] % R, R, d,
				PAT['globsolve'] % R, PAT['coefftype'] % R,
				None, None)
		else:
			write_subroutine_reconstruct_dummy(f, 0, PAT['reconstruct'] % R, R, d,
				PAT['globsolve'] % R, PAT['coefftype'] % R)

	f.write("end module %s\n" % name)
	f.close()

def write_list(f, prefix, lst):
	"""
	Utility routine for writing a comma seperated list to a file inserting
	line breaks and continuation characters in the appropriate places.

	PARAMETER

	f      -- file open for writing
	prefix -- whatever should be inserted in the beginning of the first line
	lst    -- the list to be written
	"""
	is_first = True
	s = prefix
	ofs = 1 + len(s)
	f.write(s)
	for s in lst:
		if is_first:
			is_first = False
		else:
			f.write(",")
			ofs += 1

		# len(" " + s + ", &")
		if ofs + len(s) + 4 >= F90_WIDTH:
			indent = "   &"
			f.write(" &\n" + indent)
			ofs = 1 + len(indent)
		f.write(" " + s)
		ofs += len(s) + 1
	f.write("\n")

def write_import_list(f, prefix):
	"""
	Utility for writing a very specific import statement.
	"""
	lst = []
	for R in range(1, max_rank+1):
		lst.append(PAT['coefftype'] % R)
		lst.append(PAT['reconstruct'] % R)
	write_list(f, prefix, lst)

def write_function_evaluate(f, indent, name, max_legs, max_rank, d):
	"""
	Writes to functions to a file which serve as the main interface
	for the tensorial reconstruction interface.

	PARAMETER

	f         -- file open for writing
	indent    -- indentation level
	name      -- the base name used for the function
	max_legs  -- the maximum number of legs to be considered/implemented
	max_rank  -- the maximum rank to be considered/implemented
	"""
	nin = indent + 3

	# function evaluate_s
	for line in DOC.function_evaluate(name + "_s", "set", "integer array",d-1):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent
		+ "function     %s_s(numeval, momenta, set, rank) result(amp)\n"
			% name)
	f.write(" "*nin+"! generated by: write_function_evaluate\n")
	f.write(" "*nin + "implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin + "real(ki), dimension(:,:), intent(in) :: momenta\n")
	f.write(" "*nin + "integer, dimension(:), intent(in) :: set\n")
	f.write(" "*nin + "integer, intent(in), optional :: rank\n")
	f.write(" "*nin + "type(form_factor) :: amp\n")
	f.write(" "*nin + "if (present(rank)) then\n")
	f.write(" "*(nin+3) + "amp = %s_b(numeval, momenta, packb(set), rank)\n"
			% name)
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "amp = %s_b(numeval, momenta, packb(set))\n" % name)
	f.write(" "*nin + "end if\n")
	f.write(" "*indent + "end function %s_s\n" % name)

	# function evaluate_b
	for line in DOC.function_evaluate(name + "_b", 
			"b_set", "integer (bit-set)",d-1):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent
		+ "function     %s_b(numeval, momenta, b_set, rank) result(amp)\n"
			% name)
	f.write(" "*nin + "implicit none\n")
	write_numeval_interface(f, nin, "numeval", d)
	f.write(" "*nin + "real(ki), dimension(:,:), intent(in) :: momenta\n")
	f.write(" "*nin + "integer, intent(in) :: b_set\n")
	f.write(" "*nin + "integer, intent(in), optional :: rank\n")
	f.write(" "*nin + "type(form_factor) :: amp\n")
	f.write(" "*nin + "integer :: N, r\n")
	f.write(" "*nin + "complex(ki) :: coeffs0\n")
	for r in range(1, min(max_legs+extra_rank,max_rank)+1):
		if r == 2:
			f.write(" "*nin + "type(%s) :: coeffs2, coeffs2x, coeffs2xx\n" %
				(PAT['coefftype'] % 2))
		else:
			f.write(" "*nin + "type(%s) :: coeffs%d, coeffs%dx, coeffs%dxx\n" %
				(PAT['coefftype'] % r, r, r, r))
	#f.write(" "*nin + "amp = 0.0_ki\n")
	f.write(" "*nin + "N = size(momenta,1) - countb(b_set)\n")
	f.write(" "*nin + "if (present(rank)) then\n")
	f.write(" "*(nin+3) + "r = rank\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "r = N\n")
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "select case(N)\n")
	nin += 3
	for N in range(1,max_legs+1):
		f.write(" "*(nin-3) + "case(%d)\n" % N)
		f.write(" "*nin  + "select case(r)\n")
		nin += 3
		for R in ( range(min(N+extra_rank,max_rank)+1) \
				if N <= max_leg_extra_rank else range(min(N,max_rank)+1)):
			f.write(" "*(nin-3) + "case(%d)\n" % R)

			if N <= 4:
				max_alpha = R/2
			elif N == 5 and R >= 6:
				max_alpha = R/2
			else:
				max_alpha = 0

			if R == 0:
				f.write(" "*nin + "coeffs0 = numeval(null_vec, 0.0_ki)\n")
			else:
				if max_alpha == 0:
					f.write(" "*nin + "call %s(numeval, coeffs%d)\n"
							% (PAT['reconstruct'] % R, R))
				elif max_alpha == 1:
					f.write(" "*nin + "call %s(numeval, coeffs%d, coeffs%d)\n"
							% (PAT['reconstruct'] % R, R, R-2))
				elif max_alpha == 2:
					f.write(" "*nin
							+ "call %s(numeval, coeffs%d, coeffs%d, coeffs%dx)\n"
								% (PAT['reconstruct'] % R, R, R-2, R-2))
				elif max_alpha == 3:
					f.write(" "*nin
							+ "call %s(numeval, coeffs%d, coeffs%d, coeffs%dx, coeffs%dxx)\n"
								% (PAT['reconstruct'] % R, R, R-2, R-2, R-2))
				else:
					f.write(" "*nin + "print*, \"Problem: max_alpha = %d\"\n"
							% max_alpha)
					f.write(" "*nin + "stop\n")

			if N > 1:
				if R == 0:
					f.write(" "*nin + "amp = coeffs0 * a%d0%s(b_set)\n"
							% (N, ff_suffix))
				else:
					uninitialized = True
					for alpha in range(0,max_alpha+1):
						if alpha == 0:
							proc = PAT['contract'] % (N, R)
							coeffs = "coeffs%d" % R
						elif alpha == 1:
							proc = (PAT['contract'] % (N, R)) + ("s%d" % alpha)
							coeffs = "coeffs%d" % (R - 2)
						elif alpha == 2:
							proc = (PAT['contract'] % (N, R)) + ("s%d" % alpha)
							coeffs = "coeffs%dx" % (R - 2)
						elif alpha == 3:
							proc = (PAT['contract'] % (N, R)) + ("s%d" % alpha)
							coeffs = "coeffs%dxx" % (R - 2)
						else:
							assert False, "Not yet implented: alpha=%s" % alpha

						if uninitialized:
							prev = ""
							uninitialized = False
						else:
							prev = "amp + "
						f.write(" "*nin + "amp = %s%s(%s,momenta,b_set)\n"
								% (prev, proc, coeffs))
			else:
				f.write(" "*nin + "print*, \"Tadpoles not implemented yet\"\n")
				f.write(" "*nin + "stop\n")

		f.write(" "*(nin-3) + "case default\n")
		f.write(" "*nin + "print*, \"Not yet implemented: N, r = \", %d, r\n"
				% N)
		f.write(" "*nin + "stop\n")
		nin -= 3
		f.write(" "*nin + "end select\n")
	f.write(" "*(nin-3) + "case default\n")
	f.write(" "*nin + "print*, \"Not yet implemented: N=\", N\n")
	f.write(" "*nin + "stop\n")
	nin -= 3
	f.write(" "*nin + "end select\n")
	f.write(" "*indent + "end function %s_b\n" % name)

def ff_list(N,maxrank=None):
	"""
	List all form factors for a given number of legs.
	"""
	result = []
	for r in range(N+1) if N > max_leg_extra_rank else range(N+1+extra_rank) :
		if maxrank and r > maxrank:
			break
		if not (N == 1 and r == 1):
			result.append("a%d%d" % (N, r))
		if N < 6 and r >= 2:
			result.append("b%d%d" % (N, r))
		if N < 6 and r >= 4:
			result.append("c%d%d" % (N, r))
		if N < 6 and r >= 6:
			result.append("d%d%d" % (N, r))


	return result
			
def contract_init_momenta(f, nin, N, R, d):
	"""
	Write some common definitions for both varieties of 'contractXXX'.
	"""
	f.write(" "*nin + "integer, dimension(%d) :: unpinched\n" % N)
	f.write(" "*nin+"unpinched = unpackb(pminus(b_ref, b_set), %d)\n" % N)

def write_contract_split(f, nin, N, R, d):
	"""
	Writes the 'split' variant of contracting integrals, which is used
	for N>=6.

	Implements Eq. (63) in hep-ph/0504267
	"""
	assert R > 0
	f.write(" "*nin+"! generated by: write_contract_split\n")
	f.write(" "*nin + "complex(ki), dimension(0:%d) :: C\n" % (d-1))
	if R > 1:
		f.write(" "*nin + "type(%s) :: cprime\n" % (PAT['coefftype'] % (R-1)))
	else:
		f.write(" "*nin + "complex(ki) :: cprime\n")
	f.write(" "*nin + "integer :: i, pnch, new_set\n")
	f.write(" "*nin + "integer, dimension(1) :: pnch_set\n")
	#f.write(" "*nin + "real(ki), dimension(%d,0:%d) :: mprime\n" % (N-1,d-1))
	contract_init_momenta(f, nin, N, R, d)

	f.write(" "*nin + "amp = coeffs%%c0 * a%d0%s(b_set)\n" % (N, ff_suffix))
	f.write(" "*nin + "do pnch=1,%d\n" % N)
	nin += 3

	f.write(" "*nin + "! Eq. (54) in hep-ph/0504267\n")
	f.write(" "*nin + "C(:) = 0.0_ki\n")
	if TENS_DEBUG:
		f.write(" "*nin + "write(15,'(A11)',advance='no') 'unpinched ='\n")
	f.write(" "*nin + "do i=1,%d\n" % N)
	nin += 3
	f.write(" "*nin + "C(:) = C(:) " +
		"+ inv_s(unpinched(pnch),unpinched(i),b_set) * &\n")
	f.write(" "*nin + "     & momenta(unpinched(i),:)\n")
	#f.write(" "*nin + "if (i .lt. pnch) then\n")
	#f.write(" "*(nin+3) + "mprime(i,:) = mom(i,:)\n")
	#if TENS_DEBUG:
	#	f.write(" "*(nin+3) + "write(15,'(I2)',advance='no') unpinched(i)\n")
	#f.write(" "*nin + "else if (i .gt. pnch) then\n")
	#f.write(" "*(nin+3) + "mprime(i-1,:) = mom(i,:)\n")
	#if TENS_DEBUG:
	#	f.write(" "*(nin+3) + "write(15,'(I2)',advance='no') unpinched(i)\n")
	#f.write(" "*nin + "end if\n")
	nin -= 3
	f.write(" "*nin + "end do\n")
	if TENS_DEBUG:
		f.write(" "*nin + "write(15,'(A1)') '.'\n")
	f.write(" "*nin + "! Eq. (63) in hep-ph/0504267\n")
	f.write(" "*nin + "pnch_set(1) = pnch\n")
	f.write(" "*nin + "new_set = punion(packb(pnch_set),b_set)\n")

	coeff_lsts = []
	coeff_by_rank = [{} for i in range(R+1)]

	for k in range(1, min(R,d)+1):
		lst, dic = generate_mapping(R, k)
		coeff_lsts.append(lst)
		i = 0
		for e in lst:
			i += 1
			rk = sum(e)
			coeff_by_rank[rk][e] = (k, i)

	reduced_coeffs = {}
	initialized_coeffs = {}
	for k in range(1, min(R-1,d)+1):
		lst, dic = generate_mapping(R-1, k)
		lab = 0
		for indices in select(range(d), k):
			lab += 1
			idx = 0
			for e in lst:
				idx += 1
				subscript = reduce(lambda x, y: x+y,
						map(lambda (p, q): [p]*q, zip(indices, e)), [])

				subscript = tuple(subscript)
				reduced_coeffs[subscript] = (k, lab, idx)
				initialized_coeffs[subscript] = False
	initialized_c0 = False

	for rk in range(1, R+1):
		for e, pair in coeff_by_rank[rk].iteritems():
			k, idx = pair

			tot_perm = len([q for q in permutations(e)])

			lab = 0
			for indices in select(range(d), k):
				lab += 1
				subscript = reduce(lambda x, y: x+y,
						map(lambda (p, q): [p]*q, zip(indices, e)), [])

				for i in range(k):
					new_e = list(e)
					new_e[i] -= 1
					rsubscript = reduce(lambda x, y: x+y,
						map(lambda (p, q): [p]*q, zip(indices, new_e)), [])

					if new_e[i] == 0:
						del new_e[i]
					if len(new_e) == 0:
						r_perm = 1
					else:
						r_perm = len([q for q in permutations(new_e)])

					num, den = simp_fact(r_perm, tot_perm)
					if num == den:
						RHS = ""
					elif den == 1:
						RHS = "%d.0_ki * " % num
					elif num == 1 and den == 2:
						RHS = "0.5_ki * "
					elif num == 1 and den == 4:
						RHS = "0.25_ki * "
					else:
						RHS = "%d.0_ki/%d.0_ki * " % (num, den)

					RHS = "%sC(%d) * coeffs%%c%d(%d, %d)" \
								% (RHS, indices[i], k, lab, idx)
					if sum(new_e) > 0:
						rsub = tuple(rsubscript)
						r_k, r_lab, r_idx = reduced_coeffs[rsub]
						initialized = initialized_coeffs[rsub]
						initialized_coeffs[rsub] = True
						LHS = "cprime%%c%d(%d,%d)" % (r_k, r_lab, r_idx)
					elif R > 1:
						rsubscript = []
						initialized = initialized_c0
						initialized_c0 = True
						LHS = "cprime%c0"
					else:
						rsubscript = []
						initialized = initialized_c0
						initialized_c0 = True
						LHS = "cprime"
					f.write(" "*nin + "! %s <-- %s\n"
							% (rsubscript, subscript))
					if initialized:
						f.write(" "*nin + "%s = %s + %s\n" % (LHS, LHS, RHS))
					else:
						f.write(" "*nin + "%s = %s\n" % (LHS, RHS))

	if R > 1:
		f.write(" "*nin + "amp = amp - %s(cprime, momenta, new_set)\n"
				% (PAT['contract'] % (N-1,R-1)))
	else:
		f.write(" "*nin + "amp = amp - cprime * a%d0%s(new_set)\n"
				% (N-1, ff_suffix))
	nin -= 3
	f.write(" "*nin + "end do\n")

def simp_fact(num, den):
	"""
	Simplify a rational number given by num and den.
	"""
	def gcd(a, b):
		while b > 0:
			t = b
			b = a % b
			a = t
		return a

	g = gcd(num, den)
	return (num / g, den / g)

def fact(N):
	"""
	Compute the factorial of a non-negative integer number.
	"""
	if N == 0:
		return 1
	elif N > 0:
		return reduce(lambda x, y: x*y, range(2, N), N)

def momenta_symm(lst):
	"""
	Combinatorial number, unordered selection with repetition.

	In how many ways can we assign distinct indices to not necessarily
	distinct vectors.
	"""
	count_lst = []
	p_prev = -1
	count = 0
	for p in lst:
		if p == p_prev:
			count += 1
		elif p_prev == -1:
			p_prev = p
			count = 1
		else:
			count_lst += [count]
			count = 1
			p_prev = p
	count_lst += [count]
	return fact(sum(count_lst)) / reduce(lambda x, y: x*y,
			map(fact, count_lst))

def choose_momenta(R, N, i0=1):
	"""
	Yields all possibilities of selecting R momenta from a set of
	N momenta with repetitions, however, with some ordering amongst
	the vectors.

	EXAMPLE

	>>> for m in choose_momenta(3, 2):
	...	print m
	[1, 1, 1]
	[1, 1, 2]
	[1, 2, 2]
	[2, 2, 2]
	>>> for m in choose_momenta(2, 3):
	...	print m
	[1, 1]
	[1, 2]
	[1, 3]
	[2, 2]
	[2, 3]
	[3, 3]
	"""

	if R == 0:
		yield []
	else:
		for cmidx in range(i0, N+1):
			for cmlst in choose_momenta(R - 1, N, cmidx):
				yield [cmidx] + cmlst[:]

def write_contract_simple(f, nin, N, R, d, shift):
	"""
	  N     -- number of denominators
	  R     -- overall rank
	  shift -- power of mu^2 in the numerator
	"""

	def ff_required(kind, N, R):
		if shift == 0:
			if N == 1:
				return R == 0 and kind == 0
			elif N <= 5:
				return 2*kind <= R
			else:
				# No B or C form factors for N >= 6
				return kind == 0
		elif shift == 1:
			if kind == 1:
				# B form factors
				if N == 2:
					return R >= 2
				elif N == 3:
					return R >= 2
				elif N == 5:
					return R >= 6
				else:
					return False
			elif kind == 2:
				# C form factors
				return (N == 4 and R >= 4) or \
						 (N == 3 and R >= 4) or (N == 5 and R >= 6)
			else:
				# There's no A
				return False
		elif shift == 2:
			return ( N == 4 and R >= 4 and kind == 2 ) or \
		        ( N == 3 and R >= 4 and kind == 2 ) or (N == 5 and R >= 6)
		elif shift == 3:
			return (N == 5 and R >=6)
		assert(False)

	def ff(kind, N, R, *args):
		letters = ["a", "b", "c", "d"]
		ltr = letters[kind]

		return "%s%d%d%s(%s)" % (ltr, N, R, ff_suffix,
				", ".join(["l%d" % a for a in args] + ["b_set"]))

	# if shift == 1 or shift == 2 the coeff_type is reduced by 2
	# as a concequence of the reconstruction algorithm.
	if shift == 0:
		coeff_type = R
	else:
		coeff_type = R - 2

	# residual rank = overall rank - power of mu^2
	rk_res = R - 2*shift

	f.write(" "*nin+"! generated by: write_contract_simple\n")

	# special case for rank6 pentagons (mu^2, mu^4, mu^6)
	if N==5 and R==6 and shift==1:
            f.write(" "*nin+ "amp%c = 1./64._ki * ( coeffs%c1(1,4) + coeffs%c1(2,4) + &\n")
            f.write(" "*nin+ "          & coeffs%c1(3,4) + coeffs%c1(4,4) - coeffs%c2(1,5) - coeffs%c2(2,5) - coeffs%c2(3,5) &\n")
            f.write(" "*nin+ "          &  + coeffs%c2(4,5) + coeffs%c2(5,5) + coeffs%c2(6,5) )\n")
            f.write(" "*nin+ "amp%b = 0.0_ki\n")
            f.write(" "*nin+ "amp%a = 0.0_ki\n")
            return
        elif N==5 and R==6 and shift==2:
            f.write(" "*nin+ "amp%c = - 1./48._ki * ( coeffs%c1(1,2) - coeffs%c1(2,2) - &\n")
            f.write(" "*nin+ "                      & coeffs%c1(3,2) - coeffs%c1(4,2) )\n");
            f.write(" "*nin+ "amp%b = 0.0_ki\n")
            f.write(" "*nin+ "amp%a = 0.0_ki\n")
            return
        elif N==5 and R==6 and shift==3:
            f.write(" "*nin+ "amp%c =  coeffs%c0 * 1._ki/12._ki\n")
            f.write(" "*nin+ "amp%b = 0.0_ki\n")
            f.write(" "*nin+ "amp%a = 0.0_ki\n")
            return

	for r in range(1,rk_res + 1):
		f.write(" "*nin + "real(ki), dimension(%d,0:%d) :: mom%d\n" % (r, d-1, r))

	f.write(" "*nin + "integer :: %s\n"
			% ", ".join(["l%d" % i for i in range(1,N+1)]))

	f.write(" "*nin + "complex(ki) :: tmp\n")


	contract_init_momenta(f, nin, N, R, d)

	for i in range(1, N+1):
		f.write(" "*nin+"l%d = unpinched(%d)\n" % (i, i))

	uninitialized = True
	letters = ['a', 'b', 'c', 'd']

	for r in range(R+1):
		# looping over kind:
		# - lower bound = shift because there are no
		#   A's for shift>=1 and no B's for shift == 2
		# - upper bound = r/2 because B's exist only from r>=2 and 
		#   C's only for r>=4
		for kind in range(shift,r/2+1):
			# After our preselection there are still some we do not need
			# therefore an if-statement
			if ff_required(kind, N, r):
				# How many indices are there to be summed over?
				rk_indices = r - 2 * kind

				# If ff_required works as expected, rk_indices cannot
				# become negative
				assert rk_indices >= 0, "N=%d,R=%d,r=%d,kind=%d,rk_res=%d" \
						% (N, R, r, kind, rk_res)

				for indices in choose_momenta(rk_indices, N):
					ffactor = ff(kind, N, r, *indices)
					cindex = kind - shift

					if rk_indices > 0:
						f.write(" "*nin + "mom%d = momenta((/%s/),:)\n"
								% (rk_indices, ",".join(
									map(lambda x: "l%d" % x,indices))))
					if uninitialized:
						prev = ""
						uninitialized = False
					else:
						prev = "amp + "

					if coeff_type == 0:
						assert cindex == 0
						symm = 1
						ssymm = ""
						coeff = "coeffs"
					elif cindex == 0 and rk_indices == 0:
						symm = 1
						ssymm = ""
						coeff = "coeffs%c0"
					else:
						coeff = PAT['tenscontract%s' % letters[cindex]] % coeff_type

						if rk_indices > 0:
							symm = momenta_symm(indices)
							if symm == 1:
								ssymm = ""
							else:
								ssymm = "%d.0_ki * " % symm
							coeff = "%s(coeffs, mom%d)" % (coeff, rk_indices)
						else:
							coeff = "%s(coeffs)" % (coeff)
							symm = 1
							ssymm = ""

                                        s0 = " "*nin + "tmp = %s\n" % coeff
                                        f.write(s0)
                                        if prev:
                                            s0 = " "*nin + "if (test_if_nonzero(tmp)) then\n"
                                            f.write(s0)
                                            nin=nin+1

					s1 = " "*nin + "amp = " + prev + ssymm + "tmp"
					s2 = " * " + ffactor


					if len(s1) + len(s2) <= F90_WIDTH:
						f.write(s1 + s2 + "\n")
					else:
						f.write(s1 + " &\n" + " "*nin + "& " + s2 + "\n")
                                        if prev:
                                            nin=nin-1
                                            f.write(" "*nin + "end if\n")

					if TENS_DEBUG:
						f.write(" "*nin
							+"write(15,*), \"%s =\", %s\n" % (ffactor, ffactor))
						f.write(" "*nin
							+"write(15,*), \"coeff =\", %s\n" % coeff)
						f.write(" "*nin
							+"write(15,*), \"symm =\", %d\n" % symm)
						f.write(" "*nin
							+"write(15,*), \"shift =\", %d\n" % shift)

	if uninitialized:
		f.write(" "*nin+"amp = 0.0_ki\n")
	else:
		if shift == 1:
			f.write(" "*nin+"! multiply by 2*epsilon\n")
			f.write(" "*nin+"amp%c = 2.0_ki*amp%b\n")
			f.write(" "*nin+"amp%b = 2.0_ki*amp%a\n")
			f.write(" "*nin+"amp%a = 0.0_ki\n")
		elif shift == 2:
			f.write(" "*nin+"! multiply by -4*(epsilon-epsilon^2)\n")
			f.write(" "*nin+"amp%c = -4.0_ki*(amp%b-amp%a)\n")
			f.write(" "*nin+"amp%b = -4.0_ki*amp%a\n")
			f.write(" "*nin+"amp%a = 0.0_ki\n")
		elif shift == 3:
			f.write(" "*nin+"! multiply by -8*(3 epsilon^2 - 2 epsilon)\n")
			f.write(" "*nin+"amp%c = 16._ki*amp%b - 24._ki * amp%a\n")
			f.write(" "*nin+"amp%b = 16._ki*amp%a\n")
			f.write(" "*nin+"amp%a = 0.0_ki\n")

def res_accu_line(lines, i, j, k, res_first):
	if lines[-1] == "acc = -1.0_ki":
		if res_first:
			lhs = "res = -"
		else:
			lhs = "res = res - "
		lines[-1] = lhs + "coeffs%%c%d(%d,%d)" % (k, j, i)
	elif lines[-1] == "acc = 1.0_ki":
		if res_first:
			lhs = "res = "
		else:
			lhs = "res = res + "
		lines[-1] = lhs + "coeffs%%c%d(%d,%d)" % (k, j, i)
	else:
		if res_first:
			lhs = "res = "
		else:
			lhs = "res = res + "
		lines.append(lhs + "acc * coeffs%%c%d(%d,%d)" % (k, j, i))

def write_function_a_tensor_contract(f, indent, name, R, d):
	"""
	Writes a function that contracts a A-type tensor (no g^{mu nu})
	constructed from a given set of momenta with a coefficient tensor
	"""
	coeff_lsts = []
	coeff_by_rank = [{} for i in range(R+1)]

	for k in range(1, min(R,d)+1):
		lst, dic = generate_mapping(R, k)
		coeff_lsts.append(lst)
		i = 0
		for e in lst:
			i += 1
			rk = sum(e)
			coeff_by_rank[rk][e] = (k, i)

	nin = indent + 3
	for line in DOC.function_tenscontractX(name,"A",R,PAT["coefftype"]%R,d-1):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent
		+ "pure function %s(coeffs, momenta) result(res)\n" % name)
	f.write(" "*nin+"! generated by: write_function_a_tensor_contract\n")
	f.write(" "*nin + "implicit none\n")

	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % PAT['coefftype'] % R)
	f.write(" "*nin
			+ "real(ki), dimension(:,0:), intent(in), optional :: momenta\n")
	f.write(" "*nin + "complex(ki) :: res\n")
	f.write(" "*nin + "integer :: rk\n")
	f.write(" "*nin + "real(ki) :: acc\n")

	tens_lines = {}
	registers = set()

	for r in range(1, R+1):
		res_first = True
		momenta = ["momenta(%d,%%d)" % m for m in range(1, r+1)]
		lines = []
		for key, value in coeff_by_rank[r].iteritems():
			k, i = value
			mult = coeff_lsts[k-1][i-1]
			j = 0
			for idx in select(range(d), k):
				j += 1
				terms = []
				symm = 0
				for p in permutations(mult):
					tmp = map(lambda (m, mu): m % mu,
							zip(momenta, [idx[pi] for pi in p]))
					terms.append(tmp)
					symm += 1
				calculate_tensor_terms(lines, terms, registers)

				if symm > 1:
					lines.append("acc = acc / %d.0_ki" % symm)

				res_accu_line(lines, i, j, k, res_first)
				res_first = False

		tens_lines[r] = lines

	for reg in registers:
		f.write(" "*nin + "real(ki) :: %s\n" % reg)

	f.write(" "*nin + "if (present(momenta)) then\n")
	f.write(" "*(nin+3) + "rk = size(momenta, 1)\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "rk = 0\n")
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "select case(rk)\n")
	nin += 3

	for r in range(1, R+1):
		f.write(" "*(nin-3) + "case(%d)\n" % r)

		for line in tens_lines[r]:
			f.write(" "*nin + line + "\n")

	f.write(" "*(nin-3) + "case default\n")
	f.write(" "*nin + "res = 0.0_ki\n")
	nin -= 3
	f.write(" "*nin + "end select\n")
	f.write(" "*indent + "end  function %s\n" % name)

def write_function_b_tensor_contract(f, indent, name, R, d):
	"""
	Writes a function that contracts a B-type tensor (one g^{mu nu})
	constructed from a given set of momenta with a coefficient tensor
	"""
	coeff_lsts = []
	coeff_by_rank = [{} for i in range(R+1)]

	for k in range(1, min(R,d)+1):
		lst, dic = generate_mapping(R, k)
		coeff_lsts.append(lst)
		i = 0
		for e in lst:
			i += 1
			rk = sum(e)
			coeff_by_rank[rk][e] = (k, i)

	for line in DOC.function_tenscontractX(name,"B",R,PAT["coefftype"]%R,d-1):
		f.write(" "*indent+line+"\n")
	nin = indent + 3
	f.write(" "*indent
		+ "pure function %s(coeffs, momenta) result(res)\n" % name)
	f.write(" "*nin+"! generated by: write_function_b_tensor_contract\n")
	f.write(" "*nin + "implicit none\n")

	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % PAT['coefftype'] % R)
	f.write(" "*nin
			+ "real(ki), dimension(:,0:), intent(in), optional :: momenta\n")
	f.write(" "*nin + "complex(ki) :: res\n")
	f.write(" "*nin + "integer :: rk\n")
	f.write(" "*nin + "real(ki) :: acc\n")

	tens_lines = {}
	registers = set()

	for r in range(0, (R-2)+1):
		res_first = True
		momenta = ["momenta(%d,%%d)" % m for m in range(1, r+1)]
		lines = []
		for key, value in coeff_by_rank[r+2].iteritems():
			k, i = value
			mult = coeff_lsts[k-1][i-1]
			for mi in range(len(mult)):
				if mult[mi] < 2:
					continue

				gsymm = len([X for X in select(range(mult[mi]), 2)])

				new_mult = [mm for mm in mult]
				new_mult[mi] -= 2
				j = 0
				for idx in select(range(d), k):
					j += 1
					terms = []

					if idx[mi] == 0:
						g_mi_mi = 1
					else:
						g_mi_mi = -1

					symm = 0
					for p in permutations(new_mult):
						tmp = map(lambda (m, mu): m % mu,
							zip(momenta, [idx[pi] for pi in p]))
						tmp = [g_mi_mi * gsymm] + tmp
						terms.append(tmp)
						symm += 1
					calculate_tensor_terms(lines, terms, registers)

					if symm > 1:
						lines.append("acc = acc / %d.0_ki" % symm)

					res_accu_line(lines, i, j, k, res_first)
					res_first = False

		tens_lines[r] = lines

	for reg in registers:
		f.write(" "*nin + "real(ki) :: %s\n" % reg)

	f.write(" "*nin + "if (present(momenta)) then\n")
	f.write(" "*(nin+3) + "rk = size(momenta, 1)\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "rk = 0\n")
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "select case(rk)\n")
	nin += 3

	for r in range(0, (R-2)+1):
		f.write(" "*(nin-3) + "case(%d)\n" % r)

		for line in tens_lines[r]:
			f.write(" "*nin + line + "\n")

	f.write(" "*(nin-3) + "case default\n")
	f.write(" "*nin + "res = 0.0_ki\n")
	nin -= 3
	f.write(" "*nin + "end select\n")
	f.write(" "*indent + "end  function %s\n" % name)

def select_two_pairs(mult):
	"""
	Select two pairs of indices given a list of multiplicities.
	"""
	candidates = {}
	for i in range(len(mult)):
		mi = mult[i]
		if mi >= 2:
			candidates[i] = mi

	for i, mi in candidates.iteritems():
		if mi >= 4:
			yield i, i

	for lst in select(candidates.keys(), 2):
		yield lst[0], lst[1]

def select_three_pairs(mult):
	"""
	Select three pairs of indices given a list of multiplicities.
	"""
	candidates = {}
	for i in range(len(mult)):
		mi = mult[i]
		if mi >= 2:
			candidates[i] = mi

	for i, mi in candidates.iteritems():
		if mi >= 6:
			yield i, i, i

	for i, mi in candidates.iteritems():
		if mi >= 4:
			for j, mj in candidates.iteritems():
				if mj >= 2 and i != j:
					yield i, i, j

	for lst in select(candidates.keys(), 3):
		yield lst[0], lst[1], lst[2]


def write_function_c_tensor_contract(f, indent, name, R, d):
	"""
	Writes a function that contracts a C-type tensor (two g^{mu nu})
	constructed from a given set of momenta with a coefficient tensor
	"""
	coeff_lsts = []
	coeff_by_rank = [{} for i in range(R+1)]

	for k in range(1, min(R,d)+1):
		lst, dic = generate_mapping(R, k)
		coeff_lsts.append(lst)
		i = 0
		for e in lst:
			i += 1
			rk = sum(e)
			coeff_by_rank[rk][e] = (k, i)

	for line in DOC.function_tenscontractX(name,"C",R,PAT["coefftype"]%R,d-1):
		f.write(" "*indent+line+"\n")
	nin = indent + 3
	f.write(" "*indent
		+ "pure function %s(coeffs, momenta) result(res)\n" % name)
	f.write(" "*nin+"! generated by: write_function_c_tensor_contract\n")
	f.write(" "*nin + "implicit none\n")

	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % PAT['coefftype'] % R)
	f.write(" "*nin
			+ "real(ki), dimension(:,0:), intent(in), optional :: momenta\n")
	f.write(" "*nin + "complex(ki) :: res\n")
	f.write(" "*nin + "integer :: rk\n")
	f.write(" "*nin + "real(ki) :: acc\n")

	tens_lines = {}
	registers = set()

	for r in range(0, (R-4)+1):
		res_first = True
		momenta = ["momenta(%d,%%d)" % m for m in range(1, r+1)]
		lines = []
		for key, value in coeff_by_rank[r+4].iteritems():
			k, i = value
			mult = coeff_lsts[k-1][i-1]

			for mi, mj in select_two_pairs(mult):
				new_mult = [mm for mm in mult]
				new_mult[mi] -= 2
				new_mult[mj] -= 2
				if mi == mj:
					gsymm = 3 * len([X for X in select(range(mult[mi]), 4)])
				else:
					gsymm = len([X for X in select(range(mult[mi]), 2)]) * \
					        len([X for X in select(range(mult[mj]), 2)])
				j = 0
				for idx in select(range(d), k):
					j += 1
					terms = []

					if idx[mi] == 0:
						g_mi_mi = 1
					else:
						g_mi_mi = -1

					if idx[mj] == 0:
						g_mj_mj = 1
					else:
						g_mj_mj = -1

					symm = 0
					for p in permutations(new_mult):
						tmp = map(lambda (m, mu): m % mu,
							zip(momenta, [idx[pi] for pi in p]))
						tmp = [g_mi_mi*g_mj_mj*gsymm] + tmp
						terms.append(tmp)

						symm += 1

					calculate_tensor_terms(lines, terms, registers)

					if symm > 1:
						lines.append("acc = acc / %d.0_ki" % symm)

					res_accu_line(lines, i, j, k, res_first)
					res_first = False

		tens_lines[r] = lines

	for reg in registers:
		f.write(" "*nin + "real(ki) :: %s\n" % reg)

	f.write(" "*nin + "if (present(momenta)) then\n")
	f.write(" "*(nin+3) + "rk = size(momenta, 1)\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "rk = 0\n")
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "select case(rk)\n")
	nin += 3

	for r in range(0, (R-4)+1):
		f.write(" "*(nin-3) + "case(%d)\n" % r)

		for line in tens_lines[r]:
			f.write(" "*nin + line + "\n")

	f.write(" "*(nin-3) + "case default\n")
	f.write(" "*nin + "res = 0.0_ki\n")
	nin -= 3
	f.write(" "*nin + "end select\n")
	f.write(" "*indent + "end  function %s\n" % name)

def write_function_d_tensor_contract(f, indent, name, R, d):
	"""
	Writes a function that contracts a D-type tensor (three g^{mu nu})
	constructed from a given set of momenta with a coefficient tensor
	"""
	coeff_lsts = []
	coeff_by_rank = [{} for i in range(R+1)]

	for k in range(1, min(R,d)+1):
		lst, dic = generate_mapping(R, k)
		coeff_lsts.append(lst)
		i = 0
		for e in lst:
			i += 1
			rk = sum(e)
			coeff_by_rank[rk][e] = (k, i)

	for line in DOC.function_tenscontractX(name,"D",R,PAT["coefftype"]%R,d-1):
		f.write(" "*indent+line+"\n")
	nin = indent + 3
	f.write(" "*indent
		+ "pure function %s(coeffs, momenta) result(res)\n" % name)
	f.write(" "*nin+"! generated by: write_function_d_tensor_contract\n")
	f.write(" "*nin + "implicit none\n")

	f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % PAT['coefftype'] % R)
	f.write(" "*nin
			+ "real(ki), dimension(:,0:), intent(in), optional :: momenta\n")
	f.write(" "*nin + "complex(ki) :: res\n")
	f.write(" "*nin + "integer :: rk\n")
	f.write(" "*nin + "real(ki) :: acc\n")

	tens_lines = {}
	registers = set()

	for r in range(0, (R-6)+1):
		res_first = True
		momenta = ["momenta(%d,%%d)" % m for m in range(1, r+1)]
		lines = []
		for key, value in coeff_by_rank[r+6].iteritems():
			k, i = value
			mult = coeff_lsts[k-1][i-1]

			for mi, mj, mk in select_three_pairs(mult):
				new_mult = [mm for mm in mult]
				new_mult[mi] -= 2
				new_mult[mj] -= 2
				new_mult[mk] -= 2
				if mi == mj and mi == mk:
					gsymm = 15 * len([X for X in select(range(mult[mi]), 6)])
				elif mi == mj and mi != mk:
					gsymm = 3 * len([X for X in select(range(mult[mi]), 4)]) * \
					        len([X for X in select(range(mult[mk]), 2)])
				else:
					gsymm = len([X for X in select(range(mult[mi]), 2)]) * \
					        len([X for X in select(range(mult[mj]), 2)]) * \
					        len([X for X in select(range(mult[mk]), 2)])
				j = 0
				for idx in select(range(d), k):
					j += 1
					terms = []

					if idx[mi] == 0:
						g_mi_mi = 1
					else:
						g_mi_mi = -1

					if idx[mj] == 0:
						g_mj_mj = 1
					else:
						g_mj_mj = -1

					if idx[mk] == 0:
						g_mk_mk = 1
					else:
						g_mk_mk = -1


					symm = 0
					for p in permutations(new_mult):
						tmp = map(lambda (m, mu): m % mu,
							zip(momenta, [idx[pi] for pi in p]))
						tmp = [g_mi_mi*g_mj_mj*g_mk_mk*gsymm] + tmp
						terms.append(tmp)

						symm += 1

					calculate_tensor_terms(lines, terms, registers)

					if symm > 1:
						lines.append("acc = acc / %d.0_ki" % symm)

					res_accu_line(lines, i, j, k, res_first)
					res_first = False


		tens_lines[r] = lines

	for reg in registers:
		f.write(" "*nin + "real(ki) :: %s\n" % reg)

	f.write(" "*nin + "if (present(momenta)) then\n")
	f.write(" "*(nin+3) + "rk = size(momenta, 1)\n")
	f.write(" "*nin + "else\n")
	f.write(" "*(nin+3) + "rk = 0\n")
	f.write(" "*nin + "end if\n")
	f.write(" "*nin + "select case(rk)\n")
	nin += 3

	for r in range(0, (R-6)+1):
		f.write(" "*(nin-3) + "case(%d)\n" % r)

		for line in tens_lines[r]:
			f.write(" "*nin + line + "\n")

	f.write(" "*(nin-3) + "case default\n")
	f.write(" "*nin + "res = 0.0_ki\n")
	nin -= 3
	f.write(" "*nin + "end select\n")
	f.write(" "*indent + "end  function %s\n" % name)

def calculate_tensor_terms(lines, terms, registers, level=0):
	"""
	Computes a multi-variate Horner-scheme for a set of terms,
	which represent monomials in the components of q.

	PARAMETER

	  lines     -- a list in which the generated commands are stored
	  terms     -- a list of terms (a list of lists of factors)
	  registers -- a set in which the names of intermediate expressions are
	               stored
	  level     -- internal parameter, level of recursion
	"""
	count_factors = {}
	for term in terms:
		for factor in term:
			if isinstance(factor, int):
				continue
			if factor in count_factors:
				count_factors[factor] += 1
			else:
				count_factors[factor] = 1
	if len(count_factors) == 0:
		gain = 0
	else:
		gain = max(count_factors.values())
	if gain <= 1:
		if len(terms) == 0:
			lines.append("acc = 0.0_ki")
		else:
			is_first = True
			for term in terms:
				assert len(term) > 0
				if isinstance(term[0], int):
					factor = term[0]

					if factor == -1:
						sign = -1
						ofs = 1
					elif factor == 1:
						sign = 1
						ofs = 1
					elif factor < 0:
						sign = -1
						term[0] = abs(factor)
						ofs = 0
					else:
						sign = 1
						ofs = 0
				else:
					sign = 1
					ofs = 0

				if is_first:
					if sign == -1:
						lhs = "acc = -"
					else:
						lhs = "acc = "
					is_first = False
				else:
					if sign == -1:
						lhs = "acc = acc - "
					else:
						lhs = "acc = acc + "

				if len(term[ofs:]) == 0:
					lines.append(lhs + "1.0_ki")
				else:
					lines.append(lhs + "*".join(map(str,term[ofs:])))
	else:
		pivot = max(count_factors.keys(), key=lambda x: count_factors[x])
		A0 = []
		A1 = []

		for term in terms:
			if pivot in term:
				tcp = term[:]
				tcp.remove(pivot)
				A1.append(tcp)
			else:
				A0.append(term)
		assert len(A1) > 0

		if len(A0) > 0:
			calculate_tensor_terms(lines, A1, registers, level+1)
			register = "reg%d" % level
			registers.add(register)
			lines.append(register + " = acc * " + pivot)
			calculate_tensor_terms(lines, A0, registers, level+1)
			lines.append("acc = acc + " + register)
		else:
			calculate_tensor_terms(lines, A1, registers, level+1)
			lines.append("acc = acc * " + pivot)

def write_function_contract(f, indent, name, N, R, d, shift=0):
	"""
	Writes a function to contract the coefficients with the form factors.

	shift -- 0: normal contraction
	         1: contract integrals with mu^2
				2: contract integrals with mu^4
				3: contract integrals with mu^6
	"""
	nin = indent + 3
	if shift == 0:
		coeff_type = R
	else:
		coeff_type = R - 2


	if coeff_type > 0:
		CT = PAT["coefftype"] % coeff_type
	else:
		CT = "complex(ki)"

	if shift > 0:
		SHI = "with (mu^2)^%d in the numerator " % shift
	else:
		SHI = ""

	for line in DOC.function_contract(name, N, R, CT, SHI, d-1):
		f.write(" "*indent+line+"\n")
	f.write(" "*indent
		+ "function     %s(coeffs, momenta, b_set) result(amp)\n" % name)
	f.write(" "*nin+"! generated by: write_function_contract\n")
	f.write(" "*nin + "implicit none\n")
	if coeff_type > 0:
		f.write(" "*nin + "type(%s), intent(in) :: coeffs\n" % CT)
	else:
		f.write(" "*nin + "complex(ki), intent(in) :: coeffs\n")
	f.write(" "*nin + "real(ki), dimension(:,0:), intent(in) :: momenta\n")
	f.write(" "*nin + "integer, intent(in) :: b_set\n")
	f.write(" "*nin + "type(form_factor) :: amp\n")

	if N >= 6:
		assert shift == 0
		write_contract_split(f, nin, N, R, d)
	else:
		write_contract_simple(f, nin, N, R, d, shift)

	f.write(" "*indent + "end function %s\n" % name)

def write_module_combine(name, tens_rec, max_legs, max_rank, d=4):
	f = open("%s.f90" % name, 'w')
	for line in DOC.module_combine():
		f.write(line + "\n")
	f.write("module %s\n" % name)
	f.write("! This module has been generated using a script.\n")
	f.write("! Please, refrain from modifying it directly!\n")
	f.write("use %s, only: ki\n" % precision_golem)
	f.write("use %s, only: b_ref, inv_s\n" % matrice_golem)
	f.write("use %s, only: packb, unpackb, countb, pminus, punion\n"
			% array_golem)
	f.write("use %s, only: form_factor, operator(+), operator(-), &\n"
			% ff_golem)
	f.write("    & operator(*), assignment(=)\n")
	for N in range(golem_minlegs,max_legs+1):
		write_list(f, "use %s, only:" % (ffp_golem % N),
				map(lambda s: s + ff_suffix, ff_list(N,N)))
	f.write("use form_factor_higher_ranks\n")
	lst = []

	write_import_list(f, "use %s, only:" % tens_rec)
	f.write("implicit none\n")
	f.write("private :: ki, b_ref, unpackb, packb, pminus, form_factor\n")
	write_import_list(f, "private ::")

	for N in range(golem_minlegs,max_legs+1):
		write_list(f, "private ::",
				map(lambda s: s + ff_suffix, ff_list(N,N)))


	f.write("real(ki), dimension(0:%d), parameter, private :: null_vec = &\n"
			% (d-1))
	f.write(" & (/%s/)\n" % ",".join(["0.0_ki"]*d))

	f.write("interface %s\n" % PAT['evaluate'])
	f.write("   module procedure %s_b\n" % PAT['evaluate'])
	f.write("   module procedure %s_s\n" % PAT['evaluate'])
	f.write("end interface\n")

        f.write("private :: test_if_nonzero\n")

	f.write("contains\n")

        f.write("\n\nfunction test_if_nonzero (val)\n")
        f.write("   complex(ki),intent(in):: val\n")
        f.write("   logical :: test_if_nonzero\n")
        f.write("   test_if_nonzero= (real(val,ki)*real(val,ki) + aimag(val)*aimag(val))>1E-100_ki\n")
        f.write("end function test_if_nonzero\n\n\n")

	write_function_evaluate(f, 0, PAT["evaluate"], max_legs, max_rank, d)

	for N in range(golem_minlegs, max_legs + 1):
		for R in (range(1, min(max_rank, N+extra_rank) + 1) \
				if N <= max_leg_extra_rank else range(1, min(N,max_rank)+1)):
			write_function_contract(f, 0, PAT['contract'] % (N, R), N, R, d)
			if N <= 4 and R >= 2:
				write_function_contract(f, 0, PAT['contract'] % (N, R) + "s1",
						N, R, d, 1)
			if N <= 4 and R >= 4:
				write_function_contract(f, 0, PAT['contract'] % (N, R) + "s2",
						N, R, d, 2)

			if N == 5 and R >= 6:
				write_function_contract(f, 0, PAT['contract'] % (N, R) + "s1",
						N, R, d, 1)
			if N == 5 and R >= 6:
				write_function_contract(f, 0, PAT['contract'] % (N, R) + "s2",
						N, R, d, 2)
			if N == 5 and R >= 6:
				write_function_contract(f, 0, PAT['contract'] % (N, R) + "s3",
						N, R, d, 3)


	for R in range(1, min(5,max_rank) + 1 + extra_rank):
		write_function_a_tensor_contract(f, 0, PAT['tenscontracta'] % R, R, d)
		if 2 <= R:
			write_function_b_tensor_contract(f, 0, PAT['tenscontractb'] % R, R, d)
		if 4 <= R:
			write_function_c_tensor_contract(f, 0, PAT['tenscontractc'] % R, R, d)
		if 6 <= R:
			write_function_d_tensor_contract(f, 0, PAT['tenscontractd'] % R, R, d)

	f.write("end module %s\n" % name)
	f.close()

def permutations(mult):
	"""
	Iterates over all permutations of [0]*mult[0]+...[n]*mult[n].
	"""
	N = sum(mult)

	def recselect(positions, i):
		for i_positions in select(positions, mult[i]):
			other_positions = positions[:]
			result = {}
			for p in i_positions:
				other_positions.remove(p)
				result[p] = i
			if len(other_positions) > 0:
				for dic in recselect(other_positions, i+1):
					result.update(dic)
					yield result
			else:
				yield result

	lst = [0] * N
	for dic in recselect(range(N), 0):
		for idx, value in dic.iteritems():
			lst[idx] = value
		yield lst

class RoboDocInfo:
	"""
	Class which reads text templates from a file.

	The file to be read must be in the format:

	@method_name(param1,param2,....)
	text which may contain format instructions like %(param1)s etc.

	The generated object will have methods with the same names and parameters
	iterating over the formatted lines

	EXAMPLE

	-------- file test.txt -------------
	@foo(name, age)
	Mr. %(name)s is %(age)d years old.
	@baz(name,size,cdate,mdate)
	The file %(name)s has a size of %(size)5.2f MB.
	It has been created on %(cdate)s.
	Its last modification was on %(mdate).
	---------------------------------------

	>>> DOC = RoboDoc('test.txt')
	>>> for line in DOC.baz("autoexec.bat", 3.053, "1980-01-01", "1981-03-17"):
	... 	print line
	The file autoexec.bat has a size of  3.05 MB.
	It has been created on 1980-01-01.
	Its last modification was on 1981-03-17.
	"""
	def __init__(self, fname):
		"""
		Create a new object from a file

		PARAMETER

		   fname -- the name of the file to be read.
		"""
		self._info = {}
		self._args = {}
		self._read(fname)

	def _read(self, fname):
		f = open(fname, "r")
		key = None
		for line in f.readlines():
			if line.startswith("@"):
				if "(" in line:
					tmp = line.split("(")
					line = tmp[0]
					args = map(lambda s: s.strip(),
							tmp[1].strip(" )\t\r\n").split(","))
				else:
					args = []
				key = line[1:].strip()
				self._args[key] = args
				self._info[key] = []
			elif key is not None:
				self._info[key].append(line.rstrip())
		f.close()
	
	def __getattr__(self, name):
		if name in self._info:
			lines = self._info[name]
			param = self._args[name]

			def meth(*args):
				if len(param) != len(args):
					raise TypeError, \
						"RoboDocInfo.%s(%s) expects exactly %d arguments (%d given)" \
							% (name, ", ".join(list(param)), len(param), len(args))

				opts = dict(zip(param, args))
				for line in lines:
					yield line % opts

			return meth
		else:
			raise AttributeError, \
				"This RoboDocInfo does not contain an entry for %r." % name

DOC = RoboDocInfo("robodoc.txt")

if __name__ == "__main__":
	max_legs = 6
	max_rank = 7
	max_leg_extra_rank = 5
	extra_rank = 1
	# support up to min(N+extra_rank,max_rank) rank for 1 to max_leg_extra_rank legs
	write_module_solve("tens_rec", max_rank)
	write_module_combine("tens_comb", "tens_rec", max_legs, max_rank)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
