## Prime numbers module that aggressively pre-fetches and caches primes by sieving.
# Includes prime factorization classes and methods.
# 
# License: http://www.gnu.org/licenses/gpl.html

from math import sqrt
import numpy as npy
import sys,itertools

class PrimeCache(list):
	"""The nth element of this list-like object is the n+1 th prime number, with 2 living at index 0.

	Whenever a list element out of range is requested, the prime number there is calculated as well as any others in-between.
	"""
	def isPrime(self,n):
		"""Determine if a number is prime, and get new primes in the process"""
		if n <= 1:
			return False
		i = 0
		while self[i] <= sqrt(n):
			if n%self[i] == 0:
				return False
			i += 1
		return True

	def allPrimesLessThan(self,n):
		i = 0
		p = 2
		while p < n:
			yield p
			i += 1
			p = self[i]

	def __init__(self,primes=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]):
		"""Start out with the primes below 100"""
		super(PrimeCache,self).__init__(primes)

	def __getitem__(self,i):
		"""Get the (i+1)th prime number.

		Generates new primes if necessary by aggressively sieving ahead and caching all intermediate results found."""
		if(i < self.__len__()):
			# Get a preexisting prime number
			return super(PrimeCache,self).__getitem__(i)
		else:
			## Find the next biggest prime number ##
			p1 = self[-1]
			while True:
				p0 = p1+1 # Bottom of the sieve interval
				p1 = min(p1*2,sys.maxint-1) # Top of the sieve interval
				if p0 > p1:
					raise IndexError("Ran out of space!!!")
				r = npy.ones((p1-p0,),dtype=npy.bool)
				
				for n in self:
					# Construct the sieve
					res = p0%n
					r[int(res>0)*n-res::n] = 0
				if True in r:
					# Save all the primes that were found
					self.extend(npy.flatnonzero(r)+p0)
					if i < self.__len__():
						return self[i] # If it hasn't been found, this process will be repeated.
	def export(self):
		return [p for p in self]

cachedPrime = PrimeCache()

def isPrime(n):
	global cachedPrime
	return cachedPrime.isPrime(n)


def primes():
	"""Iterator object for cachedPrime"""
	i = 0
	while i < sys.maxint - 2:
		yield cachedPrime[i]
		i += 1 
 
def decompose(n):
	"""A constructor for iterating over the list of primes that equals n when multiplied together"""
	ln = long(n)
	for p in primes():
		if p*p > ln: break
		while ln % long(p) == 0:
			yield p
			ln /= long(p)
	if ln > 1:
		yield ln

def primePi(n):
	"""Pi function; returns the number of primes less than n"""
	i = 0
	while n > cachedPrime[i]:
		i+=1;
	return i

class PrimeFactors(dict):
	"""Generate a prime factorization for a number."""
	def __init__(self,n=False):
		"""Turn an integer into a prime decomposition"""
		if(n.__class__ in [int,npy.int64,long]):
			self._n = abs(n)
			factors = list(decompose(self._n)) if self._n!=1 else [1]
			for f in factors:
				self[f] += 1
		else:
			super(PrimeFactors,self).__init__(n)
	def __getitem__(self,i):
		if i in self:
			return super(PrimeFactors,self).__getitem__(i)
		else:
			return 0
	def __str__(self):
		return ' * '.join(['*'.join([str(n)]*self[n]) for n in sorted(self.keys())])
	def __int__(self):
		return long(self.val)
	def __mul__(self,f):
		self.canOperate(f)
		return PrimeFactors(dict([(p,self[p]+f[p]) for p in self.keys()+f.keys()]))
	def __div__(self,f):
		self.canOperate(f)
		return PrimeFactors(dict([(p,self[p]-f[p]) for p in self.keys()+f.keys()]))
	def __add__(self,f):
		self.canOperate(f)
		if f.isWhole() and self.isWhole():
			gcd = self.gcd(f)
			return gcd.__mul__(PrimeFactors(self.__div__(gcd).val + f.__div__(gcd).val))
		elif self.isWhole():
			# Adding a whole number to a fraction
			num = f.numerator()
			denom = f.denominator()
			gcd = num.gcd(self)
			return ((self.__mul__(denom)).__add__(num)).__div__(denom)
		else:
			# Adding two fractions
			num1 = self.numerator()
			num2 = f.numerator()
			den1 = self.denominator()
			den2 = f.denominator()
			den = den1.lcm(den2)
			mul = den1.gcd(den2)
			num1 = num1.__mul__(den2.__div__(mul))
			num2 = num2.__mul__(den1.__div__(mul))
			return (num1.__add__(num2)).__div__(den)
	def __pow__(self,p):
		return PrimeFactors(dict([(f,self[f]*p) for f in self]))
	def isWhole(self):
		"""Returns whether the number is a natural nubmer"""
		return not False in [self[p] > 0 for p in self]
	def Phi(self):
		"""Returns the value of the Euler totient function for the number"""
		return reduce(lambda x,y:x*(1-1./y),self.keys(),self.val)
	@property
	def val(self):
		"""Evaluate the object"""
		if not hasattr(self,'_val'):
			self._val = reduce(lambda x,y: x*y**self[y],self,1)
		return self._val
	def gcd(self,f):
		"""Obtain the greatest common divisor with another prime factorization object"""
		return PrimeFactors(dict([(p,min(self[p],f[p])) for p in self]))
	def lcm(self,f):
		return PrimeFactors(dict([(p,max(self[p],f[p])) for p in self]))
	def numerator(self):
		"""Numerator of fractional representation"""
		return PrimeFactors([(pf,self[pf]) for pf in filter(lambda p: self[p]>0,self)])
	def denominator(self):
		"""Denominator of fractional representation"""
		return PrimeFactors([(pf,-1*self[pf]) for pf in filter(lambda p: self[p]<0,self)])
	def canOperate(self,f):
		if f.__class__.__name__ != "PrimeFactors":
			raise TypeError('Incompatible types.')
	def divisors(self,p=[],i=0,top=True):
		"""Constructor object for getting divisors"""
		if top:
			p = self.keys()
		for j in range(self.__getitem__(p[i])+1): # All possible powers of the i-th prime
			if i < len(p)-1:
				for f in self.divisors(p,i+1,False): # Permutations of products of possible powers of other primes beyond i
					yield f*(p[i]**j)
			else:
				yield p[i]**j
	@property
	def properDivisors(self):
		"""The set of all proper divisors of the number."""
		if not hasattr(self,'_properDivisors'):
			self._properDivisors = set(filter(lambda x: x!=self.val,self.divisors()))
			self._properDivisors.add(1)
		return self._properDivisors
	@property
	def isAbundant(self):
		"""Whether the number is abundant"""
		if not hasattr(self,'_isAbundant'):
			self._isAbundant = sum(self.properDivisors) > self.val
		return self._isAbundant
	@property
	def isPerfect(self):
		"""Whether the number is perfect"""
		if not hasattr(self,'_isPerfect'):
			self._isPerfect = sum(self.properDivisors) == self.val
	@property
	def isWeird(self):
		"""Whether the number is weird"""
		if not hasattr(self,'_isWeird'):
			self._isWeird = self.isAbundant
			if self._isWeird:
				foundCombin = False
				pd = self.properDivisors
				for n in range(len(pd)):
					for c in itertools.combinations(pd,n+1):
						if sum(c) == self.val:
							self._isWeird = False
							return self._isWeird
		return self._isWeird

def composeNum(p):
	"""Put together a number from prime factors"""
	n = 1
	for pp in p.items():
		n *= pp[0]**pp[1]
	return n

def simplifyFrac(nf):
	"""Simplify a fraction using prime factorizations. Argument should be a 2-element list or tuple (numerator and denominator)"""
	#  f[0]/f[1]
	f = [PrimeFactors(n) for n in nf]
	for p in f[1].keys():
		minPow = min(f[0][p],f[1][p]) if p in f[0].keys() and p in f[1].keys() else 0
		if minPow > 0:
			f[0][p] -= minPow
			f[1][p] -= minPow
	return [pf.val for pf in f]

def addFrac(f0,f1):
	n0 = f0[0]
	d0 = f0[1]
	n1 = f1[0]
	d1 = f1[1]
	return simplifyFrac((n0*d1+n1*d0,d0*d1))
	
def sfPf(f):
	for p in f[1].keys():
		minPow = min(f[0][p],f[1][p]) if p in f[0].keys() and p in f[1].keys() else 0
		if minPow > 0:
			f[0][p] -= minPow
			f[1][p] -= minPow
	return f
