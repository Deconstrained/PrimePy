PrimePy
=======

A Python module providing an object for finding the n-th prime number and class for prime decomposition and arithmetic with such objects.


I wrote this for problems on Project Euler involving prime numbers and factorization because I got tired of re-writing and copying/pasting prime number utilities.

## Usage:
<pre>
import prime

prime.cachedPrime[n+1] # the nth prime number (cachedPrime[0] is 2)

prime.PrimeFactors(n) # "PrimeFactors" object that behaves like a dict (prime:power pairs)
</pre>
See the methods/comments in the <tt>PrimeFactors</tt> class for more available operations. Note also that the multiplication, division and addition operators are overloaded.
