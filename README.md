PrimePy
=======

A Python module providing
- An object for finding the n-th prime number
- A class for prime decomposition / factorization, arithmetic with prime factorizations, and finding divisors
- A function for simplifying fractions


I wrote this for working on problems on Project Euler involving prime numbers and factorization because I got tired of re-writing and copying/pasting prime number utilities. Also, the prime number generators from Rosetta Code I was initially using seemed a bit inefficient, so I decided to try making something better.

## Usage:
<pre>
import prime

prime.cachedPrime[n+1] # the nth prime number (cachedPrime[0] is 2)

prime.PrimeFactors(n) # "PrimeFactors" object that behaves like a dict (prime:power pairs corresponding to natural number n)
</pre>
See the methods/comments in the <tt>PrimeFactors</tt> class for more available operations and properties. Note also that the multiplication, division and addition operators of the <tt>PrimeFactors</tt> class are overloaded, so they'll behave as integers when used in arithmetic expressions.
