#include <iostream>
#include <vector>
#include <gmp.h>
#include <cstring>
#include "math.h"

using namespace std;

/* written by Steven Hayman, summer 2011
This is an implementation of the pocklington-lehmer primality proving alg
Both n and the factorization of n-1 are input

Input:  primes = {p_1,...,p_k}, exponents = {e_1,...e_k} where
n - 1 = p_1^e_1 dots p_k^e_k

Output:  returns true if n proven to be prime otherwise returns false

Update summer 2026: Andrew rewrites this to use gmp multi-precision instead.  
Note that F isn't a parameter.  We assume the factorization of n-1 is complete, 
then the function constructs F as the smallest i such that prod_i p_1^e_1 > sqrt(n).
*/
bool isPrimePL(mpz_t n, const vector<unsigned long> primes, const vector<unsigned long> exponents, bool verbose);

/* written by Andrew Shallue, summer 2013
This prints a factorization given in prime and exponent vectors of the same length

factor is a trial division factoring algorithm.  The primes and exponents 
of the factorization are place in the input vectors.  If the input vectors 
are not empty they will be cleared.
*/
void print_factors(vector<unsigned long>& primes, vector<unsigned long>& exponents);
void trial_factor(mpz_t n, vector<unsigned long>& primes, vector<unsigned long>& exponents);

/* Given an array already allocated, and given a bound which is the size of the array,
  create a factor sieve, i.e. the array is transformed to store the smallest prime factor of index i
  Then sieve_factor fills out primes, exponents with the complete factorization of n using a factor sieve fs.
*/
void factor_sieve(unsigned long* nums, unsigned long B);
void sieve_factor(mpz_t n, vector<unsigned long>& primes, vector<unsigned long>& exponents, unsigned long* fs, unsigned long B);

/* testing functions for primality proving algs.  Written by Andrew Shallue and Gemini 3.5 flash, June 2026
*/
// this function applies isPrimePL to every prime number up to a trial division bound
bool prove_primePL_all(unsigned long trial_bound);
// this function applies isPrimePL to num_trials many random mpz_t ints of given bit length
bool prove_primePL_random(unsigned long num_trials, unsigned long bit_length);