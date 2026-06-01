#include <iostream>
#include <vector>
#include <gmp.h>
#include <cstring>

using namespace std;

/* written by Steven Hayman, summer 2011
This is an implementation of the pocklington-lehmer primality proving alg
Both n and the factorization of n-1 are input

Input:  primes = {p_1,...,p_k}, exponents = {e_1,...e_k} where
n - 1 = p_1^e_1 dots p_k^e_k

Output:  returns true if n proven to be prime otherwise returns false

Update summer 2026: Andrew rewrites this to use gmp multi-precision instead.  
*/
bool isPrimePL(mpz_t n, const vector<unsigned long> primes, const vector<unsigned long> exponents, bool verbose);

/* written by Andrew Shallue, summer 2013
This prints a factorization given in prime and exponent vectors of the same length

factor is a trial division factoring algorithm.  The primes and exponents 
of the factorization are place in the input vectors.  If the input vectors 
are not empty they will be cleared.
*/
void print_factors(vector<unsigned long>& primes, vector<unsigned long>& exponents);
void factor(mpz_t n, vector<unsigned long>& primes, vector<unsigned long>& exponents);