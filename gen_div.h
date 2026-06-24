/* Given Lambda via its prime factorization (vectors primes and exponents), generate divisors of Lambda
and keep those for which d + 1 = p is prime.  Write output to a file.

Code written by Steven Hayman, 2011.  Andrew Shallue in 2026 updated the primality proving method, 
and gemini updated ZZ to mpz_class.
*/

#include <iostream> //for writing to standard out
#include <vector>
#include <fstream> //for writing to file

#include <time.h>
#include <string.h>

#include <gmp.h>
#include <gmpxx.h>  // for mpz_class

#include "is_prime_pl.h"

using std::vector;

// outfile1 stores the primes, one per line.  outfile2 stores the product of all the primes modulo Lambda.
// n is the thread number.  Maximum value is 7, corresponding to 8 threads.
void generateDivisors(const vector<unsigned long> primes, const vector<unsigned long> exponents, int n, char* outfileName1, char* outfileName2, bool verbose);
bool nextDivisor(const vector<unsigned long> exponents, vector<unsigned long>& divisor, int n);
bool differentPrime(const mpz_class& prime, const vector<unsigned long>& primes);
void setLastThree(vector<unsigned long>& divisor, unsigned long n);