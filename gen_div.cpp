#include "gen_div.h"

/*
gen_div.cpp   written by Steven Hayman, summer 2011
Andrew Shallue in 2026 updated the primality proving method, and gemini updated ZZ to mpz_class.

Input: primes = {p_1,...,p_r} and exps = {e_1,...,e_r} where
L = p_1^{e_1} dots p_r^{e_r},
n \in [-1,7] where -1 runs the non-parallel version of gen_div
and otherwise n determines the parallel thread

This function generates all the divisors d of L and checks if d + 1
is a prime different p_1,...,p_r
*/

#include "gen_div.h"
#include <gmpxx.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

using namespace std;

// Inline helper to calculate base^exp using GMP's native exponentiation
inline mpz_class pow_mpz(unsigned long base, unsigned long exp) {
  mpz_class res;
  mpz_ui_pow_ui(res.get_mpz_t(), base, exp);
  return res;
}

/*
gen_div.cpp   written by Steven Hayman, summer 2011
Converted to GMP (mpz_class)

Input: primes = {p_1,...,p_r} and exps = {e_1,...,e_r} where
L = p_1^{e_1} dots p_r^{e_r},
n \in [-1,7] where -1 runs the non-parallel version of gen_div
and otherwise n determines the parallel thread

This function generates all the divisors d of L and checks if d + 1
is a prime different p_1,...,p_r
*/

void generateDivisors(const vector<unsigned long> primes, const vector<unsigned long> exponents, int n, char* outfileName1, char* outfileName2, bool verbose) {

  ofstream outfile1, outfile2;
  outfile1.open(outfileName1);
  outfile2.open(outfileName2);
    
  time_t start, end;
    
  time(&start);
    
  int length = primes.size();

  // Build L = p_1^{e_1} dots p_r^{e_r}
  mpz_class L = 1;
  for (size_t i = 0; i < primes.size(); i++) {
    L *= pow_mpz(primes[i], exponents[i]);
  }
  if (verbose) {
    cout << "L = " << L << "\n";
  }

  int S = 0;

  // Declare and initialize divisor    
  vector<unsigned long> divisor;
  divisor.reserve(length);
  divisor.push_back(1);
  for (int i = 1; i < length; i++) { divisor.push_back(0); }
 
  if (-1 < n && n < 8) { setLastThree(divisor, n); n = 3; }  //parallel version sets last three digits
  else if (n == -1) { n = 0; }

  mpz_class p;
  mpz_class product = 1;
  do {

    //build p = d + 1
    p = 1;
    for (int i = 0; i < length; i++) {
      p *= pow_mpz(primes[i], divisor[i]);
    }
    p += 1;

    //test p for primality and make sure d + 1 != p_i for i = 1,...,r
    if (isPrimeBLS(p.get_mpz_t(), primes, divisor, false) == 1) {
      if (differentPrime(p, primes)) {
        product = (product * p) % L;
        outfile1 << p << "\n";
        S++;
      }
    }
  } while (nextDivisor(exponents, divisor, n));
   
 
  time(&end);

  outfile2 << product << "\n";
    
  if (verbose) {
    cout << "#S: " << S << "\n";
    cout << "time: " << difftime(end, start) << "\n";
  }
    
  outfile1.close();
  outfile2.close();
}

/*
Input: exps = {e_1,...,e_r} where L = p_1^{e_1} dots p_r^{e_r}
and div = {a_1,...,a_r} representing a divisor of L

Output: If there is another divisor to generate set div to represent
the next divisor and return true otherwise return false
 */

bool nextDivisor(const vector<unsigned long> exponents, vector<unsigned long>& divisor, int n) {
    
  int length = exponents.size();    
  int current = length -  1 - n;

  while (current >= 0 && divisor[current] == exponents[current]) {
    divisor[current] = 0;
    current--;
  }
    
  if (current == -1) { return false; }
    
  divisor[current] = divisor[current] + 1;
  return true;
}

/*
Input: primes = {p_1,...p_r} where L = p_1^{e_1} dots p_r^{e_r},
and prime = p_1^{a_r} dots p_r{a_r} + 1

Output: Returns true if prime is different from p_1,...,p_r
*/

bool differentPrime(const mpz_class& prime, const vector<unsigned long>& primes) {
    
  unsigned long length = primes.size();
    
  mpz_class otherPrime; 
  otherPrime = primes[length - 1];
    
  if ( prime > otherPrime ) { return true; }  //d + 1 is greater p_r
    
  for (unsigned long i = 0; i < length; i++) {
    otherPrime = primes[i];
    if ( prime == otherPrime ) {  //check if d + 1 = p_i
      return false;
    }
  }
  return true;
}

void setLastThree(vector<unsigned long>& divisor, unsigned long n) {
    
  unsigned long length = divisor.size();
    
  if (n == 0) {
    divisor[length - 3] = 0;
    divisor[length - 2] = 0;
    divisor[length - 1] = 0;
  }
  if (n == 1) {
    divisor[length - 3] = 0;
    divisor[length - 2] = 0;
    divisor[length - 1] = 1;
  }
  if (n == 2) {
    divisor[length - 3] = 0;
    divisor[length - 2] = 1;
    divisor[length - 1] = 0;
  }
  if (n == 3) {
    divisor[length - 3] = 0;
    divisor[length - 2] = 1;
    divisor[length - 1] = 1;
  }
  if (n == 4) {
    divisor[length - 3] = 1;
    divisor[length - 2] = 0;
    divisor[length - 1] = 0;
  }
  if (n == 5) {
    divisor[length - 3] = 1;
    divisor[length - 2] = 0;
    divisor[length - 1] = 1;
  }
  if (n == 6) {
    divisor[length - 3] = 1;
    divisor[length - 2] = 1;
    divisor[length - 1] = 0;
  }
  if (n == 7) {
    divisor[length - 3] = 1;
    divisor[length - 2] = 1;
    divisor[length - 1] = 1;
  }
}