/* ModElement, based on work of Shallue and Hayman in 2011
This file written June 2026 by Andrew Shallue with the aid of Gemini Flash 3.5

CondensedInteger is as stripped-down and space-saving as possible. Each object represents 
a single prime p that satisfies p-1 | Lambda.
ModElement is designed to represent an arbitrary element in (Z/LambdaZ)*.
There is a flag which designates whether the actual residue is stored or not.  If not, 
one can recreate the residue by multiplying all the integers in the history modulo Lambda.

One should turn the flag off if there are many ModElement objects with short histories, 
and turn it on if there are few ModElement objects with long histories.
*/

#include <vector>
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include "CondensedInteger.h"

#ifndef MODELEMENT_H
#define MODELEMENT_H

class ModElement{
	public:
	    // this static information needs to be initialized before any members of the class are created
		// static means only one instance of these variables is stored for all objects in the class
	    static std::vector<long> primes; // primes in prime factorization of L
	    static std::vector<long> exponents;  // exponents in the prime factorization of L
		static mpz_class Lambda;  // the number whose factorization is given by primes and exponents

		std::vector<CondensedInteger> history;  // records which elements this one is a product of

		// flag designating whether n is stored explicitly or implicitly (via the history only)
		bool history_only;
		mpz_class n_mod_L;

		// in case you want to convert from mpz_class to mpz_t, do this:  rop = n_mod_L.get_mpz_t();

	public:
		// decided against pass-by-ref here. These are comparatively small, and copies are safer
		static void set(std::vector<long> primes_, std::vector<long> exponents_);

		ModElement(); // default constructor
		// non-default constructor assumes we already know that p-1 | Lambda
		ModElement(mpz_t p); 
		// multiply history to create n.  For times when you need n, but don't want to store it
		void to_mpz(mpz_t rop);


		// return residue of n modulo q, the prime power at index i of primes, exponents
		unsigned long residue(unsigned long index);
		// return the max index such that the residue modulo p_i^e_i is not 1.  Returns -1 if error or if n is identically 1
		long get_omega();
		// return true if n is 1 modulo L
		bool is_one();
		// flip the history_only flag to false, calculate n_mod_L and store it in the attribute
		void start_storing_n();
		// return the product of the current ModElement with another, modulo L
		ModElement product(ModElement& other, mpz_t L);
};

#endif