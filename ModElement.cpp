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

#include "ModElement.h"

std::vector<long> ModElement::primes = std::vector<long>();
std::vector<long> ModElement::exponents = std::vector<long>();

void ModElement::set(std::vector<long> primes_, std::vector<long> exponents_){
	ModElement::primes = primes_;
	ModElement::exponents = exponents_;
	CondensedInteger::set(primes_,exponents_);
}

ModElement::ModElement() {
	history = std::vector<CondensedInteger>();
}

// non-default constructor assumes we already know that p-1 | Lambda
ModElement::ModElement(mpz_t p) {
	history = std::vector<CondensedInteger>(1);
	history.at(0) = CondensedInteger(p-1); 
}

// multiply history to create n.  For times when you need n, but don't want to store it
void ModElement::to_mpz(mpz_t rop){
	// if n_mod_L holds n, simply convert from mpz_class to mpz_t
	if(!history_only){
		rop = n_mod_L.get_mpz_t();
	}else{
		mpz_set_ui(rop, 1);
		mpz_t prod; mpz_init(prod);

		// loop over history and multiply mod n
		for(unsigned long i = 0; i < history.size(); ++i){
			history.at(i).to_mpz(prod);
			mpz_mul(rop, rop, prod);
		}
		mpz_clear(prod);
	}
}


// return residue of n modulo q, the prime power at index i of primes, exponents
unsigned long ModElement::residue(unsigned long index){
	if(index >= ModElement::primes.size() || index < 0){
		std::cout << "Error in residue, index out of bounds\n";
		return -1;
	}
	
	// build q = p^e
	mpz_t q;  mpz_init(q);
	mpz_set(q, 1);
	mpz_ui_pow_ui(q, ModElement::primes[index], ModElement::exponents[index]);

	mpz_t residue; mpz_init(residue);
	mpz_t prod; mpz_init(prod);
	// if n exists, simply perform modular division
	if(!history_only){
		mpz_mod(residue, n_mod_L, q);

	// otherwise, multiply the history together modulo q
	}else{
		mpz_set_ui(residue, 1);
		for(unsigned long i = 0; i < history.size(); ++i){
			history.at(i).to_mpz(prod);
			mpz_mul(residue, residue, prod);
			mpz_mod(residue, residue, q);
		}
	}

	// either way, convert residue to unsigned long
	unsigned long output = mpz_get_ui(residue);
	mpz_clears(q, residue, prod);
	return output;
}

// return the max index such that the residue modulo p_i^e_i is not 1
unsigned long ModElement::get_omega();

// return true if n is 1 modulo L
bool ModElement::is_one();

// flip the history_only flag to true, calculate n_mod_L and store it in the attribute
void ModElement::start_storing_n();

// return the product of the current ModElement with another, modulo L
ModElement ModElement::product(ModElement& other, mpz_t L);