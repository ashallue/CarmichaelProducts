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

std::vector<unsigned long> ModElement::primes;
std::vector<unsigned long> ModElement::exponents;
mpz_class ModElement::Lambda = 1;

void ModElement::set(std::vector<unsigned long> primes_, std::vector<unsigned long> exponents_){
	ModElement::primes = primes_;
	ModElement::exponents = exponents_;
	CondensedInteger::set(primes_,exponents_);

	// compute Lambda
	ModElement::Lambda = 1;
	mpz_class power;
	for(unsigned long i = 0; i < ModElement::primes.size(); ++i){
		mpz_ui_pow_ui(power.get_mpz_t(), ModElement::primes.at(i), ModElement::exponents.at(i));
		ModElement::Lambda = ModElement::Lambda * power;
	}
}

ModElement::ModElement() {
	history = std::vector<CondensedInteger>();
	history_only = true;
}

// non-default constructor assumes we already know that p-1 | Lambda
ModElement::ModElement(mpz_t p) {
	history = std::vector<CondensedInteger>(1);

	// note that CondensedInteger constructor will subtract 1 from p
	history.at(0) = CondensedInteger(p);
	history_only = true;
}

// multiply history to create n.  For times when you need n, but don't want to store it
void ModElement::to_mpz(mpz_t rop){
	
	mpz_set_ui(rop, 1);
	mpz_t prod; mpz_init(prod);

	// loop over history and multiply mod n
	for(unsigned long i = 0; i < history.size(); ++i){
		history.at(i).to_mpz(prod);
		mpz_mul(rop, rop, prod);
	}
	mpz_clear(prod);
	
}


// return residue of n modulo q, the prime power at index i of primes, exponents
unsigned long ModElement::residue(unsigned long index){
	if(index >= ModElement::primes.size() || index < 0){
		std::cout << "Error in residue, index out of bounds\n";
		return -1;
	}
	
	// build q = p^e
	mpz_t q;  mpz_init(q);
	mpz_set_ui(q, 1);
	mpz_ui_pow_ui(q, ModElement::primes[index], ModElement::exponents[index]);

	mpz_t residue; mpz_init(residue);
	mpz_t prod; mpz_init(prod);
	// if n exists, simply perform modular division
	if(!history_only){
		mpz_mod(residue, n_mod_L.get_mpz_t(), q);

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
	mpz_clears(q, residue, prod, nullptr);
	return output;
}

// return the max index such that the residue modulo p_i^e_i is not 1
long ModElement::get_omega(){
	if (ModElement::primes.empty()) {
        return -1;
    }
	
	long index = ModElement::primes.size() - 1;
	unsigned long res = residue(index);
	while(index > -1 && res == 1){
		index--;
		if(index >= 0) res = residue(index);
	}
	
	return index;
	
}

// return true if n is 1 modulo L
bool ModElement::is_one(){
	// if we have n_mod_L, check if it is one
	if(!history_only){
		return n_mod_L == 1;

	// otherwise get n, reduce modulo L
	}else{
		mpz_t n;  mpz_init(n);
		this->to_mpz(n);

		// reduce n mod L
		mpz_mod(n, n, Lambda.get_mpz_t());
		bool output = mpz_cmp_ui(n, 1) == 0;

		// clean up and return
		mpz_clear(n);
		return output;
	}
}

// flip the history_only flag to false, calculate n_mod_L and store it in the attribute
void ModElement::start_storing_n(){
	// if we are already storing n_mod_L, do nothing
	if(!history_only){
		return;
	}else{
		history_only = false;

		// now multiply history together modulo Lambda
		n_mod_L = 1;
		mpz_t prod;  mpz_init(prod); 
		for(unsigned long i = 0; i < history.size(); ++i){
			history.at(i).to_mpz(prod);
			n_mod_L = n_mod_L * mpz_class(prod);
			n_mod_L = n_mod_L % ModElement::Lambda;
		}
		
		mpz_clear(prod);
	}
}

// return the product of the current ModElement with another, modulo L
// will update n_mod_L only if both this and other have the history_only flag turned off
ModElement ModElement::product(ModElement& other){
	ModElement result;
	
	unsigned long sizethis = this->history.size();
	unsigned long sizeother = other.history.size();
		
	// in either case, we need to concatenate the histories
	result.history.reserve(sizethis + sizeother);
	for(unsigned long i = 0; i < sizethis; ++i){
		result.history.push_back(this->history.at(i));
	}
	for(unsigned long i = sizethis; i < sizethis + sizeother; ++i){
		result.history.push_back(other.history.at(i - sizethis));
	}

	// if we are storing n_mod_L, we need to update it with the new product residue
	if(!history_only && !other.history_only){
		result.n_mod_L = this->n_mod_L * other.n_mod_L;
		result.n_mod_L = result.n_mod_L % ModElement::Lambda;
		result.history_only = false;
	}
		
	return result;
}