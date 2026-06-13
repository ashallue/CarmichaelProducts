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