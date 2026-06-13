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

#ifndef MODELEMENT_H
#define MODELEMENT_H

class ModElement{
	public:
	    // this static information needs to be initialized before any members of the class are created
		// static means only one instance of these variables is stored for all objects in the class
	    static vector<long> primes; // primes in prime factorization of L
	    static vector<long> exps;  // exponents in the prime factorization of L

	public:
		// decided against pass-by-ref here. These are comparatively small, and copies are safer
		static void set(std::vector<long> primes_, std::vector<long> exponents_);

		
};

#endif