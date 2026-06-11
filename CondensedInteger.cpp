#include "CondensedInteger.h"

std::vector<long> CondensedInteger::primes = std::vector<long>();
std::vector<long> CondensedInteger::exponents = std::vector<long>();

void CondensedInteger::set(std::vector<long> primes_, std::vector<long> exponents_){
    CondensedInteger::primes = primes_;
    CondensedInteger::exponents = exponents_;
}

CondensedInteger::CondensedInteger() {
  rep = std::vector<bool>();
}

// Accept the mpz_t by const reference/pointer (Gemini claims GMP handles this internally)
// main action here is computing p-1 as a divisor of lambda, then storing bits of the exponents in rep
CondensedInteger::CondensedInteger(const mpz_t p) {

	int length = CondensedInteger::primes.size();
	// pminus which is p-1 will be the quantity decomposed
	mpz_t pminus;  mpz_init(pminus); mpz_sub_ui(pminus, p, 1);

	std::vector<int> positions;
	long total = 0;
	for (int i = 0; i < length; i++) {
		// calculate the bit length of exponents[i]
		long exp = exponents[i];
		unsigned long bitlength = 0;
		while (exp > 0) {
		    bitlength++;
	        exp >>= 1; 
	    }

		// add to total, calculate the next position
		total += bitlength;
		positions.push_back(total - 1);
	}

	rep = std::vector<bool>(total);
	// now calculate the bit representation of each exp, push those bits onto rep
	for (int i = 0; i < length; i++) {
  
		long exponent = 0;
		long prime = CondensedInteger::primes[i];
		while (mpz_divisible_ui_p(pminus, prime)){
			mpz_divexact_ui(pminus, pminus, prime);
			exponent++;
		}
  
		int position = positions[i];
		while (exponent > 0) {
			rep[position] = exponent % 2;
			exponent /= 2;
			position--;
		}
		std::cout << "position " << position << " resulting exponent is " << rep[position] << "\n";
	}
	// testing
	std::cout << "printing rep: ";
	for(unsigned long i = 0; i < rep.size(); ++i){
		std::cout << rep[i] << " ";
	}
	std::cout << "\n";
		
	mpz_clear(pminus);
}

// convert rep back to mpz_t, set rop to that value
void CondensedInteger::to_mpz(mpz_t rop){
	//testing
	std::cout << "printing rep: ";
	for(unsigned long i = 0; i < rep.size(); ++i){
		std::cout << rep[i] << " ";
	}
	std::cout << "\n";

	// reset rop to have value 1
	mpz_set_ui(rop, 1);

	// temp var for powers
	mpz_t pow; mpz_init(pow);
	
	int length = CondensedInteger::primes.size();

	// calculate the indices from rep corresponding to each exponent
	std::vector<int> positions;
	positions.push_back(-1);
	int total = 0;
	for (int i = 0; i < length; i++) {
		// calculate the bit length of exponents[i]
		long exp = exponents[i];
		unsigned long bitlength = 0;
		while (exp > 0) {
		    bitlength++;
	        exp >>= 1; 
	    }
		
		total += bitlength;
		positions.push_back(total - 1);
	}

	// from the bit rep, calculate the exp, then multiply by p to that exp
	for (int i = 0; i < length; i++) {
	  
		long exponent = 0;
		int k = 0;
		for (int j = positions[i + 1]; j > positions[i]; j--) {
			if (rep[j]) { 
				exponent += (1 << k);  // add 2^k to the running total
			}
			k++;
		}
		mpz_ui_pow_ui(pow, CondensedInteger::primes[i], exponent);
		mpz_mul(rop, rop, pow);
	}
	mpz_clear(pow);

	// add one to get the final value
	mpz_add_ui(rop, rop, 1);
	return;
}