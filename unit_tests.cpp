/* Test functions for Carmichael products project.
Andrew Shallue, June 2026
*/

#include "unit_tests.h"

// test CondensedInteger by converting random divisors of a large Lambda
bool CI_convert_random(unsigned long num_trials) {
	// primes and exponents corresponding to 10 billion prime run
	long primes_[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
	long exponents_[] = {16,7,5,4,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	std::vector<long> primes (primes_, primes_ + sizeof(primes_)/sizeof(long));
	std::vector<long> exponents (exponents_, exponents_ + sizeof(exponents_)/sizeof(long));
	CondensedInteger::set(primes, exponents);
	
	// number to be created such that n-1 divides lambda
	mpz_t n;  mpz_init(n);
	mpz_t pow; mpz_init(pow);
	mpz_t n_after; mpz_init(n_after);
	CondensedInteger n_small;
	
	// random seed
	std::srand(static_cast<unsigned int>(std::time(nullptr)));

	// boolean accumulation on a bunch of trials
	bool all_correct = true;

	for(unsigned long i = 0; i < num_trials; ++i){
	
		// build a random divisor as a vector of exponents
		unsigned long length = primes.size();
		std::vector<unsigned long> div_exps;
		for(unsigned long i = 0; i < length; ++i){
			unsigned long e = std::rand() % (exponents.at(i) + 1);
			div_exps.push_back(e);
		}
		
		// now build the corresponding mpz.  Note adding 1 at the end
		mpz_set_ui(n, 1);
		for(unsigned long i = 0; i < length; ++i){
			mpz_set_ui(pow, 1);
			mpz_ui_pow_ui(pow, primes[i], div_exps[i]);
			mpz_mul(n, n, pow);
		}
		mpz_add_ui(n, n, 1);
	
		// convert to CondensedInteger, convert back
		n_small = CondensedInteger(n);
		n_small.to_mpz(n_after);
		all_correct = all_correct && (mpz_cmp(n, n_after) == 0);
	}
	mpz_clear(n); mpz_clear(pow); mpz_clear(n_after);
	return all_correct;
}

// test min and max divisors, make sure conversion works
bool CI_convert_extremes(){
	// primes and exponents corresponding to 10 billion prime run
	long primes_[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
	long exponents_[] = {16,7,5,4,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	std::vector<long> primes (primes_, primes_ + sizeof(primes_)/sizeof(long));
	std::vector<long> exponents (exponents_, exponents_ + sizeof(exponents_)/sizeof(long));
	CondensedInteger::set(primes, exponents);
	
	// number to be created such that n-1 divides lambda
	mpz_t n;  mpz_init(n);
	mpz_t pow; mpz_init(pow);
	mpz_t n_after; mpz_init(n_after);
	CondensedInteger n_small;

	// test n = 2, so divisor is 1
	mpz_set_ui(n, 2);
	n_small = CondensedInteger(n);
	n_small.to_mpz(n_after);
	bool compare1 = mpz_cmp(n, n_after) == 0;

	// test n which is 1 more than the maximal divisor
	mpz_set_ui(n, 1);
	// now build the corresponding mpz.  Note adding 1 at the end
	for(unsigned long i = 0; i < primes.size(); ++i){
		mpz_set_ui(pow, 1);
		mpz_ui_pow_ui(pow, primes[i], exponents[i]);
		mpz_mul(n, n, pow);
	}
	mpz_add_ui(n, n, 1);

	n_small = CondensedInteger(n);
	n_small.to_mpz(n_after);
	bool compare2 = mpz_cmp(n, n_after) == 0;

	return compare1 && compare2;
}