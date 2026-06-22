#include "is_prime_pl.h"
#include <iostream>
#include <chrono>
#include "CondensedInteger.h"
#include "unit_tests.h"
#include <memory>
#include <gmp.h>
#include <gmpxx.h>


int main(){
    std::cout << "Hello World\n";

	bool primality_proving_test_result = test_primality_proving();
	std::cout << "result of testing primality proving functions: ";
	std::cout << primality_proving_test_result << "\n";

	// primes and exponents corresponding to 10 billion prime run
	// long primes_[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
	// long exponents_[] = {16,7,5,4,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

	// Table 4, Loh/Niebuhr, entry 3
	long primes_[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43}; 
	long exponents_[] = {14,7,4,2,2,2,1,1,1,1,1,1,1,1};
	std::vector<long> primes (primes_, primes_ + sizeof(primes_)/sizeof(long));
	std::vector<long> exponents (exponents_, exponents_ + sizeof(exponents_)/sizeof(long));
	CondensedInteger::set(primes, exponents);
	
	
	bool ci = CI_testing();
	bool me1 = ME_test_set();
	std::cout << "ci result: " << ci << "\n";
	std::cout << "me result: " << me1 << "\n";

	bool me2 = ME_test_constructor();
	std::cout << "constructor: " << me2 << "\n";

	bool me3 = ME_test_history1();
	std::cout << "history1: " << me3 << "\n";
}