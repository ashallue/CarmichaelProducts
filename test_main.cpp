#include "is_prime_pl.h"
#include <iostream>

int main(){
    std::cout << "Hello World\n";

	mpz_t n, nminus;
	mpz_init(n);  mpz_init(nminus);
	mpz_set_ui(n, 3);
	mpz_set_ui(nminus, 2);

	// factor 16
	vector<unsigned long> nprimes;
	vector<unsigned long> nexps;
	trial_factor(nminus, nprimes, nexps);
	print_factors(nprimes, nexps);
	bool primeyes = isPrimePL(n, nprimes, nexps, true);
	std::cout << "is prime returns " << primeyes << "\n";

	// factor 561
	nprimes.clear();
	nexps.clear();
	mpz_set_ui(n, 561);
	trial_factor(n, nprimes, nexps);
	print_factors(nprimes, nexps);

	// prove 1009 is prime
	nprimes.clear();
	nexps.clear();
	mpz_set_ui(n, 1009);
	mpz_set_ui(nminus, 1008);
	trial_factor(nminus, nprimes, nexps);
	print_factors(nprimes, nexps);
	primeyes = isPrimePL(n, nprimes, nexps, true);
	std::cout << "is prime returns " << primeyes << "\n";

	// factor sieve work
	unsigned long B = 1000;
	unsigned long* fs = new unsigned long[B];
	factor_sieve(fs, B);

	nprimes.clear();
	nexps.clear();
	mpz_set_ui(n, 600);
	sieve_factor(n, nprimes, nexps, fs, B);
	print_factors(nprimes, nexps);
	
	mpz_clear(n); mpz_clear(nminus);
	delete[] fs;

	bool test1 = prove_primePL_all(100000);
	if(test1){
		std::cout << "test 1 passed\n";
	}else{
		std::cout << "test 1 failed\n";
	}
	
	bool test2 = prove_primePL_random(1000, 20);
	if(test2){
		std::cout << "test 2 passed\n";
	}else{
		std::cout << "test 2 failed\n";
	}
}