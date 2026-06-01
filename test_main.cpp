#include "is_prime_pl.h"
#include <iostream>

int main(){
    std::cout << "Hello World\n";

	mpz_t n, nminus;
	mpz_init(n);  mpz_init(nminus);
	mpz_set_ui(n, 101);
	mpz_set_ui(nminus, 100);

	// factor 100
	vector<unsigned long> nprimes;
	vector<unsigned long> nexps;
	factor(nminus, nprimes, nexps);
	print_factors(nprimes, nexps);
	bool primeyes = isPrimePL(n, nprimes, nexps, true);
	std::cout << "is prime returns " << primeyes << "\n";
	
	
	mpz_clear(n); mpz_clear(nminus);
}