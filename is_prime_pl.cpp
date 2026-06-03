#include "is_prime_pl.h"

using namespace std;


/*
Input: primes = {p_1,...,p_k}, exponents = {e_1,...e_k} where
n - 1 = p_1^e_1 dots p_k^e_k

Output: If n is proven to be prime then return true otherwise return false

This is a "p-1" primality test, written by Steven Hayman
Rewritten with mpz_t arithmetic by Andrew Shallue and Gemini 3.5 flash, May 2026
*/
bool isPrimePL(mpz_t n, const vector<unsigned long> primes, const vector<unsigned long> exponents, bool verbose) {

	/* Compute n - 1 = RF with F = p_0^a_0 dots p_r^a_r 
	 * where (R,F)=1 and F > sqrt{n} */
		
	int r = 0;
	
	// initialize variables
	mpz_t F, F_squared, temp;
	mpz_init(F); mpz_init(F_squared); mpz_init(temp);

	// multiply by prime powers until F is greater than sqrt(n)
	mpz_ui_pow_ui(F, primes[r], exponents[r]);
	mpz_mul(F_squared, F, F);
	while (mpz_cmp(F_squared, n) <= 0) { 
		r++; 
		mpz_ui_pow_ui(temp, primes[r], exponents[r]);
		mpz_mul(F, F, temp);
		mpz_mul(F_squared, F, F);
	}
	
	/* If verbose is on, print F */
	if (verbose) { 
		cout << "F = " << primes[0] << "^" << exponents[0]; 
		for (int i = 1; i <= r; i++) {
			cout << "*" << primes[i] << "^" << exponents[i];
		}
		cout << "\n";
	}
	
	// Search for b_i with b_i^{n-1} equiv 1 mod n
	// and gcd(b_i^{(n-1)/p_i}-1,n) = 1 for all i = 0,..., r
	mpz_t b;
	mpz_t nm1, exp, base_pow, base_pow_minus_1, gcd_val;
	
	mpz_init(b); 
	mpz_init(nm1); mpz_init(exp); mpz_init(base_pow); mpz_init(base_pow_minus_1); mpz_init(gcd_val);

	// set nm1 as n - 1, set isprime to true
	mpz_sub_ui(nm1, n, 1);
	bool is_prime = true;

	// i corresponds to a prime in the factorization of n-1
	for (int i = 0; i <= r; i++) {
		bool found_b = false;
		mpz_set_ui(b, 2);
	
		while (mpz_cmp(b, n) < 0) {  // search for b_i up to n, hope to find it faster
			if (verbose) { 
				char* b_str = mpz_get_str(NULL, 10, b);
				cout << "trying b_" << i << " = " << b_str;
				
				// Safely free the memory allocated by mpz_get_str
				void (*freefunc)(void *, size_t);
				mp_get_memory_functions(NULL, NULL, &freefunc);
				freefunc(b_str, strlen(b_str) + 1);
			} //endif

			// calculate b^{n-1} mod n
			mpz_powm(base_pow, b, nm1, n);
			if (mpz_cmp_ui(base_pow, 1) == 0) {  // if equal to 1
				
				if (verbose) { cout << ", passed fermat's test"; }
					mpz_tdiv_q_ui(exp, nm1, primes[i]);
					mpz_powm(base_pow, b, exp, n);
			
				// Equivalent to SubMod(base_pow, 1, n)
				if (mpz_cmp_ui(base_pow, 0) == 0) {
					mpz_sub_ui(base_pow_minus_1, n, 1);
				} else {
					mpz_sub_ui(base_pow_minus_1, base_pow, 1);
				}
		
				mpz_gcd(gcd_val, base_pow_minus_1, n);
				if (mpz_cmp_ui(gcd_val, 1) == 0) {  // if gcd is 1
					if (verbose) { cout << ", passed gcd test\n"; }
					found_b = true;
					break;
				} else {
				  if (verbose) { cout << ", failed gcd test\n"; }
				}
			} else {
				if (verbose) { cout << ", failed fermat's test\n"; }
				is_prime = false; // n is composite by Fermat's little theorem
				break; 
			}

			// check the next b
			mpz_add_ui(b, b, 1);
		} // end while
	
		if (!is_prime) {
		  break;
		}
	
		if (!found_b) {
		  if (verbose) { cout << "no b_i found\n"; }
		  is_prime = false; // no b_i found
		  break;
		}
	} // end for
	
	// Free all locally allocated GMP memory before exit
	mpz_clear(F); mpz_clear(F_squared); mpz_clear(temp);
	mpz_clear(b); 
	mpz_clear(nm1); mpz_clear(exp); mpz_clear(base_pow); mpz_clear(base_pow_minus_1); mpz_clear(gcd_val);
	
	return is_prime;	
}


/* Here's a trial division algorithm written by Andrew
   May 2013
   Notice that the input vectors are supposed to be empty so that they can 
   become the output.  If they are not empty they are made empty.
 */
void trial_factor(mpz_t n, vector<unsigned long>& primes, vector<unsigned long>& exponents){
	// vectors should be empty.  If not, clear them
	if(primes.size() != 0) primes.clear();
	if(exponents.size() !=0) exponents.clear();

	mpz_t p, exp, temp;
	mpz_init(p);  mpz_init(exp); mpz_init(temp);

	// continue doing trial division until n is 1
	while(mpz_cmp_ui(n, 1) != 0){
		// start with prime 2 as a potential divisor of n
		mpz_set_ui(p, 2);  mpz_set_ui(exp, 0);
		mpz_mod(temp, n, p);

		// continue checking until a divisor found; it must be prime
		while(mpz_cmp_ui(temp, 0) != 0){
			mpz_add_ui(p, p, 1);
			mpz_mod(temp, n, p);
		}

		// now find the exponent of p that divides n
		mpz_divexact(n, n, p);
		mpz_mod(temp, n, p);
		mpz_add_ui(exp, exp, 1);
		while(mpz_cmp_ui(temp, 0) == 0){
			mpz_add_ui(exp, exp, 1);
			mpz_divexact(n, n, p);
			mpz_mod(temp, n, p);
		}

		// add results to the vectors;
		primes.push_back(mpz_get_ui(p));
		exponents.push_back(mpz_get_ui(exp));

		//std::cout << "pushed p = " << mpz_get_ui(p) << " and e = " << mpz_get_ui(exp) << "\n";
	} // end while n is not 1

	// clean up mpz_t vars
	mpz_clear(p); mpz_clear(exp); mpz_clear(temp);
}

/* the following is written by Andrew Shallue, May 2013.  It is useful to be 
able to print a factorization that is in the form created above
*/
void print_factors(vector<unsigned long>& primes, vector<unsigned long>& exponents){
	if(primes.size() != exponents.size()) cout << "Error in print_factors, vectors are not the same size\n";
	for(long i =0; i < primes.size(); i++){
	    cout << primes.at(i) << "^" << exponents.at(i) << " * ";
	}
	cout << "\n";
} 

// turns nums into a factor sieve.  Assumes allocated space for B longs.  Written by Andrew in the 2021 maybe.
void factor_sieve(unsigned long* nums, unsigned long B){
	// set nums to contain i at position i for all i < B
	for(unsigned long i = 0; i < B; ++i){
		nums[i] = i;
	}
	
	// pivot refers to the current prime being sieved
	long pivot = 2;
	long index = 2;
	// outer loop is over primes up to sqrt(B)
	unsigned long prime_bound = std::ceil(std::sqrt(static_cast<double>(B)));
	
	// continue while the next prime is smaller than the bound
	while(pivot < prime_bound){
		//std::cout << "pivot = " << pivot << "\n";
	
		// while smaller than the bound replace entry in index 
		// with the pivot.
		for(index = 2 * pivot; index < B; index = index + pivot){
			if(nums[index] == index){
				nums[index] = pivot;
			}
		}
		
		// find the next prime.  It will be the first entry where 
		// the value equals the index
		++pivot;
		while(nums[pivot] != pivot){
			++pivot;
		}
	} // end outer while
}

// Then sieve_factor fills out primes, exponents with the complete factorization of n using a factor sieve fs.
void sieve_factor(mpz_t n, vector<unsigned long>& primes, vector<unsigned long>& exponents, unsigned long* fs, unsigned long B){
	// this only works if n < B.  Check bounds then convert n to unsigned long
	if(mpz_cmp_ui(n, 0) < 0 || mpz_cmp_ui(n, B) > 0){
		std::cout << "Error in sieve_factor, given n is out of bounds\n";
		return;
	}
	unsigned long n_ui = mpz_get_ui(n);
	unsigned long p; unsigned long e;
	primes.clear(); exponents.clear();

	// continue factoring until n is 1
	while(n_ui > 1){
		p = fs[n_ui];
		e = 1;
		n_ui = n_ui / p;

		// find all the factors of p
		while(fs[n_ui] == p){
			e++;
			n_ui = n_ui / p;
		}

		// now add p^e to the ouput
		primes.push_back(p);  exponents.push_back(e);
	}
}

/* testing functions for primality proving algs.  Written by Andrew Shallue and Gemini 3.5 flash, June 2026
*/
// Helper function to check if the prime factorization of n-1 is complete
bool is_factorization_complete(const mpz_t n_minus_1, const std::vector<unsigned long>& primes, const std::vector<unsigned long>& exponents) {
    mpz_t product;
    mpz_init_set_ui(product, 1);
    
    for (size_t i = 0; i < primes.size(); i++) {
        mpz_t term;
        mpz_init(term);
        mpz_ui_pow_ui(term, primes[i], exponents[i]);
        mpz_mul(product, product, term);
        mpz_clear(term);
    }
    
    int comparison = mpz_cmp(product, n_minus_1);
    mpz_clear(product);
    return comparison == 0;
}

// this function applies isPrimePL to every prime number up to a trial division bound
bool prove_primePL_all(unsigned long trial_bound) {
    if (trial_bound < 2) {
        return true; 
    }

	// [Andrew] let's use a factor sieve for the factoring
	unsigned long* fs = new unsigned long[trial_bound];
	factor_sieve(fs, trial_bound);
    
    bool all_passed = true;
    for (unsigned long p = 2; p <= trial_bound; p++) {
		// only apply the prime proving to prime numbers
        if (fs[p] == p) {
            // Pocklington-Lehmer is typically defined for n >= 3
            if (p == 2) {
                continue;
            }
            
            mpz_t n, n_minus_1;
            mpz_init_set_ui(n, p);
            mpz_init_set_ui(n_minus_1, p - 1);
            
            std::vector<unsigned long> primes;
            std::vector<unsigned long> exponents;

			// factor with the factor sieve
			sieve_factor(n_minus_1, primes, exponents, fs, trial_bound);
            
            bool result = isPrimePL(n, primes, exponents, false);
            if (!result) {
				std::cout << "failed for p = " << p << "\n";
				print_factors(primes, exponents);
                all_passed = false;
            }
            
            mpz_clear(n);
            mpz_clear(n_minus_1);
        }
    }
	// clean up and return
	delete[] fs;
    return all_passed;
}

// this function applies isPrimePL to num_trials many random mpz_t ints of given bit length
bool prove_primePL_random(unsigned long num_trials, unsigned long bit_length) {
    if (bit_length < 2) {
        return false;
    }
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    
    unsigned long successful_trials = 0;
    bool all_passed = true;
    
    mpz_t n, n_minus_1, temp_n_minus_1;
    mpz_init(n);
    mpz_init(n_minus_1);
    mpz_init(temp_n_minus_1);
    
    // To prevent an infinite loop if bit_length is too large for trial division to factor n-1
    unsigned long attempts = 0;
    unsigned long max_attempts = num_trials * 10000;
    
    while (successful_trials < num_trials && attempts < max_attempts) {
        attempts++;
        
        mpz_urandomb(n, state, bit_length);
        
        // Ensure the candidate has the specified bit length and is odd
        mpz_setbit(n, bit_length - 1);
        mpz_setbit(n, 0);
        
        // Check if the candidate is probably prime
        if (mpz_probab_prime_p(n, 25) > 0) {
            mpz_sub_ui(n_minus_1, n, 1);
            mpz_set(temp_n_minus_1, n_minus_1);
            
            std::vector<unsigned long> primes;
            std::vector<unsigned long> exponents;
            
            trial_factor(temp_n_minus_1, primes, exponents);
            
            // Pocklington-Lehmer requires a known factorization of n-1
            if (is_factorization_complete(n_minus_1, primes, exponents)) {
                bool result = isPrimePL(n, primes, exponents, false);
                if (!result) {
                    all_passed = false;
                }
                successful_trials++;
            }
        }
    }
    
    mpz_clear(n);
    mpz_clear(n_minus_1);
    mpz_clear(temp_n_minus_1);
    gmp_randclear(state);
    
    // Returns true only if we successfully performed the requested number of trials and all passed
    return (successful_trials == num_trials) && all_passed;
}