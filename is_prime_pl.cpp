#include "is_prime_pl.h"

NTL_CLIENT


/*
Input: primes = {p_1,...,p_k}, exponents = {e_1,...e_k} where
n - 1 = p_1^e_1 dots p_k^e_k

Output: If n is proven to be prime then return true otherwise return false

This is a "p-1" primality test, written by Steven Hayman
*/

bool isPrimePL(ZZ n, const vector<long> primes, const vector<long> exponents, bool verbose) {

  /* Compute n - 1 = RF with F = p_0^a_0 dots p_r^a_r 
   * where (R,F)=1 and F > sqrt{n} */
	
  int r = 0;
	
  ZZ F = power_ZZ(primes[r],exponents[r]);
	
  while (F * F <= n) { 
    r++; 
    F = F * power_ZZ(primes[r],exponents[r]); 
  }
  
  /* If verbose is on, print F */
  
  if (verbose) { 
    cout << "F = " << primes[0] << "^" << exponents[0]; 
    for (int i = 1; i <= r; i++) {
      cout << "*" << primes[i] << "^" << exponents[i];
    }
  cout << "\n";
  }
	
	//Search for b_i with b_i^{n-1} equiv 1 mod n
	//and gcd(b_i^{(n-1)/p_i}-1,n) = 1 for all i = 0,..., r
	
    ZZ b, t;
    
	for (int i = 0; i <= r; i++) {
		for (b = 2; b*b <= n; b++) {  //search for b_i from 2 up to sqrt{n} + 1
            
            if (verbose) { cout << "trying b_" << i << " = " << b; }
            
			if(IsOne(PowerMod(b, n - 1, n))) {
                if (verbose) { cout << ", passed fermat's test"; }
                
				if (IsOne(GCD(SubMod(PowerMod(b, (n - 1)/primes[i], n), 1, n),n))) {
                    if (verbose) { cout << ", passed gcd test\n"; }
					break;
				}
                else {
                    if (verbose) { cout << ", failed gcd test\n"; }
                }
            }
			else {
                if (verbose) { cout << ", failed fermat's test\n"; }
				return false; //n is composite by fermats little theorem
			}
		}
		if (b*b > n) {
            if (verbose) { cout << "no b_i found\n"; }
			return false; //no b_i found
		}
	}
	return true;	
}

/* Here's a trial division algorithm written by Andrew
   May 2013
   Notice that the input vectors are supposed to be empty so that they can 
   become the output.  If they are not empty they are made empty.
 */
void factor(ZZ n, vector<long>& primes, vector<long>& exponents){
  // vectors should be empty.  If not, clear them
  if(primes.size() != 0) primes.clear();
  if(exponents.size() !=0) exponents.clear();

  ZZ p;
  ZZ exp;  
  while(n > 1){ // we will divide n with every factor found, done when n=1
    p = to_ZZ(2);
    exp = to_ZZ(0);
    while(n % p != 0){ p++;} // increase p until it is a factor; must be prime
    while(n % p == 0){  // increase exp by all the times p divides n
       exp++;
       n = n / p;
    }
    // now add to the vector
    primes.push_back(to_long(p));
    exponents.push_back(to_long(exp));
  }
}

/* the following is written by Andrew Shallue, May 2013.  It is useful to be 
able to print a factorization that is in the form created above
*/
void print_factors(vector<long> primes, vector<long> exponents){
  if(primes.size() != exponents.size()) cout << "Error in print_factors, vectors are not the same size\n";
  for(long i =0; i < primes.size(); i++){
    cout << primes.at(i) << "^" << exponents.at(i) << " * ";
  }
  cout << "\n";
} 