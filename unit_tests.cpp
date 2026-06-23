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

// apply the above two tests in one function.  Returns true if both passed.
bool CI_testing(){
	bool test1 = CI_convert_random(100);
	std::cout << "Test Condensed Integer conversion with random divisors: " << test1 << "\n";

	bool test2 = CI_convert_extremes();
	std::cout << "Test Condensed Integer conversion with biggest, smallest divisors: " << test2 << "\n";

	return test1 && test2;
}

//////////////////// Unit tests for ModElement class
/* When I asked gemini to write me tests, it defaulted to using the Google test suite, 
which one accesses with #include <gtest/gtest.h> and has commands like TEST_F and EXPECT_EQ
I don't want to rely on that, so I'm rewriting.  So test ideas from gemini, code by Andrew
*/

// test set method for ModElement
bool ME_test_set(){
	std::vector<long> test_primes = {2, 3, 5};
    std::vector<long> test_exponents = {3, 2, 1}; // L = 2^3 * 3^2 * 5^1
    ModElement::set(test_primes, test_exponents);

	bool correct = true;
	// check that the sizes of the static vectors are correct
	if(ModElement::primes.size() != 3 || ModElement::exponents.size() != 3){
		correct = false;
	}
	// check that the values of the static vectors match what was set
	if(ModElement::primes.at(0) != 2 || ModElement::primes.at(1) != 3 || ModElement::primes.at(2) != 5){
		correct = false;
	}
	if(ModElement::exponents.at(0) != 3 || ModElement::exponents.at(1) != 2 || ModElement::exponents.at(2) != 1){
		correct = false;
	}
	// check that Lambda has value 360
	if(ModElement::Lambda != 360){
		correct = false;
	}
	
	return correct;
}

// test non-default constructor
bool ME_test_constructor(){
	std::vector<long> test_primes = {2, 3, 5};
    std::vector<long> test_exponents = {3, 2, 1}; // L = 2^3 * 3^2 * 5^1
    ModElement::set(test_primes, test_exponents);

	// create 5, 7, 31 and check that to_ZZ matches
	mpz_t p1; mpz_t p2; mpz_t p3;
	mpz_t output1; mpz_t output2; mpz_t output3;
	mpz_inits(p1, p2, p3, output1, output2, output3, nullptr);
	mpz_set_ui(p1, 5); mpz_set_ui(p2, 7); mpz_set_ui(p3, 31);

	ModElement me1(p1);  ModElement me2(p2); ModElement me3(p3);
	me1.to_mpz(output1); me2.to_mpz(output2); me3.to_mpz(output3);

	gmp_printf ("Converting back to mpz should give 5, 7, 31: %Zd, %Zd, %Zd \n", output1, output2, output3);

	// check that we get the same values back again after using to_mpz
	bool check1 = mpz_cmp_ui(output1, 5) == 0;
	bool check2 = mpz_cmp_ui(output2, 7) == 0;
	bool check3 = mpz_cmp_ui(output3, 31) == 0;

	mpz_clears(p1, p2, p3, output1, output2, output3, nullptr);
	
	// also check that history_only is set to true
	return check1 && me1.history_only && check2 && me2.history_only && check3 && me3.history_only;
}

// test history and the product function, with history_only on
bool ME_test_history1(){
	std::vector<long> test_primes = {2, 3, 5};
    std::vector<long> test_exponents = {3, 2, 1}; // L = 2^3 * 3^2 * 5^1
    ModElement::set(test_primes, test_exponents);

	// create 5, 7, 31 and set in the history
	mpz_t p1; mpz_t p2; mpz_t p3;
	mpz_t output1; mpz_t output2; mpz_t output3;
	mpz_inits(p1, p2, p3, output1, output2, output3, nullptr);
	mpz_set_ui(p1, 5); mpz_set_ui(p2, 7); mpz_set_ui(p3, 31);
	ModElement me1(p1); ModElement me2(p2); ModElement me3(p3);
	
	// check history of me1 is only p1
	bool check1 = me1.history.size() == 1;
	CondensedInteger c1 = me1.history.at(0);
	c1.to_mpz(output1);
	check1 = check1 && (mpz_cmp_ui(output1, 5) == 0);
	me1.to_mpz(output1);
	check1 = check1 && (mpz_cmp_ui(output1, 5) == 0);

	// product with me2, check history and product
	me1 = me1.product(me2);
	bool check2 = me1.history.size() == 2;
	me1.to_mpz(output2);
	check2 = check2 && (mpz_cmp_ui(output2, 35) == 0);

	// product with me3, check history and product
	me1 = me1.product(me3);
	bool check3 = me1.history.size() == 3;
	me1.to_mpz(output3);
	check3 = check3 && (mpz_cmp_ui(output3, 1085) == 0);

	mpz_clears(p1, p2, p3, output1, output2, output3, nullptr);
	
	return check1 && check2 && check3;
}

// test history and the product function, with history_only off
bool ME_test_history2(){
	std::vector<long> test_primes = {2, 3, 5};
    std::vector<long> test_exponents = {3, 2, 1}; // L = 2^3 * 3^2 * 5^1
    ModElement::set(test_primes, test_exponents);

	// create 7, 31, 41 and set in the history
	mpz_t p1; mpz_t p2; mpz_t p3;
	mpz_t output1; mpz_t output2; mpz_t output3;
	mpz_inits(p1, p2, p3, output1, output2, output3, nullptr);
	mpz_set_ui(p1, 7); mpz_set_ui(p2, 31); mpz_set_ui(p3, 41);
	ModElement me1(p1); ModElement me2(p2); ModElement me3(p3);

	// flip the flag
	me1.start_storing_n();  me2.start_storing_n();  me3.start_storing_n();
	
	// check history of me1 is only p1
	bool check1 = me1.history.size() == 1;
	CondensedInteger c1 = me1.history.at(0);
	c1.to_mpz(output1);
	check1 = check1 && (mpz_cmp_ui(output1, 7) == 0);
	me1.to_mpz(output1);
	check1 = check1 && (mpz_cmp_ui(output1, 7) == 0);
	check1 = check1 && (me1.n_mod_L == 7);

	// product with me2, check history and product
	me1 = me1.product(me2);
	bool check2 = me1.history.size() == 2;
	me1.to_mpz(output2);
	check2 = check2 && (mpz_cmp_ui(output2, 217) == 0);
	check2 = check2 && (me1.n_mod_L == 217);

	// product with me3, check history and product
	me1 = me1.product(me3);
	bool check3 = me1.history.size() == 3;
	me1.to_mpz(output3);

	// n should be 7 * 31 * 41 = 8897, n_mod_L should be 8897 % 360 = 257
	check3 = check3 && (mpz_cmp_ui(output3, 8897) == 0);
	check3 = check3 && (me1.n_mod_L == (8897 % ModElement::Lambda));

	mpz_clears(p1, p2, p3, output1, output2, output3, nullptr);
	
	return check1 && check2 && check3;
}

// test residue and omega function
bool ME_test_residue(){
	std::vector<long> test_primes = {2, 3, 5};
    std::vector<long> test_exponents = {3, 2, 1}; // L = 2^3 * 3^2 * 5^1
    ModElement::set(test_primes, test_exponents);

	// create 7, 31, 41, 361 and set in the history
	mpz_t p1; mpz_t p2; mpz_t p3; mpz_t p4;
	mpz_inits(p1, p2, p3, p4, nullptr);
	mpz_set_ui(p1, 7); mpz_set_ui(p2, 31); mpz_set_ui(p3, 41); mpz_set_ui(p4, 361);
	ModElement me1(p1); ModElement me2(p2); ModElement me3(p3); ModElement me4(p4);

	// residues of 41 for indices 0, 1, 2 are 1, 5, 1
	bool check1 = me3.residue(0) == 1 && me3.residue(1) == 5 && me3.residue(2) == 1;

	// residues of 31 for indices 0, 1, 2 are 7, 4, 1
	bool check2 = me2.residue(0) == 7 && me2.residue(1) == 4 && me2.residue(2) == 1;

	// multiply 7 by 31, residues now are 1, 1, 2
	me1 = me1.product(me2);
	bool check3 = me1.residue(0) == 1 && me1.residue(1) == 1 && me1.residue(2) == 2;

	// check values for Omega.  Should be 2, 1, 1, -1
	bool check4 = me1.get_omega() == 2 && me2.get_omega() == 1 && me3.get_omega() == 1 && me4.get_omega() == -1;

	// check that me4 is identified as 1 but me1, me2, me3 aren't
	bool check5 = me4.is_one() && !me1.is_one() && !me2.is_one() && !me3.is_one();
	
	mpz_clears(p1, p2, p3, p4, nullptr);

	return check1 && check2 && check3 && check4 && check5;
}

// combine all the above into one function call
bool ME_test_all(){
	bool me1 = ME_test_set();
	std::cout << "me result: " << me1 << "\n";

	bool me2 = ME_test_constructor();
	std::cout << "constructor: " << me2 << "\n";

	bool me3 = ME_test_history1();
	std::cout << "history1: " << me3 << "\n";

	bool me4 = ME_test_history2();
	std::cout << "history2: " << me4 << "\n";

	bool me5 = ME_test_residue();
	std::cout << "residue test: " << me5 << "\n";

	return me1 && me2 && me3 && me4 && me5;
}