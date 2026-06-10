#include "is_prime_pl.h"
#include <iostream>
#include <chrono>



int main(){
    std::cout << "Hello World\n";

	bool primality_proving_test_result = test_primality_proving();
	std::cout << "result of testing primality proving functions: ";
	std::cout << primality_proving_test_result << "\n";
}