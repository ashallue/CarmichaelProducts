/* Test functions for Carmichael products project.
Andrew Shallue, June 2026, with help from Gemini-3.5-flash
*/

#include <iostream>
#include <cstdlib> 
#include <ctime>   
#include <vector>
#include "CondensedInteger.h"
#include "ModElement.h"
#include <gmp.h>

// test CondensedInteger by converting random divisors of a large Lambda
bool CI_convert_random(unsigned long num_trials);
bool CI_convert_extremes();
// apply the above two tests in one function.  Returns true if both passed.
bool CI_testing();

// unit tests for ModElement
bool ME_test_set();
bool ME_test_constructor();
bool ME_test_history1();
bool ME_test_history2();
bool ME_test_residue();
bool ME_test_all();  // calls all ME tests above