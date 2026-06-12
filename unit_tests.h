/* Test functions for Carmichael products project.
Andrew Shallue, June 2026
*/

#include <iostream>
#include <cstdlib> 
#include <ctime>   
#include <vector>
#include "CondensedInteger.h"
#include <gmp.h>

// test CondensedInteger by converting random divisors of a large Lambda
bool CI_convert_random(unsigned long num_trials);
bool CI_convert_extremes();