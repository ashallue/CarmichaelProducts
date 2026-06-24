/* Implementation of Algorithm 1 from Alford-Grantham-Hayman-Shallue, 2013
Goal is to generate many pairwise-relatively-prime Carmichael numbers from a single Lambda

Andrew Shallue, with the aid of Gemini-3.5-flash, June 2026
*/
#include <vector>
#include <algorithm>
#include <gmpxx.h>
#include "ModElement.h"

std::vector<ModElement> many_carmichaels_subset_product(const std::vector<ModElement>& P);