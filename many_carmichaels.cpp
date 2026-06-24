/* Implementation of Algorithm 1 from Alford-Grantham-Hayman-Shallue, 2013
Goal is to generate many pairwise-relatively-prime Carmichael numbers from a single Lambda

Andrew Shallue, with the aid of Gemini-3.5-flash, June 2026
*/
#include "many_carmichaels.h"

// Note: this version written by Gemini with only light prompting.
// Context given: pseudocode from the paper, header files for ModElement, CondensedInteger
// Implementation of Algorithm 1: Many Carmichaels subset-product
std::vector<ModElement> many_carmichaels_subset_product(const std::vector<ModElement>& P) {
    // 1. S1 <- P
    std::vector<ModElement> S = P;
    size_t r = ModElement::primes.size();

    // 2. for u <- 1 to r do (0-indexed in C++)
    for (size_t u = 0; u < r; ++u) {
        
        // Calculate Q_u = primes[u]^exponents[u] using GMP
        mpz_class Q_u_mpz;
        mpz_ui_pow_ui(
            Q_u_mpz.get_mpz_t(), 
            static_cast<unsigned long>(ModElement::primes[u]), 
            static_cast<unsigned long>(ModElement::exponents[u])
        );

        // 3. sort Su by residue modulo Q_u
        std::sort(S.begin(), S.end(), [u](ModElement& x, ModElement& y) {
            return x.residue(u) < y.residue(u);
        });

        // 4. calculate bu = \prod_{a \in S_u} a mod Q_u using GMP to avoid overflow
        mpz_class b_u_mpz(1);
        for (auto& elem : S) {
            mpz_class res_mpz(elem.residue(u));
            b_u_mpz = (b_u_mpz * res_mpz) % Q_u_mpz;
        }
        unsigned long b_u = b_u_mpz.get_ui();

        // remove a from S_u that satisfies a === bu mod Q_u
        auto it = std::find_if(S.begin(), S.end(), [u, b_u](ModElement& x) {
            return x.residue(u) == b_u;
        });
        if (it != S.end()) {
            S.erase(it);
        }

        std::vector<ModElement> S_next;
        std::vector<bool> used(S.size(), false);

        // 5. for a in S_u do
        for (size_t i = 0; i < S.size(); ++i) {
            if (used[i]) continue;

            mpz_class r_a_mpz(S[i].residue(u));
            mpz_class inv_mpz;

            // 6. Find a_hat in S_u such that a * a_hat === 1 mod Q_u
            // Using GMP's native mpz_invert for modular inverse (extended gcd)
            if (mpz_invert(inv_mpz.get_mpz_t(), r_a_mpz.get_mpz_t(), Q_u_mpz.get_mpz_t())) {
                unsigned long target_residue = inv_mpz.get_ui();

                // Locate the matching element via binary search on sorted S
                auto comp_elem = [u](ModElement& elem, unsigned long val) {
                    return elem.residue(u) < val;
                };
                auto it_low = std::lower_bound(S.begin(), S.end(), target_residue, comp_elem);

                size_t match_idx = -1;
                for (auto it_match = it_low; it_match != S.end() && it_match->residue(u) == target_residue; ++it_match) {
                    size_t j = std::distance(S.begin(), it_match);
                    if (j != i && !used[j]) {
                        match_idx = j;
                        break;
                    }
                }

                if (match_idx != (size_t)-1) {
                    // 7. S_{u+1} <- a * a_hat
                    ModElement prod = S[i].product(S[match_idx]);
                    S_next.push_back(prod);

                    // 8. Remove a, a_hat from S_u (marked as used)
                    used[i] = true;
                    used[match_idx] = true;
                }
            }
        }

        // 9. product remaining a in S_u, add to S_{u+1}
        std::vector<ModElement> remaining;
        for (size_t i = 0; i < S.size(); ++i) {
            if (!used[i]) {
                remaining.push_back(S[i]);
            }
        }

        if (!remaining.empty()) {
            ModElement prod_rem = remaining[0];
            for (size_t k = 1; k < remaining.size(); ++k) {
                prod_rem = prod_rem.product(remaining[k]);
            }
            S_next.push_back(prod_rem);
        }

        // Advance to the next level
        S = std::move(S_next);
    }

    // 10. return S_{r+1}
    return S;
}