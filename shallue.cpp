
/* shallue.cpp

  Author: Steven Hayman, Andrew Shallue
  Date: July 2011

  Code for constructing Carmichael Numbers.

*/

/*
  Things to fix:
     
  1. Global M.  It would be better any way to have easier ways of 
     working with M.  UPDATE:  fixed by functor, though in this version
     it isn't even needed.
  2. Either keep j = r/2 or choose the commented out section of code.
     (currently j = length because of problem 3)
     UPDATE:  fixed
  3. Loh Niebuhr is not working for other primes.
     UPDATE:  fixed but not being used since it seems to slow things 
     down too much.

  Possible improvements: implement find_triple.  Change all the 
ModElTimes to ModElSpace, except for distinguished element.  Would need 
a new product method that can product a ModElTime and ModElSpace.
*/
#include "shallue.h"

NTL_CLIENT

/* In case we ever need to do binary searches again, this functor would 
allow us to create a compareM function that compares two ModElements
modulo M where M can change by creating a new instance of the class
Simply create as follows:  Mfunctor compareM(); compareM=Mfunctor(7);

credit to Mark Liffiton for coming up with the idea off the internet
*/
class Mfunctor{
public:
  Mfunctor() : M(to_ZZ(1)) {}
  Mfunctor(ZZ M) : M(M) {}
  bool operator()(ModElTime a, ModElTime b){ 
    if((a.to_ZZ() % M) < (b.to_ZZ() % M)) return true; else return false;
  }

  private:
    ZZ M;
};


/* Use shallue algorithm to construct a carmichael number (write more on this later) */

void shallue(vector<long>& primes, vector<long>& exponents, char* infileName1, char* infileName2, char* outfileName, bool verbose) {

  /* 3 timings.  full run, only phase2, and the algorithm except for S construction */
  time_t start, end, phase2_start, phase2_end, shallue_start, shallue_end;
  time(&start);

  int length = primes.size();  

  /* Build L = p_1^{h_1} dots p_r^{h_r} */

  ZZ L = to_ZZ(1);
  for (int i = 0; i < length; i++) {
    L *= power_ZZ(primes[i],exponents[i]);
  }
  if (verbose) {
    cout << "L = " << L << ".\n"; 
  }


  /* Compute N approx |S| */

  ZZ estimateN = getN(primes, exponents, L);
  if (verbose) {
    cout << "N = " << estimateN << ".\n";
  }


  /* Compute j smallest such that 
   *   N * prod_{i=j}^r  1 / (h_i + 1) > 2^(sqrt(log(|G|)))
   * We define overline{lambda} = prod_{i=1}^j p_i^{e_i} 
   * and G = (Z / overline{lambda} Z)*.
   */
  int j = length-1;

  //construct G
  ZZ G = L;
  long l = to_long(sqrt(NumBits(G)));
  l = l+1; // the computation of l is a lower bound, so +1 makes it the ceiling

  if(verbose) cout << "l = " << l << "\n";

  if(estimateN < power_ZZ(4,l)){
    cout << "Error in shallue.cpp!  We don't anticipate the algorithm will work since num primes is " << estimateN << " and the number needed is " << power_ZZ(4,l) << "\n";
  }

  // shrink j until the inequality flips, experiments revealed j 
  // was a little too small, so multiplied by log(N) 
  while ( estimateN > power_ZZ(4,l) ) {
    estimateN /= (exponents[j] + 1);
    G /= power_ZZ(primes[j],exponents[j]);
    j--;
    l = to_long(sqrt(NumBits(G)));
    l = l+1;

    if(verbose){
      cout << "j = " << j << " estimateN is " << estimateN << "\n";
      cout << " and number of primes needed is " << power_ZZ(4,l) << "\n";
    }
  }
  if(j < length-1) j++;


  if (verbose) {
    cout << "j = " << j << ".\n";
  }
 
  /* Compute l = sqrt{ log{ overline{ lambda } }
   * approx sqrt{ log { |G| } } */
  
  /* Start with o = overline{ lambda } */
  ZZ o = to_ZZ(1);
  for (int i = 0; i <= j; i++) {
    o *= power_ZZ(primes[i], exponents[i]);
  }
    // note: the +1 in the next line gives the ceiling
    // this l doesn't get used much.  Maybe I should just delete it
  l = SqrRoot(NumBits(o))+1; // l = log { overline{ lambda } }  + 1
  if (verbose) {
    cout << "l = " << l << ".\n"; }
  
  /* Before using ModElTime elements we must set primes and exponents */
  ModElTime::set(primes, exponents);

  /* Compute b = a^{-1}, where a is the product of all primes in infileName2, ie product file */
  ModElTime b = getb(L, infileName2);

  /* Declare and initialize S to have length + 1 vectors.
   * Then read in S from a file */
  vector< vector<ModElTime> > S(length + 1, vector<ModElTime>()); // Figure out the vector initalization statements sometime
  getS(S, infileName1, power_long(3,l), j+1);

  // print out number of elements in each bin of S
  if(verbose){
  cout << "\n";
    for(int i=0; i < S.size(); i++){
      cout << "S[" << i << "] has " << S.at(i).size() << " many elements\n";
    }
  }
  cout.flush();

  time(&shallue_start); //S is constructed, actual alg starts here

  /* Initial M value is the product of prime powers of Lambda we have skipped 
     over */
  ZZ M = to_ZZ(1);
  for (int i = j+1; i < length; i++) {
    M *= power_ZZ(primes[i], exponents[i]);
  }
  
  if(verbose) cout << "M = " << M << "\n"; 
 
  /* will contain the bins for residues of the current working M
     note that bin i contains the residue mod i, so that all bins i
     with gcd(i,workM) != 1 will be empty
  */
  vector< list<ModElTime> > products  = vector< list<ModElTime> >();

  if (verbose) {
    cout << "Phase 1.\n";
  }

  // based off of the Loh/Niebuhr paper, does the hard work of phase 1
  // including backtracking, modifying b until it has omega value j+1
  // note omega value is one more than the index of the corresponding prime
  if (!LohNiebuhr(primes, exponents, L, S, j+1, b, true, true)) {
    cout << "Error: Could not find a suitable product to get b in G.\n";
    exit(1);
  }

  if (verbose) {
    cout << "b = " << b.to_ZZ() << " now has omega value " << b.getOmega() << ".\n";
  }

  /* Used for storage*/
  ModElTime temp;
  
  /* */
  quad_float levelbits = sqrt(to_quad_float(NumBits(L/M)));
  ZZ nextM = getNextM(primes, exponents, M, levelbits);
  ZZ nextnextM = getNextM(primes, exponents, nextM, levelbits);
  long workM = getWorkM(primes, exponents, M,nextM);
  long nextworkM = getWorkM(primes, exponents, nextM, nextnextM);

  // the initial number of lists needed in products is workM
  for(int i=0; i < workM; i++){
    products.push_back( list<ModElTime>() );
  }

  if(verbose) cout << "workM is " << workM << "\n";

  /*
   * Add more elements to products so that our probably of 
   * success is higher.  In particular,  use Loh-Niebuhr
   * method to find products for primes in S that are in G
   */
/*
  for (int i = length; i > j; i--) {
    while (!S.at(i).empty()) {
      // Remove the first element in S, say p. Then try to
       // find a suitable product that has p in G. 
      temp = S.at(i).back();
      S.at(i).pop_back();
      if (LohNiebuhr(primes, exponents, L, S, j, temp, false, false)) {
        products.at(temp.to_ZZ() % workM).push_front(temp);
      }
    }  
  }
  if (verbose) {
    cout << "Pushed all primes with omega value greater than " << j <<  " into G.\n";
  }
 */
  /* Copy each S[i] where 0 LEQ i LEQ j+1 and clear S 
     Note that we want primes with omega value at most j+1, since omega = j+1
     corresponds to the prime at position j
  */
  for (int i = 0; i <= j+1; i++) { 
    for (int k = 0; k < S.at(i).size(); k++) {
      temp = S.at(i).at(k);
      products.at(temp.to_ZZ() % workM).push_front(temp);
    }
  }
  S.clear();

  /* Pair elements in product by increasing modulii until we
   * reach the identity */
  if (verbose) {
    cout << "Phase 2.\n";
  }

  time(&phase2_start);

  if (!phase2(primes, exponents, b, products, M, L)) {
    exit(1); //Fail
  }
  cout << "Solution found and printed to file.\n";

  time(&phase2_end);
  time(&shallue_end);

  ofstream outfile;
  outfile.open("solution");
  for (int i = 0; i < b.history.size(); i++) {
    outfile << b.history.at(i).toZZ()+1 << "\n";
  }
  outfile.close();

  time(&end);
  cout << "\n";
  cout << "Shallue took " << difftime(end,start) << " seconds.\n";
  cout << "phase2 took " << difftime(phase2_end,phase2_start) << " seconds.\n";
  cout << "Shallue without S construction took " << difftime(shallue_end,shallue_start) << " seconds.\n";
}

/*
 * Find a product bp_1p_2 dots p_u such that bp_1p_2 dot p_u has omega value
 * less than or equal to j, p_1,..,p_u in S (u LEQ r-j).  
 * Uses the basic structure of the algorithm given by Loh and Niebuhr in 
 * A New Algorithm for Constructin Carmichael Numbers.
 * Note: This is urrently is written to be mostly space efficient.  There
 * many small changes we could make to increase time efficiency.  
 */

// Added a parameter called dist (distinguished) which changes how we remove the primes in a solution from s

bool LohNiebuhr(vector<long>& primes, vector<long>& exponents, ZZ& L, vector< vector<ModElTime> >& S, int j, ModElTime& b, bool dist, bool verbose) {

  int length = primes.size();

  int u = 0; // u is the current level

  /* We use position as a vector of iterators that are pointing
   * to possible primes to add to the product */
  vector< vector<ModElTime>::iterator > position = vector< vector<ModElTime>::iterator>();

  /* We use a as a temporary variable which always
   * is the target, i.e. b = b^{-1} */
  ModElTime a;

  /* We use omega as a temporary variable which
   * always stores omega values */
  int omega;

  /* Do a backtracking search for suitable product
   * that has b in G. */

  omega = b.getOmega();
  if (omega <= j) return true;
 
  if(verbose) cout << "inital LohNie: j = " << j << " and omega = " << omega << "\n";
 
 /* note use of alternate constructor so a does not have a history*/
  a = ModElTime(InvMod(b.to_ZZ(), L),-1);
  position.push_back(std::lower_bound(S[omega].begin(), S[omega].end(), a, compareResidue));

  while ( u >= 0) {

    while (nextPrime(position.at(u), b, S, L)) {

      if (verbose) {
        cout << "\nu = " << u << ".\n";
        cout << "b = " << b.to_ZZ() << ".\n";
        cout << "Select: " << (*(position.at(u))).to_ZZ() << ".\n";
      }

      /* Set b := b * p where p = *position.back() and
       * it is safe because nextPrime gaurentees it */
      b = product(b, *(position.at(u)), L);


      /* Test to see if b is now in G */
      omega = b.getOmega();
      if ( omega <= j ) {
        /* Remove primes in solution product from S*/
        if (dist) {
          for (int i = 0; i < b.history.size(); i++) {
            a = ModElTime(b.history.at(i).toZZ() + 1, -1); //Don't need -1 but we also dont need a history
            S[a.getOmega()].erase(position.at(i));
          }
        }
        else {
          for (int i = 1; i < b.history.size(); i++) {
            a = ModElTime(b.history.at(i).toZZ() + 1, -1);
            S[a.getOmega()].erase(position.at(i-1));
          }
        }
        return true;
      }

      /* Advance to next level */
      u++;

      /* Set position for next level */
      a = ModElTime(InvMod(b.to_ZZ(), L), -1);

      position.push_back(std::lower_bound(S[omega].begin(), S[omega].end(), a, compareResidue));    
    }

    /* Backtrack */

    if (verbose) {
      cout << "Backtrack.\n";
    }
      
    if (u > 0) {
      removeLast(b);
    }
    
    position.pop_back();
    u--;
    if (u >= 0) {
      position.at(u)++;
    }
  }

  /* No suitable product was found in backtracking search */
  return false; 
}

ModElTime getb(ZZ& L, char* infileName) {
  
  ifstream infile;
  infile.open(infileName);

  /* Build product equal to the product of the numbers in a file */

  ZZ product = to_ZZ(1);
  ZZ temp;  

  while (SkipWhiteSpace(infile)) {
    infile >> temp;
    product = MulMod(product, temp, L);
  }

  infile.close();

  return ModElTime(InvMod(product, L), -1); // Return b^{-1} 
}

void getS(vector< vector<ModElTime> >& S, char* infileName, long limit, long j) {

  /* Open a file with name fileName */
  ifstream infile;
  infile.open(infileName);

  /* The following are used for storage */
  ZZ temp;
  ModElTime p;
  int omega;

  /* Read in each primes and put it into the 
   * appropriate vector in S */
  while (SkipWhiteSpace(infile)) {
    infile >> temp;
    p = ModElTime(temp);
    omega = p.getOmega();

    // only apply the limit to primes with large omega values
    if (omega <= j || S[omega].size() < limit) {
      S[omega].push_back(p);
    }
  }

  /* Sort each vector in S by residue */
  for (int i = 0; i < S.size(); i++) {
    sort (S[i].begin(), S[i].end(), compareResidue);
  }

  infile.close();
}

bool nextPrime(vector<ModElTime>::iterator& it, ModElTime& b, vector< vector<ModElTime> >& S, ZZ& L) {
  
  int omega = b.getOmega();
  ModElTime a = ModElTime(InvMod(b.to_ZZ(), L),-1);
 
  if (!(it < S[omega].end())) return false;
  if ((*it).residue(omega) == a.residue(omega)) return true;
  return false;
}

bool compareResidue(ModElTime a, ModElTime b) {

  int omega = a.getOmega();

  /* Error Checking */
  if (b.getOmega() != omega) {
    cout << "Error, compare residue was called for elements with different omega values.\n"; 
  }
  
  if (a.residue(omega) < b.residue(omega)) return true;
  else return false;
}

/* Written by Andrew Shallue */
ZZ getN(vector<long>& primes, vector<long>& exponents, ZZ& L) {

  int length = primes.size();

  ZZ phi = to_ZZ(1);
  for (int i = 0; i < length; i++) {
    phi *= power_long(primes[i], exponents[i]) - power_long(primes[i], exponents[i] - 1);
  }

  quad_float correction = to_quad_float(1);
  for (int i = 0; i < length; i++) {
    correction *= (exponents[i] + ((double) primes[i] - 2) / (primes[i] - 1));
  }

  quad_float probability = log(sqrt(to_quad_float(2*L)));

  return FloorToZZ((to_RR(L) * to_RR(correction)) / (to_RR(phi)* to_RR(probability)));
}

/*
The M that gets passed in is the initial M that is computed with respect to j
*/

bool phase2(vector<long>& primes, vector<long>& exponents, ModElTime& b, vector< list<ModElTime> >& products, ZZ& M, ZZ& L) {
 
  bool verbose = true; // Move this later
 
  /* levelbits is l, the number of bits the group shrinks by each step*/
  quad_float levelbits = sqrt(to_quad_float(NumBits(L/M))); //would be better to pass this in maybe?

  /* We use temp to store next level of products */
  vector< list<ModElTime> > temp;

  /* Iterators are used for storage */
  list<ModElTime>::iterator pos1, pos2, end1, end2;

  /* We use residue to store residue mod M */
  long residue;

  /* We use a for temporary storage */
  ModElTime a;

  /* Used for storing sizes of vectors when we update them v*/
  long nextsize, tempsize;

  /* Used for determining wether to continue searching or not*/
  bool proceed;

  /* Set up the family of M's (we recompute these to simplying
  * parameters passed into phase 2. */
  ZZ nextM = getNextM(primes, exponents, M, levelbits);
  ZZ nextnextM = getNextM(primes, exponents, nextM, levelbits);
  long workM = getWorkM(primes, exponents, M, nextM);
  long nextworkM = getWorkM(primes, exponents, nextM, nextnextM);
  
  /* A counter for levels */
  int level = 0;

  while (M < L) {  
  
    level++;
    if (verbose) {
      long products_total = 0;
      for(long i = 0; i < products.size(); i++){
        products_total += products.at(i).size();
      }

      cout << "Pairing mod " << workM << " at level " << level << ", we have " << products_total << " elements to work with.\n";
   /* 
      cout << "\n";
      for(long i = 0; i < products.size(); i++){
        cout << products.at(i).size() << " ";
      }
      cout << "\n";
    */
    }

    /* Set proceed to true for this level */
    proceed = true;

    /* Reconcile size of temp with nextworkM */
    nextsize = nextworkM;
    tempsize = temp.size();
    if (tempsize < nextsize) {
      for (int i = 0; i < nextsize - tempsize; i++) {
        temp.push_back( list<ModElTime>() );
      }
    }
    else if (tempsize > nextsize) {
      for (int i = 0; i < tempsize - nextsize; i++) {
        temp.pop_back();
      }
    }
 
    /* Find b's inverse (mod workM) */
    
    residue = InvMod(b.to_ZZ() % workM, workM); //
    
    if (b.to_ZZ() % workM == 1) {} //Do nothing
    else {
      /* Check to see if a pair exists to put b into the next level*/
      residue = InvMod(b.to_ZZ() % workM, workM); //
      if (!products.at(residue).empty()) {
        b = product(b, *(products.at(residue).begin()), L); //Ugly
        products.at(residue).erase(products.at(residue).begin());
      }
      else {
        /* Search for a triple */
        for (int i = 1; i < products.size() && proceed; i++) { //Ugly
          if (!products.at(i).empty()) {
            residue = InvMod(b.to_ZZ() * i % workM, workM);
            if (!products.at(residue).empty()) {
              proceed = false;
              b = product(b, *(products.at(i).begin()), L);
              b = product(b, *(products.at(residue).begin()), L);
              products.at(i).erase(products.at(i).begin()); //Ugly
              products.at(residue).erase(products.at(residue).begin());
            }
          }
          /* If there were no triples then we have failed */
          if (proceed) {
            cout << "Error: could not get b into G\n";
            return false;
          }
        }
      }
    }

    /* Push every product in products[1] (products
     * with residue 1 (mod M))  into next level */
 
    /* Set up iterators */
    pos1 = products.at(1).begin();
    end1 = products.at(1).end();
    while (pos1 != end1) {
      /* Add product to temp */
      a = *pos1;
      temp.at(a.to_ZZ() % nextworkM).push_back(a); //Ugly
      /* Remove *pos1 from products */
      products.at(1).erase(pos1);
      pos1 = products.at(1).begin();
    }

    /* Find pairs that product to 1 (mod M) */

    /* Pair each product in products[i] with
     * an inverse for it (mod M). If we let
     * residue be it's inverse, then inverses
     * for elements in products[i] can be found
     * in products[residue]. 
     */
    for (long i = 2; i < workM - 1; i++) {
      if (products.at(i).size() > 0) {
        residue = InvMod(i, workM);
        if (i < residue) {
          /* Set up iterators */
          pos1 = products.at(i).begin();
          end1 = products.at(i).end();
          pos2 = products.at(residue).begin();
          end2 = products.at(residue).end();

          while (pos1 != end1 && pos2 != end2) {
            /* Add product of *pos1 and *pos2 to temp */
            a = product(*pos1, *pos2, L);

            temp.at(a.to_ZZ() % nextworkM).push_back(a);
            /* Remove *pos1 and *pos2 from products */
            products.at(i).erase(pos1);
            products.at(residue).erase(pos2);
            /* Update iterators to point to the beginning of the list again. */
            pos1 = products.at(i).begin();
            pos2 = products.at(residue).begin();
          } //endwhile
        }//endif
      }
    }

    /* Pair every product in products[M-1] with another
     * products from products[M-1]. This is because
     * (M-1)^2 = 1 (mod M). */

    /* Set iterators */
    pos1 = products.at(workM-1).begin();
    pos2 = products.at(workM-1).begin();
    pos2++; //Ugly
    end1 = products.at(workM-1).end();
    while (pos1 != end1 && pos2 != end1) {
      /* Add produt of *pos1 and *(pos1 + 1) to temp */
      a = product(*pos1, *pos2, L);
      temp.at(a.to_ZZ() % nextworkM).push_back(a);
      products.at(workM-1).erase(pos1);
      products.at(workM-1).erase(pos2);
      /* Update iterators*/
      pos1 = products.at(workM-1).begin();
      pos2 = products.at(workM-1).begin();
      pos2++; //Ugly
    }

    /* Find triples that product to 1 (mod M) */
    // finish this 

    /* After all the matching there may still be elements left over.  They 
       need to be discarded */
    for(long j = 0; j < products.size(); j++){
      products.at(j).clear();
    }

    // after swapping products is ready for the next level, temp is empty
    // but with the wrong number of bins.  This is reconciled at the top
    // of the next level
    swap(products, temp); 

    /* Update M's */
    M = nextM;
    nextM = nextnextM;
    nextnextM = getNextM(primes, exponents, nextM, levelbits);
    workM = getWorkM(primes, exponents, M, nextM);
    nextworkM = getWorkM(primes, exponents, nextM, nextnextM);

  }//end big while
  
  return true;
}// end phase 2

/*
  Written by Andrew Shallue
 
  There's a built in fudge factor here: I pretend that going
  from no restriction (mod p) to restricting to 1 (mod p) is
  a subgroup of index p, but in reality its a subgroup of index
  p - 1.

  Also, currently each set is guaranteed to be too big, though
  hopefully not by too much.  Testing might reveal that I should
  undershoot at each step instead of overshooting.
*/
ZZ getNextM(vector<long>& primes, vector<long>& exponents, ZZ& M, quad_float& levelbits) {

  // M will not be divisible by small primes.  Find first prime dividing it
  long currentprime = 2;
  int primepos = 0;
  while(M % primes.at(primepos) != 0 && primepos < primes.size()-1){  
    primepos++;
  }
  currentprime = primes.at(primepos);

  // now find the power of that prime that divides M
  long currentexp = 0;
  ZZ tempM = M;
  while(tempM % primes.at(primepos) == 0){
    tempM = tempM / primes.at(primepos);
    currentexp++;
  }

  quad_float temp_bits = to_quad_float(0);

  tempM = M;

  // special case; when M = 1 just start at the beginning
  if(tempM == 1){
     currentexp = exponents.at(primepos);
     temp_bits = temp_bits + log(to_quad_float(currentprime))/log(2);
     tempM = tempM * currentprime;
  }

  /* Temp_bits records how many bits have been added to M at each step.
     each step of the while loop adds a single prime to M
     When we have added at least level_bits number of bits we stop
  */
  while(temp_bits < levelbits){
    // first check if we're already done, i.e. M = Lambda
    if(currentprime == 2 && currentexp == exponents.at(0)){
        //cout << "M shouldn't have to be changed, already max power of 2\n";
        // in real program maybe just output M? Yes!
        break;
    }

    // next deal with the case that we've maxed out current prime
    if(currentexp == exponents.at(primepos)){
      primepos--;
      currentprime = primes.at(primepos);
      currentexp = 0;
    }
 
    //increment exponent.  In the 2 case an exponent of 1 isn't helpful
    currentexp++;
    if(currentprime ==2 && currentexp == 1){
      currentexp++;
      tempM = tempM * 2;
    }
    //number of bits added is log_2(currentprime)
    temp_bits = temp_bits + log(to_quad_float(currentprime))/log(2);

    //M is enlarged by the multiple of the new prime
    tempM = tempM * currentprime;

  } //endwhile

  return tempM;
}

/* Requirement:  M < N and M needs to divide N

  The idea of this subroutine is to return N / M, i.e. the extra factors that 
separate the M from the previous level (M) and the M from the current level (N)
  However, an annoying thing is that if L is divisible by 7^2, one workM has 
a factor of 7, the next needs to be divisible by 7^2, not 7.  For if elements 
are congruent to 1 mod 7, workM needs to be divisible by 7^2 so that we are 
looking for products congruent to 1 mod 7^2.
*/

long getWorkM(vector<long>& primes, vector<long>& exponents, ZZ& M, ZZ& N){
  // make sure N and M fit the correct input profile
  if(N < M || N % M != 0) cout << "Error in getWorkM! Integer input not correct\n";

  long workM = 1;
  long divM = to_long(N/M);

  int length = primes.size();
  for(int i=0; i < length; i++){
    if(divM % primes[i] == 0){
      workM = workM * to_long(GCD(power_ZZ(primes[i],exponents[i]), N));
    }
  }

  return workM;

}