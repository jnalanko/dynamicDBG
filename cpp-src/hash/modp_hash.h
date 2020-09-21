#ifndef GEN_HASH
#define GEN_HASH

#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;
#include "HashUtil.cpp"

using namespace std;

typedef boomphf::SingleHashFunctor<u_int64_t> hasher_t;
typedef boomphf::mphf<u_int64_t, hasher_t> boophf_t;

//For now, this is our big int class. If we need more than 1024 bits, we can increase
//LARGE_BITS should be a power of 2
#define LARGE_BITS 1024
typedef number<cpp_int_backend<LARGE_BITS, LARGE_BITS, unsigned_magnitude, checked, void>> largeUnsigned;

typedef u_int64_t HashInt;
/**
 * Take in a set of k-mers, generate a hash function
 */
class generate_hash
{

public:
    HashInt n_kmer;  //the number of k-mers
    unsigned k_kmer; //the lengths of the k-mers (max 32)

    vector<string> kmer_data; // pointer to kmer_data TODO: should read from file

    std::unordered_set<u_int64_t> KRHash; //image of k-mers through Karp-Rabin hash function
                                          //std::vector<HashInt> KRHash_vec; //vector form

    //the image of our k-mers under our Karp-Rabin hash function
    //vector<HashInt> KR_hash_val;

    /* 
	 * Will store the precomputed powers of 
	 * Needs to be a vector will have k distinct powers
	 * Using 128 bit type to prevent overflow
	 * Stores powers in order r^1, r^2, ..., r^k
	 */
    vector<largeUnsigned> powersOfR;
    vector<largeUnsigned> powersOfRModP;

    boophf_t *bphf; //MPHF we will generate

    HashInt r;     // the base for our Karp-Rabin Hash function
    HashInt rinv;  //inverse of r modulo Prime
    HashInt Prime; // the prime for our Karp-Rabin Hash function

    const static short sigma = 4; // alphabet size

    /**
         * Create hash function out of n k-mers of length k
         */
    generate_hash(unordered_set<kmer_t> &kmers, HashInt n, unsigned k)
    {
        //	  std::srand(std::time(NULL));
        std::srand(0);

        construct_hash_function(kmers, n, k);
    }

    generate_hash()
    {
        //default constructor
        std::srand(std::time(NULL));
    }

    void construct_hash_function(unordered_set<kmer_t> &kmers, HashInt n, unsigned k)
    {
        this->n_kmer = n; // number of k-mers
        this->k_kmer = k; // length of each k-mer

        BOOST_LOG_TRIVIAL(info) << "Constructing the hash function ...";
        build_KRHash(kmers);        // build KR hash function
        build_minimalPerfectHash(); // build minimal perfect hash function
        int256_t kr = findInverse(static_cast<int256_t>(r), static_cast<int256_t>(Prime));

        while (kr < 0)
            kr = kr + Prime;

        if (kr > Prime)
            kr = kr % Prime;

        if ((kr * r) % Prime == 1)
        {
            BOOST_LOG_TRIVIAL(info) << "r-inverse correctly computed modulo prime.";
        }
        else
        {
            BOOST_LOG_TRIVIAL(fatal) << "r-inverse incorrectly computed modulo prime!";
            exit(1);
        }

        rinv = static_cast<HashInt>(kr);
    }

    /*
	 * Once we know k (k_kmer) and r
	 * we can precompute the powers of r
	 */
    void precomputePowers()
    {
        largeUnsigned ri;
        powersOfR.clear();
        //NEED 1 to k. Not 0 to (k - 1)
        for (unsigned i = 1; i <= k_kmer; ++i)
        {
            ri = mypower(r, i);
            powersOfR.push_back(ri);
        }
    }

    /*
	 * Once we know k (k_kmer) and r
	 * we can precompute the powers of r
	 */
    void precomputePowers_mod()
    {
        largeUnsigned ri;
        powersOfRModP.clear();
        //NEED 1 to k. Not 0 to (k - 1)
        for (unsigned i = 1; i <= k_kmer; ++i)
        {
            ri = mypower_mod(r, i);
            powersOfRModP.push_back(ri);
        }
    }

    /**
         * Find the hash value of a k-mer
         */
    u_int64_t get_hash_value(const kmer_t &seq)
    {
        u_int64_t krv = generate_KRHash_val(seq, k_kmer);
        u_int64_t res = this->bphf->lookup(krv); // still need only 64 bits for kmer_t
        return res;
    }

    /**
         * Find the hash value of a k-mer
	 * Allows f( v ) notation
         */
    HashInt operator()(const kmer_t &seq)
    {
        return get_hash_value(seq);
    }

    // Task4: generate_KRHash_val
    // data is a k-mer
    // k is the length of the k-mer
    // r is the base
    // P is the prime
    void build_KRHash(unordered_set<kmer_t> &kmers)
    {

        BOOST_LOG_TRIVIAL(info) << "Constructing Karp-Rabin hash function ...";

        HashInt v; // holder for KRH value

        // prime we will mod out by
        const HashInt tau = 1;
        BOOST_LOG_TRIVIAL(info) << "Theoretical prime lower bound: " << tau * k_kmer * n_kmer * n_kmer;

        //Prime = getPrime(max((HashInt)this->sigma, (HashInt)tau*k_kmer*n_kmer*n_kmer));

        //Having problem with overflows because of large primes
        //Even though the theoretical bound is above, let's try smaller ones.
        double smallerPrime = n_kmer * n_kmer / 5.0;
        Prime = getPrime((HashInt)smallerPrime);

        BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;

        // keep generating new base until we find one that is injective over our k-mers
        bool f_injective;

        unsigned n_failures = 0; //After a few failures, double the prime
        do
        {
            if (n_failures == 5)
            {
                BOOST_LOG_TRIVIAL(info) << "Trying a larger prime... ";
                Prime = getPrime(Prime * 2);
                n_failures = 0;
                BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;
            }

            f_injective = true; //assume f is injective until evidence otherwise
            this->r = randomNumber((HashInt)1, Prime - 1);
            //Once we have a candidate base r
            //we should avoid recomputing its powers all the time
            precomputePowers();

            for (unordered_set<kmer_t>::iterator
                     it1 = kmers.begin();
                 it1 != kmers.end();
                 ++it1)
            {
                v = generate_KRHash_val(*it1, k_kmer);
                //		BOOST_LOG_TRIVIAL(trace) << "hash of kmer: " << v;
                if (this->KRHash.find(v) == this->KRHash.end())
                {
                    // this is a new value
                    this->KRHash.insert(v);
                }
                else // not injective
                {
                    BOOST_LOG_TRIVIAL(trace) << "Base " << this->r << " with prime "
                                             << Prime << " failed injectivity.";
                    this->KRHash.clear(); // clear it out and start over
                    f_injective = false;
                    ++n_failures;
                    break;
                }
            }
        } while (!f_injective);

        BOOST_LOG_TRIVIAL(info) << "Base " << this->r << " with prime " << Prime << " is injective.";

        return;
    }

    /*
	 * Computes powers with HashInt and integer exponents
	 */
    largeUnsigned mypower(const HashInt &base, unsigned exponent)
    {
        largeUnsigned rvalue(1);
        while (exponent > 0)
        {
            rvalue *= static_cast<largeUnsigned>(base);
            --exponent;
        }

        return rvalue;
    }

    /*
	 * Computes powers modulo P
	 */
    largeUnsigned mypower_mod(const HashInt &base, unsigned exponent)
    {
        HashInt rvalue(1);
        while (exponent > 0)
        {
            rvalue = (rvalue * base) % Prime;
            --exponent;
        }

        return rvalue;
    }

    /**
         * Given a kmer, find out its KRH using base r and prime P
         */
    HashInt generate_KRHash_val(const kmer_t &kmer,
                                const unsigned &k)
    {

        //	  BOOST_LOG_TRIVIAL(trace) << "Generating KRHash val";
        //use 128 bits to prevent overflow
        largeUnsigned val = 0; // what will be the KRH value

        // go through each bp and add value
        for (unsigned i = 0;
             i < k;
             ++i)
        {
            // val += baseNum(kmer.at(i)) * pow(r, i);
            val +=
                static_cast<largeUnsigned>(access_kmer(kmer, k, static_cast<unsigned>(i))) *
                powersOfR[i]; //powersOfR[i] = r^{i + 1}
        }

        val = val % Prime;

        return static_cast<HashInt>(val);
    }

    /**
         * Given a kmer, find out its KRH
	 * MODULO Prime
	 * Need to use 128 bits since 4 * x might overflow
         */
    u_int64_t generate_KRHash_val_mod(const kmer_t &kmer,
                                      const unsigned &k)
    {
        uint128_t val = 0; // what will be the KRH value

        // go through each bp and add value
        for (unsigned i = 0;
             i < k;
             ++i)
        {
            // val += baseNum(kmer.at(i)) * pow(r, i);
            val = val + ((access_kmer(kmer, k, i) *
                          static_cast<uint128_t>(powersOfRModP[i])));
            val = val % Prime;
        }
        val = val % Prime;

        return static_cast<u_int64_t>(val);
    }

    /**
         * Given a kmer, find out its KRH using base r and prime P
	 * HOWEVER: does not mod out by P. So returns large unsigned
         */
    largeUnsigned generate_KRHash_raw(const kmer_t &kmer,
                                      const unsigned &k)
    {
        //	  BOOST_LOG_TRIVIAL(trace) << "Generating KRHash val";
        //use 128 bits to prevent overflow
        largeUnsigned val = 0; // what will be the KRH value

        // go through each bp and add value
        for (unsigned i = 0;
             i < k;
             ++i)
        {
            // val += baseNum(kmer.at(i)) * pow(r, i);
            val +=
                static_cast<largeUnsigned>(access_kmer(kmer, k, static_cast<unsigned>(i))) *
                powersOfR[i]; //powersOfR[i] = r^{i + 1}
        }

        return val;
    }

    /*
	 * This function takes as input a Karp-Rabin value (KR_val)
	 * Then updates it by subtracting the value from 'first' character source kmer
	 * Then dividing by r (at this point, it has shifted last k-1 characters up)
	 * Finally adding the last term corresponding to the 'last' character
	 *
	 * target k-mer is OUT neighbor of source k-mer
	 *
	 */
    void update_KRHash_val_OUT(largeUnsigned &KR_val, //KR hash of source kmer
                               const unsigned &first, //character at front of source k-mer
                               const unsigned &last)
    { //last character in target k-mer
        //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by OUT...";
        //	   BOOST_LOG_TRIVIAL(debug) << "First of source: " << first;
        //	   BOOST_LOG_TRIVIAL(debug) << "Last of target: " << last;

        //	   largeUnsigned before_div = KR_val;
        //	   BOOST_LOG_TRIVIAL(debug) << "Value before division: " << before_div;

        //	   BOOST_LOG_TRIVIAL(debug) << "Division check: " << (KR_val * static_cast< largeUnsigned >( r ) == before_div );

        //	   KR_val = before_div;
        //largeUnsigned q;
        //	   largeUnsigned rem;
        //	   divide_qr( KR_val, static_cast< largeUnsigned >( r ),  q, rem );

        //	   BOOST_LOG_TRIVIAL(debug) << "Division check 2: " << (q * static_cast< largeUnsigned >( r ) == before_div );

        //	   BOOST_LOG_TRIVIAL(debug) << "Remainder: " << rem;
        KR_val = KR_val / static_cast<largeUnsigned>(r);
        KR_val = KR_val - static_cast<largeUnsigned>(first);
        KR_val = KR_val + static_cast<largeUnsigned>(last) * powersOfR[k_kmer - 1]; // last * r^k
    }

    /*
	 * This function takes as input a Karp-Rabin value (KR_val)
	 * Then updates it by subtracting the value from 'first' character source kmer
	 * Then dividing by r (at this point, it has shifted last k-1 characters up)
	 * Finally adding the last term corresponding to the 'last' character
	 *
	 * target k-mer is OUT neighbor of source k-mer
	 *
	 */

    typedef number<cpp_int_backend<1024, 1024, unsigned_magnitude, checked, void>> moduloInt;
    largeUnsigned update_KRHash_val_OUT_mod(largeUnsigned &KR_val, //KR hash of source kmer (mod P)
                                            const unsigned &first, //character at front of source k-mer
                                            const unsigned &last)
    { //last character in target k-mer
        //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by OUT(mod)...";
        //	   BOOST_LOG_TRIVIAL(debug) << "First of source: " << first;
        //	   BOOST_LOG_TRIVIAL(debug) << "Last of target: " << last;
        moduloInt rinv = this->rinv;
        moduloInt Prime = this->Prime;
        moduloInt kr = KR_val;
        moduloInt llast = last;
        moduloInt ffirst = first;
        moduloInt rk = powersOfRModP[k_kmer - 1];
        moduloInt four = 4;

        moduloInt sub_val = four * Prime - ffirst * r;
        kr = (kr + sub_val);
        kr = (kr * rinv);
        kr = (kr + llast * rk);

        //	   moduloInt q, rem;
        //	   divide_qr( kr, Prime, q, rem );
        //	   HashInt kr2 = integer_modulus( kr, this->Prime );
        kr = kr % Prime;

        //	   while (kr > Prime)
        //	      kr = kr - Prime;

        //	   kr = kr - (kr / Prime)*Prime;

        return static_cast<largeUnsigned>(kr);
    }

    /*
	 * This function takes as input a Karp-Rabin value (KR_val) modulo Prime
	 *
	 * target k-mer is IN neighbor of source k-mer
	 */
    largeUnsigned update_KRHash_val_IN_mod(largeUnsigned &KR_val, //KR hash of source kmer
                                           const unsigned &first, //character at front of target k-mer
                                           const unsigned &last)
    { //last character in source k-mer
        //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by IN(mod)...";
        //	   BOOST_LOG_TRIVIAL(debug) << "First of target: " << first;
        //	   BOOST_LOG_TRIVIAL(debug) << "Last of source: " << last;

        moduloInt r = this->r;
        moduloInt Prime = this->Prime;
        moduloInt kr = KR_val;
        moduloInt llast = last;
        moduloInt ffirst = first;
        moduloInt rk = powersOfRModP[k_kmer - 1];
        moduloInt four = 4;

        moduloInt sub_val = four * Prime - llast * rk;

        kr = (kr + sub_val); // last * r^k
        kr = (kr * r);
        kr = (kr + ffirst * r);

        //	   moduloInt q, rem;
        //	   divide_qr( kr, Prime, q, rem );
        //	   HashInt kr2 = integer_modulus( kr, this->Prime );
        kr = kr % Prime;

        //	   while (kr > Prime)
        //	      kr = kr - Prime;

        //	   kr = kr - (kr / Prime)*Prime;

        return static_cast<largeUnsigned>(kr);
    }

    /*
	 * This function takes as input a Karp-Rabin value (KR_val)
	 *
	 * target k-mer is IN neighbor of source k-mer
	 */
    void update_KRHash_val_IN(largeUnsigned &KR_val, //KR hash of source kmer
                              const unsigned &first, //character at front of target k-mer
                              const unsigned &last)
    { //last character in source k-mer
        //	   BOOST_LOG_TRIVIAL(debug) << "Updating a KR value by IN...";
        //	   BOOST_LOG_TRIVIAL(debug) << "First of target: " << first;
        //	   BOOST_LOG_TRIVIAL(debug) << "Last of source: " << last;

        KR_val = KR_val - last * powersOfR[k_kmer - 1]; // last * r^k
        KR_val = KR_val * r;
        KR_val = KR_val + first * r;
    }

    /*
	 * Looks up the minimal perfect hash value, given the Karp-Rabin (raw value)
	 * KR raw value means not modded out by the prime yet.
	 */
    HashInt perfect_from_KR(const largeUnsigned &KR_val)
    {
        largeUnsigned KR2 = KR_val % Prime;
        return this->bphf->lookup(static_cast<u_int64_t>(KR2));
    }

    /**
         * Build a minimal perfect hash function on the set of integers that our kmers are
         * mapped to via KRH
         */
    void build_minimalPerfectHash()
    {

        //std::sort(KR_hash_val, KR_hash_val+n_kmer);
        //HashInt jj = 0;
        //for (int ii = 1; ii < n_kmer; ii++) {
        //    if (KR_hash_val[ii] != KR_hash_val[jj])
        //        KR_hash_val[++jj] = KR_hash_val[ii];
        //}
        //printf("Found %lli duplicated items from KR_hash_val.  \n", n_kmer-(jj + 1) );

        //auto data_iterator = boomphf::range(static_cast<const HashInt*>(KR_hash_val), static_cast<const HashInt*>(KR_hash_val+n_kmer));

        //            bphf = new boomphf::mphf<HashInt, hasher_t>(n_kmer, data_iterator, nthreads, gammaFactor);

        BOOST_LOG_TRIVIAL(info) << "Building minimal perfect hash function ...";

        std::vector<HashInt> KRHash_vec = std::vector<HashInt>(this->KRHash.begin(),
                                                               this->KRHash.end());

        // MPHF for our KRHash function values
        this->bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, KRHash_vec, 4, 2.0, true, false);

        BOOST_LOG_TRIVIAL(info) << "Minimal perfect hash function created with "
                                << (float)(bphf->totalBitSize()) / n_kmer << " bits per element.";
    }
};

#endif
