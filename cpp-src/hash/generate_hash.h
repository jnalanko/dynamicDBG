#ifndef GEN_HASH
#define GEN_HASH
//#define NDEBUG

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <boost/multiprecision/cpp_int.hpp>
#include "sux/function/RecSplit.hpp"

using namespace boost::multiprecision;
#include "HashUtil.cpp"

using namespace std;

typedef uint128_t kr_hash_t;
typedef sux::function::RecSplit<8> RecSplit;

/**
 * Take in a set of k-mers, generate a hash function
 */
class generate_hash
{

public:
    static const int BUCKET_SIZE = 100;
    uint64_t n_kmer_orig; // the number of k-mers not considering deletions/insertions
    uint64_t max_hash;    // the max value this hash function has ever been able to return
                          // used to produce and check for invalid hash values
    unsigned k_kmer;      //the lengths of the k-mers (max 32)

    /* 
    * Will store the precomputed powers of 
    * Needs to be a vector will have k distinct powers
    * Using 256 bit type to prevent overflow
    * Stores powers in order r^1, r^2, ..., r^k
    */
    vector<uint256_t> powersOfRModP;

    //image of k-mers through Karp-Rabin hash function
    //discarded after use
    std::unordered_set<kr_hash_t> KRHash;

    // A map of kmers to their hash values for newly added
    // nodes. A hack until we have a dynamic hash function
    // implementation
    std::map<kmer_t, u_int64_t> new_nodes;

    RecSplit *recsplit;

    uint256_t r;     // the base for our Karp-Rabin Hash function
    uint256_t rinv;  // inverse of r modulo Prime
    uint256_t Prime; // the prime for our Karp-Rabin Hash function

    const static short sigma = 4; // alphabet size

    void save(ostream &of)
    {
        of.write((char *)&n_kmer_orig, sizeof(n_kmer_orig));
        of.write((char *)&max_hash, sizeof(max_hash));
        of.write((char *)&k_kmer, sizeof(k_kmer));
        of.write((char *)&r, sizeof(r));
        of.write((char *)&rinv, sizeof(rinv));
        of.write((char *)&Prime, sizeof(Prime));

        // Write recsplit
        of << (*recsplit);
    }

    void load(istream &of)
    {
        of.read((char *)&n_kmer_orig, sizeof(n_kmer_orig));
        of.read((char *)&max_hash, sizeof(max_hash));
        of.read((char *)&k_kmer, sizeof(k_kmer));
        of.read((char *)&r, sizeof(r));
        of.read((char *)&rinv, sizeof(rinv));
        of.read((char *)&Prime, sizeof(Prime));
        recsplit = new RecSplit(); // Todo: leaks memory
        of >> (*recsplit);
        precomputePowers_mod();
    }

    /**
    * Create hash function out of n k-mers of length k
    */
    generate_hash(unordered_set<kmer_t> &kmers, u_int64_t n, unsigned k)
    {
        std::srand(0);

        construct_hash_function(kmers, n, k);
    }

    generate_hash()
    {
        //default constructor
        std::srand(0);
    }

    generate_hash(istream &of)
    {
        std::srand(0);
        load(of);
    }

    void construct_hash_function(unordered_set<kmer_t> &kmers, u_int64_t n, unsigned k)
    {
        this->n_kmer_orig = n;  // number of k-mers
        this->max_hash = n - 1; // max possible hash value
        this->k_kmer = k;       // length of each k-mer

        BOOST_LOG_TRIVIAL(info) << "Constructing the hash function ...";
        build_KRHash(kmers);        // build KR hash function
        build_minimalPerfectHash(); // build minimal perfect hash function

        //compute inverse of r modulo prime. The inverse is r^(p-2) mod Prime.  Proof:
        // r^(p-2) * r = r^(p-1) = 1, by Fermat's little theorem.
        uint256_t kr = modexp(r,Prime-2,Prime);

        if ((kr * r) % Prime == 1)
        {
            BOOST_LOG_TRIVIAL(info) << "r-inverse correctly computed modulo prime.";
        }
        else
        {
            BOOST_LOG_TRIVIAL(fatal) << "r-inverse incorrectly computed modulo prime!";
            exit(1);
        }

        rinv = kr;
    }

    void precomputePowers_mod()
    {
        uint256_t ri = 1;
        powersOfRModP.clear();
        //NEED 1 to k. Not 0 to (k - 1)
        for (unsigned i = 1; i <= k_kmer; ++i)
        {
            ri = (ri * r) % Prime;
            powersOfRModP.push_back(ri);
        }
    }


    /**
	 * Add a node to the hash function
	 * HACK: For now, we keep the new node in a map from kmer to its hash val
	 * Don't actually have a dynamic hash function
	 * Returns new hash value
	 */
    u_int64_t add_node(const kmer_t &new_kmer)
    {

        assert(this->new_nodes.find(new_kmer) == this->new_nodes.end());

        // Add it to the new nodes map with the next available hash value
        u_int64_t new_hash = this->max_hash + 1;
        this->max_hash++;

        this->new_nodes.insert(make_pair(new_kmer, new_hash));

        assert(this->new_nodes.find(new_kmer) != this->new_nodes.end());

        return new_hash;
    }

    //void remove_node(const kmer_t& node, const u_int64_t& hash) {
    void remove_node(const kmer_t &node)
    {

        auto it = this->new_nodes.find(node);

        if (it != this->new_nodes.end())
        {

            // this was a new node
            this->new_nodes.erase(it);
        }
        assert(this->new_nodes.find(node) == this->new_nodes.end());
    }

    /**
	 * Checks whether this kmer has been stored as one of the new nodes
	 * And returns the hash if it has
	 * o/w returns this->max_hash + 1
	 * which is never a valid hash
	 */
    u_int64_t new_node_hash(const kmer_t &kmer)
    {

        if (!this->new_nodes.empty())
        {

            map<kmer_t, u_int64_t>::iterator found_kmer = this->new_nodes.find(kmer);

            if (found_kmer != this->new_nodes.end())
            {
                return found_kmer->second;
            }
        }

        return this->max_hash + 1;
    }

    /**
	 * Returns whether the hash value is in the range of possible hash values
	 * Between 0 and this->max_hash
	 */
    bool hash_in_range(const u_int64_t &hash)
    {

        if ((hash < 0) || (hash > this->max_hash))
        {

            return false;
        }
        else
        {

            return true;
        }
    }

    /**
	 * Find the hash value of a k-mer
	 * Checks added and original nodes
	 * Will still return hash value of deleted nodes, that must be
	 * checked separately.
	 */
    u_int64_t get_hash_value(const kmer_t &seq)
    {
        u_int64_t res = this->new_node_hash(seq);

        if (res == this->max_hash + 1)
        {
            // not an added node
            kr_hash_t krv = generate_KRHash_val_mod(seq, this->k_kmer);
            res = (*recsplit)(KR_to_hash128_t(krv));

            //if (!this->isValidOriginalHash(res)) {
            // can't be valid
            //	res = -1;
            //}
        }

        return res;
    }

    /**
    * Find the hash value of a k-mer
    * Allows f( v ) notation
    */
    u_int64_t operator()(const kmer_t &seq)
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

        BOOST_LOG_TRIVIAL(info) << "Theoretical prime lower bound: "
                                << this->k_kmer * this->n_kmer_orig * this->n_kmer_orig;

        //Having problem with overflows because of large primes
        //Even though the theoretical bound is above, let's try smaller ones.
        //double smallerPrime = this->n_kmer_orig * this->n_kmer_orig;
        //Prime = getPrime((uint256_t)smallerPrime);
        Prime = 1;
        for(int i = 0; i < 128; i++) Prime *= 2;
        Prime -= 159; // 2^128 - 159 is a prime

        BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;

        unsigned n_failures = 0;
        // keep generating new base until we find one that is injective over our k-mers
        bool f_injective;

        do
        {
            ++n_failures;
            /*if (n_failures == 5)
            {

                BOOST_LOG_TRIVIAL(info) << "Trying a larger prime... ";
                Prime = getPrime(Prime * 2);
                n_failures = 0;
                BOOST_LOG_TRIVIAL(info) << "Trying prime: " << Prime;
            }*/

            f_injective = true; //assume f is injective until evidence otherwise
            this->r = randomNumber((uint256_t)1, Prime - 1);
            //Once we have a candidate base r
            //we should avoid recomputing its powers all the time

            precomputePowers_mod();

            for (unordered_set<kmer_t>::iterator
                     it1 = kmers.begin();
                 it1 != kmers.end();
                 ++it1)
            {
                kr_hash_t v1 = generate_KRHash_val_mod(*it1, k_kmer);
                if (KRHash.find(v1) == KRHash.end())
                {
                    // this is a new value
                    KRHash.insert(v1);
                }
                else // not injective
                {
                    BOOST_LOG_TRIVIAL(trace) << "Base " << this->r << " with prime "
                                             << Prime << " failed injectivity.";
                    KRHash.clear(); // clear it out and start over
                    f_injective = false;
                    break;
                }
            }
        } while (!f_injective);

        BOOST_LOG_TRIVIAL(info) << "Base " << this->r << " with prime " << Prime << " is injective.";

        return;
    }

    /**
    * Given a kmer, find out its KRH
    * MODULO Prime
    * Need to use 256 bits since 4 * x might overflow
	 * Note that this shouldn't be used for newly added nodes
	 * that are not part of the main hash function
    */
    kr_hash_t generate_KRHash_val_mod(const kmer_t &kmer,
                                      const unsigned &k)
    {

        uint256_t val = 0; // what will be the KRH value

        // go through each bp and add value
        for (unsigned i = 0; i < k; ++i)
        {
            // val += baseNum(kmer.at(i)) * pow(r, i);
            val = val + ((access_kmer(kmer, k, i) *
                          static_cast<uint256_t>(powersOfRModP[i])));
            val = val % Prime;
        }

        return kr_hash_t(val);
    }

    typedef number<cpp_int_backend<256, 256, unsigned_magnitude, checked, void>> moduloInt;
    kr_hash_t update_KRHash_val_OUT_mod(const kr_hash_t &KR_val, //KR hash of source kmer (mod P)
                                        const unsigned &first,   //character at front of source k-mer
                                        const unsigned &last)
    {
        uint256_t rinv = this->rinv;
        uint256_t Prime = this->Prime;
        uint256_t kr = KR_val;
        uint256_t llast = last;
        uint256_t ffirst = first;
        uint256_t rk = powersOfRModP[k_kmer - 1];
        uint256_t four = 4;

        uint256_t sub_val = four * Prime - ffirst * r;
        kr = (kr + sub_val) % Prime;
        kr = (kr * rinv)  % Prime;
        kr = (kr + llast * rk)  % Prime;

        return static_cast<kr_hash_t>(kr);
    }

    /*
    * This function takes as input a Karp-Rabin value (KR_val) modulo Prime
    *
    * target k-mer is IN neighbor of source k-mer
      */

    kr_hash_t update_KRHash_val_IN_mod(const kr_hash_t &KR_val, //KR hash of source kmer
                                       const unsigned &first,   //character at front of target k-mer
                                       const unsigned &last)
    {

        uint256_t r = this->r;
        uint256_t Prime = this->Prime;
        uint256_t kr = KR_val;
        uint256_t llast = last;
        uint256_t ffirst = first;
        uint256_t rk = powersOfRModP[k_kmer - 1];
        uint256_t four = 4;

        uint256_t sub_val = four * Prime - llast * rk;

        kr = (kr + sub_val) % Prime; // last * r^k
        kr = (kr * r) % Prime;
        kr = (kr + ffirst * r) % Prime;

        //	   while (kr > Prime)
        //	      kr = kr - Prime;

        //	   kr = kr - (kr / Prime)*Prime;

        return static_cast<kr_hash_t>(kr);
    }

    /*
    * Looks up the minimal perfect hash value, given the Karp-Rabin (value modulo Prime)
	 * If in added nodes, just returns the hash for that kmer. Doesn't use KR.
    */
   
    u_int64_t perfect_from_KR_mod(const kmer_t &kmer, const kr_hash_t &KR_val)
    {

        u_int64_t hash = this->new_node_hash(kmer);

        if (hash == this->max_hash + 1)
        {
            // not a new node
            hash = (*recsplit)(KR_to_hash128_t(KR_val));
        }

        return hash;
    }
    

    /**
    * Build a minimal perfect hash function on the set of integers that our kmers are
    * mapped to via KRH
    */
    void build_minimalPerfectHash()
    {

        BOOST_LOG_TRIVIAL(info) << "Building minimal perfect hash function ...";

        std::vector<sux::function::hash128_t> hash128_vec;
        for(kr_hash_t x : KRHash){
            hash128_vec.push_back(KR_to_hash128_t(x));
        }
        

        // MPHF for our KRHash function values
        this->recsplit = new RecSplit(hash128_vec, BUCKET_SIZE);

        BOOST_LOG_TRIVIAL(info) << "Minimal perfect hash function created.";
    }

    // Added to Jarno N. Alanko
    sux::function::hash128_t KR_to_hash128_t(kr_hash_t x){
        uint256_t mask = ((uint256_t)1 << 64) - 1;
        uint64_t low_bits = (uint64_t)(x & mask);
        uint64_t high_bits = (uint64_t)(x >> 64);
        return sux::function::hash128_t(high_bits, low_bits);
    }

};

#endif
