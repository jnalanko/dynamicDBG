#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <iostream>
#include <string>
#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <chrono>

#include "FDBG.cpp"
#include "TestUtil.cpp"
#include "formatutil.cpp"
#include "input_reading.hh"
#include "throwing_streams.hh"
#include <chrono>

long long cur_time_millis() {
    return (std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()))
        .count();
}


int main(int argc, char* argv[]) {

    string graphfile = argv[1];
    string queryfile = argv[2];
    string outfile = argv[3];
    

    FDBG Graph;
    ifstream ifile_ds(graphfile, ios::in | ios::binary);
    Graph.load(ifile_ds);
    long long nodemer_k = Graph.k; // k is node length in fdbg
    long long edgemer_k = nodemer_k + 1;

    BOOST_LOG_TRIVIAL(info) << "Data structure loaded ";
    BOOST_LOG_TRIVIAL(info) << "k value or node size is: " << Graph.k
                            << " graph has " << Graph.n << " nodes ";
    BOOST_LOG_TRIVIAL(info) << "Original Bits per element:"
                            << Graph.bitSize() / static_cast<double>(Graph.n);

    // Query edges

    long long query_start_jarno = cur_time_millis();
    throwing_ofstream ofs(outfile);
    Sequence_Reader sr(queryfile, FASTA_MODE);
    while (!sr.done()) {
        string seq = sr.get_next_query_stream().get_all();
        vector<bool> hits(std::max(0LL, (long long)seq.size() - edgemer_k + 1));
        for (long long i = 0; i < (long long)seq.size() - edgemer_k + 1; i++) {
            kmer_t prefix_bin = mer_string_to_binary(seq, 0, nodemer_k);
            kmer_t suffix_bin = mer_string_to_binary(seq, 1, nodemer_k);
            bool present = Graph.IsEdgeInGraph(prefix_bin, suffix_bin);
            hits[i] = present;
        }
        // Write out
        for(long long i = 0; i < (long long)hits.size(); i++)
            ofs << hits[i];
        ofs << "\n";
    }
    long long elapsed_jarno = cur_time_millis() - query_start_jarno;
    std::cerr << "Time for all queries: " << (double)elapsed_jarno / 1e3
                << " seconds" << std::endl;
    
    //for (size_t i = 0; i < nodes.size(); i++) {
        //if (!Graph.IsEdgeInGraph(nodes[i].first, nodes[i].second))

 
}
