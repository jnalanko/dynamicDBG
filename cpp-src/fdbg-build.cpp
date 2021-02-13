/**

Copyright 2017 Alan Kuhnle, Victoria Crawford

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <chrono>
#include "FDBG.cpp"
#include "TestUtil.cpp"
#include "formatutil.cpp"

using namespace std;


/**
 * Run like "./a.out kmers.fasta 30 index.bin"
 */
int main(int argc, char* argv[]) {
   // Set debug level
   boost::log::core::get()->set_filter(boost::log::trivial::severity
      >= boost::log::trivial::info);

   BOOST_LOG_TRIVIAL(info) << "Beginning to build data structure ...";

   // Check if the user put in the correct command line arguments
   if (argc < 3) {
         BOOST_LOG_TRIVIAL(fatal) << "Missing required arguments. Usage:"
				  << argv[0] << " <reads.fasta> <k>";
         exit(1);
   }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);
   string dsfile = argv[3];
   FDBG Graph;

   // get k-mers and edgemers from file
   unordered_set<kmer_t> kmers;
   unordered_set<kmer_t> edgemers;
   auto start = std::chrono::system_clock::now();
   handle_mers( filename, k, kmers, edgemers );
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   BOOST_LOG_TRIVIAL(info) << "Getting " << kmers.size() + edgemers.size() << " mers took " << elapsed_seconds.count() << " s";

   BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
   Graph.build( kmers, edgemers, kmers.size(), k);

   BOOST_LOG_TRIVIAL(info) << "Data structure built in " << Graph.construction_time << " s";
   BOOST_LOG_TRIVIAL(info) << "Writing data structure to file " << dsfile;
   ofstream ofile( dsfile.c_str(), ios::out | ios::binary );
   Graph.save( ofile );
   ofile.close();

   return 0;
}

