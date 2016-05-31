////////////////////////////////////////////////////////////////////// 
// genome.h 
//////////////////////////////////////////////////////////////////////////////
//              COPYRIGHT NOTICE FOR GENOME CODE
//
// Copyright (C) 2006 - 2009, Liming Liang and Goncalo Abecasis,
// All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <string>

using namespace std;

void genome(string popSize, 						// effective population size of each population
		  int nSubPOP, 							// number of populations
		  vector<int> & nSubSample, 			// number of samples (haplotypes) draw from each populations
		  int numPieces, 						// number of fragments for each sample (chromosome)
		  int pieceLen, 						// length in base pair of each fragment
		  int numIndepRegion, 					// number of independent regions (independent chromosome)
		  int s,								// fixed number of SNPs want to simulate, randomly place s SNPs on the genealogy 
		  string rec,							// recombination rate between consecutive fragments per generation
		  double mut, 							// mutation rate per generation per base pair
		  double mig, 							// migration rate per generation
		  vector< vector<bool> > &return_chromosome, // the return chromosome by connecting all independent chromosome into one long chromosome
		  long t,								// random seed 
		  int drawtree,						// drawtree=1 output the genealogy trees for each fragment and chromosome; drawtree=0 do not output the trees	
		  bool printparameters =false 	 		// =1 print out all input parameters, =0 do not print input parameters
		  );



 
