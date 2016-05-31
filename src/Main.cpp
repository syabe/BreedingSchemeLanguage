//'@useDynLib BreedingSchemeLanguage
//'@importFrom Rcpp evalCpp

//////////////////////////////////////////////////////////////////////
// Main.cpp
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

#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <time.h>
#include "genome.h"
#include "Random.h"
#include <stdio.h>
#include <iomanip>
#include <fstream>

using namespace std;

RcppExport SEXP mainR(SEXP par1, SEXP par2, SEXP par3, SEXP par4, SEXP par5, SEXP par6, SEXP par7, SEXP par8, SEXP par9, SEXP par10)
{
  string popSize = Rcpp::as<std::string>(par3);
  string rec = Rcpp::as<std::string>(par6);
  int numPieces = Rcpp::as<int>(par5);
  int pieceLen = 10000;
  int numIndepRegion = Rcpp::as<int>(par4);
  int nSubPOP = Rcpp::as<int>(par1);
  vector<int> nSubSample = Rcpp::as< std::vector<int> >(par2);
  double mut=1e-8, mig=2.5e-4;
  double maf = Rcpp::as<double>(par7);
  double prop = 1.0;
  int SNP = Rcpp::as<int>(par10);

  long seed = Rcpp::as<long>(par8);
  int drawtree = Rcpp::as<int>(par9);

  vector< vector<bool> > data;

//  string popSize("10000"), rec("0.0001");

//  int numPieces=100, pieceLen=10000, numIndepRegion=1, SNP=-1;
//  int nSubPOP=2;
//  vector<int> nSubSample(nSubPOP);
//  for(int i=0;i<nSubPOP;i++)nSubSample[i]=10;

//  double mut=1e-8, mig=2.5e-4;
//  double maf=0.0,prop=1.0;

//  vector< vector<bool> > data;

//  long seed=-1;

//  int drawtree=0;

//  if(argc==1){
//    cout << endl  << "GENOME-0.2: Whole Genome Coalescent Simulator. (2009.2) Liming Liang, Goncalo Abecasis"  << endl << endl;
//    cout << "     Parameters and default values:" << endl;
//    cout << "     -pop     " << "number of subpopulations and size of subsamples [ " <<nSubPOP<<" ";
//    for(int i=0;i<nSubPOP;i++)cout<<nSubSample[i]<<" ";cout<<"]"<< endl;
//    cout << "     -N       " << "effective size of each subpopulation or the filename for population profile [" <<popSize.c_str()<<"]"<< endl;
//    cout << "     -c       " << "number of independent regions (chromosomes) to simulate [" <<numIndepRegion<<"]"<< endl;
//    cout << "     -pieces  " << "number of fragments per independent region [" <<numPieces<<"]"<< endl;
//    cout << "     -len     " << "length in base of each fragment [" <<pieceLen<<"]"<< endl;
//    cout << "     -s       " << "fixed number of SNPs per independent region (chromosome), -1 = number of SNPs follows Poisson distribution [" <<SNP<<"]"<< endl;
//    cout << "     -rec     " << "recombination rate bewteen consecutive fragments"<<" or the filename for recombination rate distribution ["<<rec.c_str()<<"]"<< endl;
//    cout << "     -mut     " << "mutation rate per generation per base pair [" <<mut<<"]"<< endl;
//    cout << "     -mig     " << "migration rate per generation per individual [" <<mig<<"]"<< endl;
//    cout << "     -seed    " << "random seed, -1 = use time as the seed [" <<seed<<"]"<< endl;
//    cout << "     -tree    " << "1=draw the genealogy trees, 0=do not output [" <<drawtree<<"]"<< endl;
//    cout << "     -maf     " << "Output SNPs with minor allele frequency greater than [" <<maf<<"]"<< endl;
//    cout << "     -prop    " << "To keep this proportion of SNPs with MAF < the value of -maf parameter [" <<prop<<"]"<< endl;
//    exit(0);
//  }


//  int k=0;
//  for(int i=1;i<argc;i++){
//    if(strcmp(argv[i],"-pop")==0){
//      nSubPOP=atoi(argv[++i]); k++;
//      nSubSample.clear();
//      nSubSample.resize(nSubPOP);
//      for(int iter=0;iter<nSubPOP;iter++)nSubSample[iter]=atoi(argv[++i]);
//    }
//    else if(strcmp(argv[i],"-N")==0){ popSize.assign(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-c")==0){ numIndepRegion=atoi(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-pieces")==0){ numPieces=atoi(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-len")==0){ pieceLen=atoi(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-s")==0){ SNP=atoi(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-rec")==0){ rec.assign(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-mut")==0){ mut=atof(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-mig")==0){ mig=atof(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-seed")==0){ seed=atol(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-tree")==0){ drawtree=atoi(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-maf")==0){ maf=atof(argv[++i]); k++;}
//    else if(strcmp(argv[i],"-prop")==0){ prop=atof(argv[++i]); k++;}
//  }

//  if(k!=(argc-1)/2 && k!=(argc-1-nSubPOP)/2){
//    cout << "Parameters do not match! " << argc<< " " << k << endl;
//    exit(1);
//  }

  Rcpp::Rcout << endl
              << "GENOME-0.2: Whole Genome Coalescent Simulator. (2009.2) Liming Liang, Goncalo Abecasis"  << endl << endl;
  Rcpp::Rcout << "     Parameters and effective values:" << endl;
  Rcpp::Rcout << "     -pop     " << "number of subpopulations and size of subsamples [ " <<nSubPOP<<" ";
  for(int i=0;i<nSubPOP;i++)Rcpp::Rcout<<nSubSample[i]<<" ";Rcpp::Rcout<<"]"<< endl;
  Rcpp::Rcout << "     -N       " << "effective size of each subpopulation or the filename for population profile [" <<popSize.c_str()<<"]"<< endl;
  Rcpp::Rcout << "     -c       " << "number of independent regions (chromosomes) to simulate [" <<numIndepRegion<<"]"<< endl;
  Rcpp::Rcout << "     -pieces  " << "number of fragments per independent region [" <<numPieces<<"]"<< endl;
  Rcpp::Rcout << "     -len     " << "length in base of each fragment [" <<pieceLen<<"]"<< endl;
  Rcpp::Rcout << "     -s       " << "fixed number of SNPs per independent region (chromosome), -1 = number of SNPs follows Poisson distribution [" <<SNP<<"]"<< endl;
  Rcpp::Rcout << "     -rec     " << "recombination rate bewteen consecutive fragments or the filename for recombination rate distribution ["<<rec.c_str()<<"]"<< endl;
  Rcpp::Rcout << "     -mut     " << "mutation rate per generation per base pair [" <<mut<<"]"<< endl;
  Rcpp::Rcout << "     -mig     " << "migration rate per generation per individual [" <<mig<<"]"<< endl;
  Rcpp::Rcout << "     -seed    " << "random seed, -1 = use time as the seed [" <<seed<<"]"<< endl;
  Rcpp::Rcout << "     -tree    " << "1=draw the genealogy trees, 0=do not output [" <<drawtree<<"]"<< endl;
  Rcpp::Rcout << "     -maf     " << "Output SNPs with minor allele frequency greater than [" <<maf<<"]"<< endl;
  Rcpp::Rcout << "     -prop    " << "To keep this proportion of SNPs with MAF < the value of -maf parameter [" <<prop<<"]"<< endl<<endl;

  if(atoi(popSize.c_str())>0 && atof(rec.c_str())>0){
    Rprintf("     Scaled recombination rate = %.4e\n", 2 * atoi(popSize.c_str()) * (numPieces - 1) * atof(rec.c_str()));
    Rprintf("     Scaled migration rate = %.4e\n", 2 * atoi(popSize.c_str()) * mig);
    Rprintf("     Scaled mutation rate = %.4e\n\n", 2 * atoi(popSize.c_str()) * mut * numPieces * pieceLen);
  }

  genome(popSize, nSubPOP, nSubSample, numPieces, pieceLen, numIndepRegion, SNP, rec, mut, mig, data, seed,drawtree);

  // output format for haploview

  // 	for(int i=0;i<data.size();i++){
    // 		printf("FAMILY1 MEMBER%d ",i/2+1);
    // 		for(int j=0;j<data[i].size();j++){
      // 			cout<<data[i][j]+1<<" ";
      // 		}
    // 		cout<<endl;
    // 	}


  // 0 = ancestral state, 1 = mutated state

  Rcpp::Rcout<<"Total number of SNPs = "<<data[0].size()<<endl;
  Rcpp::Rcout<<"Samples:"<<endl;

  if(maf>0.0 && prop < 1.0){ // filter out the SNPs with MAF < maf

     Rcpp::Rcout<<"Rare SNPs filtered according to parameters."<<endl;

     int n=0;
     for(int k=0;k<nSubPOP;k++)n+=nSubSample[k];

     int m=(int)(data[0].size());
     vector<bool> keep(m,true);

     float p;
     for(int i=0;i<m;i++){
       p=0.0;
       for(int j=0;j<n;j++){
         p += data[j][i];
       }
       p /= n;

       if( (p<maf || (1-p)<maf) && globalRandom.Next()>prop)keep[i]=0;
     }


     int index=0;
     for(int k=0;k<nSubPOP;k++){
       for(int i=0;i<nSubSample[k];i++){
         Rcpp::Rcout<<"POP"<<k+1<<": ";
         for(unsigned int j=0;j<data[index].size();j++){
           if(keep[j])Rcpp::Rcout<<data[index][j];
         }
         Rcpp::Rcout<<endl;
         index++;
       }
     }

  }
  else{
    int index=0;
    for(int k=0;k<nSubPOP;k++){
      for(int i=0;i<nSubSample[k];i++){
        Rcpp::Rcout<<"POP"<<k+1<<": ";
        for(unsigned int j=0;j<data[index].size();j++){
          Rcpp::Rcout<<data[index][j];
        }
        Rcpp::Rcout<<endl;
        index++;
      }
    }
  }
  return R_NilValue;
}
