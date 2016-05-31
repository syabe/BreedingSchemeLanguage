//////////////////////////////////////////////////////////////////////
// genome.cpp
//////////////////////////////////////////////////////////////////////////////
//              COPYRIGHT NOTICE FOR GENOME CODE
//
// Copyright (C) 2006 - 2009, Liming Liang and Goncalo Abecasis,
// Version 0.2
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


#define _CRT_SECURE_NO_DEPRECATE

#include <Rcpp.h>
#include "Random.h"
#include "Error.h"
#include "stochastic.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <algorithm>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string.h>

using namespace std;

//static bool debug = false;
static bool fixedPOP = true;
static bool noChange = true;


static int N = 10000, n = 100, L = 1, poolsize = 1024 * 1024 * 64, length=10000, numChr=1, SNP=-1; // length = number of base per fragment
static int maxN;
static unsigned int currentProfile;
static int nPOP=1;
static vector<int> nSample(nPOP,3);
static double recombination = 1e-8, migration = 0.0, mutation=1e-8;

static int * genepool, * genetimes, * fragments, * MRCA, * MRCAtimes;
static int *** children, *** parents, * cblocks, * pblocks;
static int * parent_First, * parent_Last, * child_First, * child_Last;
static bool * child_flags, * parent_flags;
static int nextFree = 0, activeSegments = 0, pblockcounter = 0;
static int numBLOCKs;
static int BUFFER=1500;
static int INCREMENT=1000;

#define BLOCKSIZE   128
#define SWAP(tmp,a,b)  tmp=a; a=b; b=tmp;

static vector< vector< vector<bool> > > chromosome;
static int * fragmentLength;
static int * branchFragment;
static int * geneIndex;
static int * branchMutation;

static float * recombVector=NULL;

struct popProfilesStruct{
	int generation;
	vector<int> popsize;
	vector<int> popStart;
	vector< vector<float> > selectProb;
	vector< vector<int> > popChange;
};

static vector<struct popProfilesStruct> popProfile;

void print(){
  Rcpp::Rcout<<"\n*****************\nParents:\t";
	for(int i=0;i<N;i++){
		if(parent_flags[i]){
			for(int j=0;j<L;j++) Rcpp::Rcout<<"["<<((parents[i][j / BLOCKSIZE]==NULL)?999:(parents[i][j / BLOCKSIZE][j % BLOCKSIZE]))<<"]";
		  Rcpp::Rcout<<"\t";
		}
	}
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"Children:\t";
	for(int i=0;i<N;i++){
		if(child_flags[i]){
			for(int j=0;j<L;j++) Rcpp::Rcout<<"["<<((children[i][j / BLOCKSIZE]==NULL)?999:(children[i][j / BLOCKSIZE][j % BLOCKSIZE]))<<"]";
		  Rcpp::Rcout<<"\t";
		}
	}
	Rcpp::Rcout<<endl;



// 	cout<<"parent_flag:\t";
// 	for(int i=0;i<N;i++){
// 		cout<<parent_flags[i]<<"\t";
// 	}
// 	cout<<endl;
//
// 	cout<<"children_flag:\t";
// 	for(int i=0;i<N;i++){
// 		cout<<child_flags[i]<<"\t";
// 	}
// 	cout<<endl;

	Rcpp::Rcout<<"genepool:\t";
	for(int i=0;i<2*n*L-L;i++)Rcpp::Rcout<<genepool[i]<<" ";
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"genetimes:\t";
	for(int i=0;i<2*n*L-L;i++)Rcpp::Rcout<<genetimes[i]<<" ";
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"MCRA:\t";
	for(int i=0;i<L;i++)Rcpp::Rcout<<MRCA[i]<<" ";
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"MCRAtime:\t";
	for(int i=0;i<L;i++)Rcpp::Rcout<<MRCAtimes[i]<<" ";
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"fragments:\t";
	for(int i=0;i<L;i++)Rcpp::Rcout<<fragments[i]<<" ";
	Rcpp::Rcout<<endl;

	Rcpp::Rcout<<"NextFree="<<nextFree<<endl;
}

int *** AllocateIntPtrMatrix(int rows, int columns){
   int *** result = new int ** [rows];

   for (int i = 0; i < rows; i++)
      result[i] = new int * [columns];

   return result;
}

void FreeIntPtrMatrix(int *** & matrix, int rows){

   for (int i = 0; i < rows; i++)
      delete [] matrix[i];

   delete [] matrix;

   matrix = NULL;
}

void NewGeneration(){

   int *** ptrmatrix, * vector;
   bool * boolvector;

//    for(int i=0;i<maxN;i++){
// 	if(parent_flags[i])printf("numBlock=%d parent=%d first=%d last=%d\n",numBLOCKs,i,parent_First[i],parent_Last[i]);
//    }

   SWAP(ptrmatrix, children, parents);
   SWAP(vector, child_First, parent_First);
   SWAP(vector, child_Last, parent_Last);
   SWAP(boolvector, child_flags, parent_flags);
   SWAP(vector, cblocks, pblocks);

   if(!noChange){
	   currentProfile++;
	   nPOP = (int)(popProfile[currentProfile].popsize.size());
	   N = popProfile[currentProfile].popStart.back()+popProfile[currentProfile].popsize.back();
	   noChange = true;
   }

   for (int i = 0; i < maxN; i++)
      parent_flags[i] = 0;

//    for (int i = 0; i < maxN; i++)
//        for (int j = 0; j < numBLOCKs; j++)
//            parents[i][j] = NULL;

   pblockcounter = 0;

}


void AllocateMemory(vector< vector<bool> > &return_chromosome){


//    if(SNP<=0){
//  	   	cout<<"Reserve size for return_chromosome="<<max(1,(int)(mutation*(double)numChr*(double)L*(double)length*40*(double)N/(double)nPOP))<<endl;
//    }
//    else{
//   		cout<<"Reserve size for return_chromosome="<<SNP*numChr<<endl;
//    }
//    cout<<"Reserve size for working chromosome="<<max(1,(int)(mutation*(double)length*40*(double)N/(double)nPOP))<<endl;


   return_chromosome.clear();
   return_chromosome.resize(n);
//    for(int i=0;i<n;i++){
// 	   if(SNP<=0){
//  	   		//return_chromosome[i].reserve(max(1,40*(int)(mutation*(double)numChr*(double)L*(double)length*(double)maxN/(double)nPOP)));
// 	   }
//  	   else{
//  	   		//return_chromosome[i].reserve(SNP*numChr);
// 	   }
// 	}

   chromosome.resize(n);
   for(int i=0;i<n;i++){
 		chromosome[i].resize(L);
// 		for(int j=0;j<L;j++){
//  		   if(SNP<=0){
// 	 		   //chromosome[i][j].reserve(max(1,40*(int)(mutation*(double)length*(double)maxN/(double)nPOP)));
// 		   }
// 	 	   else{
// 	 	   	   //chromosome[i][j].reserve(SNP/L+1);
// 		   }
// 		}
   }

   genepool = new int [poolsize];

   numBLOCKs= (int)ceil((double)L / (double)BLOCKSIZE);

   pblockcounter = 0;

}


void AllocateMutationMemory(){

	fragmentLength = new int [L];

	branchFragment = new int [poolsize];

	geneIndex = new int [poolsize];

}



void FreeMemory(){

   delete [] genepool;

   if(recombVector!=NULL)delete [] recombVector;

}


void FreeMutationMemory(){

   delete [] fragmentLength;
   delete [] geneIndex;
   delete [] branchMutation;


}

void FreeCoalescentMemory(){

   delete [] MRCA;
   delete [] MRCAtimes;

   delete [] fragments;

   FreeIntPtrMatrix(children,maxN);

   FreeIntPtrMatrix(parents,maxN);

   delete [] child_flags;
   delete [] parent_flags;

   delete [] cblocks;
   delete [] pblocks;

   delete [] parent_First;
   delete [] parent_Last;

   delete [] child_First;
   delete [] child_Last;

}

void reallocBlock(){

	pblocks = (int *)realloc(pblocks, (unsigned)(sizeof(int)*(n*(numBLOCKs)*BLOCKSIZE+BUFFER*BLOCKSIZE+INCREMENT*BLOCKSIZE)) );
	cblocks = (int *)realloc(cblocks, (unsigned)(sizeof(int)*(n*(numBLOCKs)*BLOCKSIZE+BUFFER*BLOCKSIZE+INCREMENT*BLOCKSIZE)) );
	BUFFER += INCREMENT;

}

void setnull(int parent, int position){

  if (parent_flags[parent] == 0){
		parent_First[parent]=position / BLOCKSIZE;
		parent_Last[parent]=position / BLOCKSIZE;
		parents[parent][position / BLOCKSIZE] = NULL;
  }
  else if(parent_First[parent]>position / BLOCKSIZE){
	  for(int i=position / BLOCKSIZE;i<parent_First[parent];i++){
		  parents[parent][i] = NULL;
	  }
	  parent_First[parent]=position / BLOCKSIZE;
  }
  else if(parent_Last[parent]<position / BLOCKSIZE){
	  for(int i=position / BLOCKSIZE;i>parent_Last[parent];i--){
		  parents[parent][i] = NULL;
	  }
	  parent_Last[parent]=position / BLOCKSIZE;
  }

}



void TouchParent(int parent, int position){


  setnull(parent,position);

  if (parents[parent][position / BLOCKSIZE] != NULL)return;


  if (pblockcounter >= n*numBLOCKs+BUFFER){

	    int * pblocks_old=pblocks;
	    int * cblocks_old=cblocks;

	    reallocBlock();

	    for(int i=0;i<N;i++)
	      if(parent_flags[i])
	    	for(int j=parent_First[i];j<=parent_Last[i];j++)
	    		if(parents[i][j]!=NULL)parents[i][j] += pblocks - pblocks_old;

	    for(int i=0;i<N;i++)
	      if(child_flags[i])
	    	for(int j=child_First[i];j<=child_Last[i];j++)
	    		if(children[i][j]!=NULL)children[i][j] += cblocks - cblocks_old;

 		Rprintf("pblock and cblock enlarged, BUFFER=%d.!\n",BUFFER);

		if(pblocks==NULL || cblocks==NULL)error("Blocks memory reallocation error.\n");
  }

  parents[parent][position / BLOCKSIZE] = pblocks + BLOCKSIZE * pblockcounter++;

  for (int i = 0; i < BLOCKSIZE; i++)
      parents[parent] [position / BLOCKSIZE] [i] = -1;

  parent_flags[parent] = 1;

}

void Initialize(){

   genetimes = new int [poolsize];

   MRCA = new int [L];
   MRCAtimes = new int [L];

   fragments = new int [L];

   children = AllocateIntPtrMatrix(maxN, numBLOCKs);
   parents = AllocateIntPtrMatrix(maxN, numBLOCKs);

   parent_First = new int [maxN];
   parent_Last = new int [maxN];

   child_First = new int [maxN];
   child_Last = new int [maxN];

   cblocks = (int *)malloc((unsigned)(sizeof(int)*(n*(numBLOCKs)*BLOCKSIZE+BUFFER*BLOCKSIZE)));
   pblocks = (int *)malloc((unsigned)(sizeof(int)*(n*(numBLOCKs)*BLOCKSIZE+BUFFER*BLOCKSIZE)));

   child_flags = new bool [maxN];
   parent_flags = new bool [maxN];


   noChange = true;

   nPOP = (int)(popProfile[0].popsize.size());

   N = popProfile[0].popStart.back()+popProfile[0].popsize.back();

   currentProfile = 0;

   nextFree = 0;

   for (int i = 0; i < N; i++)
      child_flags[i] = 0;

   int count=0;
   for(int k=0;k<nPOP;k++){
	   for (int i = 0; i < nSample[k]; i++)
	      {
		  child_First[popProfile[0].popStart[k]+i]=0;
		  child_Last[popProfile[0].popStart[k]+i]=numBLOCKs-1;

	      child_flags[popProfile[0].popStart[k]+i] = 1;

	      for (int j = 0; j < numBLOCKs; j++)
         		children[popProfile[0].popStart[k]+i][j] = cblocks + (count * L + j * BLOCKSIZE);

	      for (int j = 0; j < L; j++)
	         {
	         genepool[nextFree] = 1;
	         genetimes[nextFree] = 0;
	         children[popProfile[0].popStart[k]+i][j / BLOCKSIZE][j % BLOCKSIZE] = nextFree++;
	         }

	         count++;
	      }
   }

   for (int i = 0; i < L; i++)
      fragments[i] = n;

   for (int i = 0; i < N; i++)
      parent_flags[i] = 0;

   activeSegments = L;

   //for (int i = 0; i < N; i++)
   //  for (int j = 0; j < numBLOCKs; j++)
   //     parents[i][j] = NULL;

   pblockcounter = 0;

}

void InitializeMutation(vector< vector<bool> > &return_chromosome){

// 	   return_chromosome.clear();
//    return_chromosome.resize(n);
//    for(int i=0;i<n;i++){
//  	   return_chromosome[i].reserve(max(1,(int)(mutation*(double)numChr*(double)L*(double)length*40*(double)N/(double)nPOP)));
// 	}
//
//    chromosome.resize(n);
//    for(int i=0;i<n;i++){
//  		chromosome[i].resize(L);
// 		for(int j=0;j<L;j++){
//  		   chromosome[i][j].reserve(max(1,(int)(mutation*(double)length*40*(double)N/(double)nPOP)));
// 		}
//    }

//    cout<<"WGCS-- chromosomes capacity="<<chromosome[0][0].capacity()<<" size="<<chromosome[0][0].size()<<endl;

   for(int i=0;i<n;i++){
	   for(int j=0;j<L;j++){
		   chromosome[i][j].clear();
	   }
   }

//    cout<<"WGCS-- after clear() chromosomes capacity="<<chromosome[0][0].capacity()<<" size="<<chromosome[0][0].size()<<endl;

   branchMutation = new int [nextFree];
   for(int i=0;i<nextFree;i++)branchMutation[i]=0;

}

void PrepareMutation(){

   for(int i=0;i<L;i++)fragmentLength[i]=0;

   for(int i=0;i<nextFree;i++){
	   geneIndex[i] = -1;
	   branchFragment[i] = -1;
   }

}

int SampleParent(int individual, int previous){

	if(fixedPOP){
	   int population = individual/popProfile[0].popsize[0];

	   if (previous == -1 && nPOP>1){
	      if (globalRandom.Next() < migration){
		      int random = globalRandom.NextInt()%(nPOP-1);
	          if(population<=random)population=random+1;
	          else population=random;
	      }
	   }

	   return globalRandom.NextInt() % (popProfile[0].popsize[0]) + popProfile[0].popStart[population];
   }
   else{

	   if(noChange){
		   int population=0;

		   for(int i=0;i<nPOP;i++){
			   if(individual>=popProfile[currentProfile].popStart[i])population=i;
			   else break;
		   }

//		   cout<<"No change: ind="<<individual<<" pop="<<population<<endl;

		   if (previous == -1 && nPOP>1){
		      if (globalRandom.Next() < migration){
			      int random = globalRandom.NextInt()%(nPOP-1);
		          if(population<=random)population=random+1;
		          else population=random;
		      }
		   }

//		   cout<<"new pop="<<population<<endl;

		   return globalRandom.NextInt() % (popProfile[currentProfile].popsize[population]) + popProfile[currentProfile].popStart[population];
	   }
	   else{
		   int population=0;

		   for(int i=0;i<nPOP;i++){
			   if(individual>=popProfile[currentProfile].popStart[i])population=i;
			   else break;
		   }
//		   cout<<"Change: ind="<<individual<<" pop="<<population<<endl;

		   if (previous == -1 && nPOP>1){
		      if (globalRandom.Next() < migration){
			      int random = globalRandom.NextInt()%(nPOP-1);
		          if(population<=random)population=random+1;
		          else population=random;
		      }
		   }

//		   cout<<"new pop="<<population<<endl;

		   if(popProfile[currentProfile].popChange[population].size()==1){
//			   cout<<"To next pop="<<popProfile[currentProfile].popChange[population][0]<<endl;

				return globalRandom.NextInt() % (popProfile[currentProfile+1].popsize[popProfile[currentProfile].popChange[population][0]])
						+ popProfile[currentProfile+1].popStart[popProfile[currentProfile].popChange[population][0]];
		   }
		   else{
		   		double selectRandom = globalRandom.Next();
		   		for(unsigned int i=0;i<popProfile[currentProfile].popChange[population].size();i++){
			   		if(selectRandom < popProfile[currentProfile].selectProb[population][i]){

//				   		cout<<"To next pop="<<popProfile[currentProfile].popChange[population][i]<<endl;

				   		return globalRandom.NextInt() % (popProfile[currentProfile+1].popsize[popProfile[currentProfile].popChange[population][i]])
						+ popProfile[currentProfile+1].popStart[popProfile[currentProfile].popChange[population][i]];
			   		}
				}
		   	Rcpp::Rcout<<"Random number generator error!"<<endl;
		   	Rcpp::stop("genome error! \n");
				//exit(1);
				return 0;
	   	   }


	   }

   }
}

void readPopProfile(string filename, unsigned int numPopulation){

	   int popProfile_index=0;

	   popProfile.clear();
	   struct popProfilesStruct temp_popProfile;
	   popProfile.push_back(temp_popProfile);

	   ifstream popFile;
	   string fileline;
	   char *term;
	   char *term1;

	   popFile.open(filename.c_str());
	   if(!popFile)error("Population profile: %s cannot be opened!",filename.c_str());

	   getline(popFile,fileline);					// the first line is the population sizes at generation 0

	   term = strtok ((char*)fileline.c_str()," \t");  // the first term is the generation
	   popProfile[popProfile_index].generation=atoi(term);

	   if(popProfile[popProfile_index].generation!=0)error("The first line in population profile should be generation 0!");

	   term = strtok (NULL, " \t");					   // the second term and after are the sizes of that population
	   while(term!=NULL){
			   popProfile[popProfile_index].popsize.push_back(atoi(term));
			   term = strtok (NULL, " \t");
	   }

	   popProfile[popProfile_index].popStart.push_back(0); // the first population starts from 0
	   for(unsigned int i=0;i<popProfile[popProfile_index].popsize.size()-1;i++){
	   	   popProfile[popProfile_index].popStart.push_back(popProfile[popProfile_index].popsize[i]+popProfile[popProfile_index].popStart.back());
	   }

	   maxN = popProfile[popProfile_index].popStart.back()+popProfile[popProfile_index].popsize.back();

	   if(popProfile[0].popsize.size()!=numPopulation)error("The number of populations specified by -pop is not consistent with generation 0 in population profile:%s",filename.c_str());

	   // read in the population correspondence to the next population profile


	   while(popFile.peek()!=EOF){
		   getline(popFile,fileline);						// the even number lines are the correspondence of populations of different generation
		   popProfile[popProfile_index].popChange.resize(popProfile[popProfile_index].popsize.size());

		   term = strtok ((char*)fileline.c_str(),"- \t");  // the FROM population
		   while(term!=NULL){

			   term1 = strtok (NULL, "- \t");				// the TO population

			   popProfile[popProfile_index].popChange[atoi(term)-1].push_back(atoi(term1)-1);

			   term = strtok (NULL, "- \t");
		   }

		   popProfile_index++;
		   popProfile.push_back(temp_popProfile);

		   if(popFile.peek()==EOF)error("Error in population profile: the last line should specify the population sizes.");

		   // read the next population profile

	   	   getline(popFile,fileline);					// the population sizes

		   term = strtok ((char*)fileline.c_str()," \t");  // the first term is the generation
		   popProfile[popProfile_index].generation=atoi(term);

		   term = strtok (NULL, " \t");					   // the second term and after are the sizes of that population
		   while(term!=NULL){
				   popProfile[popProfile_index].popsize.push_back(atoi(term));
				   term = strtok (NULL, " \t");
		   }

		   popProfile[popProfile_index].popStart.push_back(0); // the first population starts from 0
		   for(unsigned int i=0;i<popProfile[popProfile_index].popsize.size()-1;i++){
		   	   popProfile[popProfile_index].popStart.push_back(popProfile[popProfile_index].popsize[i]+popProfile[popProfile_index].popStart.back());
		   }

		   if(maxN < popProfile[popProfile_index].popStart.back()+popProfile[popProfile_index].popsize.back())
		   			maxN = popProfile[popProfile_index].popStart.back()+popProfile[popProfile_index].popsize.back();
	   }

	   for(unsigned int i=0;i<popProfile.size();i++){
		   popProfile[i].selectProb.resize(popProfile[i].popChange.size());
		   for(unsigned int j=0;j<popProfile[i].popChange.size();j++){
			   float sum=0.0;
			   for(unsigned int k=0;k<popProfile[i].popChange[j].size();k++)sum += popProfile[i+1].popsize[popProfile[i].popChange[j][k]];

			   popProfile[i].selectProb[j].push_back(((float)popProfile[i+1].popsize[popProfile[i].popChange[j][0]])/sum);
			   for(unsigned int k=1;k<popProfile[i].popChange[j].size();k++)
			   		popProfile[i].selectProb[j].push_back(popProfile[i].selectProb[j].back()+popProfile[i+1].popsize[popProfile[i].popChange[j][k]]/sum);
		   }
	   }
	   // check for consistence

	   int min, max;

	   for(unsigned int i=0;i<popProfile.size()-1;i++){
			min=99999;max=-1;
			for(unsigned int j=0;j<popProfile[i].popChange.size();j++){
				if(popProfile[i].popChange[j].size()==0)error("One population does not have its parent population: profile file line %d, population %d\n",i+1,j+1);
				for(unsigned int k=0;k<popProfile[i].popChange[j].size();k++){
					if(min>popProfile[i].popChange[j][k])min=popProfile[i].popChange[j][k];
					if(max<popProfile[i].popChange[j][k])max=popProfile[i].popChange[j][k];
				}
			}
			if(min!=0 || max!=(int)(popProfile[i+1].popsize.size()-1))
				error("The population changes does not consistent with the next population profile. line=%d, min=%d, max=%d, num POP in next line=%d\n",i+1,min+1,max+1,popProfile[i+1].popsize.size());
	   }

}


void readRecombination(string filename, int pieces){

	   ifstream recombFile;
	   string fileline;
	   char *term;

	   vector<float> recDistribution;
	   vector<float> recRate;

	   float sum=0.0;
	   float random;

	   recombFile.open(filename.c_str());
	   if(!recombFile)error("Recombination profile: %s cannot be opened!",filename.c_str());

	   getline(recombFile,fileline);		// the first line is head of distribution table
	   term = strtok ((char*)fileline.c_str()," \t");
	   term = strtok (NULL, " \t");

	   if(strcmp(term,"Freq")==0 || strcmp(term,"freq")==0){
		   while(recombFile.peek()!=EOF){
			   getline(recombFile,fileline);

			   term = strtok ((char*)fileline.c_str()," \t"); // the first number is the recombination rate
			   recRate.push_back((float)(atof(term)));

			   if(atof(term)<0 || atof(term)>1)error("Recombination rate should be between 0 and 1. input=%f",atof(term));

			   term = strtok (NULL, " \t");					  // the second number is the frequency of recombination rate
			   recDistribution.push_back((float)(atof(term)));

			   if(atof(term)<=0)error("Frequency should be positive. input=%f",atof(term));

			   sum += (float)(atof(term));
		   }

		   if(recRate.size()!=recDistribution.size() || recRate.size()==0)error("Errors in the recombination profile: %s",filename.c_str());

		   // normalize the frequency
		   for(unsigned int i=0;i<recDistribution.size();i++) recDistribution[i] /= sum;

		   // compute the cumulative distribution
		   for(unsigned int i=1;i<recDistribution.size();i++) recDistribution[i] += recDistribution[i-1];

		   recombVector = new float [pieces-1];

		   for(int i=0;i<pieces-1;i++){
			   random = (float)(globalRandom.Next());
			   for(unsigned int j=0;j<recDistribution.size();j++){
					if(random < recDistribution[j]){
						recombVector[i]=recRate[j];
						break;
					}
			   }
		   }

			//check the input
			Rcpp::Rcout<<"*** RECOMBINATION RATES PROFILE ***"<<endl;
		  Rcpp::Rcout<<"Recombination_Rate\tCumulative_frequency"<<endl;
			for(unsigned int i=0;i<recRate.size();i++){
			  Rcpp::Rcout<<recRate[i]<<"\t\t\t"<<recDistribution[i]<<endl;
			}

   	   }
   	   else if(strcmp(term,"Pos")==0 || strcmp(term,"pos")==0){
	   	   float currentRec=0.0,nextRec=0.0;
	   	   int nextPos=-1;

	   	   if(recombFile.peek()!=EOF){
		   	   getline(recombFile,fileline);
		   	   term = strtok ((char*)fileline.c_str()," \t");
	   		   currentRec=(float)(atof(term));

	   		   if(recombFile.peek()!=EOF){
		   		   getline(recombFile,fileline);
		   	   	   term = strtok ((char*)fileline.c_str()," \t");
	   		       nextRec=(float)(atof(term));
	   		       term = strtok (NULL, " \t");
	   		       nextPos=atoi(term);
	   		   }
	   		   else{
		   		   nextPos=-1;
	   		   }

   		   }
   		   else{
	   		   error("Recombination file does not contain recombination rate information.\n");
   		   }

		   recombVector = new float [pieces-1];

   		   for(int i=0;i<pieces-1;i++){
	   		   if(i+1==nextPos){
		   		    currentRec=nextRec;
		   		    if(recombFile.peek()!=EOF){
		   		   		getline(recombFile,fileline);
		   	   	   		term = strtok ((char*)fileline.c_str()," \t");
	   		       		nextRec=(float)(atof(term));
	   		       		term = strtok (NULL, " \t");
	   		       		nextPos=atoi(term);
	   		   		}
	   		   		else{
		   		   		nextPos=-1;
	   		   		}
		   		}

		   		recombVector[i]=currentRec;
		   }



   	   }
   	   else{
	   	   error("The second column is not \"freq\" or \"pos\".\n");
   	   }




	   Rcpp::Rcout<<"The recombination rates between fragments are:"<<endl<<"\t";
	   for(int i=0;i<pieces-1;i++)Rcpp::Rcout<<"("<<i+1<<") "<<recombVector[i]<<" ";
	   Rcpp::Rcout<<"("<<pieces<<")"<<endl;
	   Rcpp::Rcout<<"*** END OF RECOMBINATION RATES PROFILE ***"<<endl;

}

void drawtrees(){

	map<int,string > tree;
	ostringstream  temp;


	for(int i=0;i<L;i++){
		tree.clear();

		for(int j=0;j<n;j++){
			if(tree[genepool[j*L+i]]!=""){
				temp.str("");
				temp <<tree[genepool[j*L+i]]<<","<<j+1<<":"<<genetimes[genepool[j*L+i]]-genetimes[j*L+i];
				tree[genepool[j*L+i]]=	temp.str();
			}
			else{
				temp.str("");
				temp <<j+1<<":"<<genetimes[genepool[j*L+i]]-genetimes[j*L+i];
				tree[genepool[j*L+i]]=	temp.str();
			}
		}


		for(int j=n*L;j<nextFree;j++){
			if(tree[j]!=""){
				if(genepool[j]==0){
				  Rcpp::Rcout<<"The tree for fragment "<<i+1<<"= ("<<tree[j]<<");"<<endl;
				}
				else{
					if(tree[genepool[j]]!=""){
						temp.str("");
						temp <<tree[genepool[j]]<<",("<<tree[j]<<"):"<<genetimes[genepool[j]]-genetimes[j];
						tree[genepool[j]]=	temp.str();
					}
					else{
						temp.str("");
						temp <<"("<<tree[j]<<"):"<<genetimes[genepool[j]]-genetimes[j];
						tree[genepool[j]]=temp.str();
					}
				}
			}
		}
	}

	tree.clear();
}




void genome(string popSize, 						// effective population size of each population
		  int nSubPOP, 							// number of populations
		  vector<int> & nSubSample, 			// number of samples (haplotypes) draw from each populations
		  int numPieces, 						// number of fragments for each sample (chromosome)
		  int pieceLen, 						// length in base pair of each fragment
		  int numIndepRegion, 					// number of independent regions (independent chromosome)
		  int s,								// fixed number of SNPs want to simulate, randomly place s SNPs on the genealogy, -1 means the number of mutation follows a Poisson distribution
		  string rec,							// recombination rate between consecutive fragments per generation
		  double mut, 							// mutation rate per generation per base pair
		  double mig, 							// migration rate per generation
		  vector< vector<bool> > &return_chromosome, // the return chromosome by connecting all independent chromosome into one long chromosome
		  long t,								// random seed
		  int drawtree,						// drawtree=1 output the genealogy trees for each fragment and chromosome; drawtree=0 do not output the trees
		  bool printparameters		   	 		// =1 print out all input parameters, =0 do not print input parameters
		  )
{
   if(atoi(popSize.c_str())>0){            	// if the population size is fixed during evoluation
	   N = atoi(popSize.c_str())*nSubPOP;
	   maxN = N;

	   fixedPOP = true;

	   popProfile.clear();
	   struct popProfilesStruct temp_popProfile;
	   popProfile.push_back(temp_popProfile);

	   popProfile[0].generation = 0;

	   for(int i=0;i<nSubPOP;i++)popProfile[0].popsize.push_back(atoi(popSize.c_str()));
	   for(int i=0;i<nSubPOP;i++)popProfile[0].popStart.push_back(i*atoi(popSize.c_str()));

	   Rcpp::Rcout<<"*** POPULATION PROFILE ***"<<endl;
	   Rcpp::Rcout<<"The size of each population is fixed to "<<atoi(popSize.c_str())<<endl;
	   Rcpp::Rcout<<"Sizes of populations = ";
	   for(unsigned int j=0;j<popProfile[0].popsize.size();j++)Rcpp::Rcout<<popProfile[0].popsize[j]<<" ";
	   Rcpp::Rcout<<endl;
// 	   cout<<"POP Start index = ";
// 	   for(int j=0;j<popProfile[0].popStart.size();j++)cout<<popProfile[0].popStart[j]<<" ";
// 	   cout<<endl;
	   Rcpp::Rcout<<"*** END OF POPULATION PROFILE ***"<<endl<<endl;

   }
   else{									// popSize stores the filename of the population profile during evoluation
 	   fixedPOP = false;

	   currentProfile = 0;

	   readPopProfile(popSize,nSubPOP);

	   N = popProfile[0].popStart.back()+popProfile[0].popsize.back();

	   // print out the population profile information.

	   Rcpp::Rcout<<"*** POPULATION PROFILE ***"<<endl;
	   Rcpp::Rcout<<"Maximum total size of populations at all generations: N = "<<maxN<<endl;
	   Rcpp::Rcout<<"Total size of populations at generation 0: N = "<<N<<endl<<endl;

	   for(unsigned int i=0;i<popProfile.size();i++){
	    Rcpp::Rcout<<"Generation = "<<popProfile[i].generation<<endl;
	    Rcpp::Rcout<<"Sizes of populations= ";
			for(unsigned int j=0;j<popProfile[i].popsize.size();j++)Rcpp::Rcout<<popProfile[i].popsize[j]<<" ";
			Rcpp::Rcout<<endl;
// 			cout<<"POP Start index = ";
// 			for(int j=0;j<popProfile[i].popStart.size();j++)cout<<popProfile[i].popStart[j]<<" ";
// 			cout<<endl;
			if(popProfile[i].popChange.size()>0){
			  Rcpp::Rcout<<"Populations change backward in time:"<<endl;
				for(unsigned int j=0;j<popProfile[i].popChange.size();j++){
				  Rcpp::Rcout<<" From "<<j+1<<" to ";
					for(unsigned int k=0;k<popProfile[i].popChange[j].size();k++)Rcpp::Rcout<<popProfile[i].popChange[j][k]+1<<" ";
					Rcpp::Rcout<<endl;
				}

// 				cout<<"Selection cumulative probability:"<<endl;
// 				for(int j=0;j<popProfile[i].selectProb.size();j++){
// 					cout<<" From "<<j+1<<" with prob = ";
// 					for(int k=0;k<popProfile[i].selectProb[j].size();k++)cout<<popProfile[i].selectProb[j][k]<<" ";
// 					cout<<endl;
// 				}
			}
			if(i==popProfile.size()-1)Rcpp::Rcout<<"*** END OF POPULATION PROFILE ***"<<endl;
			Rcpp::Rcout<<endl;
	   }
   }

   n=0;
   nPOP=nSubPOP;
   nSample.resize(nPOP);
   for(int i=0;i<nPOP;i++){
	   if(popProfile[0].popsize[i]<nSubSample[i])error("Population size should be larger than sample size!");
	   n += nSubSample[i];
	   nSample[i]=nSubSample[i];
   }

   L = numPieces;
   length = pieceLen;
   numChr = numIndepRegion;
   SNP=s;

   if(atof(rec.c_str())>0){				// if the recombination rate is fixed
	   recombination = atof(rec.c_str());
   }
   else{						// the distribution of recombination rate is described in the file (filename stored in rec)
   	   recombination = -1;
	   readRecombination(rec,L);
   }



   mutation = mut;
   migration = mig;

   if(t<=0)t=(long)(time(0)); // if t >0 we set the seed to t
   Rcpp::Rcout<<endl<<"random seed="<<t<<endl;
   globalRandom.Reset(t);

   if(printparameters){
	   Rprintf("popSize=%s, n=%d, nPOP=%d, L=%d, length=%d, numChr=%d, SNP=%d, rec=%s, mut=%.4e, mig=%.4e, drawtree=%d, printparameters=%d\n",popSize.c_str(),n,nPOP,L,length,numChr,SNP,rec.c_str(),mut,mig,drawtree,printparameters);
     Rcpp::Rcout<<"subSample = ";
	   for(int i=0;i<nPOP;i++){
	     Rcpp::Rcout<<nSample[i]<<" ";
	   }
	   Rcpp::Rcout<<endl;
   }

   if(N<=0 || n<=0 || L<=0 || length <=0 || numChr <=0 || atof(rec.c_str()) <0 || mut <0 || mig<0){
     Rcpp::Rcout<<"Input parameters have to be positive."<<endl;
	   Rprintf("N=%d, n=%d, L=%d, length=%d, numChr=%d, recombination=%s, mutation=%lf, migration=%lf\n",N,n,L,length,numChr,rec.c_str(),mut,mig);
	   Rcpp::stop("Normal exit.\n");
	   //exit(0);
   }

   poolsize = 2 * n * L - L;

   AllocateMemory(return_chromosome);

//    cout<<"memory OK!"<<endl<<endl;

	for(int chr=0;chr<numChr;chr++){

	   Initialize();

// 	   cout<<"Initialization OK!"<<endl;

	   bool done = false;
	   int  generation = 0, units = n * L;

	   while (!done)
	      {
	//       printf("\n\nGeneration %d -- Tracking %d units of sequence [pool used %.2f%%]   %c",
	//              generation++, units, nextFree * 100. / poolsize,
	//              !debug ? '\r' : '\n');

		  generation++;

		  if(!fixedPOP && currentProfile<popProfile.size()-1)
				if(generation == popProfile[currentProfile+1].generation)noChange=false;

	      //if (debug)
	      //   {
	      //   printf("F  : ");
	      //   for (int i = 0; i < L; i++)
	      //      printf("%3d ", fragments[i]);
	      //   printf("\n");
	      //   }

	      // Position of the first segment allocated to the current generation
	      int generationStart = nextFree;

	      // For each individual in the population
	      for (int i = 0; i < N && !done; i++)
	         if (child_flags[i])
	            {

	            /*if (debug) printf("%2d : ", i);*/

	            // The position we are currently tracking and its distance
	            // from the previous position
	            int distance = 0, parent = -1;
	            int startPiece = 0;
	            float recombProb;

			    for (int block = child_First[i]; block <= child_Last[i]; block++)
                if (children[i][block] == NULL)
                   distance += BLOCKSIZE;
                else
                   {
                   int pos = BLOCKSIZE * block;
                   int pos_in_block;
				   int pos_bound = (L < BLOCKSIZE * (block + 1) ? L : BLOCKSIZE * (block + 1) );

                   while (pos <  pos_bound )
                       {

	                   pos_in_block = pos % BLOCKSIZE;
                       // Skip through uninformative positions along each chromosome
                       if (children[i][block][pos_in_block] == -1)
                          {
                          //if (debug) printf("  . \n");
                          distance++;
                          pos++;
                          continue;
                          }

                       // Pick an ancestor at random for the current piece of chromosome

                       if(recombination>0){			// if the recombination rate is fixed over the genome
	                       recombProb = (float)(pow(1.0 - recombination, distance));
                       }
                       else{						// if the recombination rate is specified by the profile
	                       recombProb = 1;
	                       for(int recV_iter=startPiece;recV_iter<pos;recV_iter++)recombProb *= (1-recombVector[recV_iter]);
                       }


                       if (parent == -1 ||
                           (distance > 0 && globalRandom.Next() > recombProb)){
	                            parent = SampleParent(i, parent);
                        }


                       TouchParent(parent, pos);


                       // Reset distance between segments
                       distance = 1;
                       startPiece=pos;

                       //if (debug) printf("parent=%3d pos=%d address=%d address2=%d\n", parent,pos,parents[parent][block],parents[parent][pos/BLOCKSIZE]);

                       // Check if we sampled a new ancestral segment ...
                       if (parents[parent][block][pos_in_block] == -1)
                          {
                          // We delay sampling a new gene from the pool until we know
                          // a coalescent event occured ...
                          parents[parent][block][pos_in_block] = children[i][block][pos_in_block];

                          }
                       // Or if a coalescent event occured
                       else
                          {
                          // If we previously delayed sampling this parent
                          if (parents[parent][block][pos_in_block] < generationStart)
                            {
	                            genetimes[nextFree] = generation;
	                            genepool[parents[parent][block][pos_in_block]] = nextFree;
	                            parents[parent][block][pos_in_block] = nextFree++;

	                            if (nextFree == poolsize && units != 2)
	                               error("Genepool exhausted\n");
                            }

                          genepool[children[i][block][pos_in_block]] = parents[parent][block][pos_in_block];

                          fragments[pos]--;
                          units--;

                         if (fragments[pos] == 1)
                            {
                            // Reached the MRCA for this fragment
                            MRCAtimes[pos] = generation;
                            MRCA[pos] = parents[parent][block][pos_in_block];

                            genepool[parents[parent][block][pos_in_block]]=0;

                            parents[parent][block][pos_in_block] = -1;
                            units--;

                            // One less segment to track
                            if (--activeSegments == 0)
                               {
                               done = true;
                               break;
                               }
                            }
                         }
                      pos++;
                      }
                   }
	            }

	      NewGeneration();

	      }

	   if(drawtree){
	     Rcpp::Rcout<<endl<<"Genealogy trees for chromosome:"<<chr+1<<endl;
		   drawtrees();
	   }

	   // assign mutations

	   Rcpp::Rcout<<"Simulate mutation for Chromosome "<<chr+1<<endl;

	   FreeCoalescentMemory();

	   Rcpp::Rcout<<"Coalescent Memory free"<<endl;

	   InitializeMutation(return_chromosome);

	   Rcpp::Rcout<<"Mutation process initialized"<<endl;

// 	    cout<<"genepool:\t";
// 		for(int i=0;i<2*n*L-L;i++)cout<<genepool[i]<<" ";
// 		cout<<endl;
//
// 		cout<<"genetimes:\t";
// 		for(int i=0;i<2*n*L-L;i++)cout<<genetimes[i]<<" ";
// 		cout<<endl;
//
// 		cout<<"nextFree="<<nextFree<<endl;

// 		cout<<"mutation:\t";
// 		for(int i=0;i<nextFree;i++)cout<<branchMutation[i]<<" ";
// 		cout<<endl;

	   if(SNP<=0){
		   for(int i=0;i<nextFree;i++){
				if(genepool[i]!=0)branchMutation[i]=poissonInt((double)length*(double)(genetimes[genepool[i]]-genetimes[i])*mutation);
		   }
	   }
	   else{
		   unsigned long * cumulativeCount = new unsigned long [nextFree];

		   if(genepool[0]!=0)cumulativeCount[0] = genetimes[genepool[0]]-genetimes[0];
		   else cumulativeCount[0] = 0;

		   for(int i=1;i<nextFree;i++){
				if(genepool[i]!=0)cumulativeCount[i] = cumulativeCount[i-1] + genetimes[genepool[i]]-genetimes[i];
				else cumulativeCount[i] = cumulativeCount[i-1];
		   }

		   unsigned long totalGeneration = cumulativeCount[nextFree - 1];

// 		   cout<<"Total Generation = "<<totalGeneration<<endl;
//
// 		   if((unsigned long)SNP>totalGeneration){
// 				cout<<"Required number of SNPs("<<SNP<<") is more than total number of generations("<<totalGeneration<<")!"<<endl;
// 		   }

		   for(int i=0;i<SNP;i++){

			   unsigned long randomLong = globalRandom.NextInt()%totalGeneration;

			   int begin = 0;
			   int end = nextFree - 1;
			   int middle = (begin+end)/2;

			   while(1){
				   if(cumulativeCount[middle]<randomLong){
						begin = middle;
				   }
				   else{
						end = middle;
				   }

				   if(end == begin+1)break;

				   middle = (begin+end)/2;
			   }

			   if(cumulativeCount[end] < randomLong || (cumulativeCount[begin] >= randomLong && begin != 0)){
			    Rcpp::Rcout<<"searching error: begin="<<begin<<" end="<<end<<" random="<<randomLong<<" c1="<<cumulativeCount[begin]<<" c2="<<cumulativeCount[end]<<endl;
					error("see above message!\n");
			   }

			   int sampleIndex;

			   if(cumulativeCount[begin]>=randomLong)sampleIndex = begin;
			   else sampleIndex = end;

			   if(genepool[sampleIndex]==0){
			     Rcpp::Rcout<<"genepool error: begin="<<begin<<" end="<<end<<" random="<<randomLong<<" c1="<<cumulativeCount[begin]<<" c2="<<cumulativeCount[end]<<endl;
				   error("see above message!\n");
			   }

			   branchMutation[sampleIndex]++;
		   }

		   delete [] cumulativeCount;

	   }

	   delete [] genetimes;

	   Rcpp::Rcout<<"Mutation assigned"<<endl;

//  	cout<<"mutation:\t";
// 		for(int i=0;i<nextFree;i++)cout<<branchMutation[i]<<" ";
// 		cout<<endl;

	   AllocateMutationMemory();

	   Rcpp::Rcout<<"Mutation Memory allocated"<<endl;

	   PrepareMutation();

	   for(int i=0;i<n;i++){
	   		for(int j=0;j<L;j++){
		   		if(branchMutation[i*L+j]>0){
			   		geneIndex[i*L+j] = fragmentLength[j];
			   		fragmentLength[j] += branchMutation[i*L+j];
		   		}
		   		branchFragment[i*L+j] = j;
	   		}
	   }


	   for(int i=0;i<nextFree;i++){
			if(genepool[genepool[i]]!=0 && genepool[i]!=0){
			   if(branchMutation[genepool[i]]>0 && branchFragment[genepool[i]]== -1 ){
				   geneIndex[genepool[i]] = fragmentLength[branchFragment[i]];
				   fragmentLength[branchFragment[i]] += branchMutation[genepool[i]];
			   }
			   branchFragment[genepool[i]]=branchFragment[i];
			}
	   }

	   //cout<<"branchfragment done"<<endl;

//    	cout<<"branchFragment:\t";
// 		for(int i=0;i<nextFree;i++)cout<<branchFragment[i]<<" ";
// 		cout<<endl;

//    	cout<<"geneIndex:\t";
// 		for(int i=0;i<nextFree;i++)cout<<geneIndex[i]<<" ";
// 		cout<<endl;
//
// 		cout<<"fragmentLength:\t";
// 		for(int i=0;i<L;i++)cout<<fragmentLength[i]<<" ";
// 		cout<<endl;


// 	    cout<<"genepool genetime mutation frag index\n";
// 		for(int i=0;i<nextFree;i++)cout<<genepool[i]<<" "<<genetimes[i]<<" "<<branchMutation[i]<<" "<<branchFragment[i]<<" "<<geneIndex[i]<<endl;
// 		cout<<endl;



       delete [] branchFragment;

       Rcpp::Rcout<<"Internal memory free"<<endl;

	   // print the fragment index for each SNP
	   Rcpp::Rcout<<"SNP genetic position scaled to [0,1]:"<<endl;
	   if(L>1){
		   for(int i=0;i<L;i++){
				for(int j=0;j<fragmentLength[i];j++)Rcpp::Rcout<<(float)i/(float)(L-1)<<" ";
		   }
	   }
	   else{
	     Rcpp::Rcout<<"Only one fragment to be simulated. No recombination between SNPs.";
	   }
	   Rcpp::Rcout<<endl;

	   int totalLength=0;
	   if(SNP <=0){
		   for(int i=0;i<L;i++)totalLength += fragmentLength[i];
	   }


	   for(int i=0;i<n;i++){
	   		for(int j=0;j<L;j++){
	   			chromosome[i][j].resize(fragmentLength[j],0);

	   			for(int k=0;k<branchMutation[i*L+j];k++)chromosome[i][j][geneIndex[i*L+j]+k] = 1;

		   		int current=i*L+j;

		   		while(genepool[genepool[current]]!=0){
			   		for(int k=0;k<branchMutation[genepool[current]];k++)chromosome[i][j][geneIndex[genepool[current]]+k]=1;
			   		current=genepool[current];
		   		}
	   		}
	   }

	   Rcpp::Rcout<<"Chromosome done"<<endl;

	   FreeMutationMemory();

	   if(SNP <=0){
		   for(int i=0;i<n;i++){
				return_chromosome[i].reserve(return_chromosome[i].size()+totalLength);
		   }
	   }

	   for(int i=0;i<n;i++){
			for(int j=0;j<L;j++){
				return_chromosome[i].insert(return_chromosome[i].end(),chromosome[i][j].begin(),chromosome[i][j].end());
			}
	   }

	   Rcpp::Rcout<<"Return chromosome done"<<endl;

	   // print the chromosomes
// 	   for(int i=0;i<n;i++){
// 			printf("Chr %d: ",i);
// 			for(int j=0;j<L;j++){
// 				for(int k=0;k<chromosome[i][j].size();k++)cout<<chromosome[i][j][k];
// 				cout<<" ";
// 			}
// 			cout<<endl;
// 	   }
//
//
// 	  for(int i=0;i<n;i++){
// 			printf("return Chr %d: ",i);
// 			for(int j=0;j<return_chromosome[i].size();j++){
// 				cout<<return_chromosome[i][j];
// 			}
// 			cout<<endl;
// 	   }



   }

   FreeMemory();

   Rprintf("\nDone simulating ARG ...\n");



}


