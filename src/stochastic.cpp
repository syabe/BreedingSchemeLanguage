////////////////////////////////////////////////////////////////////// 
// stochastic.cpp 
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

#include "Random.h"
#include <math.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;

#define Pi 3.1415926535897932

double lngamma(double z)
{
	double result,sum;
	static double c[7]={1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
	
	sum=c[0];
	for(int i=1;i<=6;i++)sum += c[i]/(z+i);
	
	result = (z+0.5)*log(z+5.5)-(z+5.5)+log(2.5066282746310005*sum/z);
	
	return result;
}



int poissonInt(double lambda) 
{
	static double mean=-1, waitTime,s,logmean,c;
	double count, eventTime,ratio,y;
	
	static double maxInt=pow(2.0,(double)(8*sizeof(int)-1))-1;
	
	if(lambda <0){
	  Rcpp::Rcout<<"The mean of Poisson should be >=0, input="<<lambda<<endl;
	  Rcpp::stop("poissonInt error! \n");
		//exit(1);	
	}
	else if(lambda < 12){
		if(lambda != mean){
			mean = lambda;
			waitTime = exp(-lambda);	
		}
		
		count = 0;
		eventTime = globalRandom.Next();
		
		while(eventTime > waitTime){
			count++;
			eventTime *= globalRandom.Next();	
		}
	}
	else{
		if(lambda != mean){
			mean = lambda;
			s=sqrt(2*lambda);
			logmean=log(lambda);
			c=lambda*log(lambda)-lngamma(lambda+1);
		}
		
		do{
			do{
				y=tan(Pi*globalRandom.Next());
				count = floor(s*y+lambda);
			}while(count<0);
			
			ratio = 0.9*(1+y*y)*exp(count*logmean-lngamma(count+1)-c);
			
		}while(globalRandom.Next()>ratio);
		
	}
	
	if(count>maxInt)count=maxInt;
	
	return (int)count;		
}




