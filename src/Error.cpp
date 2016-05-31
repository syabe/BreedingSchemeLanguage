////////////////////////////////////////////////////////////////////// 
// Error.cpp 
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

#include "Error.h"
#include <Rcpp.h>

#include "stdlib.h"
#include "stdarg.h"
#include "stdio.h"

// Declare a dummy class to ensure that compilers recognize this as C++ code
class String;

void error ( const char * msg, ... )
   {
   va_list  ap;

   va_start(ap, msg);

   Rprintf("\nFATAL ERROR - \n");
   Rprintf(msg);
   Rprintf("\n\n");

   va_end(ap);
   
   Rcpp::stop("Error!\n");

//   exit(EXIT_FAILURE);
   }

void warning ( const char * msg, ... )
   {
   va_list  ap;

   va_start(ap, msg);

   Rprintf("\n\aWARNING - \n");
   Rprintf(msg);
   Rprintf("\n");

   va_end(ap);
   }

void numerror ( const char * msg , ... )
   {
   va_list  ap;

   va_start(ap, msg);

   Rprintf("\nFATAL NUMERIC ERROR - ");
   Rprintf(msg);
   Rprintf("\n\n");

   va_end(ap);
   
   Rcpp::stop("Error!\n");

//   exit(EXIT_FAILURE);
   }
