////////////////////////////////////////////////////////////////////// 
// MathConstant.h 
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



#ifndef __MATHCONSTANT_H__
#define __MATHCONSTANT_H__

#ifdef  _MSC_VER
#define _USE_MATH_DEFINES 
#endif

#include "math.h"
#include "stdlib.h"

// Constants for numerical routines
//

#define TINY    1.0e-30        // A small number
#define ITMAX   200            // Maximum number of iterations
#define EPS     3.0e-7         // Relative accuracy
#define ZEPS    3.0e-10        // Precision around zero
#define FPMIN   1.0e-30        // Number near the smallest representable number
#define FPMAX   1.0e+100       // Number near the largest representable number
#define TOL     1.0e-6         // Zero SVD values below this
#define GOLD    0.61803399     // Golden ratio
#define CGOLD   0.38196601     // Complement of golden ratio

inline double square(double a)         { return a * a; }
inline double sign(double a, double b) { return b >= 0 ? fabs(a) : -fabs(a); }
inline double min(double a, double b)  { return a < b ? a : b; }
inline double max(double a, double b)  { return a > b ? a : b; }

inline int square(int a)      { return a * a; }
inline int sign(int a, int b) { return b >= 0 ? abs(a) : -abs(a); }
inline int min(int a, int b)  { return a < b ? a : b; }
inline int max(int a, int b)  { return a > b ? a : b; }

#endif
