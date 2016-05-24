/*

 Copyright (c)  2016 Kornev Andrey A.
                     Ozeritsky  Alexey V.

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:
 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 3. The name of the author may not be used to endorse or promote products
    derived from this software without specific prior written permission

 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FIX_UNUSED(X) (void)(X)

typedef struct
{
  int Nx;
  int Ny;
  int N;
  int NA;

  double Lx;
  double Ly;

  double Hx;
  double Hy;

  double *diag;

  double mu;

} UserDataCurr_struct;



void initparam_UserDataCurr_struct (UserDataCurr_struct * udc);

int  L_op (double * Lapl, const double *u, const UserDataCurr_struct * udc);

void print_2dfun_double (FILE* f, const char * name, const double  * u,
                        const int nx, const int ny);

int convert_u_to_au (double * au, const double  * u,
                     const UserDataCurr_struct * udc);

int convert_au_to_u (double *u, const double  * au,
                     const UserDataCurr_struct * udc);

void A_op (double *Aau, const double *au, int n, void * udc);

double * make_vector_double (int n,const char *info_1, const char *info_2);

int numsds_spectral_problem (
    double *eigen_values,
    double *eigen_functions,
    int dim,
    int eigenvalues_number,
    int max_iterations,
    double tolerance,
    const char *spectralSubSet,
    void (*)(double *, const double *, int, void *),
    void * user_data
    );
