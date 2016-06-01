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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#define FIX_UNUSED(X) (void)(X)
#define FREE_ARRAY(X) if ((X)) free ((X))

#define SQUARE 1

typedef struct
{
  int Nx;
  int Ny;
  int Nx_0;
  int Ny_0;
  int N;
  int NA;

  double Lx;
  double Ly;
  double Lx_0;
  double Ly_0;

  double Hx;
  double Hy;

  double p_ro;
  double p_2ro;
  double mu;

} UserDataCurr_struct;



void initparam_UserDataCurr_struct (UserDataCurr_struct * udc);

int L_op (double *Lu, const double *u, const UserDataCurr_struct * udc,
          const double *G, const double *V1, const double *V2, const int *st, const int *M0L, const int *M0R);

void fill_node_phys_prop (int m, // number of mesh node
                          double *p00, double *pL0, double *pR0,
                          double *p0L, double *p0R, double *pLL,
                          double *pRL, double *pLR, double *pRR,
                          const double *const p, const int *const M0L, const int *const M0R);

void calc_mesh_params (int *st, double *X, double *Y, int *M0L,
                       int *M0R, const UserDataCurr_struct * udc);

void print_2dfun_double (FILE* f, const char * name, const double  * u,
                        const int nx, const int ny);

int convert_u_to_au (double * au, const double  * u,
                     const UserDataCurr_struct * udc, const int *st);

int convert_au_to_u (double *u, const double  * au,
                     const UserDataCurr_struct * udc, const int *st);

void A_op (double *Aau, const double *au, int n, void * udc, const double *G, const double *V1, const double *V2, const int *st, const int *M0L, const int *M0R);

double *make_vector_double (int n, const char *info_1, const char *info_2);

int *make_vector_int (int n, const char *info_1, const char *info_2);

int *sp_alloc_i_vector (int n, const char *info_1, const char *info_2);

void sp_free_i_vector (int *);

double *sp_alloc_d_vector (int n, const char *info_1, const char *info_2);

void sp_free_d_vector (double *);

int read_stationary_solution (const char *fname, int N, double *G, double *V1, double *V2);

void dnaupd_(
    int *ido, char *bmat, int *n, char *which, int*nev, double *tol,
    double *resid, int *ncv, double *v, int *ldv, int *iparam,
    int *ipntr, double *workd, double *workl, int *lworkl, int *info);

void dneupd_(
    int *vec, char *c, int *select, double *d, double * /*d(1,2)*/,
    double *v, int *ldv, double *sigmar, double *sigmai, double *workev,
    char *bmat, int *n, char *which, int *nev, double *tol,
    double *resid, int *ncv, double *vv, int *ldvv, int *iparam,
    int *ipntr, double *workd, double *workl, int *lworkl, int *ierr);

int nummsds_check_dnaupd_status (int info);

int nummsds_check_dneupd_status (int info);

int numsds_spectral_problem (
    double *eigen_values,
    double *eigen_functions,
    int dim,
    int eigenvalues_number,
    int max_iterations,
    double tolerance,
    const char spectralSubSet[3],
    void (*)(double *, const double *, int, void *, const double *, const double *,  const double *, const int *, const int *, const int *),
    void * user_data,
    const double *G,
    const double *V1,
    const double *V2,
    const int *st,
    const int *M0L,
    const int *M0R
    );

void recreate_coefficients (
        double *a_J_0L, double *a_W1_0L, double *a_W2_0L,
        double *a_J_L0, double *a_W1_L0, double *a_W2_L0,
        double *a_J_00, double *a_W1_00, double *a_W2_00,
        double *a_J_R0, double *a_W1_R0, double *a_W2_R0,
        double *a_J_0R, double *a_W1_0R, double *a_W2_0R,
        double *a_J_LL, double *a_W1_LL, double *a_W2_LL,
        double *a_J_LR, double *a_W1_LR, double *a_W2_LR,
        double *a_J_RL, double *a_W1_RL, double *a_W2_RL,
        double *a_J_RR, double *a_W1_RR, double *a_W2_RR);
