/*
 Copyright (c)  2016 Kornev Andrey A.
                     Ozeritsky  Alexey V.
                     and
                     Afanasieva A.
                     Bekbulatov R.
                     Nazarov V.
                     Halikov P.
*/

#include "f.h"

int numsds_spectral_problem (
    double *eigen_values,
    double *eigen_functions,
    int dim,
    int eigenvalues_number,
    int max_iterations,
    double tolerance,
    const char spectralSubSet[3],
    void (*A_op)(double *, const double *, int, void *),
    void * user_data
    )

{
  int n;
  int enx;

  //	int nx[1];

  int i, k, count;
  double dl;

  int info = 0;

  /* Number of eigenvalues of OP to be computed. 0 < NEV < N-1.*/
  int nev;

  /* Number of columns of the matrix V. NCV must satisfy the two
                inequalities 2 <= NCV-NEV and NCV <= N.*/
  int ncv;

  int lworkl;
  int ldv;

  int *c;

  double * resid;
  double * v;
  double * workd;
  double * workev;
  double * workl;

  double *wr;
  double *wi;

  /* BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x*/
  /* BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x*/
  char bmat[] = "I";

  /* 'LM' -> want the NEV eigenvalues of largest magnitude. */
  /* 'SM' -> want the NEV eigenvalues of smallest magnitude.*/
  /* 'LR' -> want the NEV eigenvalues of largest real part. */
  /* 'SR' -> want the NEV eigenvalues of smallest real part.*/
  /* 'LI' -> want the NEV eigenvalues of largest imaginary part. */
  /* 'SI' -> want the NEV eigenvalues of smallest imaginary part.*/
  char which[] = "LM";
  int  iparam[11], ipntr[14];
  int  ishfts = 1;
  int  mode   = 1;
  int  ido    = 0;
  int it1 = 0;
  //	double * tmp=NULL;
  double tol;

  char ch[] = "A";
  int rvec   = 1;
  int *select;
  double sigmar = 0, sigmai = 0;
  int ierr = 0;
  int nconv;


  n = dim;
  enx = eigenvalues_number;
  //	nx[0] = enx;
  nev = enx;
  ncv  = ((2 * nev + 2) < n) ? (2 * nev + 2) : n;
  lworkl= 3 * ncv * ncv + 6 * ncv;
  ldv   = n;

  c      = sp_alloc_i_vector (ncv, __FUNCTION__, __FILE__);
  resid  = sp_alloc_d_vector (n, __FUNCTION__, __FILE__);
  v      = sp_alloc_d_vector (ldv * ncv, __FUNCTION__, __FILE__);
  workd  = sp_alloc_d_vector (3 * n, __FUNCTION__, __FILE__);
  workev = sp_alloc_d_vector (3 * ncv, __FUNCTION__, __FILE__);
  workl  = sp_alloc_d_vector (lworkl, __FUNCTION__, __FILE__);

  wr   = sp_alloc_d_vector(n, __FUNCTION__, __FILE__);
  wi   = sp_alloc_d_vector(n, __FUNCTION__, __FILE__);


  if (spectralSubSet != 0 && strlen (spectralSubSet) > 1)
    {
      which[0] = spectralSubSet[0];
      which[1] = spectralSubSet[1];
    }
  else
    {
      printf ("Arnoldi: incorrect spectralSubSet parametr\n");
      exit (1);
    }

  // tmp = sp_alloc_d_vector(n, __FUNCTION__, __FILE__);
  tol = tolerance;

  iparam[1-1] = ishfts;
  iparam[3-1] = max_iterations;
  iparam[7-1] = mode;

  do {
      dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
              &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl,
              &info);

      if (nummsds_check_dnaupd_status (info) != 0)
        {
          printf ("Arnoldi: dnaupd_ call error\n");
          exit (1);
        }


      if (ido == -1 || ido == 1)
        {
          int addr2 = ipntr[2-1]-1;
          int addr1 = ipntr[1-1]-1;
          double *w1 = &workd[addr2];
          double *v1 = &workd[addr1];

          A_op (w1, v1, n, user_data);

          ++it1;
          printf ("Arnoldi iter -> %d\n",it1);
        }
    } while (ido == -1 || ido == 1);


  select = sp_alloc_i_vector (ncv, __FUNCTION__, __FILE__);

  dneupd_ ( &rvec, ch, select, wr, wi, v, &ldv,
            &sigmar, &sigmai, workev, bmat, &n, which, &nev, &tol,
            resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, &ierr );

  if (nummsds_check_dneupd_status (ierr) != 0)
    {
      printf ("Arnoldi: dneupd_ call error\n");
      exit (1);
    }

  nconv = iparam[5 - 1];

  if (nconv < enx)
    {
      printf("Arnoldi: warning: found %d vectors - less then requested %d\n",
             nconv, enx);
    }


  {
    printf("\n\nArnoldi: found %d vectors. Accuracy=%.2e\n", nconv, tol);
  }



  //------------------------------------------------------
  if((1 == 0) && (nconv == 3) )
    {//!!!!!!!!!!!!!!!!!!!!!!!!
      //
      int ko;
      double tmp;

      if( fabs(wr[0]-wr[1])+fabs(wi[0]-wi[1]) >2.*tol*10.)
        {
          if( fabs(wr[2]-wr[1])+fabs(wi[2]-wi[1]) >2.*tol*10.)
            ko=1;
          else
            ko=0;

          for(k=0;k<n;k++)
            {	tmp = v[ko * n + k];
              v[ko * n + k]=v[2*n+k];
              v[2*n+k] = tmp;
            }

        }


    }
  //------------------------------------------------------
  if( (1==0)&&(nconv==3) )
    {//!!!!!!!!!!!!!!!!!!!!!!!!
      //
      int ko;
      double tmp;
      double dl0;
      double dl1;
      double dl2;

      i=0;
      dl0=sqrt(wr[i]*wr[i]+ wi[i]*wi[i]);
      i=1;
      dl1=sqrt(wr[i]*wr[i]+ wi[i]*wi[i]);
      i=2;
      dl2=sqrt(wr[i]*wr[i]+ wi[i]*wi[i]);

      ko=2;
      if(dl2>dl0)ko=0;
      if(dl2>dl1)ko=1;
      if(ko!=2){

          for(k=0;k<n;k++)
            {	tmp = v[ko * n + k];
              v[ko * n + k]=v[2*n+k];
              v[2*n+k] = tmp;
            }
          tmp = wr[ko];
          wr[ko]=wr[2];
          wr[2] = tmp;

          tmp = wi[ko];
          wi[ko]=wi[2];
          wi[2] = tmp;
        }

    }
  //------------------------------------------------------



  for (i = 0; i < nconv; ++i) {
      dl = sqrt (wr[i] * wr[i] + wi[i] * wi[i]);

      printf("         lambda[%2d]=(%24.16e,%24.16e), |lambda| = %24.16e\n",
             i+1, wr[i], wi[i],dl);

    }


  for (count = 0; count < eigenvalues_number; count++) {

      for (k = 0; k < n; k++) {
          eigen_functions[count * n + k] = v[count * n + k];
        }

      eigen_values[2*count+0] = wr[count];
      eigen_values[2*count+1] = wi[count];
    }


  sp_free_d_vector(wr);
  sp_free_d_vector(wi);
  sp_free_i_vector(c);

  return nconv;
}

int nummsds_check_dnaupd_status (int info)
{
  switch (info) {
    case 0:
      //          =  0: Normal exit.
      break;
    case 1:
      printf("1: Maximum number of iterations taken.\n");
      printf("All possible eigenvalues of OP has been found. IPARAM(5)\n");
      printf("returns the number of wanted converged Ritz values.\n");
      break;
    case 2:
      printf("2: No inter an informational error. Deprecated starting\n");
      printf("with release 2 of ARPACK.\n");
      break;
    case 3:
      printf("3: No shifts could be applied during a cycle of the\n");
      printf("Implicitly restarted Arnoldi iteration. One possibility\n");
      printf("is to increase the size of NCV relative to NEV.\n");
    case -1:
      printf("-1: N must be positive.\n");
      break;
    case -2:
      printf("-2: NEV must be positive.\n");
      break;
    case -3:
      printf("-3: NCV-NEV >= 2 and less than or equal to N.\n");
      break;
    case -4:
      printf("-4: The maximum number of Arnoldi update iteration\n");
      printf("must be greater than zero.\n");
      break;
    case -5:
      printf("5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
      break;
    case -6:
      printf("6: BMAT must be one of 'I' or 'G'.\n");
      break;
    case -7:
      printf("-7: Length of private work array is not sufficient.\n");
      break;
    case -8:
      printf("-8: Error return from LAPACK eigenvalue calculation;\n");
      break;
    case -9:
      printf("-9: Starting vector is zero.\n");
      break;
    case -10:
      printf("-10: IPARAM(7) must be 1,2,3,4.\n");
      break;
    case -11:
      printf("-11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.\n");
      break;
    case -12:
      printf("-12: IPARAM(1) must be equal to 0 or 1.\n");
      break;
    case -9999:
      printf("-9999: Could not build an Arnoldi factorization.\n");
      printf("IPARAM(5) returns the size of the current Arnoldi\n");
      printf("factorization.\n");
      break;
    default:
      printf("unknown error\n");
      return -100000;
      break;
    }
  return info;
}

int nummsds_check_dneupd_status (int info)
{
  switch(info) {
    case 0:
      //c          =  0: Normal exit.
      break;
    case 1:
      printf("1: The Schur form computed by LAPACK routine dlahqr\n");
      printf("could not be reordered by LAPACK routine dtrsen .\n");
      printf("Re-enter subroutine dneupd  with IPARAM(5)=NCV and \n");
      printf("increase the size of the arrays DR and DI to have \n");
      printf("dimension at least dimension NCV and allocate at least NCV \n");
      printf("columns for Z. NOTE: Not necessary if Z and V share \n");
      printf("the same space. Please notify the authors if this error\n");
      printf("occurs.\n");
      break;
    case -1:
      printf("-1: N must be positive.\n");
      break;
    case -2:
      printf("= -2: NEV must be positive.\n");
      break;
    case -3:
      printf("-3: NCV-NEV >= 2 and less than or equal to N.\n");
      break;
    case -5:
      printf("-5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
      break;
    case -6:
      printf("-6: BMAT must be one of 'I' or 'G'.\n");
      break;
    case -7:
      printf("-7: Length of private work WORKL array is not sufficient.\n");
      break;
    case -8:
      printf("-8: Error return from calculation of a real Schur form.\n");
      printf("Informational error from LAPACK routine dlahqr .\n");
      break;
    case -9:
      printf("-9: Error return from calculation of eigenvectors.\n");
      printf("Informational error from LAPACK routine dtrevc .\n");
      break;
    case -10:
      printf("-10: IPARAM(7) must be 1,2,3,4.\n");
      break;
    case -11:
      printf("-11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
      break;
    case -12:
      printf("-12: HOWMNY = 'S' not yet implemented\n");
      break;
    case -13:
      printf("-13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
      break;
    case -14:
      printf("-14: DNAUPD  did not find any eigenvalues to sufficient\n");
      printf("accuracy.\n");
      break;
    case -15:
      printf("-15: DNEUPD got a different count of the number of converged\n");
      printf("Ritz values than DNAUPD got.  This indicates the user\n");
      printf("probably made an error in passing data from DNAUPD to\n");
      printf("DNEUPD or that the data was modified before entering\n");
      printf("DNEUPD\n");
      break;
    default:
      printf("unknown error\n");
      return -10000;
    }
  return info;
}
