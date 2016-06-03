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

int find_eigen_values (
  double *eigen_values,
  double *eigen_functions,
  const int dim,
  const int n_eigen_values,
  const int max_iterations,
  const double tolerance,
  const char *spectralSubSet,
  const char *bmat,
  const void *user_data,
  const double *G,
  const double *V1,
  const double *V2,
  const int *st,
  const int *M0L,
  const int *M0R,
  const double *X,
  const double *Y
)

{
  int n = dim;
  int ido = 0;
  int iters = 0;
  int rvec = 1;
  int info = 0;
  double tol = tolerance;
  double sigma_real = 0, sigma_imag = 0;
  int n_found_eigen_values;

  /* Number of eigenvalues of OP to be computed. 0 < NEV < N-1. */
  int nev = n_eigen_values;

  /* Number of columns of the matrix V. NCV must satisfy the two inequalities 2 <= NCV-NEV and NCV <= N. */
  int ncv    = ((2 * nev + 2) < n) ? (2 * nev + 2) : n;
  int lworkl = 3 * ncv * ncv + 6 * ncv;
  int ldv    = n;

  double *resid  = make_vector_double (n, __FUNCTION__, __FILE__);
  double *v      = make_vector_double (ldv * ncv, __FUNCTION__, __FILE__);
  double *workd  = make_vector_double (3 * n, __FUNCTION__, __FILE__);
  double *workev = make_vector_double (3 * ncv, __FUNCTION__, __FILE__);
  double *workl  = make_vector_double (lworkl, __FUNCTION__, __FILE__);

  double *w_real = make_vector_double (n, __FUNCTION__, __FILE__);
  double *w_imag = make_vector_double (n, __FUNCTION__, __FILE__);

  int *select = make_vector_int (ncv, __FUNCTION__, __FILE__);

  int iparam[11], ipntr[14];
  iparam[1 - 1] = 1; // ishfts
  iparam[3 - 1] = max_iterations;
  iparam[7 - 1] = 1; // mode

  /// Calculations run

  printf ("Arnoldi iter = %4.d", iters);

  do
    {
      dnaupd_ (&ido, bmat, &n, spectralSubSet, &nev, &tol, resid,
               &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl,
               &info);

      if (check_dnaupd_status (info) != 0)
        {
          printf ("Arnoldi: dnaupd_ call error\n");
          exit (1);
        }

      if (ido == -1 || ido == 1)
        {
          int addr2 = ipntr[2 - 1] - 1;
          int addr1 = ipntr[1 - 1] - 1;
          double *w1 = &workd[addr2];
          double *v1 = &workd[addr1];

          A_op (w1, v1, user_data, G, V1, V2, st, M0L, M0R, X, Y);

          iters++;

          printf ("\b\b\b\b");
          printf ("%4.d", iters);
          fflush (stdout);
        }
    }
  while (ido == -1 || ido == 1);

  printf ("\n");

  dneupd_ ( &rvec, "A", select, w_real, w_imag, v, &ldv,
            &sigma_real, &sigma_imag, workev, bmat, &n, spectralSubSet, &nev, &tol,
            resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
            &lworkl, &info );

  if (check_dneupd_status (info) != 0)
    {
      printf ("Arnoldi: dneupd_ call error\n");
      exit (1);
    }

  n_found_eigen_values = iparam[5 - 1];

  if (n_found_eigen_values < n_eigen_values)
    {
      printf ("Arnoldi: warning: found %d vectors - less then requested %d\n", n_found_eigen_values, n_eigen_values);
    }

  printf ("\nArnoldi: found %d vectors. Accuracy = %.2e\n", n_found_eigen_values, tol);

  for (int i = 0; i < n_found_eigen_values; ++i)
    {
      double dl = sqrt (w_real[i] * w_real[i] + w_imag[i] * w_imag[i]);

      printf ("\tlambda_%d = %11.5e %c %11.5e * i, |lambda| = %11.5e\n", i + 1, w_real[i], w_imag[i] < 0. ? '-' : '+', fabs (w_imag[i]), dl);
    }

  for (int count = 0; count < n_eigen_values; count++)
    {
      for (int k = 0; k < n; k++)
        {
          eigen_functions[count * n + k] = v[count * n + k];
        }

      eigen_values[2 * count + 0] = w_real[count];
      eigen_values[2 * count + 1] = w_imag[count];
    }

  FREE_ARRAY (resid);
  FREE_ARRAY (v);
  FREE_ARRAY (workd);
  FREE_ARRAY (workev);
  FREE_ARRAY (workl);

  FREE_ARRAY (w_real);
  FREE_ARRAY (w_imag);

  FREE_ARRAY (select);

  return n_found_eigen_values;
}

int check_dnaupd_status (int info)
{
  if (info != 0)
    {
      printf ("\n");
    }

  switch (info)
    {
      case 0:
        //          =  0: Normal exit.
        break;

      case 1:
        printf ("1: Maximum number of iterations taken.\n");
        printf ("All possible eigenvalues of OP has been found. IPARAM(5)\n");
        printf ("returns the number of wanted converged Ritz values.\n");
        break;

      case 2:
        printf ("2: No inter an informational error. Deprecated starting\n");
        printf ("with release 2 of ARPACK.\n");
        break;

      case 3:
        printf ("3: No shifts could be applied during a cycle of the\n");
        printf ("Implicitly restarted Arnoldi iteration. One possibility\n");
        printf ("is to increase the size of NCV relative to NEV.\n");

      case -1:
        printf ("-1: N must be positive.\n");
        break;

      case -2:
        printf ("-2: NEV must be positive.\n");
        break;

      case -3:
        printf ("-3: NCV-NEV >= 2 and less than or equal to N.\n");
        break;

      case -4:
        printf ("-4: The maximum number of Arnoldi update iteration\n");
        printf ("must be greater than zero.\n");
        break;

      case -5:
        printf ("5: spectralSubSet must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
        break;

      case -6:
        printf ("6: BMAT must be one of 'I' or 'G'.\n");
        break;

      case -7:
        printf ("-7: Length of private work array is not sufficient.\n");
        break;

      case -8:
        printf ("-8: Error return from LAPACK eigenvalue calculation;\n");
        break;

      case -9:
        printf ("-9: Starting vector is zero.\n");
        break;

      case -10:
        printf ("-10: IPARAM(7) must be 1,2,3,4.\n");
        break;

      case -11:
        printf ("-11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.\n");
        break;

      case -12:
        printf ("-12: IPARAM(1) must be equal to 0 or 1.\n");
        break;

      case -9999:
        printf ("-9999: Could not build an Arnoldi factorization.\n");
        printf ("IPARAM(5) returns the size of the current Arnoldi\n");
        printf ("factorization.\n");
        break;

      default:
        printf ("unknown error\n");
        return -100000;
        break;
    }

  return info;
}

int check_dneupd_status (int info)
{
  switch (info)
    {
      case 0:
        //c          =  0: Normal exit.
        break;

      case 1:
        printf ("1: The Schur form computed by LAPACK routine dlahqr\n");
        printf ("could not be reordered by LAPACK routine dtrsen .\n");
        printf ("Re-enter subroutine dneupd  with IPARAM(5)=NCV and \n");
        printf ("increase the size of the arrays DR and DI to have \n");
        printf ("dimension at least dimension NCV and allocate at least NCV \n");
        printf ("columns for Z. NOTE: Not necessary if Z and V share \n");
        printf ("the same space. Please notify the authors if this error\n");
        printf ("occurs.\n");
        break;

      case -1:
        printf ("-1: N must be positive.\n");
        break;

      case -2:
        printf ("= -2: NEV must be positive.\n");
        break;

      case -3:
        printf ("-3: NCV-NEV >= 2 and less than or equal to N.\n");
        break;

      case -5:
        printf ("-5: spectralSubSet must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'\n");
        break;

      case -6:
        printf ("-6: BMAT must be one of 'I' or 'G'.\n");
        break;

      case -7:
        printf ("-7: Length of private work WORKL array is not sufficient.\n");
        break;

      case -8:
        printf ("-8: Error return from calculation of a real Schur form.\n");
        printf ("Informational error from LAPACK routine dlahqr .\n");
        break;

      case -9:
        printf ("-9: Error return from calculation of eigenvectors.\n");
        printf ("Informational error from LAPACK routine dtrevc .\n");
        break;

      case -10:
        printf ("-10: IPARAM(7) must be 1,2,3,4.\n");
        break;

      case -11:
        printf ("-11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n");
        break;

      case -12:
        printf ("-12: HOWMNY = 'S' not yet implemented\n");
        break;

      case -13:
        printf ("-13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.\n");
        break;

      case -14:
        printf ("-14: DNAUPD  did not find any eigenvalues to sufficient\n");
        printf ("accuracy.\n");
        break;

      case -15:
        printf ("-15: DNEUPD got a different count of the number of converged\n");
        printf ("Ritz values than DNAUPD got.  This indicates the user\n");
        printf ("probably made an error in passing data from DNAUPD to\n");
        printf ("DNEUPD or that the data was modified before entering\n");
        printf ("DNEUPD\n");
        break;

      default:
        printf ("unknown error\n");
        return -10000;
    }

  return info;
}
