/*

 Copyright (c)  2016 Kornev Andrey A.
                     Ozeritsky  Alexey V.
                     and
                     Afanasieva A.
                     Bekbulatov R.
                     Nazarov V.
                     Halikov P.
*/

#include <stdio.h>
#include "f.h"

#define ST_SOL_FILE "stat_sol.txt"

int main (void)
{
  UserDataCurr_struct  udc;

  double *eigen_values;
  double *eigen_functions_A_op;
  double *eigen_functions;

  // aux vectors for mesh parameters
  int *st;
  int *M0L;
  int *M0R;
  double *X;
  double *Y;

  // Values from file
  double *G;
  double *V1;
  double *V2;

  int eigenvalues_number;
  char spectralSubSet[3];
  int max_iterations;
  double tolerance;
  int eignum=0;
  char fn[1024];
  int i,len;

  FILE *out;

  initparam_UserDataCurr_struct (&udc);

  // aux vectors for all mesh elements
  st  = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  M0L = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  M0R = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  X  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  Y  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  G  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  V1 = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  V2 = make_vector_double (udc.N, __FILE__, __FUNCTION__);

  int ret = read_stationary_solution (ST_SOL_FILE, udc.N, G, V1, V2);
  if (ret < 0)
    return 0;

  // it equals to number of non-trivial equations
  eigenvalues_number = 6;

  /* 'LM' -> eigenvalues of largest magnitude. */
  /* 'SM' -> eigenvalues of smallest magnitude.*/
  /* 'LR' -> eigenvalues of largest real part. */
  /* 'SR' -> eigenvalues of smallest real part.*/
  /* 'LI' -> eigenvalues of largest imaginary part. */
  /* 'SI' -> eigenvalues of smallest imaginary part.*/

  strcpy (spectralSubSet, "SM");
  max_iterations = 1000;
  tolerance = 1.e-12;


  // TODO: change amount of memory to correct value
  eigen_values = make_vector_double (2 * eigenvalues_number,
                                     __FILE__, __FUNCTION__);

  // TODO: change amount of memory to correct value
  eigen_functions_A_op = make_vector_double (eigenvalues_number * udc.NA,
                                             __FILE__, __FUNCTION__);

  // TODO: change amount of memory to correct value
  eigen_functions = make_vector_double (eigenvalues_number * udc.N,
                                        __FILE__, __FUNCTION__);


  // TODO: look through this function to check correctness
  eignum = numsds_spectral_problem (eigen_values, eigen_functions_A_op,
                                    udc.NA, eigenvalues_number,
                                    max_iterations, tolerance,
                                    spectralSubSet, A_op, (void *) (&udc));

  // TODO: change to correct values of []
  for (i = 0; i < eignum; i++)
    {
      convert_au_to_u (&eigen_functions[i * udc.N],
                       &eigen_functions_A_op[i * udc.NA],
                       &udc, st);
    }


  // TODO: change to correct values of []
  for (i = 0; i < eignum; i++)
    {
      len = snprintf (fn, sizeof(fn)-1, "./Numres/eigenfun_%02d", i);
      fn[len] = '\0';
      out = fopen (fn,"w");
      if (out == NULL)
        {
          printf ("%s ERR!\n",fn);
          exit (1);
        }

      print_2dfun_double (out, "egenfun",
                          &eigen_functions[i * udc.N],
                          udc.Nx, udc.Ny);
      fclose(out);
    }

  FREE_ARRAY (eigen_values);
  FREE_ARRAY (eigen_functions_A_op);
  FREE_ARRAY (eigen_functions);

  FREE_ARRAY (st);
  FREE_ARRAY (M0L);
  FREE_ARRAY (M0R);

  FREE_ARRAY (X);
  FREE_ARRAY (Y);

  FREE_ARRAY (G);
  FREE_ARRAY (V1);
  FREE_ARRAY (V2);

  return 0;
}
