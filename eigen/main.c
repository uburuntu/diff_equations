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

#define ST_SOL_FILE "stat_sol.txt"

int main (void)
{
  user_data ud;

  int n_found_eigen_values;

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

  init_user_data (&ud);

  // aux vectors for all mesh elements
  st  = make_vector_int (ud.N, __FILE__, __FUNCTION__);
  M0L = make_vector_int (ud.N, __FILE__, __FUNCTION__);
  M0R = make_vector_int (ud.N, __FILE__, __FUNCTION__);
  X  = make_vector_double (ud.N, __FILE__, __FUNCTION__);
  Y  = make_vector_double (ud.N, __FILE__, __FUNCTION__);
  G  = make_vector_double (ud.N, __FILE__, __FUNCTION__);
  V1 = make_vector_double (ud.N, __FILE__, __FUNCTION__);
  V2 = make_vector_double (ud.N, __FILE__, __FUNCTION__);

  int ret = read_stationary_solution (ST_SOL_FILE, ud.N, G, V1, V2);

  if (ret < 0)
    {
      return 0;
    }

  /* 'LM' -> eigenvalues of largest magnitude. */
  /* 'SM' -> eigenvalues of smallest magnitude.*/
  /* 'LR' -> eigenvalues of largest real part. */
  /* 'SR' -> eigenvalues of smallest real part.*/
  /* 'LI' -> eigenvalues of largest imaginary part. */
  /* 'SI' -> eigenvalues of smallest imaginary part.*/
  const char *spectralSubSet = "LR";
  /* BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x*/
  /* BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x*/
  const char *bmat = "I";

  const int n_eigen_values = 6;
  const int max_iterations = 1000;
  const double tolerance = 1.e-12;

  // eigen_values[2 * i + 0] = LAMBDA_I_REAL_PART
  // eigen_values[2 * i + 1] = LAMBDA_I_IMAG_PART
  eigen_values = make_vector_double (2 * n_eigen_values, __FILE__, __FUNCTION__);

  eigen_functions_A_op = make_vector_double (n_eigen_values * ud.NA, __FILE__, __FUNCTION__);

  eigen_functions = make_vector_double (n_eigen_values * (3 * ud.N), __FILE__, __FUNCTION__);

  calc_mesh_params (st, X, Y, M0L, M0R, (void *) (&ud));

  // TODO: look through this function to check correctness
  // i have checked, but i`m not sure it`s correct...
  n_found_eigen_values = find_eigen_values (eigen_values, eigen_functions_A_op,
                         ud.NA, n_eigen_values,
                         max_iterations, tolerance,
                         spectralSubSet, bmat, (void *) (&ud),
                         G, V1, V2, st, M0L, M0R, X, Y);

  for (int i = 0; i < n_found_eigen_values; i++)
    {
      convert_au_to_u (&eigen_functions[i * (3 * ud.N)],
                       &eigen_functions_A_op[i * ud.NA],
                       &ud, st, X, Y);
    }

  for (int i = 0; i < n_found_eigen_values; i++)
    {
      char fn[1024];

      int len = snprintf (fn, sizeof (fn) - 1, "./results/eigenfun_%02d.txt", i);
      fn[len] = 0;

      FILE *out = fopen (fn, "w");

      if (out == NULL)
        {
          printf ("%s open error!\n", fn);
          exit (1);
        }

      print_2dfun_double (out, &eigen_functions[i * (3 * ud.N)], 3 * ud.N);
      fclose (out);
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
