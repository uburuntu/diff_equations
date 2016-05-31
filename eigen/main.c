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

#define STSOLFILE "stat_sol.txt"

int main (void)
{
  UserDataCurr_struct  udc;

  double *eigen_values;
  double *eigen_functions_A_op;
  double *eigen_functions;

  // aux vectors for mesh parameters
  int *st;
  int *MOL;
  int *MOR;
  double *X;
  double *Y;

  int eigenvalues_number;
  char spectralSubSet[3];
  int max_iterations;
  double tolerance;
  int eignum=0;
  char fn[1024];
  int i,len;
  int readedN;

  FILE *out;

  // TODO: read G, V1, V2 from file stat_sol.txt,
  // which has following syntax:
  // N -- number of vector elements
  // {N elements of G}
  // {N elements of V1}
  // {N elements of V2}
  // Make a check that N == udc->N
  // and you read all 3 * N elements
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FILE *input_stationary_solution;
  double *G;
  double *V1;
  double *V2;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  initparam_UserDataCurr_struct (&udc);

  // aux vectors for all mesh elements
  st  = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  MOL = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  MOR = make_vector_int (udc.N, __FILE__, __FUNCTION__);
  X  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  Y  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  G  = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  V1 = make_vector_double (udc.N, __FILE__, __FUNCTION__);
  V2 = make_vector_double (udc.N, __FILE__, __FUNCTION__);

  input_stationary_solution = fopen(STSOLFILE, "r");
  if(!input_stationary_solution)
  {
      printf("fopen error: cannot open %s", STSOLFILE);
      return -1;
  }
  if(!fscanf (input_stationary_solution, "%d", &readedN))
  {
      printf("fread error: incorrect N by file %s", STSOLFILE);
      fclose(input_stationary_solution);
      return -1;
  }
  if(readedN != udc.N)
  {
      printf("fread error: incorrect N=%d by file %s", readedN, STSOLFILE);
      fclose(input_stationary_solution);
      return -1;
  }
  readedN = 0;
  for(i = 0; i < 3 * udc.N; i++)
  {
      if(!fscanf (input_stationary_solution, "%lf", i%udc.N + (i/udc.N == 0 ? G : (i/udc.N == 1 ? V1 : V2))))
      {
          printf("fread error: incorrect data by %d in line %d by file %s", i%udc.N, i/udc.N, STSOLFILE);
          fclose(input_stationary_solution);
          return -1;
      }
      readedN++;
  }
  fclose(input_stationary_solution);
  printf("Totally readed %d numbers from file %s", readedN, STSOLFILE);

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
                       &udc);
    }


  // TODO: change to correct values of []
  for(i = 0; i < eignum; i++)
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

  // TODO: where are memory free function calls for *eigen_values?
  // Kornev really hates memory...

  if (st) free (st);
  if (MOL) free (MOL);
  if (MOR) free (MOR);
  if (X) free (X);
  if (Y) free (Y);

  return 0;
}
