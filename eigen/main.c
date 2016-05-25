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



int main (void)
{

  UserDataCurr_struct  udc;

  double *eigen_values;
  double *eigen_functions_A_op;
  double *eigen_functions;

  int eigenvalues_number;
  char spectralSubSet[3];
  int max_iterations;
  double tolerance;
  int eignum=0;
  char fn[1024];
  int i,len;

  FILE*out;

  initparam_UserDataCurr_struct (&udc);

  eigenvalues_number=6;

  /* 'LM' -> eigenvalues of largest magnitude. */
  /* 'SM' -> eigenvalues of smallest magnitude.*/
  /* 'LR' -> eigenvalues of largest real part. */
  /* 'SR' -> eigenvalues of smallest real part.*/
  /* 'LI' -> eigenvalues of largest imaginary part. */
  /* 'SI' -> eigenvalues of smallest imaginary part.*/

  strcpy(spectralSubSet,"SM");
  max_iterations = 1000;
  tolerance =1.e-12;


  eigen_values = make_vector_double (2 * eigenvalues_number,
                                     __FILE__, __FUNCTION__);

  // TODO: now we don`t use NA, because we calc for full area
  eigen_functions_A_op = make_vector_double (eigenvalues_number * udc.NA,
                                             __FILE__, __FUNCTION__);

  eigen_functions = make_vector_double (eigenvalues_number * udc.N,
                                        __FILE__, __FUNCTION__);


  // TODO: now we don`t use NA, because we calc for full area
  eignum = numsds_spectral_problem (eigen_values, eigen_functions_A_op,
                                    udc.NA, eigenvalues_number,
                                    max_iterations, tolerance,
                                    spectralSubSet, A_op, (void *) (&udc));
  
  // TODO: now we don`t use NA, because we calc for full area
  for (i = 0; i < eignum; i++)
    {
      convert_au_to_u (&eigen_functions[i * udc.N],
                       &eigen_functions_A_op[i * udc.NA],
                       &udc);
    }


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

  return 0;

}
