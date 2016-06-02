#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "tabtex.h"
#include "func.h"
#include "gnuplot.h"

/*
 Macros of Daenerys Stormborn Targaryen,
 Mother of Dragons,
 Breaker of Chains,
 Queen of the Andels and the First Men,
 the rightful ruler of Westeros
*/
#define FREE_ALL() {                                               \
  FREE_ARRAY(nc_g);     FREE_ARRAY(nl2_g);   FREE_ARRAY(nc_v1);    \
  FREE_ARRAY(nc_v2);    FREE_ARRAY(nl2_v1);  FREE_ARRAY(nl2_v2);   \
  FREE_ARRAY(time);     FREE_ARRAY(tauit);   FREE_ARRAY(st);       \
  FREE_ARRAY(M0L);      FREE_ARRAY(M0R);     FREE_ARRAY(X);        \
  FREE_ARRAY(Y);        FREE_ARRAY(G);       FREE_ARRAY(V1);       \
  FREE_ARRAY(V2);       FREE_ARRAY(G_prev);  FREE_ARRAY(V1_prev);  \
  FREE_ARRAY(V2_prev);  FREE_ARRAY(G_eigen); FREE_ARRAY(V1_eigen); \
  FREE_ARRAY(V2_eigen); FREE_ARRAY(G_stat);  FREE_ARRAY(V1_stat);  \
  FREE_ARRAY(V2_stat);}

#define CLOSE_ALL() {       \
  CLOSE_FILE(fout);         \
  CLOSE_FILE(stat_sol_out); \
  CLOSE_FILE(stat_sol_in);  \
  CLOSE_FILE(eig_func_in);}

int main ()
{
  int ret;
  int i;
  FILE *stat_sol_out = NULL;
  FILE *stat_sol_in  = NULL;
  FILE *eig_func_in  = NULL;

  P_dif p_d;
  param_dif (&p_d);

  int it_t, it_t_max, it_sp, it_sp_max, n_ver, it;
  it_t_max = 0;
  it_sp_max = 0;
  it_t = 0;
  it_sp = 0;

  n_ver = (it_t_max + 1) * (it_sp_max + 1);

  double *nc_g, *nl2_g, *nc_v1, *nc_v2, *nl2_v1, *nl2_v2;
  double *time, *tauit;

  if (!NEW_INIT)
    {
      nc_g = (double *) malloc ((n_ver) * sizeof (double));
      nc_v1 = (double *) malloc ((n_ver) * sizeof (double));
      nc_v2 = (double *) malloc ((n_ver) * sizeof (double));
      nl2_g = (double *) malloc ((n_ver) * sizeof (double));
      nl2_v1 = (double *) malloc ((n_ver) * sizeof (double));
      nl2_v2 = (double *) malloc ((n_ver) * sizeof (double));
    }
  else
    {
      nc_g   = NULL;
      nc_v1  = NULL;
      nc_v2  = NULL;
      nl2_g  = NULL;
      nl2_v1 = NULL;
      nl2_v2 = NULL;
    }

  time = (double *) malloc ((n_ver) * sizeof (double));
  tauit = (double *) malloc ((it_t_max + 1) * sizeof (double));

  P_she p_s;

  int *st, *M0L, *M0R;
  double *X, *Y, *G, *V1, *V2;
  double *G_prev, *V1_prev, *V2_prev;
  double *G_eigen, *V1_eigen, *V2_eigen;
  double *G_stat, *V1_stat, *V2_stat;

  st = NULL;
  M0L = NULL;
  M0R = NULL;
  X = NULL;
  Y = NULL;
  G = NULL;
  V1 = NULL;
  V2 = NULL;
  G_prev = NULL;
  V1_prev = NULL;
  V2_prev = NULL;
  G_eigen = NULL;
  V1_eigen = NULL;
  V2_eigen = NULL;
  G_stat = NULL;
  V1_stat = NULL;
  V2_stat = NULL;

  clock_t BegClock, EndClock;
  //-------------------------------------------------
  FILE *fout;

  if (!NEW_INIT)
    {
      fout = fopen (OUTTEX_SMOOTH, "w");
    }
  else
    {
      fout = fopen (OUTTEX_ABRUPT, "w");
    }

  if (fout)
    {
      printhead (fout);
      fprintf (fout, "\\section* {Графики:}  \n");
      fclose (fout);
    }
  else
    {
      printf ("Cannot open OUTTEX file\n");
      FREE_ALL();
      return -1;
    }

  if (EIG_FUNC_INIT)
    {
      eig_func_in = fopen ("./eigen/results/eigenfun_00.txt", "r");
      if (eig_func_in == NULL)
        {
          printf ("Cannot open ./eigen/results/eigenfun_00.txt");
          printf (" in EIG_FUNC_INIT mode.\n");
          FREE_ALL();
          return -1;
        }
      stat_sol_in = fopen ("./eigen/stat_sol.txt", "r");
      if (stat_sol_in == NULL)
        {
          printf ("Cannot open ./eigen/stat_sol.txt");
          printf (" in EIG_FUNC_INIT mode.\n");
          FREE_ALL();
          return -1;
        }
    }


  //-------------------------------------------------

  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          BegClock = clock ();

          // Define area params
          param_she_step (&p_s, &p_d, it_t, it_sp);
          printf ("N = %3.d, M1 = %3.d, M2 = %3.d, Dim = %6.d. \n", p_s.N, p_s.M_x, p_s.M_y, p_s.Dim);
          fflush (stdout);

          if (it_sp == 0)
            {
              tauit[it_t] = p_s.tau;
            }

          X = (double *) malloc ((p_s.Dim) * sizeof (double)); // x-coord array of nodes
          Y = (double *) malloc ((p_s.Dim) * sizeof (double)); // y-coord array of nodes

          G  = (double *) malloc ((p_s.Dim) * sizeof (double)); // press array of nodes
          V1 = (double *) malloc ((p_s.Dim) * sizeof (double)); // v1 array of nodes
          V2 = (double *) malloc ((p_s.Dim) * sizeof (double)); // v2 array of nodes

          if (STAT_SOL_SRCH && it_sp == it_t && it_sp == 0)
            {
              G_prev  = (double *) malloc ((p_s.Dim) * sizeof (double)); // press array of nodes
              V1_prev = (double *) malloc ((p_s.Dim) * sizeof (double)); // v1 array of nodes
              V2_prev = (double *) malloc ((p_s.Dim) * sizeof (double)); // v2 array of nodes
            }

          if (EIG_FUNC_INIT && it_sp == it_t && it_sp == 0)
            {
              ret = 0;
              G_eigen  = (double *) malloc ((p_s.Dim) * sizeof (double)); // press array of nodes
              V1_eigen = (double *) malloc ((p_s.Dim) * sizeof (double)); // v1 array of nodes
              V2_eigen = (double *) malloc ((p_s.Dim) * sizeof (double)); // v2 array of nodes
              G_stat   = (double *) malloc ((p_s.Dim) * sizeof (double)); // press array of nodes
              V1_stat  = (double *) malloc ((p_s.Dim) * sizeof (double)); // v1 array of nodes
              V2_stat  = (double *) malloc ((p_s.Dim) * sizeof (double)); // v2 array of nodes
              for (i = 0; i < p_s.Dim; i++)
                {
                  if (!fscanf (eig_func_in, "%lf ", G_eigen + i))
                    {
                      ret = -999;
                      break;
                    }
                  if (!fscanf (eig_func_in, "%lf ", V1_eigen + i))
                    {
                      ret = -999;
                      break;
                    }
                  if (!fscanf (eig_func_in, "%lf ", V2_eigen + i))
                    {
                      ret = -999;
                      break;
                    }
                }
              if (ret < 0)
                {
                  printf ("fread error: incorrect G_eigen, V1_eigen, V2_eigen filling from %s\n",
                          "./eigen/results/eigenfun_00.txt");
                  FREE_ALL();
                  CLOSE_ALL();
                  return -1;
                }
              CLOSE_FILE(eig_func_in);

              G_stat   = (double *) malloc ((p_s.Dim) * sizeof (double)); // press array of nodes
              V1_stat  = (double *) malloc ((p_s.Dim) * sizeof (double)); // v1 array of nodes
              V2_stat  = (double *) malloc ((p_s.Dim) * sizeof (double)); // v2 array of nodes

              if (!fscanf (eig_func_in, "%lf ", &i) || i != p_s.Dim)
                {
                  ret = -999;
                }

              for (i = 0; i < p_s.Dim; i++)
                {
                  if (!fscanf (eig_func_in, "%lf ", G_stat + i))
                    {
                      ret = -999;
                      break;
                    }
                  if (!fscanf (eig_func_in, "%lf ", V1_stat + i))
                    {
                      ret = -999;
                      break;
                    }
                  if (!fscanf (eig_func_in, "%lf ", V2_stat + i))
                    {
                      ret = -999;
                      break;
                    }
                }
              if (ret < 0)
                {
                  printf ("fread error: incorrect G_stat, V1_stat, V2_stat filling from %s\n",
                          "./eigen/results/eigenfun_00.txt \n");
                  FREE_ALL();
                  CLOSE_ALL();
                  return -1;
                }
              CLOSE_FILE(stat_sol_in);
            }

          st  = (int *) malloc ((p_s.Dim) * sizeof (int));     // status of nodes
          M0L = (int *) malloc ((p_s.Dim) * sizeof (int));
          M0R = (int *) malloc ((p_s.Dim) * sizeof (int));

          // Define properties of nodes
          Setka (st, X, Y, M0L, M0R, &p_s, &p_d);
          // Run calculations
          ret = Sxema (G, V1, V2,
                       G_prev, V1_prev, V2_prev,
                       G_eigen, V1_eigen, V2_eigen,
                       G_stat, V1_stat, V2_stat,
                       st, X, Y, M0L, M0R, &p_s, &p_d);

          EndClock = clock();
          time[it] = (double) (EndClock - BegClock) / CLOCKS_PER_SEC;
          printf ("Elapsed time: %.2f sec.\n", time[it]);

          if (!NEW_INIT)
            {
              nc_g[it] = Norm_c (G, p_s.Dim, X, Y, p_d.Segm_T, gg);
              nl2_g[it] = Norm_l2 (G, p_s.Dim, X, Y, p_d.Segm_T, gg);
              nc_v1[it] = Norm_c (V1, p_s.Dim, X, Y, p_d.Segm_T, u1);
              nl2_v1[it] = Norm_l2 (V1, p_s.Dim, X, Y, p_d.Segm_T, u1);
              nc_v2[it] = Norm_c (V2, p_s.Dim, X, Y, p_d.Segm_T, u2);
              nl2_v2[it] = Norm_l2 (V2, p_s.Dim, X, Y, p_d.Segm_T, u2);

              // debug
              printf (" %lf %lf %lf \n", nl2_g[it], nl2_v1[it], nl2_v2[it]);
            }

          if (STAT_SOL_SRCH && it_sp == it_t && it_sp == 0 &&
              ret == 1 /* stat_sol was found in SRCH mode*/)
            {
              stat_sol_out = fopen ("./eigen/stat_sol.txt", "w");
              fprintf (stat_sol_out, "%d \n", p_s.Dim);

              for (i = 0; i < p_s.Dim; i++)
                {
                  fprintf (stat_sol_out, "%lf ", G[i]);
                }

              fprintf (stat_sol_out, "\n");

              for (i = 0; i < p_s.Dim; i++)
                {
                  fprintf (stat_sol_out, "%lf ", V1[i]);
                }

              fprintf (stat_sol_out, "\n");

              for (i = 0; i < p_s.Dim; i++)
                {
                  fprintf (stat_sol_out, "%lf ", V2[i]);
                }

              fprintf (stat_sol_out, "\n");

              FREE_ALL();
              CLOSE_ALL();
              return 0;
            }

          if (EIG_FUNC_INIT && it_sp == it_t && it_sp == 0 &&
              ret == 2 /* stat_sol was found in EIG mode*/)
            {
              FREE_ALL();
              CLOSE_ALL();
              return 0;
            }

          it++;

          free (X);
          free (Y);
          free (G);
          free (V1);
          free (V2);
          free (st);
          free (M0L);
          free (M0R);
        }
    }

  //-------------------------------------------------
  if (!NEW_INIT)
    {
      fout = fopen (OUTTEX_SMOOTH, "a+");
    }
  else
    {
      fout = fopen (OUTTEX_ABRUPT, "a+");
    }

  if (fout)
    {
      printtail (fout);
    }
  else
    {
      printf ("Can't open file OUTTEX\n");
      return -1;
    }

  fclose (fout);
  //-------------------------------------------------

  if (!NEW_INIT)
    {
      tabtex_nc_g (it_t_max, it_sp_max, nc_g, tauit, p_d.p_ro, p_d.mu);
      tabtex_nc_v1 (it_t_max, it_sp_max, nc_v1, tauit, p_d.p_ro, p_d.mu);
      tabtex_nc_v2 (it_t_max, it_sp_max, nc_v2, tauit, p_d.p_ro, p_d.mu);
      tabtex_nl2_g (it_t_max, it_sp_max, nl2_g, tauit, p_d.p_ro, p_d.mu);
      tabtex_nl2_v1 (it_t_max, it_sp_max, nl2_v1, tauit, p_d.p_ro, p_d.mu);
      tabtex_nl2_v2 (it_t_max, it_sp_max, nl2_v2, tauit, p_d.p_ro, p_d.mu);
      tabtex_time (it_t_max, it_sp_max, time, tauit, p_d.p_ro, p_d.mu);

      make_tabletex();
    }

  FREE_ALL();
  CLOSE_ALL();
  return 0;
}

