#include <stdio.h>
#include <math.h>
#include "laspack/getopts.h"
#include "laspack/vector.h"
#include "laspack/errhandl.h"
#include "laspack/qmatrix.h"
#include "laspack/itersolv.h"
#include "laspack/rtc.h"
#include "laspack/operats.h"
#include "laspack/version.h"
#include "laspack/copyrght.h"
#include "func.h"
#include "gnuplot.h"

#define DEBUG_VARIANT_I   0
#define DEBUG_VARIANT_II  0
#define DEBUG_VARIANT_III 0
#define EIG_USAGE_TIME    5

int  Sxema (double *G, double *V1, double *V2,
            double *G_prev, double *V1_prev, double *V2_prev,
            const double *G_eigen, const double *V1_eigen, const double *V2_eigen,
            const double *G_stat, const double *V1_stat, const double *V2_stat,
            int *st, double *X, double *Y, int *M0L,
            int *M0R, P_she *p_s, P_dif *p_d)
{
  // Arguments
  int k;
  int N, Dim;
  double norm;
  double hx, hy, tau, mu, p_ro;

  N   = p_s->N;
  Dim = p_s->Dim;
  hx  = p_s->h_x;
  hy  = p_s->h_y;
  tau = p_s->tau;

  mu   = p_d->mu;
  p_ro = p_d->p_ro;

  // Local variables
  int nn, m, mm;
  int mmg0R, mmgL0, mmv1L0, mmv2L0, mmgR0, mmv1R0, mmv2R0, mmg0L, mmv10L, mmv20L;
  int mmv10R, mmv20R;

  char char_A[2] = "A";
  char char_B[2] = "B";
  char char_D[2] = "D";

  char plotname [50];
  char texname [50];
  int nameiter = 1;

  double tt, xx, yy;
  double tmp, tmp1;
  double g00;
  double v1L0, v100, v1R0, v10L, v10R, v1LL, v1LR, v1RL, v1RR;
  double v2L0, v200, v2R0, v20L, v20R, v2LL, v2LR, v2RL, v2RR;
  double MUM, MU43x, MU43y, MUx, MUy, MUv1, MUv2;

  double thx, thy, thx_05, thy_05;
  double thxy_1_12, thxp_05, thyp_05;
  double thxx_4_3, thyy_4_3, thxx, thyy;

  thx       = tau / hx;
  thy       = tau / hy;
  thx_05    = 0.5 * thx;
  thy_05    = 0.5 * thy;
  thxy_1_12 = tau / (12 * hx * hy);
  thxp_05   = 0.5 * tau * p_ro / (hx);
  thyp_05   = 0.5 * tau * p_ro / (hy);
  thxx_4_3  = 4. * tau / (3. * hx * hx);
  thyy_4_3  = 4. * tau / (3. * hy * hy);
  thxx      = tau / (hx * hx);
  thyy      = tau / (hy * hy);

  // A -- sparse matrix of the system, D -- solution vector, B -- rhs vector
  QMatrix A;
  Vector D, B;
  Q_Constr (&A, char_A, 3 * Dim, False, Rowws, Normal, True);
  V_Constr (&B, char_B, 3 * Dim, Normal, True);
  V_Constr (&D, char_D, 3 * Dim, Normal, True);
  SetRTCAccuracy (1e-9);

  // Initial values
  tt = 0.;
  mm = 1;

  for (m = 0; m < Dim; m++)
    {
      xx = X[m];
      yy = Y[m];

      tmp = gg (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      G[m] = tmp;

      if (EIG_FUNC_INIT && G_eigen)
        {
          if (DEBUG_VARIANT_III)
            {
              G[m] += 0;
            }
          else
            {
              if (DEBUG_VARIANT_I)
                {
                  G[m] += G_eigen[m];
                }
              else if (DEBUG_VARIANT_II)
                {
                  G[m] -= G_eigen[m];
                }
            }
        }

      if (STAT_SOL_SRCH && G_prev)
        {
          G_prev[m] = tmp;
        }

      mm++;

      tmp = u1 (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      V1[m] = tmp;

      if (EIG_FUNC_INIT && G_eigen)
        {
          if (DEBUG_VARIANT_III)
            {
              V1[m] += 0;
            }
          else
            {
              if (DEBUG_VARIANT_I)
                {
                  V1[m] += V1_eigen[m];
                }
              else if (DEBUG_VARIANT_II)
                {
                  V1[m] -= V1_eigen[m];
                }
            }
        }

      if (STAT_SOL_SRCH && V1_prev)
        {
          V1_prev[m] = tmp;
        }

      mm++;

      tmp = u2 (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      V2[m] = tmp;

      if (EIG_FUNC_INIT && G_eigen)
        {
          if (DEBUG_VARIANT_III)
            {
              V2[m] += 0;
            }
          else
            {
              if (DEBUG_VARIANT_I)
                {
                  V2[m] += V2_eigen[m];
                }
              else if (DEBUG_VARIANT_II)
                {
                  V2[m] -= V2_eigen[m];
                }
            }
        }

      if (STAT_SOL_SRCH && V2_prev)
        {
          V2_prev[m] = tmp;
        }

      mm++;
    }


  for (nn = 1; nn <= N; nn++)
    {
      tt = nn * tau;

      //-///////////////// MUM //////////////////////////////////
      MUM = 0.0;

      for (m = 0; m < Dim; m++)
        {
          if (st[m] == 0)
            {
              tmp = exp (-G[m]);

              if (MUM < tmp)
                {
                  MUM = tmp;
                }
            }
        }

      MUM *= mu;
      MU43x = MUM * tau * 4. / (3. * hx * hx);
      MUx   = MUM * tau / (hx * hx);
      MU43y = MUM * tau * 4. / (3. * hy * hy);
      MUy   = MUM * tau / (hy * hy);
      MUv1  = (8. / (3. * hx * hx) + 2. / (hy * hy)) * MUM * tau;
      MUv2  = (8. / (3. * hy * hy) + 2. / (hx * hx)) * MUM * tau;
      //-////////////////////////////////////////////////////////

      mm = 1;

      for (m = 0; m < Dim; m++)
        {
          // setting A and B elements for every node
          xx = X[m];
          yy = Y[m];

          switch (st[m])
            {
              case 0:
                {
                  mmgL0  = mm - 3;
                  mmv1L0 = mm - 2;
                  mmv2L0 = mm - 1;
                  mmgR0  = mm + 3;
                  mmv1R0 = mm + 4;
                  mmv2R0 = mm + 5;
                  mmg0L  = 3 * M0L[m] + 1;
                  mmv10L = mmg0L + 1;
                  mmv20L = mmv10L + 1;
                  mmg0R  = 3 * M0R[m] + 1;
                  mmv10R = mmg0R + 1;
                  mmv20R = mmv10R + 1;
                  g00  = G[m];
                  v100 = V1[m];
                  v1L0 = V1[m - 1];
                  v1R0 = V1[m + 1];
                  v10L = V1[M0L[m]];
                  v10R = V1[M0R[m]];
                  v1LL = V1[M0L[m] - 1];
                  v1RL = V1[M0L[m] + 1];
                  v1LR = V1[M0R[m] - 1];
                  v1RR = V1[M0R[m] + 1];
                  v200 = V2[m];
                  v2L0 = V2[m - 1];
                  v2R0 = V2[m + 1];
                  v20L = V2[M0L[m]];
                  v20R = V2[M0R[m]];
                  v2LL = V2[M0L[m] - 1];
                  v2RL = V2[M0L[m] + 1];
                  v2LR = V2[M0R[m] - 1];
                  v2RR = V2[M0R[m] + 1];

                  // g(mx, my)--------------------------------------------
                  Q_SetLen (&A, mm, 9);
                  tmp = 1. + thx * fabs (v100) + thy * fabs (v200);   // +
                  Q_SetEntry (&A, mm, 0, mm, tmp);                    // +
                  tmp = thx_05 * (v100 - fabs (v100));                // +
                  Q_SetEntry (&A, mm, 1, mmgR0, tmp);                 // +
                  tmp = -thx_05 * (v100 + fabs (v100));               // +
                  Q_SetEntry (&A, mm, 2, mmgL0, tmp);                 // +
                  tmp = thy_05 * (v200 - fabs (v200));                // +
                  Q_SetEntry (&A, mm, 3, mmg0R, tmp);                 // +
                  tmp = -thy_05 * (v200 + fabs (v200));               // +
                  Q_SetEntry (&A, mm, 4, mmg0L, tmp);                 // +
                  Q_SetEntry (&A, mm, 5, mmv1R0, thx_05);             // +
                  Q_SetEntry (&A, mm, 6, mmv20R, thy_05);             // +
                  Q_SetEntry (&A, mm, 7, mmv1L0, -thx_05);            // +
                  Q_SetEntry (&A, mm, 8, mmv20L, -thy_05);            // +
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy)); // +
                  mm++;

                  // v1(mx, my)-----------------------------------------
                  Q_SetLen (&A, mm, 7);
                  tmp = 1. + thx * fabs (v100) + thy * fabs (v200) + MUv1;    // +
                  Q_SetEntry (&A, mm, 0, mm, tmp);                            // +
                  tmp = thx_05 * (v100 - fabs (v100)) - MU43x;                // +
                  Q_SetEntry (&A, mm, 1, mmv1R0, tmp);                        // +
                  tmp = thy_05 * (v200 - fabs (v200)) - MUy;                  // +
                  Q_SetEntry (&A, mm, 2, mmv10R, tmp);                        // +
                  tmp = thxp_05;                                              // +
                  Q_SetEntry (&A, mm, 3, mmgL0, -tmp);                        // +
                  Q_SetEntry (&A, mm, 4, mmgR0, tmp);                         // +
                  tmp = -thx_05 * (v100 + fabs (v100)) - MU43x;               // +
                  Q_SetEntry (&A, mm, 5, mmv1L0, tmp);                        // +
                  tmp = -thy_05 * (v200 + fabs (v200)) - MUy;                 // +
                  Q_SetEntry (&A, mm, 6, mmv10L, tmp);                        // +
                  tmp1 = mu * exp (-g00);
                  tmp = v100 + (tmp1 - MUM)
                        * (thxx_4_3 * (v1R0 - 2. * v100 + v1L0) + thyy * (v10R - 2. * v100 + v10L))
                        + tmp1 * thxy_1_12 * (v2RR - v2RL - v2LR + v2LL) + tau * Func_v1 (tt, xx, yy, p_ro, mu);
                  V_SetCmp (&B, mm, tmp);
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 7);
                  tmp = 1. + thy * fabs (v200) + thx * fabs (v100) + MUv2;
                  Q_SetEntry (&A, mm, 0, mm, tmp);
                  tmp = thx_05 * (v100 - fabs (v100)) - MUx;
                  Q_SetEntry (&A, mm, 1, mmv2R0, tmp);
                  tmp = thy_05 * (v200 - fabs (v200)) - MU43y;
                  Q_SetEntry (&A, mm, 2, mmv20R, tmp);
                  tmp = thyp_05;
                  Q_SetEntry (&A, mm, 3, mmg0L, -tmp);
                  Q_SetEntry (&A, mm, 4, mmg0R, tmp);
                  tmp = -thx_05 * (v100 + fabs (v100)) - MUx;
                  Q_SetEntry (&A, mm, 5, mmv2L0, tmp );
                  tmp = -thy_05 * (v200 + fabs (v200)) - MU43y;
                  Q_SetEntry (&A, mm, 6, mmv20L, tmp );
                  tmp1 = mu * exp (-g00);
                  tmp = v200 + (tmp1 - MUM)
                        * (thyy_4_3 * (v20R - 2. * v200 + v20L) + thxx * (v2R0 - 2. * v200 + v2L0))
                        + tmp1 * thxy_1_12 * (v1RR - v1RL - v1LR + v1LL) + tau * Func_v2 (tt, xx, yy, p_ro, mu);
                  V_SetCmp (&B, mm, tmp);
                  mm++;

                  break;
                }

              case 1:
                {
                  // On boundary where the velocity vector is directed inside
                  // the value of g is known according to your test conditions.
                  g00 = (NO_SMOOTH && is_equal (xx, 0.0) ? RHO_G : G[m]);
                  mmv1R0 = mm + 4;

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 3);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  Q_SetEntry (&A, mm, 1, mm + 1, -thx);
                  Q_SetEntry (&A, mm, 2, mmv1R0, thx);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 2:
                {
                  g00 = G[m];
                  v100 = (NO_SMOOTH && SQUARE ? V1[m - 1] : u1 (tt, xx, yy));
                  v200 = u2 (tt, xx, yy);
                  mmv1L0 = mm - 2;

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 3);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  Q_SetEntry (&A, mm, 1, mm + 1, thx);
                  Q_SetEntry (&A, mm, 2, mmv1L0, -thx);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, v100);
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 3:
                {
                  g00 = G[m];
                  v100 = u1 (tt, xx, yy);
                  v200 = (NO_SMOOTH && !SQUARE && is_equal (yy, 0.) ? V2[M0R[m]] : u2 (tt, xx, yy));
                  mmv20R = (3 * M0R[m] + 1) + 2;

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 3);
                  Q_SetEntry (&A, mm, 0, mm, 1);
                  Q_SetEntry (&A, mm, 1, mm + 2, -thy);
                  Q_SetEntry (&A, mm, 2, mmv20R, thy);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, v100);
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, v200);
                  mm++;

                  break;
                }

              case 4:
                {
                  g00 = G[m];
                  mmv20L  = (3 * M0L[m] + 1) + 2;

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 3);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  Q_SetEntry (&A, mm, 1, mm + 2, thy);
                  Q_SetEntry (&A, mm, 2, mmv20L, -thy);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 5:         //   (x,y) in (0,0)
                {
                  g00 = G[m];

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 6:         //   (x,y) in (X,0)
                {
                  g00 = G[m];

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 7:         //   (x,y) in (0,Y)
                {
                  g00 = G[m];

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              case 8:      //   (x,y) in (X,Y)
                {
                  g00 = G[m];

                  // g(mx,my)--------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, g00 + tau * Func_g (tt, xx, yy));
                  mm++;

                  // v1(mx,my)-----------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u1 (tt, xx, yy));
                  mm++;

                  // v2(mx,my)------------------------------------------------------
                  Q_SetLen (&A, mm, 1);
                  Q_SetEntry (&A, mm, 0, mm, 1.);
                  V_SetCmp (&B, mm, u2 (tt, xx, yy));
                  mm++;

                  break;
                }

              default:
                {
                  printf ("default case section... %d \n", m);
                }
            }
        }

      // Solver
      static const int max_iters = 2000;
      static const double precond_omega = 1.;

      switch (4)
        {
          case 1:
            {
              BiCGIter (&A, &D, &B, max_iters, JacobiPrecond, precond_omega);
              break;
            }

          case 2:
            {
              BiCGIter (&A, &D, &B, max_iters, SSORPrecond, precond_omega);
              break;
            }

          case 3:
            {
              BiCGIter (&A, &D, &B, max_iters, NULL, precond_omega);
              break;
            }

          case 4:
            {
              BiCGSTABIter (&A, &D, &B, max_iters, NULL, precond_omega);
              break;
            }

          default:
            {
              BiCGIter (&A, &D, &B, max_iters, NULL, precond_omega);
              break;
            }
        }

      // Copy solutiong to G, V1, V2 arrays
      mm = 1;

      for (m = 0; m < Dim; m++)
        {
          xx = X[m];
          yy = Y[m];

          G[m] = V_GetCmp (&D, mm);
          mm++;
          //printf ("%lf ", fabs (G[m] - gg (tt, xx, yy)));

          V1[m] = V_GetCmp (&D, mm);
          mm++;
          //printf ("%lf ", fabs (V1[m] - u1 (tt, xx, yy)));

          V2[m] = V_GetCmp (&D, mm);
          mm++;
          //printf ("%lf, ", fabs (V2[m] - u2 (tt, xx, yy)));

          //printf ("%lf %lf %lf \n", tt, xx, yy);
        }


      if (nn == nameiter || nn == N)
        {
          nameiter += N / DIVISOR;
          //printf("%d, %d \n", nn, N);
          print_plot (plot_name (plotname, tau, hx, hy, nn), X, Y, G, V1, V2, Dim, tt);
          make_graph (tex_name (texname, tau, hx, hy, nn), plotname, hx, hy, tau, tt);
        }

      if (STAT_SOL_SRCH && G_prev && V1_prev && V2_prev)
        {
          norm = calc_sol_residual_norm (Dim, G, V1, V2, G_prev, V1_prev, V2_prev);

          if (norm < STAT_SOL_EPS)
            {
              printf ("Stationary solution has been found at T = %d. \n", nn);
              printf ("Accuracy = %E. \n", norm);
              Q_Destr (&A);
              V_Destr (&D);
              V_Destr (&B);
              return 1;
            }
          else if (nn == 1 || nn % 10 == 0)
            {
              printf ("t = %3.d, norm = %E \n", nn, norm);
            }

          init_prev_with_curr (Dim, G, V1, V2, G_prev, V1_prev, V2_prev);
        }

      if (EIG_FUNC_INIT && G_stat && V1_stat && V2_stat)
        {
          norm = calc_sol_residual_norm (Dim, G, V1, V2, G_stat, V1_stat, V2_stat);

          if (norm < STAT_SOL_EPS)
            {
              printf ("Stationary solution has been found at T = %d. \n", nn);
              printf ("Accuracy = %E. \n", norm);
              Q_Destr (&A);
              V_Destr (&D);
              V_Destr (&B);
              return 2;
            }
          else if (nn == 1 || nn % 10 == 0)
            {
              printf ("t = %3.d, norm = %E \n", nn, norm);
            }

          if (nn == EIG_USAGE_TIME)
            {
              if (DEBUG_VARIANT_III)
                {
                  if (DEBUG_VARIANT_I)
                    {
                      for (k = 0; k < Dim; k++)
                        {
                          G[m]  += G_eigen[k];
                          V1[m] += V1_eigen[k];
                          V2[m] += V2_eigen[k];
                        }
                    }
                  else if (DEBUG_VARIANT_I)
                    {
                      for (k = 0; k < Dim; k++)
                        {
                          G[m]  += G_eigen[k];
                          V1[m] += V1_eigen[k];
                          V2[m] += V2_eigen[k];
                        }
                    }
                }
            }
        }
    }

  Q_Destr (&A);
  V_Destr (&D);
  V_Destr (&B);

  if (STAT_SOL_SRCH && G_prev && V1_prev && V2_prev)
    {
      printf ("Stationary solution has not been found at T = %d. \n", nn);
      printf ("Accuracy = %E. \n", STAT_SOL_EPS);
      return -1;
    }

  return 0;
}

