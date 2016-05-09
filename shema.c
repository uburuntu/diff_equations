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

void Sxema (double *G, double *V1, double *V2, int *st, double *X, double *Y, int *M0L, int *M0R, P_she *p_s, P_dif *p_d)
{
  // local variables ///////////////////////////////////////////////////////////////////////////////
  int M1, M2, N, Dim;
  FIX_UNUSED (M1); FIX_UNUSED (M2);

  double hx, hy, tau, mu, p_ro;
  M1 = p_s->M_x;
  M2 = p_s->M_y;
  N = p_s->N;
  Dim = p_s->Dim;
  hx = p_s->h_x;
  hy = p_s->h_y;
  tau = p_s->tau;
  mu = p_d->mu;
  p_ro = p_d->p_ro;

  int nn, m, mm, mx, my;
  FIX_UNUSED (mx); FIX_UNUSED (my);

  int mmg0R, mmgL0, mmv1L0, mmv2L0, mmgR0, mmv1R0, mmv2R0, mmg0L, mmv10L, mmv20L;
  int mmv10R, mmv20R;

  double tt, xx, yy;
  double tmp, tmp1;

  double gL0, g00, gR0, g0L, g0R;
  FIX_UNUSED (gL0); FIX_UNUSED (gR0); FIX_UNUSED (g0L); FIX_UNUSED (g0R);

  double v1L0, v100, v1R0, v10L, v10R, v1LL, v1LR, v1RL, v1RR;
  double v2L0, v200, v2R0, v20L, v20R, v2LL, v2LR, v2RL, v2RR;

  char char_A[2] = "A";
  char char_B[2] = "B";
  char char_D[2] = "D";

  char plotname [50];
  char texname [50];

  int nameiter = 1;


  // local variable /////////////////////////////////////////////////////////////////////////////////

  double thx, thy, thx_05, thy_05, thx_2, thy_2, thx_4, thy_4, tau_2, tau_4, tau_6, thx_3_2, thy_3_2;
  double thxx_6, thxx_8, thyy_6, thyy_8, thxy_1_12, thxp_05, thyp_05;
  double thxx_4_3, thyy_4_3, thxx, thyy;

  thx = tau / hx;
  thy = tau / hy;
  thx_05 = 0.5 * thx;
  thy_05 = 0.5 * thy;
  thx_2 = 2 * thx; FIX_UNUSED (thx_2);
  thy_2 = 2 * thy; FIX_UNUSED (thy_2);
  thx_4 = 4 * thx; FIX_UNUSED (thx_4);
  thy_4 = 4 * thy; FIX_UNUSED (thy_4);
  thx_3_2 = 1.5 * thx; FIX_UNUSED (thx_3_2);
  thy_3_2 = 1.5 * thy; FIX_UNUSED (thy_3_2);
  tau_2 = 2 * tau; FIX_UNUSED (tau_2);
  tau_4 = 4 * tau; FIX_UNUSED (tau_4);
  tau_6 = 6 * tau; FIX_UNUSED (tau_6);
  thxx_8 = tau * 8. / (hx * hx); FIX_UNUSED (thxx_8);
  thxx_6 = 6 * tau / (hx * hx); FIX_UNUSED (thxx_6);
  thyy_8 = tau * 8. / (hy * hy); FIX_UNUSED (thyy_8);
  thyy_6 = 6 * tau / (hy * hy); FIX_UNUSED (thyy_6);
  thxy_1_12 = tau / (12 * hx * hy);
  thxp_05 = 0.5 * tau * p_ro / (hx);
  thyp_05 = 0.5 * tau * p_ro / (hy);
  thxx_4_3 = 4. * tau / (3. * hx * hx);
  thyy_4_3 = 4. * tau / (3. * hy * hy);
  thxx  = tau / (hx * hx);
  thyy  = tau / (hy * hy);
  //------------------------------------------------------------
  double MUM, MU43x, MU43y, MUx, MUy, MUv1, MUv2;

  // A -- sparse matrix of the system, D -- solution vector, B -- rhs vector
  QMatrix A;
  Vector D, B;
  Q_Constr (&A, char_A, 3 * Dim, False, Rowws, Normal, True);
  V_Constr (&B, char_B, 3 * Dim, Normal, True);
  V_Constr (&D, char_D, 3 * Dim, Normal, True);
  SetRTCAccuracy (1e-9);

  // initial values
  tt = 0.;
  mm = 1;
  for (m = 0; m < Dim; m++)
    {
      xx = X[m];
      yy = Y[m];

      tmp = gg (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      G[m] = tmp;
      mm++;

      tmp = u1 (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      V1[m] = tmp;
      mm++;

      tmp = u2 (tt, xx, yy);
      V_SetCmp (&D, mm, tmp);
      V2[m] = tmp;
      mm++;
    }


  for (nn = 1; nn <= N; nn++)
    {
      tt = nn * tau;
      //-///////////////// MUM ///////////////////////////////
      MUM = 0.0;
      for(m = 0; m < Dim; m++)
      {
        if (st[m] == 0)
          {
            tmp = exp(-G[m]);
            if (MUM < tmp)
              MUM = tmp;
          }
      }
      MUM *= mu;
      MU43x = MUM * tau * 4. / (3. * hx * hx);
      MUx   = MUM * tau / (hx * hx);
      MU43y = MUM * tau * 4. / (3. * hy * hy);
      MUy   = MUM * tau / (hy * hy);
      MUv1  = (8. / (3. * hx * hx) + 2. / (hy * hy)) * MUM * tau;
      MUv2  = (8. / (3. * hy * hy) + 2. / (hx * hx)) * MUM * tau;
      //-/////////////////////////////////////////////////////////
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
                V_SetCmp(&B, mm, g00 + tau * Func_g (tt, xx, yy));  // +
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
                Q_SetEntry(&A, mm, 6, mmv10L, tmp);                         // +
                tmp1 = mu * exp(-g00);
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
                Q_SetEntry (&A, mm, 5, mmv2L0,tmp );
                tmp = -thy_05 * (v200 + fabs (v200)) - MU43y;
                Q_SetEntry (&A, mm, 6, mmv20L, tmp );
                tmp1 = mu * exp(-g00);
                tmp = v200 + (tmp1 - MUM)
                      * (thyy_4_3 * (v20R - 2. * v200 + v20L) + thxx * (v2R0 - 2. * v200 + v2L0))
                      + tmp1 * thxy_1_12 * (v1RR - v1RL - v1LR + v1LL) + tau * Func_v2 (tt, xx, yy, p_ro, mu);
                V_SetCmp (&B, mm, tmp);
                mm++;

                break;
              }

            case 1:
              {
                g00 = G[m];
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

                // v2(mx,my)------------------------------------------------------
                Q_SetLen (&A, mm, 1);
                Q_SetEntry (&A, mm, 0, mm, 1.);
                V_SetCmp (&B, mm, u2 (tt, xx, yy));
                mm++;

                break;
              }

            case 2:
              {
                g00 = G[m];
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
                V_SetCmp (&B, mm, u1 (tt, xx, yy));
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
                V_SetCmp (&B, mm, u1 (tt, xx, yy));
                mm++;

                // v2(mx,my)------------------------------------------------------
                Q_SetLen (&A, mm, 1);
                Q_SetEntry (&A, mm, 0, mm, 1.);
                V_SetCmp (&B, mm, u2 (tt, xx, yy));
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
                Q_SetEntry (&A, mm, 0, mm,1.);
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
                Q_SetEntry (&A, mm, 0, mm,1.);
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
                V_SetCmp(&B, mm, g00 + tau * Func_g (tt, xx, yy));
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

      // BiCGIter(&A, &D, &B, 2000, SSORPrecond, 1.);
      BiCGIter(&A, &D, &B, 2000, JacobiPrecond, 1.);
      // BiCGIter(&A, &D, &B, 2000, NULL, 1.);

      mm = 1;
      // Copy solutiong to G, V1, V2 arrays
      for(m = 0; m < Dim; m++)
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
          //scanf ("%lf \n", &xx);
        }
      if (nn == nameiter || nn == N)
        {
            nameiter += N / DIVISOR;
            printf("%d, %d \n", nn, N);
            print_plot(plot_name(plotname, tau, hx, hy, nn), X, Y, G, V1, V2, Dim, tt);
            make_graph(tex_name(texname, tau, hx, hy, nn), plotname, hx, hy, tau, tt);
        }
    }

  Q_Destr (&A);
  V_Destr (&D);
  V_Destr (&B);
}

