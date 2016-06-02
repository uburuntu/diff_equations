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

#define RECREATE_COEFFS       \
  a_J_0L = a_W1_0L = a_W2_0L =\
  a_J_L0 = a_W1_L0 = a_W2_L0 =\
  a_J_00 = a_W1_00 = a_W2_00 =\
  a_J_R0 = a_W1_R0 = a_W2_R0 =\
  a_J_0R = a_W1_0R = a_W2_0R =\
  a_J_LL = a_W1_LL = a_W2_LL =\
  a_J_LR = a_W1_LR = a_W2_LR =\
  a_J_RL = a_W1_RL = a_W2_RL =\
  a_J_RR = a_W1_RR = a_W2_RR = 0.

void init_user_data (
  user_data *ud)
{
  // Number of mesh nodes
  ud->Nx   =  61;
  ud->Ny   =  61;
  ud->Nx_0 =  41;
  ud->Ny_0 =  21;

#if SQUARE
  ud->N  = ud->Nx * ud->Ny;
  ud->NA = 3 * (ud->Nx - 2) * (ud->Ny - 2) + // inner nodes
           2 * (ud->Ny - 2) +                 // right boundary
           1 * (ud->Nx - 2) +                 // down boundary
           1 * (ud->Nx - 2) +                 // top boundary
           0 * (ud->Ny - 2) +                 // left boundary
           0 * 4;                              // vertices of square
#else
  ud->N = ud->Nx * ud->Ny - (ud->Nx_0 - 1) * (ud->Ny_0 - 1);
  ud->NA = (ud->Nx - 2) * (ud->Ny - 2) - (ud->Nx_0 - 1) * (ud->Ny_0 - 1);
  ud->NA = 3 * ((ud->Nx - 2) * (ud->Ny - 2) - (ud->Nx_0 - 1) * (ud->Ny_0 - 1) + 1) +
           2 * (ud->Nx - ud->Nx_0 - 1) +     // II-part of down boundary
           1 * (ud->Ny_0 - 2) +               // I-part of left boundary
           1 * (ud->Ny - 2) +                 // right boundary
           1 * (ud->Nx - 2) +                 // top boundary
           1 * (ud->Nx_0 - 2) +               // I-part of down boundary
           0 * (ud->Ny - ud->Ny_0 - 1) +     // II-part of left boundary
           0 * 5;                              // vertices
#endif

  ud->Lx = 3.;
  ud->Ly = 3.;
  ud->Hx = ud->Lx / (ud->Nx - 1);
  ud->Hy = ud->Ly / (ud->Ny - 1);

  ud->p_ro = 1.0;
  ud->p_2ro = 0.0;
  ud->mu = 0.1;
  return;
}

int L_op (double *Lu, const double *u, const user_data *ud,
          const double *G, const double *V1,  const double *V2, const int *st,
          const int *M0L, const int *M0R, const double *X, const double *Y)
{
  //
  //  Lapl[i+j*Nx]= mu *  \delta u[i+j*Nx] + diag[i+j*Nx]*u[i+j*Nx] + u_x u[i+j*Nx]
  //

  int m;
  int N = ud->N;
  double Hx = ud->Hx;
  double Hy = ud->Hy;
  double mu = ud->mu;

  double p_ro = ud->p_ro;
  double p_2ro = ud->p_2ro;

  // DON`T DELETE THIS PARAMETER
  // IT MAY BE USED FOR ANOTHER NO-SQUARE MESH
  // !!!!
  FIX_UNUSED(X);

  double *J   = make_vector_double (N, __FILE__, __FUNCTION__);
  double *W1  = make_vector_double (N, __FILE__, __FUNCTION__);
  double *W2  = make_vector_double (N, __FILE__, __FUNCTION__);

  for (m = 0; m < N; m++)
    {
      J[m]  = u[3 * m + 0];
      W1[m] = u[3 * m + 1];
      W2[m] = u[3 * m + 2];
    }

  /*
   * TODO -- now this loop is for laplacian equation:
   * mu * (u_{x \bar{x}} + u_{y \bar{y}}) + du + 10 * u_x + 7 * u_y = lambda * u
   * We should rewrite it for our linearized scheme equation (23.5 - 23.6)
   * on page 48-49.
  */

  /*

  for (i = 1; i < Nx - 1; i++)
    for(j = 1; j < Ny - 1; j++)
      Lu[i + j * Nx] =
          - mu*(2.*u[i+j*Nx] -u[(i+1)+j*Nx] - u[(i-1)+j*Nx])/Hx/Hx
          - mu*(2.*u[i+j*Nx] -u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy/Hy
          + diag[i+j*Nx]*u[i+j*Nx] + 10.*(u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy
          + 7.*(u[(i+1)+(j)*Nx] - u[(i-1)+(j)*Nx])/Hx;
  */

  // Local variables of all properties for current mesh node
  double g00, gL0, gR0, g0L, g0R, gLL, gRL, gLR, gRR;
  double v100, v1L0, v1R0, v10L, v10R, v1LL, v1RL, v1LR, v1RR;
  double v200, v2L0, v2R0, v20L, v20R, v2LL, v2RL, v2LR, v2RR;
  double j00, jL0, jR0, j0L, j0R, jLL, jRL, jLR, jRR;
  double w100, w1L0, w1R0, w10L, w10R, w1LL, w1RL, w1LR, w1RR;
  double w200, w2L0, w2R0, w20L, w20R, w2LL, w2RL, w2LR, w2RR;

  // Сoefficients of linearization
  double a_J_0L, a_W1_0L, a_W2_0L;
  double a_J_L0, a_W1_L0, a_W2_L0;
  double a_J_00, a_W1_00, a_W2_00;
  double a_J_R0, a_W1_R0, a_W2_R0;
  double a_J_0R, a_W1_0R, a_W2_0R;
  double a_J_LL, a_W1_LL, a_W2_LL;
  double a_J_LR, a_W1_LR, a_W2_LR;
  double a_J_RL, a_W1_RL, a_W2_RL;
  double a_J_RR, a_W1_RR, a_W2_RR;

  double Hx2 = (1. / (2  / Hx));
  double Hy2 = (1. / (2  / Hy));

  // Чтобы избежать багов заполняем здесь ВСЕ уравнения,
  // ненужные потом в конверте не копируем в укороченный вектор au.
  int mm = 0;

  for (m = 0; m < N; m++)
    {
      // Fill all phys properties of current mesh node
      fill_node_phys_prop (m,
                           &g00, &gL0, &gR0, &g0L, &g0R, &gLL, &gRL, &gLR, &gRR,
                           G, M0L, M0R);

      fill_node_phys_prop (m,
                           &v100, &v1L0, &v1R0, &v10L, &v10R, &v1LL, &v1RL, &v1LR, &v1RR,
                           V1, M0L, M0R);

      fill_node_phys_prop (m,
                           &v200, &v2L0, &v2R0, &v20L, &v20R, &v2LL, &v2RL, &v2LR, &v2RR,
                           V2, M0L, M0R);

      fill_node_phys_prop (m,
                           &j00, &jL0, &jR0, &j0L, &j0R, &jLL, &jRL, &jLR, &jRR,
                           J, M0L, M0R);

      fill_node_phys_prop (m,
                           &w100, &w1L0, &w1R0, &w10L, &w10R, &w1LL, &w1RL, &w1LR, &w1RR,
                           W1, M0L, M0R);


      fill_node_phys_prop (m,
                           &w200, &w2L0, &w2R0, &w20L, &w20R, &w2LL, &w2RL, &w2LR, &w2RR,
                           W2, M0L, M0R);

      switch (st[m])
        {
          case 0:
            {
              // first equation
              RECREATE_COEFFS;

              a_J_0L  = - (v200 + fabs (v200)) * Hy2;
              a_W2_0L = -Hy2;

              a_J_L0  = - (v100 + fabs (v100)) * Hx2;
              a_W1_L0 = -Hx2;

              a_J_00  = fabs (v100) / Hx + fabs (v200) / Hy;
              a_W1_00 = (v100 > 0. ? (g00 - gL0) / Hx : (gR0 - g00) / Hx);
              a_W2_00 = (v200 > 0. ? (g00 - g0L) / Hy : (g0R - g00) / Hy);

              a_J_R0  = (v100 - fabs (v100)) * Hx2;
              a_W1_R0 = Hx2;

              a_J_0R  = (v200 - fabs (v200)) * Hy2;
              a_W2_0R = Hy2;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation
              RECREATE_COEFFS;
              a_W1_0L = - (v200 + fabs (v200)) * Hy2 - mu * exp (-g00)  / Hy  / Hy;

              a_J_L0  = -p_ro * Hx2; // p_ro = p_ro(g00)
              a_W1_L0 = - (v100 + fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.)  / Hx  / Hx;

              a_J_00  = p_2ro * (gR0 - gL0) * Hx2 + mu * exp (-g00) * (
                          (4. / 3.) * (v1R0 - 2. * v100 + v1L0)  / Hx  / Hx +
                          (v10R - 2. * v100 + v10L)  / Hy  / Hy +
                          (1. / 3.) * (v2RR - v2RL - v2LR + v2LL) * Hx2 * Hy2);
              a_W1_00 = fabs (v100)  / Hx + fabs (v200)  / Hy +
                        (v100 > 0. ? (v100 - v1L0) / Hx : (v1R0 - v100) / Hx) +
                        mu * exp (-g00) * 2. * ((4. / 3.)  / Hx  / Hx + Hy  / Hy);
              a_W2_00 = (v200 > 0. ? (v100 - v10L) / Hy : (v10R - v100) / Hy);

              a_J_R0  = p_ro * Hx2;
              a_W1_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.)  / Hx  / Hx;

              a_W1_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00)  / Hy  / Hy;

              a_W2_RR = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W2_RL = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W2_LR = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W2_LL = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation
              RECREATE_COEFFS;
              a_J_0L  = - p_ro * Hy2;
              a_W2_0L = - (v200 + fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.)  / Hy  / Hy;

              a_W2_L0 = - (v100 + fabs (v100)) * Hx2 - mu * exp (-g00)  / Hx  / Hx;

              a_J_00  = p_2ro * (g0R - g0L) * Hx2 + mu * exp (-g00) * (
                          (4. / 3.) * (v20R - 2. * v200 + v20L)  / Hy  / Hy +
                          (v2R0 - 2. * v200 + v2L0)  / Hx  / Hx +
                          (1. / 3.) * (v1RR - v1RL - v1LR + v1LL) * Hx2 * Hy2);
              a_W1_00 = (v100 > 0. ? (v200 - v2L0) / Hx : (v2R0 - v200) / Hx);
              a_W2_00 = fabs (v100)  / Hx + fabs (v200)  / Hy +
                        (v200 > 0. ? (v200 - v20L) / Hy : (v20R - v200) / Hy) +
                        mu * exp (-g00) * 2. * ((4. / 3.)  / Hy  / Hy + Hx  / Hx);

              a_W2_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00)  / Hx  / Hx;

              a_J_0R  = p_ro * Hy2;
              a_W2_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.)  / Hy  / Hy;;

              a_W1_RR = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W1_RL = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W1_LR = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
              a_W1_LL = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 1:
            {
              // first equation

              RECREATE_COEFFS;
              a_W1_R0 = 1.;
              a_W1_00 = -1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 2:
            {
              // first equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;
              a_W1_L0 = -1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;
              a_W1_L0 = (SQUARE ? -1. : 0.);

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 3:
            {
              // first equation

              RECREATE_COEFFS;
              a_W2_0R = 1.;
              a_W2_00 = -1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;
              a_W2_0R = (!SQUARE && is_equal(Y[m], 0.) ? -1. : 0.);

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 4:
            {
              // first equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;
              a_W2_0L = -1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 5:
            {
              // first equation

              RECREATE_COEFFS;
              a_J_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }
          case 6:
            {
              // first equation

              RECREATE_COEFFS;
              a_J_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 7:
            {
              // first equation

              RECREATE_COEFFS;
              a_J_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          case 8:
            {
              // first equation

              RECREATE_COEFFS;
              a_J_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // second equation

              RECREATE_COEFFS;
              a_W1_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;

              // third equation

              RECREATE_COEFFS;
              a_W2_00 = 1.;

              // New Lu elements:
              Lu[mm] =
                a_J_00 * j00 + a_J_L0 * jL0 + a_J_R0 * jR0 + a_J_0L * j0L +
                a_J_0R * j0R + a_J_LL * jLL + a_J_RL * jRL + a_J_LR * jLR +
                a_J_RR * jRR +
                a_W1_00 * w100 + a_W1_L0 * w1L0 + a_W1_R0 * w1R0 + a_W1_0L * w10L +
                a_W1_0R * w10R + a_W1_LL * w1LL + a_W1_RL * w1RL + a_W1_LR * w1LR +
                a_W1_RR * w1RR +
                a_W2_00 * w200 + a_W2_L0 * w2L0 + a_W2_R0 * w2R0 + a_W2_0L * w20L +
                a_W2_0R * w20R + a_W2_LL * w2LL + a_W2_RL * w2RL + a_W2_LR * w2LR +
                a_W2_RR * w2RR;
              mm++;
              break;
            }

          default:
            break;
        }
    }

  // Check that we fill all elements of Lu
  assert (mm == 3 * N);
  FREE_ARRAY (J);
  FREE_ARRAY (W1);
  FREE_ARRAY (W2);
  return 0;
}

int convert_u_to_au (double *au, const double  *u, const user_data *ud,
                     const int *st, const double *X, const double *Y)
{
  int m;
  int N  = ud->N;
  int NA = ud->NA;
  int m1 = 0; // u-index
  int m2 = 0; // au-index

  // DON`T DELETE THIS PARAMETER
  // IT MAY BE USED FOR ANOTHER NO-SQUARE MESH
  // !!!!
  FIX_UNUSED(X);

  for (m = 0; m < N; m++)
    {
      switch (st[m])
        {
          case 0:
            {
              // first equation
              au[m2] = u[m1];
              m1++;
              m2++;
              // second equation
              au[m2] = u[m1];
              m1++;
              m2++;
              // third equation
              au[m2] = u[m1];
              m1++;
              m2++;
              break;
            }

          case 1: // left boundary
            {
              // 0 non-trivial equations
              // first equation
              m1++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          case 2: // right boundary
            {
              // 1 or 2 non-trivial equation
              // first equation
              au[m2] = u[m1];
              m1++;
              m2++;
              // second equation
              if (SQUARE)
                {
                  au[m2] = u[m1];
                  m2++;
                }
              m1++;
              // third equation
              m1++;
              break;
            }

          case 3: // down boundary
            {
              // 1 or 2 non-trivial equation
              // first equation
              au[m2] = u[m1];
              m1++;
              m2++;
              // second equation
              m1++;
              // third equation
              if (!SQUARE && is_equal (Y[m1], 0.))
                {
                  au[m2] = u[m1];
                  m2++;
                }
              m1++;
              break;
            }

          case 4: // top boundary
            {
              // 1 non-trivial equation
              // first equation
              au[m2] = u[m1];
              m1++;
              m2++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          case 5:
            {
              // 0 non-trivial equations
              // first equation
              m1++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          case 6:
            {
              // 0 non-trivial equations
              // first equation
              m1++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          case 7:
            {
              // 0 non-trivial equations
              // first equation
              m1++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          case 8:
            {
              // 0 non-trivial equations
              // first equation
              m1++;
              // second equation
              m1++;
              // third equation
              m1++;
              break;
            }

          default:
            break;
        }
    }

  // Check that we fill all elements correctly
  assert (m1 == 3 * N);
  assert (m2 == NA);
  return 0;
}

int convert_au_to_u (double *u, const double  *au, const user_data *ud,
                     const int *st, const double *X, const double *Y)
{
  int m;
  int N  = ud->N;
  int NA = ud->NA;
  int m1 = 0; // u-index
  int m2 = 0; // au-index

  // DON`T DELETE THIS PARAMETER
  // IT MAY BE USED FOR ANOTHER NO-SQUARE MESH
  // !!!!
  FIX_UNUSED(X);

  for (m = 0; m < N; m++)
    {
      switch (st[m])
        {
          case 0:
            {
              // first equation
              u[m1] = au[m2];
              m1++;
              m2++;
              // second equation
              u[m1] = au[m2];
              m1++;
              m2++;
              // third equation
              u[m1] = au[m2];
              m1++;
              m2++;
              break;
            }

          case 1: // left boundary
            {
              // 0 non-trivial equations
              // first equation
              u[m1] =  0.;
              m1++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 2: // right boundary
            {
              // 1 or 2 non-trivial equation
              // first equation
              u[m1] = au[m2];
              m1++;
              m2++;
              // second equation
              if (SQUARE)
                {
                  u[m1] = au[m2];
                  m2++;
                }
              else
                {
                  u[m1] = 0.;
                }
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 3: // down boundary
            {
              // 1 or 2 non-trivial equation
              // first equation
              u[m1] = au[m2];
              m1++;
              m2++;
              // second equation
              u[m1] = 0.;
              m1++;
              // third equation
              if (!SQUARE && is_equal (Y[m1], 0.))
                {
                  u[m1] = au[m2];
                  m2++;
                }
              else
                {
                  u[m1] =  0.;
                }
              m1++;
              break;
            }

          case 4: // top boundary
            {
              // 1 non-trivial equation
              // first equation
              u[m1] = au[m2];
              m1++;
              m2++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 5:
            {
              // 0 non-trivial equations
              // first equation
              u[m1] =  0.;
              m1++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 6:
            {
              // 0 non-trivial equations
              // first equation
              u[m1] =  0.;
              m1++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 7:
            {
              // 0 non-trivial equations
              // first equation
              u[m1] =  0.;
              m1++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          case 8:
            {
              // 0 non-trivial equations
              // first equation
              u[m1] =  0.;
              m1++;
              // second equation
              u[m1] =  0.;
              m1++;
              // third equation
              u[m1] =  0.;
              m1++;
              break;
            }

          default:
            break;
        }
    }

  // Check that we fill all elements correctly
  assert (m1 == 3 * N);
  assert (m2 == NA);
  return 0;
}

void A_op (double *Aau, const double *au, const user_data *ud,
           const double *G, const double *V1,  const double *V2, const int *st,
           const int *M0L, const int *M0R, const double *X, const double *Y)
{
  double *u = NULL, *Lu = NULL;

  u  = make_vector_double (3 * ud->N, __FILE__, __FUNCTION__);
  Lu = make_vector_double (3 * ud->N, __FILE__, __FUNCTION__);

  // convert solution vector au (dim = NA < 3 * N) into solution vector u
  convert_au_to_u (u, au, ud, st, X, Y);

  // multiplication L * u (dim(L) = (3 * N) ^ 2, dim(u) = 3 * N
  L_op (Lu, u, ud, G, V1, V2, st, M0L, M0R, X, Y);

  // convert solution vector Lu (dim = 3 * N) into solution vector Aau
  convert_u_to_au (Aau, Lu, ud, st, X, Y);

  FREE_ARRAY (u);
  FREE_ARRAY (Lu);

  return;
}

void fill_node_phys_prop (int m, // number of mesh node
                          double *p00, double *pL0, double *pR0,
                          double *p0L, double *p0R, double *pLL,
                          double *pRL, double *pLR, double *pRR,
                          const double *const p,
                          const    int *const M0L,
                          const    int *const M0R)
{
  *p00 = p[m];
  *pL0 = p[m - 1];
  *pR0 = p[m + 1];
  *p0L = p[M0L[m]];
  *p0R = p[M0R[m]];
  *pLL = p[M0L[m] - 1];
  *pRL = p[M0L[m] + 1];
  *pLR = p[M0R[m] - 1];
  *pRR = p[M0R[m] + 1];
}
