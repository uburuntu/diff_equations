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


void initparam_UserDataCurr_struct (
    UserDataCurr_struct * udc)
{
  int i,j;

  udc->Nx    =  61;
  udc->Ny    =  61;
  udc->Nx_0  =  41;
  udc->Ny_0  =  21;

  // NA is not used in our method, because we don`t cut area

#if SQUARE
  udc->N = udc->Nx * udc->Ny;
  udc->NA=(udc->Nx - 2) * (udc->Ny - 2);
#else
  udc->N = udc->Nx * udc->Ny - (udc->Nx_0 - 1) * (udc->Ny_0 - 1);
  udc->NA=(udc->Nx - 2) * (udc->Ny - 2) - (udc->Nx_0 - 1) * (udc->Ny_0 - 1);
#endif

  udc->matrix_dim = 3 * udc->N;

  udc->Lx = 3.;
  udc->Ly = 3.;
  udc->Hx = udc->Lx / (udc->Nx - 1);
  udc->Hy = udc->Ly / (udc->Ny - 1);

  udc->mu = 1.;

  udc->diag = make_vector_double (udc->N, __FILE__, __FUNCTION__);

  for(i = 0; i < udc->Nx; i++)
    {
      for(j = 0; j < udc->Ny; j++)
        {
          // TODO: change these for our problem
          udc->diag[i + j * udc->Nx] = 1.;
          // * i * udc->Hx * sin(PI * i * udc->Hx / udc->Lx) * sin(PI * j * udc->Hy / udc->Ly);
        }
    }
  return ;
}



int L_op (double *Lu, const double *u, const UserDataCurr_struct *udc,
          const double *G, const double *V1, const double *V2)
{
  //
  //  Lapl[i+j*Nx]= mu *  \delta u[i+j*Nx] + diag[i+j*Nx]*u[i+j*Nx] + u_x u[i+j*Nx]
  //

  int m;
  int N = udc->N;
  double Hx = udc->Hx;
  double Hy = udc->Hy;
  double mu = udc->mu;
  double p_ro, p_2ro;

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

  double Hx2 = (1. / (2 * Hx));
  double Hy2 = (1. / (2 * Hy));

  // Чтобы избежать багов заполняем здесь ВСЕ уравнения,
  // ненужные потом в конверте не копируем в укороченный вектор au.
  int mm = 0;
  for (m = 0; m < N; m++)
    {
      // Fill all phys properties of current mesh node
      fill_node_phys_prop (m, g00, gL0, gR0, g0L, g0R, gLL,
                           gRL, gLR, gRR, G, M0L, M0R);

      fill_node_phys_prop (m, v100, v1L0, v1R0, v10L, v10R, v1LL,
                           v1RL, v1LR, v1RR, V1, M0L, M0R);

      fill_node_phys_prop (m, v200, v2L0, v2R0, v20L, v20R, v2LL,
                           v2RL, v2LR, v2RR, V2, M0L, M0R);

      fill_node_phys_prop (3 * m + 0, j00, jL0, jR0, j0L, j0R, jLL,
                           jRL, jLR, jRR, u, M0L, M0R);

      fill_node_phys_prop (3 * m + 1, w100, w1L0, w1R0, w10L, w10R, w1LL,
                           w1RL, w1LR, w1RR, u, M0L, M0R);


      fill_node_phys_prop (3 * m + 2, w200, w2L0, w2R0, w20L, w20R, w2LL,
                           w2RL, w2LR, w2RR, u, M0L, M0R);

      switch (st[m])
        {
        case 0:
          {
            // first equation
            a_J_0L  = -(v200 + fabs (v200)) * Hy2;
            a_W1_0L = 0.;
            a_W2_0L = -Hy2;

            a_J_L0  = -(v100 + fabs (v100)) * Hx2;
            a_W1_L0 = -Hx2;
            a_W2_L0 = 0.;

            a_J_00  = fabs (v100) / Hx + fabs (v200) / Hy;
            a_W1_00 = (v100 > 0. ? (g00 - gL0) / Hx : (gR0 - g00) / Hx);
            a_W2_00 = (v200 > 0. ? (g00 - g0L) / Hy : (g0R - g00) / Hy);

            a_J_R0  = (v100 - fabs (v100)) * Hx2;
            a_W1_R0 = Hx2;
            a_W2_R0 = 0.;

            a_J_0R  = (v200 - fabs (v200)) * Hy2;
            a_W1_0R = 0.;
            a_W2_0R = Hy2;

            a_J_RR = 0.;
            a_J_RL = 0.;
            a_J_LR = 0.;
            a_J_LL = 0.;

            a_W1_RR = 0.;
            a_W1_RL = 0.;
            a_W1_LR = 0.;
            a_W1_LL = 0.;

            a_W2_RR = 0.;
            a_W2_RL = 0.;
            a_W2_LR = 0.;
            a_W2_LL = 0.;

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
            a_J_0L  = 0.;
            a_W1_0L = -(v200 + fabs (v200)) * Hy2 - mu * exp (-g00) * Hy * Hy;
            a_W2_0L = 0.;

            a_J_L0  = - p_ro * Hx2; // p_ro = p_ro(g00)
            a_W1_L0 = -(v100 + fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.) * Hx * Hx;
            a_W2_L0 = 0.;

            a_J_00  = p_2ro * (gR0 - gL0) * Hx2 + mu * exp (-g00) * (
                        (4. / 3.) * (v1R0 - 2. * v100 + v1L0) * Hx * Hx +
                        (v10R - 2. * v100 + v10L) * Hy * Hy +
                        (1. / 3.) * (v2RR - v2RL - v2LR + v2LL) * Hx2 * Hy2);
            a_W1_00 = fabs (v100) * Hx + fabs (v200) * Hy +
                    (v100 > 0. ? (v100 - v1L0) / Hx : (v1R0 - v100) / Hx) +
                    mu * exp (-g00) * 2. * ((4. / 3.) * Hx * Hx + Hy * Hy);
            a_W2_00 = (v200 > 0. ? (v100 - v10L) / Hy : (v10R - v100) / Hy);

            a_J_R0  = p_ro * Hx2;
            a_W1_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.) * Hx * Hx;
            a_W2_R0 = 0.;

            a_J_0R  = 0.;
            a_W1_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00) * Hy * Hy;
            a_W2_0R = 0.;

            a_J_RR = 0.;
            a_J_RL = 0.;
            a_J_LR = 0.;
            a_J_LL = 0.;

            a_W1_RR = 0.;
            a_W1_RL = 0.;
            a_W1_LR = 0.;
            a_W1_LL = 0.;

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

            a_J_0L  = - p_ro * Hy2;
            a_W1_0L = 0.;
            a_W2_0L = -(v200 + fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.) * Hy * Hy;

            a_J_L0  = 0.;
            a_W1_L0 = 0.;
            a_W2_L0 = -(v100 + fabs (v100)) * Hx2 - mu * exp (-g00) * Hx * Hx;

            a_J_00  = p_2ro * (g0R - g0L) * Hx2 + mu * exp (-g00) * (
                        (4. / 3.) * (v20R - 2. * v200 + v20L) * Hy * Hy +
                        (v2R0 - 2. * v200 + v2L0) * Hx * Hx +
                        (1. / 3.) * (v1RR - v1RL - v1LR + v1LL) * Hx2 * Hy2);
            a_W1_00 = (v100 > 0. ? (v200 - v2L0) / Hx : (v2R0 - v200) / Hx);
            a_W2_00 = fabs (v100) * Hx + fabs (v200) * Hy +
                    (v200 > 0. ? (v200 - v20L) / Hy : (v20R - v200) / Hy) +
                    mu * exp (-g00) * 2. * ((4. / 3.) * Hy * Hy + Hx * Hx);

            a_J_R0  = 0.;
            a_W1_R0 = 0.;
            a_W2_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00) * Hx * Hx;

            a_J_0R  = p_ro * Hy2;
            a_W1_0R = 0.;
            a_W2_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.) * Hy * Hy;;

            a_J_RR = 0.;
            a_J_RL = 0.;
            a_J_LR = 0.;
            a_J_LL = 0.;

            a_W1_RR = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            a_W1_RL = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            a_W1_LR = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            a_W1_LL = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;

            a_W2_RR = 0.;
            a_W2_RL = 0.;
            a_W2_LR = 0.;
            a_W2_LL = 0.;

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
          }
        case 1:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 2:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 3:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 4:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 5:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 6:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        case 7:
          {
            // first equation
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
            // TODO: fill coefficients
            a_J_0L = a_W1_0L = a_W2_0L = 0.;
            a_J_L0 = a_W1_L0 = a_W2_L0 = 0.;
            a_J_00 = a_W1_00 = a_W2_00 = 0.;
            a_J_R0 = a_W1_R0 = a_W2_R0 = 0.;
            a_J_0R = a_W1_0R = a_W2_0R = 0.;
            a_J_LL = a_W1_LL = a_W2_LL = 0.;
            a_J_LR = a_W1_LR = a_W2_LR = 0.;
            a_J_RL = a_W1_RL = a_W2_RL = 0.;
            a_J_RR = a_W1_RR = a_W2_RR = 0.;

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
          }
        default:
          break;
        }
    }

  // Check that we fill all elements of Lu
  assert (mm == 3 * N);
  return 0;
}



void print_2dfun_double(FILE* f, const char * name, const double  * u,
                        const int nx, const int ny)
{
  int i;
  int j;

  fprintf(f,"#%s\n",name);

  for(j=ny-1;j>=0;j--){
      for(i=0;i<nx;i++){
          fprintf(f,"%9.2e ",u[i+j*nx]);
        }
      fprintf(f,"\n");
    }
  return ;
}


int convert_u_to_au(double * au, const double  * u,
                    const UserDataCurr_struct * udc)
{

  int i,j;

  int Nx,Ny;

  Nx=udc->Nx;
  Ny=udc->Ny;

  for(i=0;i<Nx-2;i++)
    for(j=0;j<Ny-2;j++)
      au[i+j*(Nx-2)]=u[(i+1)+(j+1)*Nx];


  return 0;
}

int convert_au_to_u(double *u, const double  * au,
                    const UserDataCurr_struct * udc)
{

  int i,j;
  int Nx,Ny;

  Nx=udc->Nx;
  Ny=udc->Ny;

  i=0;
  for(j=0;j<Ny;j++)
    u[i+j*Nx]=0.;

  i=Nx-1;
  for(j=0;j<Ny;j++)
    u[i+j*Nx]=0.;

  j=0;
  for(i=0;i<Nx;i++)
    u[i+j*Nx]=0.;

  j=Ny-1;
  for(i=0;i<Nx;i++)
    u[i+j*Nx]=0.;



  for(i=1;i<Nx-1;i++)
    for(j=1;j<Ny-1;j++)
      u[i+j*Nx]=au[(i-1)+(j-1)*(Nx-2)];

  return 0;

}


void A_op (double *Aau, const double *au, int n, void * ud)
{
  double *u = NULL;
  double *Lu = NULL;
  UserDataCurr_struct * udc = (UserDataCurr_struct *)ud;

  FIX_UNUSED (n);

  u  = make_vector_double (3 * udc->N, __FILE__, __FUNCTION__);
  /*
   *  TODO -- rename Lu, L_op, Lapl names, because
   *  they suppose Laplace operator case, which was initial
   *  for this propgram.
  */
  Lu = make_vector_double (3 * udc->N, __FILE__, __FUNCTION__);


  // convert solution vector au (dim = NA < 3 * N) into solution vector u
  convert_au_to_u (u, au, udc);

  // multiplication L * u (dim(L) = (3 * N) ^ 2, dim(u) = 3 * N
  L_op (Lu, u, udc);

  // convert solution vector Lu (dim = 3 * N) into solution vector Aau
  convert_u_to_au (Aau, Lu, udc);

  free(u);
  free(Lu);

  return;
}



double *make_vector_double (int n, const char *info_1, const char *info_2)
{
  double *u;
  int i;

  u = (double*) malloc (sizeof(double) * n);
  if (u == NULL)
    {
      printf ("Error in %s %s: Fail to allocate %lu bytes\n",
              info_1, info_2, sizeof(double) * n);
      exit(1);
    }

  for(i = 0; i < n; i++)
      u[i] = 0.;
  return u;
}


int *make_vector_int (int n, const char *info_1, const char *info_2)
{
  int *u;
  int i;

  u = (int*) malloc (sizeof(int) * n);
  if (u == NULL)
    {
      printf ("Error in %s %s: Fail to allocate %lu bytes\n",
              info_1, info_2, sizeof(int) * n);
      exit(1);
    }

  for(i = 0; i < n; i++)
      u[i] = 0.;
  return u;
}

void fill_node_phys_prop (int m, // number of mesh node
                          double &p00, double &pL0, double &pR0,
                          double &p0L, double &p0R, double &pLL,
                          double &pRL, double &pLR, double &pRR,
                          const double *p, const int *M0L, const int *M0R)
{
  p00 = p[m];
  pL0 = p[m - 1];
  pR0 = p[m + 1];
  p0L = p[M0L[m]];
  p0R = p[M0R[m]];
  pLL = p[M0L[m] - 1];
  pRL = p[M0L[m] + 1];
  pLR = p[M0R[m] - 1];
  pRR = p[M0R[m] + 1];
}
