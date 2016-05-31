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



int L_op(double * Lapl, const double *u, const UserDataCurr_struct * udc)
{
  //
  //  Lapl[i+j*Nx]= mu *  \delta u[i+j*Nx] + diag[i+j*Nx]*u[i+j*Nx] + u_x u[i+j*Nx]
  //

  int i,j;
  int m;
  int N = udc->N;

  int Nx = udc->Nx;
  int Ny = udc->Ny;

  double Hx = udc->Hx;
  double Hy = udc->Hy;

  double mu = udc->mu;
  const double * diag = udc->diag;
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
      Lapl[i + j * Nx] =
          - mu*(2.*u[i+j*Nx] -u[(i+1)+j*Nx] - u[(i-1)+j*Nx])/Hx/Hx
          - mu*(2.*u[i+j*Nx] -u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy/Hy
          + diag[i+j*Nx]*u[i+j*Nx] + 10.*(u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy
          + 7.*(u[(i+1)+(j)*Nx] - u[(i-1)+(j)*Nx])/Hx;
  */


  // Нужно обратное отображение для выбора i и j по номеру узла m,
  // чтобы заполнить матрицу

  double mmgL0, mmv1L0, mmv2L0, mmgR0, mmv1R0, mmv2R0;
  double mmg0L, mmv10L, mmv20L, mmg0R, mmv10R, mmv20R;
  double g00, gL0, gR0, g0L, g0R, gLL, gRL, gLR, gRR;
  double v100, v1L0, v1R0, v10L, v10R, v1LL, v1RL, v1LR, v1RR;
  double v200, v2L0, v2R0, v20L, v20R, v2LL, v2RL, v2LR, v2RR;

  double J_0L, W1_0L, W2_0L;
  double J_L0, W1_L0, W2_L0;
  double J_00, W1_00, W2_00;
  double J_R0, W1_R0, W2_R0;
  double J_0R, W1_0R, W2_0R;

  double W1_LL, W1_RL, W1_LR, W1_RR;
  double W2_LL, W2_RL, W2_LR, W2_RR;

  double Hx2 = (1. / (2 * Hx));
  double Hy2 = (1. / (2 * Hy));

  int mm = 1;
  for (m = 0; m < N; m++)
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
      gL0 = G[m - 1];
      gR0 = G[m + 1];
      g0L = G[M0L[m]];
      g0R = G[M0R[m]];
      gLL = G[M0L[m] - 1];
      gRL = G[M0L[m] + 1];
      gLR = G[M0R[m] - 1];
      gRR = G[M0R[m] + 1];
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
      switch (st[m])
        {
        case 0:
          {
            // first equation
            J_0L  = -(v200 + fabs (v200)) * Hy2;
            W1_0L = 0.;
            W2_0L = -Hy2;
            J_L0  = -(v100 + fabs (v100)) * Hx2;
            W1_L0 = -Hx2;
            W2_L0 = 0.;
            J_00  = fabs (v100) / Hx + fabs (v200) / Hy;
            W1_00 = (v100 > 0. ? (g00 - gL0) / Hx : (gR0 - g00) / Hx);
            W2_00 = (v200 > 0. ? (g00 - g0L) / Hy : (g0R - g00) / Hy);
            J_R0  = (v100 - fabs (v100)) * Hx2;
            W1_R0 = Hx2;
            W2_R0 = 0.;
            J_0R  = (v200 - fabs (v200)) * Hy2;
            W1_0R = 0.;
            W2_0R = Hy2;

            W1_RR = 0.;
            W1_RL = 0.;
            W1_LR = 0.;
            W1_LL = 0.;

            W2_RR = 0.;
            W2_RL = 0.;
            W2_LR = 0.;
            W2_LL = 0.;

            // second equation
            J_0L  = 0.;
            W1_0L = -(v200 + fabs (v200)) * Hy2 - mu * exp (-g00) * Hy * Hy;
            W2_0L = 0.;
            J_L0  = - p_ro * Hx2; // p_ro = p_ro(g00)
            W1_L0 = -(v100 + fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.) * Hx * Hx;
            W2_L0 = 0.;
            J_00  = p_2ro * (gR0 - gL0) * Hx2 + mu * exp (-g00) * (
                        (4. / 3.) * (v1R0 - 2. * v100 + v1L0) * Hx * Hx +
                        (v10R - 2. * v100 + v10L) * Hy * Hy +
                        (1. / 3.) * (v2RR - v2RL - v2LR + v2LL) * Hx2 * Hy2);
            W1_00 = fabs (v100) * Hx + fabs (v200) * Hy +
                    (v100 > 0. ? (v100 - v1L0) / Hx : (v1R0 - v100) / Hx) +
                    mu * exp (-g00) * 2. * ((4. / 3.) * Hx * Hx + Hy * Hy);
            W2_00 = (v200 > 0. ? (v100 - v10L) / Hy : (v10R - v100) / Hy);
            J_R0  = p_ro * Hx2;
            W1_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00) * (4. / 3.) * Hx * Hx;
            W2_R0 = 0.;
            J_0R  = 0.;
            W1_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00) * Hy * Hy;
            W2_0R = 0.;

            W1_RR = 0.;
            W1_RL = 0.;
            W1_LR = 0.;
            W1_LL = 0.;

            W2_RR = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W2_RL = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W2_LR = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W2_LL = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;

            // third equation

            J_0L  = - p_ro * Hy2;
            W1_0L = 0.;
            W2_0L = -(v200 + fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.) * Hy * Hy;
            J_L0  = 0.;
            W1_L0 = 0.;
            W2_L0 = -(v100 + fabs (v100)) * Hx2 - mu * exp (-g00) * Hx * Hx;
            J_00  = p_2ro * (g0R - g0L) * Hx2 + mu * exp (-g00) * (
                        (4. / 3.) * (v20R - 2. * v200 + v20L) * Hy * Hy +
                        (v2R0 - 2. * v200 + v2L0) * Hx * Hx +
                        (1. / 3.) * (v1RR - v1RL - v1LR + v1LL) * Hx2 * Hy2);
            W1_00 = (v100 > 0. ? (v200 - v2L0) / Hx : (v2R0 - v200) / Hx);
            W2_00 = fabs (v100) * Hx + fabs (v200) * Hy +
                    (v200 > 0. ? (v200 - v20L) / Hy : (v20R - v200) / Hy) +
                    mu * exp (-g00) * 2. * ((4. / 3.) * Hy * Hy + Hx * Hx);
            J_R0  = 0.;
            W1_R0 = 0.;
            W2_R0 = (v100 - fabs (v100)) * Hx2 - mu * exp (-g00) * Hx * Hx;
            J_0R  = p_ro * Hy2;
            W1_0R = 0.;
            W2_0R = (v200 - fabs (v200)) * Hy2 - mu * exp (-g00) * (4. / 3.) * Hy * Hy;;

            W1_RR = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W1_RL = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W1_LR = mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;
            W1_LL = - mu * exp (-g00) * (1. / 3.) * Hx2 * Hy2;

            W2_RR = 0.;
            W2_RL = 0.;
            W2_LR = 0.;
            W2_LL = 0.;
          }
        case 1:
          {
            ;
          }
        case 2:
          {
            ;
          }
        case 3:
          {
            ;
          }
        case 4:
          {
            ;
          }
        case 5:
          {
            ;
          }
        case 6:
          {
            ;
          }
        case 7:
          {
            ;
          }
        case 8:
          {
            ;
          }
        default:
          break;
        }
    }

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
  double *u=NULL;
  double *Lu=NULL;
  UserDataCurr_struct * udc = (UserDataCurr_struct *)ud;

  FIX_UNUSED (n);

  u  = make_vector_double(udc->N, __FILE__, __FUNCTION__);
  /*
   *  TODO -- rename Lu, L_op, Lapl names, because
   *  they suppose Laplace operator case, which was initial
   *  for this propgram.
  */
  Lu = make_vector_double(udc->N, __FILE__, __FUNCTION__);


  /*
   *  TODO -- in our case it`s easier not to use convert function,
   *  because we have unknown functions on boudary. We should modify only
   *  L_op method. However, main eigen code uses convert function, so this
   *  solution propbably could not work, even though Popov said it should...
  */
  convert_au_to_u (u,au,udc);

  L_op (Lu, u, udc);

  // TODO -- simular to previus TODO case.
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
