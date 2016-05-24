/*

 Copyright (c)  2016 Kornev Andrey A.
                     Ozeritsky  Alexey V.
*/


#include "f.h"

void initparam_UserDataCurr_struct (
    UserDataCurr_struct * udc)
{

  double PI=3.14159265358979323846;
  int i,j;

  udc->Nx    =  51;
  udc->Ny    =  51;
  udc->N=udc->Nx*udc->Ny;
  udc->NA=(udc->Nx-2)*(udc->Ny-2);

  udc->Lx =PI;
  udc->Ly =PI;
  udc->Hx = udc->Lx/(udc->Nx-1);
  udc->Hy = udc->Ly/(udc->Ny-1);


  udc->mu = 1.;

  udc->diag = make_vector_double (udc->N, __FILE__, __FUNCTION__);

  for(i=0;i<udc->Nx;i++)
    {
      for(j=0;j<udc->Ny;j++)
        {
          udc->diag[i+j*udc->Nx] = 1.*i*udc->Hx* sin(PI*i*udc->Hx/udc->Lx)
              * sin(PI*j*udc->Hy/udc->Ly);
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

  int Nx=udc->Nx;
  int Ny=udc->Ny;

  double Hx=udc->Hx;
  double Hy=udc->Hy;

  double mu=udc->mu;
  const double * diag=udc->diag;

  /*
   * TODO -- now this loop is for laplacian equation:
   * mu * (u_{x \bar{x}} + u_{y \bar{y}}) + du + 10 * u_x + 7 * u_y = lambda * u
   * We should rewrite it for our linearized scheme equation (23.5 - 23.6)
   * on page 48-49.
  */

  for (i = 1; i < Nx - 1; i++)
    for(j = 1; j < Ny - 1; j++)
      Lapl[i + j * Nx] =
          - mu*(2.*u[i+j*Nx] -u[(i+1)+j*Nx] - u[(i-1)+j*Nx])/Hx/Hx
          - mu*(2.*u[i+j*Nx] -u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy/Hy
          + diag[i+j*Nx]*u[i+j*Nx] + 10.*(u[i+(j+1)*Nx] - u[i+(j-1)*Nx])/Hy
          + 7.*(u[(i+1)+(j)*Nx] - u[(i-1)+(j)*Nx])/Hx;

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


double *  make_vector_double(int n,const char *info_1, const char *info_2)
{
  double *u;
  int i;

  u = (double*)malloc( sizeof(double)*n);
  if(u==NULL){
      printf("Error in %s %s: Fail to allocate %lu bytes\n",
             info_1,info_2,sizeof(double)*n);
      exit(1);
    }

  for(i=0;i<n ;i++){
      u[i] = 0.;

    }
  return u;
}
