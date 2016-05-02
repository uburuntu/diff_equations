#include <math.h>
#include "func.h"

#define cave 0
#define square 1
#define light_funcs 1
#define light_g 1
#define light_u1 1
#define light_u2 1

void param_dif(P_dif *p_d)
{
  p_d->Segm_T = 1.;
  p_d->Segm_X = 3.;
  p_d->Segm_Y = 3.;
  p_d->Segm_X_0 = 2.;
  p_d->Segm_Y_0 = 1.;
  p_d->p_ro = 1.0;
  p_d->mu = 0.1;
}

void param_she_step(P_she *p_s, P_dif *p_d, int it_t, int it_sp)
{
  // main rectangle
  p_s->M_x = 60;
  p_s->M_y = 60;
  p_s->N   = 20;

  // cut rectangle
  p_s->M_x_0 = 40;
  p_s->M_y_0 = 20;

  int i, k_sp, k_t;
  k_sp = 1;
  k_t = 1;

  if (it_sp > 0)
    {
      for (i = 1; i <= it_sp; i++)
        k_sp *= 2;
    }
  if (it_t > 0)
    {
      for (i = 1; i <= it_t; i++)
        k_t *= 2;
    }
  p_s->M_x *= k_sp;
  p_s->M_y *= k_sp;
  p_s->N   *= k_t;

  p_s->M_x_0 *= k_sp;
  p_s->M_y_0 *= k_sp;

#if square
  p_s->Dim = (p_s->M_x + 1) * (p_s->M_y + 1);
#else
  p_s->Dim = (p_s->M_x + 1) * (p_s->M_y + 1) - (p_s->M_x_0 + 0) * (p_s->M_y_0 + 0);
#endif
  p_s->h_x = p_d->Segm_X / p_s->M_x;
  p_s->h_y = p_d->Segm_Y / p_s->M_y;
  p_s->tau = p_d->Segm_T / p_s->N;
  p_s->eta = 0.;
}

//*******************************************************

double Norm_c (double *a, int Dim, double *X, double *Y, double t, double (*f) (double tt,double x1,double x2))
{
  int m;
  double norma = 0.;
  double tmp;

  for(m = 0; m < Dim; m++)
    {
      tmp = fabs (a[m] - (*f) (t, X[m], Y[m]));
      if(tmp > norma)
        {
          norma = tmp;
        }
    }
  return norma;
}

double Norm_l2 (double *a, int Dim, double *X, double *Y, double t, double (*f) (double tt,double x1,double x2))
{
  int m;
  double norma = 0.;
  double tmp;

  for(m = 0; m < Dim; m++)
    {
      tmp = a[m] - (*f) (t, X[m], Y[m]);
      norma += tmp * tmp;
    }
  return sqrt (norma / Dim);
}

//*******************************************************

#if cave

double ro(double t, double x, double y)
{
  return (double)((cos(2*M_PI*x)+1.5)*(sin(2*M_PI*y)+1.5)*exp(t));
}

double gg(double t, double x, double y)
{
  return log(ro(t,x,y));
}


double p(double t, double x, double y, double p_ro)
{
  return p_ro*ro(t,x,y);
}


double u1(double t, double x, double y)
{
  return (double)(sin(2*M_PI*x)*sin(2*M_PI*y)*exp(t));
}

double u2(double t, double x, double y)
{
  return (double)(sin(2*M_PI*x)*sin(2*M_PI*y)*exp(-t));
}

double Func_g(double t, double x, double y)
{
  return (double)(1+u1(t,x,y)*
                  (-2*M_PI*sin(2*M_PI*x))/(cos(2*M_PI*x)+1.5)
                  +u2(t,x,y)*(2*M_PI*cos(2*M_PI*y))
                  /(sin(2*M_PI*y)+1.5)+2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y)*exp(t)
                  +2*M_PI*cos(2*M_PI*y)*sin(2*M_PI*x)*exp(-t));
}

double Func_v1(double t, double x, double y, double p_ro, double mu)
{
  return (double)( u1(t,x,y)*(1 +
                              2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y)*exp(t))
                   + u2(t,x,y)*2*M_PI*sin(2*M_PI*x)*
                   cos(2*M_PI*y)*exp(t) + p_ro*(-2*M_PI*sin(2*M_PI*x))/(cos(2*M_PI*x)+1.5)
                   -mu*exp(-gg(t,x,y))*(4./3.*
                                        (-4)*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*exp(t)
                                        - 4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*exp(t) +
                                        4./3.*M_PI*M_PI*cos(2*M_PI*x)*cos(2*M_PI*y)*exp(-t) ));
}

double Func_v2(double t, double x, double y, double p_ro, double mu)
{
  return (double)( u2(t,x,y)*
                   (-1 + 2*M_PI*cos(2*M_PI*y)*sin(2*M_PI*x)*exp(-t) )
                   + u1(t,x,y)*2*M_PI*sin(2*M_PI*y)*
                   cos(2*M_PI*x)*exp(-t) +p_ro*( (2*M_PI*cos(2*M_PI*y))/
                                                 (sin(2*M_PI*y)+1.5) ) -mu*exp(-gg(t,x,y))*
                   (-4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*exp(-t)
                    -16./3.*M_PI*M_PI*sin(2*M_PI*y)*sin(2*M_PI*x)*exp(-t) +
                    4./3.*M_PI*M_PI*cos(2*M_PI*x)*cos(2*M_PI*y)*exp(t)));
}

#else
inline double ro (double t, double x, double y)
{
#if light_g
  return (double) ((cos (M_PI * x) + 1.5) * (sin (M_PI * y) + 1.5));
#else
  return (double) ((cos (M_PI * x) + 1.5) * (sin (M_PI * y) + 1.5) * exp (t));
#endif
}

inline double gg (double t, double x, double y)
{
  return log (ro (t, x, y));
}

#define d_gg_dt 0.
#define d_gg_dx (-M_PI) * sin (M_PI * x) / (cos (M_PI * x) + 1.5)
#define d_gg_dy (+M_PI) * cos (M_PI * y) / (sin (M_PI * y) + 1.5)

inline double p (double t, double x, double y, double p_ro)
{
  return p_ro * ro (t, x, y);
}

inline double u1 (double t, double x, double y)
{
#if light_u1
  return (double) (sin (M_PI * x));
#else
  return (double) (sin (M_PI * x) * sin (M_PI * y) * exp (t));
#endif
}

#define d_u1_dt 0.
#define d_u1_dx M_PI * cos (M_PI * x)
#define d_u1_dy 0.
#define d_u1_dxdx -M_PI * M_PI * sin (M_PI * x)
#define d_u1_dydy 0.
#define d_u1_dxdy 0.

inline double u2 (double t, double x, double y)
{
#if light_u2
  return 0.;
#else
  return (double) (sin (M_PI * x) * sin (M_PI * y) * exp (-t));
#endif
}

#define d_u2_dt 0.
#define d_u2_dx 0.
#define d_u2_dy 0.
#define d_u2_dxdx 0.
#define d_u2_dydy 0.
#define d_u2_dxdy 0.

inline double Func_g(double t, double x, double y)
{
#if light_funcs
  return (double) (d_gg_dt + u1 (t, x, y) * d_gg_dx + u2 (t, x, y) * d_gg_dy + d_u1_dx + d_u2_dy);

#else
  return (double) (1
                   + u1 (t, x, y) * (-1 * M_PI * sin (M_PI * x)) / (cos (M_PI * x) + 1.5)
                   + u2 (t, x, y) * (+1 * M_PI * cos (M_PI * y)) / (sin (M_PI * y) + 1.5)
                   + M_PI * cos (M_PI * x) * sin (M_PI * y) * exp (t)
                   + M_PI * cos (M_PI * y) * sin (M_PI * x) * exp (-t));
#endif
}

inline double Func_v1(double t, double x, double y, double p_ro, double mu)
{
#if light_funcs
  return (double) (d_u1_dt + u1 (t, x, y) * d_u1_dx + u2 (t, x, y) * d_u1_dy + p_ro * d_gg_dx
                   - mu * exp (-gg (t, x, y)) * ((4. / 3) * d_u1_dxdx));
#else
  return (double) (u1 (t, x, y)
                   + u1 (t, x, y) * (M_PI * cos (M_PI * x) * sin (M_PI * y) * exp (t))
                   + u2 (t, x, y) * (M_PI * sin (M_PI * x) * cos (M_PI * y) * exp (t))
                   + p_ro * (-1 * M_PI * sin (M_PI * x)) / (cos (M_PI*  x) + 1.5)
                   - mu * exp (-gg (t, x, y))
                   * ((4. / 3.) * (-1) * M_PI * M_PI * sin (M_PI * x) * sin (M_PI * y) * exp (t)
                      + (-1) * M_PI * M_PI * sin (M_PI * x) * sin (M_PI * y) * exp (t)
                      + (1. / 3.) * M_PI * M_PI * cos (M_PI * x) * cos (M_PI * y) * exp (-t)));
#endif
}

inline double Func_v2(double t, double x, double y, double p_ro, double mu)
{
#if light_funcs
  return (double) (d_u2_dt + u1 (t, x, y) * d_u2_dx + u2 (t, x, y) * d_u2_dy + p_ro * d_gg_dy);
#else
  return (double) ((-1) * u2 (t, x, y)
                   + u2 (t, x, y) * (M_PI * cos (M_PI * y) * sin (M_PI * x) * exp (-t))
                   + u1 (t, x, y) * (M_PI * sin (M_PI * y) * cos (M_PI * x) * exp (-t))
                   + p_ro * ((M_PI * cos (M_PI * y)) / (sin (M_PI * y) + 1.5))
                   - mu * exp (-gg (t, x, y))
                   * ((4. / 3.) * (-1) * M_PI * M_PI * sin (M_PI * x) * sin (M_PI * y) * exp (-t)
                      + (-1) * M_PI * M_PI * sin (M_PI * x) * sin (M_PI * y) * exp (-t)
                      + (1. / 3.) * M_PI * M_PI * cos (M_PI * y) * cos (M_PI * x) * exp (t)));
#endif
}
#endif
