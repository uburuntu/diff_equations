#include <stdio.h>
#include <assert.h>
#include "f.h"

/*
 *
 *     Area:
 *     ________________________________
 *     |         |          |          |
 *     |         |          |          |
 *     |   02    |    12    |    22    |
 *     |         |          |          |
 *     |_________|__________|__________|
 *     |         |          |          |
 *     |         |          |          |
 *     |   01    |    11    |    21    |
 *     |         |          |          |
 *     |_________|__________|__________|
 *                          |          |
 *                          |          |
 *         00         10    |    20    |
 *                          |          |
 *                          |__________|
 *
 *     <-------------------><---------->
 *               c               a
 *     <------------------------------->
 *                     b
 *
 *
 *  Statuses:
 *  0 - inner                                  +
 *  1 - x = 0, y = (a, b) && x = c, y = (0, a) modified
 *  2 - x = b, y = (0, b)                      +
 *  3 - y = 0, x = (c, b) && y = a, x = (0, c) modified
 *  4 - y = b, x = (0, b)                      +
 *  5 - x = 0, y = a && x = c, y = 0           modified
 *  6 - x = b, y = 0                           +
 *  7 - x = 0. y = b                           +
 *  8 - x = b, y = b                           +
 *
 */

#define a 1.
#define c 2.
#define b 3.

void calc_mesh_params (int *st, double *X, double *Y, int *M0L, 
                       int *M0R, const UserDataCurr_struct * udc)
{

#if SQUARE
    int M1, M2;
    double hx,hy;
    // Minus 1 because its number of line-segments
    M1 = udc->Nx - 1;
    M2 = udc->Ny - 1;
    hx = udc->Hx;
    hy = udc->Hy;
    int j, j1, j2;

    st[0] = 5;
    M0L[0] = -1;
    M0R[0] = M1 + 1;
    X[0] = 0.;
    Y[0] = 0.;

    for(j1 = 1; j1 < M1; j1++)
      {
        st[j1] = 3;
        M0L[j1] = -1;
        M0R[j1] = M1 + j1 + 1;
        X[j1] = j1 * hx;
        Y[j1] = 0.;
      }

    st[M1] = 6;
    M0L[M1] = -1;
    M0R[M1] = 2 * M1 + 1;
    X[M1] = M1 * hx;
    Y[M1] = 0.;
    j = M1 + 1;

    for(j2 = 1; j2 < M2; j2++)
      {
        st[j] = 1;
        M0L[j] = j - M1 - 1;
        M0R[j] = j + M1 + 1;
        X[j] = 0.;
        Y[j] = j2 * hy;
        j++;

        for(j1 = 1; j1 < M1; j1++)
          {
            st[j] = 0;
            M0L[j] = j - M1 - 1;
            M0R[j] = j + M1+1;
            X[j] = j1 * hx;
            Y[j] = j2 * hy;
            j++;
          }

        st[j] = 2;
        M0L[j] = j - M1 - 1;
        M0R[j] = j + M1 + 1;
        X[j] = M1 * hx;
        Y[j] = j2 * hy;
        j++;
      }

    st[j] = 7;
    M0L[j] = j - M1 - 1;
    M0R[j] = -1;
    X[j] = 0;
    Y[j] = M2 * hy;
    j++;

    for(j1 = 1; j1 < M1; j1++)
      {
        st[j] = 4;
        M0L[j] = j - M1 - 1;
        M0R[j] = -1;
        X[j] = j1 * hx;
        Y[j] = M2 * hy;
        j++;
      }

    st[j] = 8;
    M0L[j] = j - M1 - 1;
    M0R[j] = -1;
    X[j] = M1 * hx;
    Y[j] = M2 * hy;
#else
  int M1, M2, M1_0, M2_0;
  double hx, hy;
  M1 = (udc->Nx - 1);
  M2 = (udc->Ny - 1);
  M1_0 = M1 - (udc->Nx_0 - 1);
  M2_0 = (udc->Ny_0 - 1);
  hx = udc->Hx;
  hy = udc->Hy;
  int j, j1, j2;

  // (x, y) = (c, 0)
  st[0] = 5;
  M0L[0] = -1;
  M0R[0] = M1_0 + 1;
  X[0] = c;
  Y[0] = 0.;

  // (X, Y) = {X = (c, b) & Y = {0}}
  for(j1 = 1; j1 < M1_0; j1++)
    {
      st[j1] = 3;
      M0L[j1] = -1;
      M0R[j1] = M1_0 + j1 + 1;
      X[j1] = c + j1 * hx;
      Y[j1] = 0.;
    }

  // (x, y) = (b, 0)
  st[M1_0] = 6;
  M0L[M1_0] = -1;
  M0R[M1_0] = 2 * M1_0 + 1;
  X[M1_0] = c + M1_0 * hx;
  Y[M1_0] = 0.;
  j = M1_0 + 1;

  for(j2 = 1; j2 < M2_0 - 1; j2++)
    {
      // (X, Y) = {X = {c} & Y = (0, a - 1)}
      st[j] = 1;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1_0 + 1;
      X[j] = c;
      Y[j] = j2 * hy;
      j++;
      // internal nodes
      for(j1 = 1; j1 < M1_0; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M1_0 - 1;
          M0R[j] = j + M1_0 + 1;
          X[j] = c + j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }
      // (X, Y) = {X = {b} & Y = (0, a - 1)}
      st[j] = 2;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1_0 + 1;
      X[j] = c + M1_0 * hx;
      Y[j] = j2 * hy;
      j++;
    }

  assert (j == (M1_0 + 1) * (M2_0 - 1));

  // last non-boundary layer of smaller rectangle
  st[j] = 1;
  M0L[j] = j - M1_0 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = c;
  Y[j] = (M2_0 - 1) * hy;
  j++;
  // internal nodes of last non-boundary layer of smaller rectangle
  for(j1 = 1; j1 < M1_0; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = c + j1 * hx;
      Y[j] = (M2_0 - 1) * hy;
      j++;
    }
  st[j] = 2;
  M0L[j] = j - M1_0 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = c + M1_0 * hx;
  Y[j] = (M2_0 - 1) * hy;
  j++;

  assert (j == (M1_0 + 1) * (M2_0));

  // (x, y) = (0, a)
  st[j] = 5;
  M0L[j] = -1;
  M0R[j] = j + M1 + 1;
  X[j] = 0.;
  Y[j] = a;
  j++;

  // (X, Y) = {X = (0, c) & Y = {a}}
  for(j1 = 1; j1 < (udc->Nx_0 - 1); j1++)
    {
      st[j] = 3;
      M0L[j] = -1;
      M0R[j] = j + M1 + 1;
      X[j] = j1 * hx;
      Y[j] = a;
      j++;
    }

  // (X, Y) = {X = [c, b) & Y = a}
  for(j1 = 0; j1 < M1_0; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = c + j1 * hx;
      Y[j] = a;
      j++;
    }

  // (x, y) = (b, a)
  st[j] = 2;
  M0L[j] = j - M1 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = b;
  Y[j] = a;
  j++;

  assert (j == (M1_0 + 1) * (M2_0) + (M1 + 1));

  for(j2 = M2_0 + 1; j2 < M2; j2++)
    {
      // (X, Y) = {X = {0} & Y = (a, b)}
      st[j] = 1;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = 0.;
      Y[j] = j2 * hy;
      j++;
      // internal nodes
      for(j1 = 1; j1 < M1; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M1 - 1;
          M0R[j] = j + M1 + 1;
          X[j] = j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }
      // (X, Y) = {X = {b} & Y = (0, a - 1)}
      st[j] = 2;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = M1 * hx;
      Y[j] = j2 * hy;
      j++;
    }

  assert (j == (M1_0 + 1) * (M2_0) + (M1 + 1) * (M2 - M2_0));

  // (x, y) = (0, b)
  st[j] = 7;
  M0L[j] = j - M1 - 1;
  M0R[j] = -1;
  X[j] = 0.;
  Y[j] = M2 * hy;
  j++;

  for(j1 = 1; j1 < M1; j1++)
    {
      st[j] = 4;
      M0L[j] = j - M1 - 1;
      M0R[j] = -1;
      X[j] = j1 * hx;
      Y[j] = M2 * hy;
      j++;
    }

  // (x, y) = (b, b)
  st[j] = 8;
  M0L[j] = j - M1 - 1;
  M0R[j] = -1;
  X[j] = M1 * hx;
  Y[j] = M2 * hy;
  j++;
  assert (j == (M1_0 + 1) * (M2_0) + (M1 + 1) * (M2 + 1 - M2_0));
  assert (j == udc->N);
#endif
}
#undef a
#undef c
#undef b
