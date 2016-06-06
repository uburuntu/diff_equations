#include "f.h"
#include "../config.h"

void calc_mesh_params_square (int *st, double *X, double *Y, int *M0L, int *M0R, const user_data *ud)
{
  int M1, M2;
  double hx, hy;
  int j, j1, j2;

  // Minus 1 because it`s a number of line-segments
  M1 = ud->Nx - 1;
  M2 = ud->Ny - 1;
  hx = ud->Hx;
  hy = ud->Hy;

  st[0] = 5;
  M0L[0] = -1;
  M0R[0] = M1 + 1;
  X[0] = 0.;
  Y[0] = 0.;

  for (j1 = 1; j1 < M1; j1++)
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

  for (j2 = 1; j2 < M2; j2++)
    {
      st[j] = 1;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = 0.;
      Y[j] = j2 * hy;
      j++;

      for (j1 = 1; j1 < M1; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M1 - 1;
          M0R[j] = j + M1 + 1;
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

  for (j1 = 1; j1 < M1; j1++)
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
}

void calc_mesh_params_volodya_9 (int *st, double *X, double *Y, int *M0L, int *M0R, const user_data *ud)
{
  int M1, M2, M1_0, M2_0;
  int j, j1, j2;
  double hx, hy;

  M1 = (ud->Nx - 1);
  M2 = (ud->Ny - 1);
  M1_0 = M1 - (ud->Nx_0 - 1);
  M2_0 = (ud->Ny_0 - 1);
  hx = ud->Hx;
  hy = ud->Hy;

  // (x, y) = (c, 0)
  st[0] = 5;
  M0L[0] = -1;
  M0R[0] = M1_0 + 1;
  X[0] = C_LENGHT;
  Y[0] = 0.;

  // (X, Y) = {X = (c, b) & Y = {0}}
  for (j1 = 1; j1 < M1_0; j1++)
    {
      st[j1] = 3;
      M0L[j1] = -1;
      M0R[j1] = M1_0 + j1 + 1;
      X[j1] = C_LENGHT + j1 * hx;
      Y[j1] = 0.;
    }

  // (x, y) = (b, 0)
  st[M1_0] = 6;
  M0L[M1_0] = -1;
  M0R[M1_0] = 2 * M1_0 + 1;
  X[M1_0] = C_LENGHT + M1_0 * hx;
  Y[M1_0] = 0.;
  j = M1_0 + 1;

  for (j2 = 1; j2 < M2_0 - 1; j2++)
    {
      // (X, Y) = {X = {c} & Y = (0, a - 1)}
      st[j] = 1;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1_0 + 1;
      X[j] = C_LENGHT;
      Y[j] = j2 * hy;
      j++;

      // internal nodes
      for (j1 = 1; j1 < M1_0; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M1_0 - 1;
          M0R[j] = j + M1_0 + 1;
          X[j] = C_LENGHT + j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }

      // (X, Y) = {X = {b} & Y = (0, a - 1)}
      st[j] = 2;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1_0 + 1;
      X[j] = C_LENGHT + M1_0 * hx;
      Y[j] = j2 * hy;
      j++;
    }

  assert (j == (M1_0 + 1) * (M2_0 - 1));

  // last non-boundary layer of smaller rectangle
  st[j] = 1;
  M0L[j] = j - M1_0 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = C_LENGHT;
  Y[j] = (M2_0 - 1) * hy;
  j++;

  // internal nodes of last non-boundary layer of smaller rectangle
  for (j1 = 1; j1 < M1_0; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M1_0 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = C_LENGHT + j1 * hx;
      Y[j] = (M2_0 - 1) * hy;
      j++;
    }

  st[j] = 2;
  M0L[j] = j - M1_0 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = C_LENGHT + M1_0 * hx;
  Y[j] = (M2_0 - 1) * hy;
  j++;

  assert (j == (M1_0 + 1) * (M2_0));

  // (x, y) = (0, a)
  st[j] = 5;
  M0L[j] = -1;
  M0R[j] = j + M1 + 1;
  X[j] = 0.;
  Y[j] = A_LENGHT;
  j++;

  // (X, Y) = {X = (0, c) & Y = {a}}
  for (j1 = 1; j1 < (ud->Nx_0 - 1); j1++)
    {
      st[j] = 3;
      M0L[j] = -1;
      M0R[j] = j + M1 + 1;
      X[j] = j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  // (X, Y) = {X = [c, b) & Y = a}
  for (j1 = 0; j1 < M1_0; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = C_LENGHT + j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  // (x, y) = (b, a)
  st[j] = 2;
  M0L[j] = j - M1 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = B_LENGHT;
  Y[j] = A_LENGHT;
  j++;

  assert (j == (M1_0 + 1) * (M2_0) + (M1 + 1));

  for (j2 = M2_0 + 1; j2 < M2; j2++)
    {
      // (X, Y) = {X = {0} & Y = (a, b)}
      st[j] = 1;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = 0.;
      Y[j] = j2 * hy;
      j++;

      // internal nodes
      for (j1 = 1; j1 < M1; j1++)
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

  for (j1 = 1; j1 < M1; j1++)
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
  assert (j == ud->N);
}


void calc_mesh_params_ramzan_10 (int *st, double *X, double *Y, int *M0L, int *M0R, const user_data *ud)
{
  int j = 0, j1, j2;
  double hx, hy;

  int M_block = ud->Nx_0 - 1;
  int M_block_wide = 2 * M_block;

  int M_x = ud->Nx - 1;
  int M_y = ud->Ny - 1;

  hx = ud->Hx;
  hy = ud->Hy;

  // 1 step
  st[0] = 5;
  M0L[0] = -1;
  M0R[0] = j + M_block_wide + 1;
  X[0] = A_LENGHT;
  Y[0] = 0.;
  j++;

  for (j1 = 1; j1 < M_block_wide; j1++)
    {
      st[j1] = 3;
      M0L[j1] = -1;
      M0R[j1] = j + M_block_wide + 1;
      X[j1] = A_LENGHT + j1 * hx;
      Y[j1] = 0.;
      j++;
    }

  st[M_block_wide] = 6;
  M0L[M_block_wide] = -1;
  M0R[M_block_wide] = j + M_block_wide + 1;
  X[M_block_wide] = B_LENGHT;
  Y[M_block_wide] = 0.;
  j++;

  // 2 step
  for (j2 = 1; j2 < M_block - 1; j2++)
    {
      st[j] = 1;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = j + M_block_wide + 1;
      X[j] = A_LENGHT;
      Y[j] = j2 * hy;
      j++;

      for (j1 = 1; j1 < M_block_wide; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M_block_wide - 1;
          M0R[j] = j + M_block_wide + 1;
          X[j] = A_LENGHT + j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }

      st[j] = 2;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = j + M_block_wide + 1;
      X[j] = B_LENGHT;
      Y[j] = j2 * hy;
      j++;
    }

  st[j] = 1;
  M0L[j] = j - M_block_wide - 1;
  M0R[j] = j + M_x + 1;
  X[j] = A_LENGHT;
  Y[j] = (M_block - 1) * hy;
  j++;

  for (j1 = 1; j1 < M_block_wide; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = j + M_x + 1;
      X[j] = A_LENGHT + j1 * hx;
      Y[j] = (M_block - 1) * hy;
      j++;
    }

  st[j] = 2;
  M0L[j] = j - M_block_wide - 1;
  M0R[j] = j + M_x + 1;
  X[j] = B_LENGHT;
  Y[j] = (M_block - 1) * hy;
  j++;

  assert (j == (M_block_wide + 1) * M_block);

  // 3 step

  st[j] = 5;
  M0L[j] = -1;
  M0R[j] = j + M_x + 1;
  X[j] = 0.;
  Y[j] = A_LENGHT;
  j++;

  for (j1 = 1; j1 < M_block; j1++)
    {
      st[j] = 3;
      M0L[j] = -1;
      M0R[j] = j + M_x + 1;
      X[j] = j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  for (j1 = 0; j1 <= M_block; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M_x - 1;
      M0R[j] = j + M_x + 1;
      X[j] = A_LENGHT + j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  for (j1 = 1; j1 < M_block; j1++)
    {
      st[j] = 4;
      M0L[j] = j - M_x - 1;
      M0R[j] = -1;
      X[j] = C_LENGHT + j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  st[j] = 8;
  M0L[j] = j - M_x - 1;
  M0R[j] = -1;
  X[j] = B_LENGHT;
  Y[j] = A_LENGHT;
  j++;

  assert (j == ((M_block_wide + 1) * M_block + (M_x + 1)));

  // 4 step

  st[j] = 1;
  M0L[j] = j - M_x - 1;
  M0R[j] = j + M_block_wide + 1;
  X[j] = 0.;
  Y[j] = (M_block + 1) * hy;
  j++;

  for (j1 = 1; j1 < M_block_wide; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M_x - 1;
      M0R[j] = j + M_block_wide + 1;
      X[j] = j1 * hx;
      Y[j] = (M_block + 1) * hy;
      j++;
    }

  st[j] = 2;
  M0L[j] = j - M_x - 1;
  M0R[j] = j + M_block_wide + 1;
  X[j] = M_block_wide * hx;
  Y[j] = (M_block + 1) * hy;
  j++;

  for (j2 = M_block + 2; j2 < M_y; j2++)
    {
      st[j] = 1;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = j + M_block_wide + 1;
      X[j] = 0.;
      Y[j] = j2 * hy;
      j++;

      for (j1 = 1; j1 < M_block_wide; j1++)
        {
          st[j] = 0;
          M0L[j] = j - M_block_wide - 1;
          M0R[j] = j + M_block_wide + 1;
          X[j] = j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }

      st[j] = 2;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = j + M_block_wide + 1;
      X[j] = M_block_wide * hx;
      Y[j] = j2 * hy;
      j++;
    }

  // 5 step
  st[j] = 7;
  M0L[j] = j - M_block_wide - 1;
  M0R[j] = -1;
  X[j] = 0.;
  Y[j] = M_y * hy;
  j++;

  for (j1 = 1; j1 < M_block_wide; j1++)
    {
      st[j] = 4;
      M0L[j] = j - M_block_wide - 1;
      M0R[j] = -1;
      X[j] = j1 * hx;
      Y[j] = M_y * hy;
      j++;
    }

  st[j] = 8;
  M0L[j] = j - M_block_wide - 1;
  M0R[j] = -1;
  X[j] = M_block_wide * hx;
  Y[j] = M_y * hy;
  j++;

  assert (j == ud->N);
}


void calc_mesh_params_nastya_11 (int *st, double *X, double *Y, int *M0L, int *M0R, const user_data *ud)
{
  int M1, M2;
  double hx, hy;
  int j, j1, j2;

  int MA = ud->Nx_0 - 1; // = M1 / 3;
  int MC = 2 * MA;

  M1 = ud->Nx - 1;
  M2 = ud->Ny - 1 - MA; // without upper layers
  hx = ud->Hx;
  hy = ud->Hy;

  // (x, y) = (0, 0)
  st[0] = 5;
  M0L[0] = -1;
  M0R[0] = MA + MA + 1;
  X[0] = 0.;
  Y[0] = 0.;

  // (X, Y) = {X = (0, a) & Y = {0}}
  for (j1 = 1; j1 < MA; j1++)
    {
      st[j1] = 3;
      M0L[j1] = -1;
      M0R[j1] = MA + MA + j1 + 1;
      X[j1] = j1 * hx;
      Y[j1] = 0.;
    }

  // (x, y) = (a, 0)
  st[MA] = 6;
  M0L[MA] = -1;
  M0R[MA] = MA + MA + 1;
  X[MA] = A_LENGHT;
  Y[MA] = 0.;

  j = MA + 1;
  // (x, y) = (c, 0)
  st[j] = 5;
  M0L[j] = -1;
  M0R[j] = MC + 1;
  X[j] = C_LENGHT;
  Y[j] = 0.;

  // (X, Y) = {X = (c, b) & Y = {0}}
  for (j1 = 1; j1 < MA; j1++)
    {
      st[j + j1] = 3;
      M0L[j + j1] = -1;
      M0R[j + j1] = MC + j1 + 1;
      X[j + j1] = C_LENGHT + j1 * hx;
      Y[j + j1] = 0.;
    }

  // (x, y) = (b, 0)
  st[MC + 1] = 6;
  M0L[MC + 1] = -1;
  M0R[MC + 1] = MC + 1;
  X[MC + 1] = C_LENGHT + MA * hx;
  Y[MC + 1] = 0.;
  j = MC + 2;

  //__________________________________________

  // Y = (0, a - 1)
  for (j2 = 1; j2 < MA - 1; j2++)
    {
      // (X, Y) = {X = {0} & Y = (0, a - 1)}
      st[j] = 1;
      M0L[j] = j - MC - 1;
      M0R[j] = j + MC + 1;
      X[j] = 0;
      Y[j] = j2 * hy;
      j++;

      // internal nodes
      for (j1 = 1; j1 < MA; j1++)
        {
          st[j] = 0;
          M0L[j] = j - MC - 1;
          M0R[j] = j + MC + 1;
          X[j] = j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }

      // (X, Y) = {X = {a} & Y = (0, a - 1)}
      st[j] = 2;
      M0L[j] = j - MC - 1;
      M0R[j] = j + MC + 1;
      X[j] = MA * hx;
      Y[j] = j2 * hy;
      j++;
      // -------------------------------------
      // (X, Y) = {X = {c} & Y = (0, a - 1)}
      st[j] = 1;
      M0L[j] = j - MC - 1;
      M0R[j] = j + MC + 1;
      X[j] = C_LENGHT;
      Y[j] = j2 * hy;
      j++;

      // internal nodes
      for (j1 = 1; j1 < MA; j1++)
        {
          st[j] = 0;
          M0L[j] = j - MC - 1;
          M0R[j] = j + MC + 1;
          X[j] = C_LENGHT + j1 * hx;
          Y[j] = j2 * hy;
          j++;
        }

      // (X, Y) = {X = {b} & Y = (0, a - 1)}
      st[j] = 2;
      M0L[j] = j - MC - 1;
      M0R[j] = j + MC + 1;
      X[j] = C_LENGHT + MA * hx;
      Y[j] = j2 * hy;
      j++;
    }

  assert (j == 2 * (MA + 1) * (MA - 1));
  // ----------------------- 2:10 -------------------------
  // last non-boundary layer of smaller rectangle x = 0
  st[j] = 1;
  M0L[j] = j - MC - 1;
  M0R[j] = j + M1 + 1;
  X[j] = 0;
  Y[j] = (MA - 1) * hy;
  j++;

  // internal nodes of last non-boundary layer of smaller rectangle
  for (j1 = 1; j1 < MA; j1++)
    {
      st[j] = 0;
      M0L[j] = j - MC - 1;
      M0R[j] = j + M1 + 1;
      X[j] = j1 * hx;
      Y[j] = (MA - 1) * hy;
      j++;
    }

  st[j] = 2;
  M0L[j] = j - MC - 1;
  M0R[j] = j + M1 + 1;
  X[j] = MA * hx;
  Y[j] = (MA - 1) * hy;
  j++;

  // last non-boundary layer of smaller rectangle x = C
  st[j] = 1;
  M0L[j] = j - MC - 1;
  M0R[j] = j + M1 + 1;
  X[j] = C_LENGHT;
  Y[j] = (MA - 1) * hy;
  j++;

  // internal nodes of last non-boundary layer of smaller rectangle
  for (j1 = 1; j1 < MA; j1++)
    {
      st[j] = 0;
      M0L[j] = j - MC - 1;
      M0R[j] = j + M1 + 1;
      X[j] = C_LENGHT + j1 * hx;
      Y[j] = (MA - 1) * hy;
      j++;
    }

  st[j] = 2;
  M0L[j] = j - MC - 1;
  M0R[j] = j + M1 + 1;
  X[j] = C_LENGHT + MA * hx;
  Y[j] = (MA - 1) * hy;
  j++;

  assert (j == 2 * (MA + 1) * (MA));
  // ----------------------- 2:58 -------------------------

  // (x, y) = (0, a)
  st[j] = 1;
  M0L[j] = j - MC - 1;
  M0R[j] = j + M1 + 1;
  X[j] = 0.;
  Y[j] = A_LENGHT;
  j++;

  // (X, Y) = {X = (0, a] & Y = {a}}
  for (j1 = 1; j1 <= MA; j1++)
    {
      st[j] = 0;
      M0L[j] = j - MC - 1;
      M0R[j] = j + M1 + 1;
      X[j] = j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  // (X, Y) = {X = (a, c) & Y = {a}}
  for (j1 = 1; j1 < MA; j1++)
    {
      st[j] = 3;
      M0L[j] = -1;
      M0R[j] = j + M1 + 1;
      X[j] = A_LENGHT + j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  // (X, Y) = {X = [c, b) & Y = a}
  for (j1 = 0; j1 < MA; j1++)
    {
      st[j] = 0;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = C_LENGHT + j1 * hx;
      Y[j] = A_LENGHT;
      j++;
    }

  // (x, y) = (b, a)
  st[j] = 2;
  M0L[j] = j - M1 - 1;
  M0R[j] = j + M1 + 1;
  X[j] = B_LENGHT;
  Y[j] = A_LENGHT;
  j++;

  assert (j == 2 * (MA + 1) * (MA) + (M1 + 1));

  // ----------------------- 2:30 -------------------------
  for (j2 = MA + 1; j2 < M2; j2++)
    {
      // (X, Y) = {X = {0} & Y = (a, b)}
      st[j] = 1;
      M0L[j] = j - M1 - 1;
      M0R[j] = j + M1 + 1;
      X[j] = 0.;
      Y[j] = j2 * hy;
      j++;

      // internal nodes
      for (j1 = 1; j1 < M1; j1++)
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

  assert (j == 2 * (MA + 1) * (MA) + (M1 + 1) * (M2 - MA));
  assert (j == 2 * (MA + 1) * (MA) + (M1 + 1) * MA);
  // ----------------------- 3:00 -------------------------

  // (x, y) = (0, b)
  st[j] = 7;
  M0L[j] = j - M1 - 1;
  M0R[j] = -1;
  X[j] = 0.;
  Y[j] = M2 * hy;
  j++;

  for (j1 = 1; j1 < M1; j1++)
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
  assert (j == 2 * (MA + 1) * (MA) + (M1 + 1) * (M2 + 1 - MA));
  assert (j == 2 * (MA + 1) * (MA) + (M1 + 1) * (MA + 1));
  assert (j == ud->N);
}
