#include "f.h"

int read_stationary_solution (const char *fname, int N, double *G, double *V1, double *V2)
{
    int i, readedN;
    FILE *input_file = fopen(fname, "r");

    if (!input_file)
      {
        printf ("fopen error: cannot open %s\n", fname);
        return -1;
      }

    if (!fscanf (input_file, "%d", &readedN))
      {
        printf ("fread error: incorrect N by file %s\n", fname);
        fclose (input_file);
        return -1;
      }

    if (readedN != N)
      {
        printf ("fread error: incorrect N=%d by file %s\n", readedN, fname);
        fclose (input_file);
        return -1;
      }

    for (i = 0; i < N; i++)
      {
        fscanf(input_file, "%lf ", G + i);
        readedN++;
      }

    for (i = 0; i < N; i++)
      {
        fscanf(input_file, "%lf ", V1 + i);
        readedN++;
      }

    for (i = 0; i < N; i++)
      {
        fscanf(input_file, "%lf ", V2 + i);
        readedN++;
      }

  //  for (i = 0; i < 3 * N; i++)
  //    {
  //      if (!fscanf (input_file, "%lf", i % N + (i / N == 0 ? G : (i / N == 1 ? V1 : V2))))
  //        {
  //          printf ("fread error: incorrect data by %d in line %d by file %s\n", i % N, i / N, fname);
  //          fclose (input_file);
  //          return -1;
  //        }
  //      readedN++;
  //    }
    fclose (input_file);
    printf ("Totally readed %d numbers from file %s\n", readedN, fname);
    printf ("%d data missing.\n", 3 * N - readedN);

  return 0;
}

void print_2dfun_double (FILE* f, const char * name, const double  * u,
                         const int nx, const int ny)
{
  int i;
  int j;

  fprintf(f,"#%s\n",name);

  for(j = ny - 1; j >= 0; j--)
    {
      for(i = 0; i < nx; i++)
        {
          fprintf(f, "%9.2e ", u[i + j * nx]);
        }
      fprintf(f, "\n");
    }

  return;
}

int *sp_alloc_i_vector (
    int n,
    const char *info_1,
    const char *info_2)
{
  int *tmp = (int *)malloc(n * sizeof(int ));

  if (tmp == NULL)
    {
      printf("Error in %s %s : Fail to allocate %lu bytes\n",
             info_1, info_2, n * sizeof(int));
      exit(1);
    }

  memset (tmp, 0, n * sizeof(int));

  return tmp;
}

void sp_free_i_vector (int * d)
{
  if (d)
    free(d);
}

double *sp_alloc_d_vector (
    int n,
    const char *info_1,
    const char *info_2)
{
  double *tmp = (double *)malloc(n * sizeof(double));

  if (tmp == NULL)
    {
      printf ("Error in %s %s : Fail to allocate %lu bytes\n",
              info_1, info_2, n * sizeof(double));
      exit (1);
    }

  memset (tmp, 0, n * sizeof(double));

  return tmp;
}

void sp_free_d_vector (double * d)
{
  if (d)
    free(d);
}


double *make_vector_double (int n, const char *info_1, const char *info_2)
{
  double *u;
  int i;

  u = (double *) malloc (sizeof(double) * n);
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

  u = (int *) malloc (sizeof(int) * n);
  if (u == NULL)
    {
      printf ("Error in %s %s: Fail to allocate %lu bytes\n",
              info_1, info_2, sizeof(int) * n);
      exit(1);
    }

  for (i = 0; i < n; i++)
    u[i] = 0.;

  return u;
}
