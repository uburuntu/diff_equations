#include "f.h"

int read_stationary_solution (const char *fname, int N, double *G, double *V1, double *V2)
{
  int i, readedN;
  FILE *input_file = fopen (fname, "r");

  if (!input_file)
    {
      printf ("fopen error: cannot open %s\n", fname);
      return -1;
    }

  if (!fscanf (input_file, "%d", &readedN) || readedN != N)
    {
      printf ("fread error: incorrect N by file %s\n", fname);
      fclose (input_file);
      return -1;
    }

  readedN = 0;

  for (i = 0; i < N; i++)
    {
      if (!fscanf (input_file, "%lf ", G + i))
        {
          printf ("fread error: incorrect G filling from %s\n", fname);
        }

      readedN++;
    }

  for (i = 0; i < N; i++)
    {
      if (!fscanf (input_file, "%lf ", V1 + i))
        {
          printf ("fread error: incorrect V1 filling from %s\n", fname);
        }

      readedN++;
    }

  for (i = 0; i < N; i++)
    {
      if (!fscanf (input_file, "%lf ", V2 + i))
        {
          printf ("fread error: incorrect V2 filling from %s\n", fname);
        }

      readedN++;
    }

  fclose (input_file);
  printf ("Totally readed %d numbers from file %s\n", readedN, fname);

  return 0;
}

void print_2dfun_double (FILE *f, const double *u, const int n)
{
  int i;

  for (i = 0; i < n; i++)
    {
      fprintf (f, "%9.2e ", u[i]);
    }

  fprintf (f, "\n");
  return;
}

int *make_vector_int (
  int n,
  const char *info_1,
  const char *info_2)
{
  unsigned int size = n * sizeof (int);

  int *tmp = (int *)malloc (size);

  if (tmp == NULL)
    {
      printf ("Error in %s %s : Fail to allocate %u bytes\n",
              info_1, info_2, size);
      exit (1);
    }

  memset (tmp, 0, size);

  return tmp;
}

double *make_vector_double (
  int n,
  const char *info_1,
  const char *info_2)
{
  unsigned int size = n * sizeof (double);

  double *tmp = (double *)malloc (size);

  if (tmp == NULL)
    {
      printf ("Error in %s %s : Fail to allocate %u bytes\n",
              info_1, info_2, size);
      exit (1);
    }

  memset (tmp, 0, size);

  return tmp;
}
