#include "tabtex.h"
#include <stdio.h>
#include "gnuplot.h"
#include "func.h"

void tabtex_nc_g (int it_t_max, int it_sp_max, double *nc_g, double *tauit,
                  double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("gc.tex", "a");
  fprintf (fi1, "\\documentstyle{article}\n");
  fprintf (fi1, "\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;C\\;for\\;g: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf \\\\ $\n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nc_g[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  fprintf (fi1, "\\end{document}\n");
  fclose (fi1);
}

void tabtex_nc_v1 (int it_t_max, int it_sp_max, double *nc_v1, double *tauit,
                   double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("v1c.tex", "a");
  // fprintf(fi1,"\\documentstyle{article}\n");
  // fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;C\\;for\\;v1: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nc_v1[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  //  fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void tabtex_nc_v2 (int it_t_max, int it_sp_max, double *nc_v2, double *tauit,
                   double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("v2c.tex", "a");
  //  fprintf(fi1,"\\documentstyle{article}\n");
  // fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;C\\;for\\;v2: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nc_v2[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  // fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void tabtex_nl2_g (int it_t_max, int it_sp_max, double *nl2_g, double *tauit,
                   double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("gl2.tex", "a");
  // fprintf(fi1,"\\documentstyle{article}\n");
  //fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;L_2\\;for\\;g: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nl2_g[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  //  fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void tabtex_nl2_v1 (int it_t_max, int it_sp_max, double *nl2_v1, double *tauit,
                    double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("v1l2.tex", "a");
  //  fprintf(fi1,"\\documentstyle{article}\n");
  // fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;L_2\\;for\\;v1: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nl2_v1[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  //  fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void tabtex_nl2_v2 (int it_t_max, int it_sp_max, double *nl2_v2, double *tauit,
                    double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("v2l2.tex", "a");
  //  fprintf(fi1,"\\documentstyle{article}\n");
  // fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $norm\\;of\\;the\\;error\\;in\\;L_2\\;for\\;v2: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", nl2_v2[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  //  fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void tabtex_time (int it_t_max, int it_sp_max, double *time, double *tauit,
                  double p_ro, double mu)
{
  int it_t, it_sp, it;
  FILE *fi1 = fopen ("time.tex", "a");
  // fprintf(fi1,"\\documentstyle{article}\n");
  //fprintf(fi1,"\\begin{document}\n");
  fprintf (fi1, "\n $ time, p_{\\rho}=%.3lf, \\mu = %.3lf $ \\\\ \n", p_ro, mu);
  fprintf (fi1, "\\begin{tabular}{|p{0.6in}|p{0.7in}|p{0.7in}|p{0.7in}|p{0.7in}|} \\hline\n");
  fprintf (fi1, "$\\tau\\setminus h$ & $0.05$ & 0.025& 0.0125 & 0.00625 \\\\ \\hline");
  it = 0;

  for (it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf (fi1, "\n");
      fprintf (fi1, "$%.5lf$ & ", tauit[it_t]);

      for (it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf (fi1, "$%.3le$ ", time[it]);
          it++;

          if (it_sp < it_sp_max)
            {
              fprintf (fi1, "&");
            }
        }

      fprintf (fi1, " \\\\ \\hline");
    }

  fprintf (fi1, "\n\\end{tabular}\\\\[20pt]\n");
  //  fprintf(fi1,"\\end{document}\n");
  fclose (fi1);
}

void make_tabletex()
{
  FILE *fout;

  fout = fopen (TABLETEX_SMOOTH, "a+");

  if (fout != NULL)
    {
      printhead (fout);
      fprintf (fout, "\\section* {Таблицы:}  \n");
      fprintf (fout, "\\input{v1c.tex}  \n");
      fprintf (fout, "\\input{v1l2.tex} \n");
      fprintf (fout, "\\input{v2c.tex}  \n");
      fprintf (fout, "\\input{v2l2.tex} \n");
      fprintf (fout, "\\input{time.tex} \n");
      printtail (fout);
      fclose (fout);
    }
  else
    {
      printf ("Can't open texfile\n");
    }
}
