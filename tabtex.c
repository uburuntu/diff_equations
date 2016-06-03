#include "tabtex.h"
#include <stdio.h>
#include "gnuplot.h"
#include "func.h"

void tabtex(const char *filename, const char *message, int it_t_max,int it_sp_max,double *  data,double *tauit,
                 double p_ro, double mu, int standalone)
{
  int it_t,it_sp,it;
  FILE *fi1 = fopen(filename,"w");
  if(standalone)
    {
      fprintf(fi1,"\\documentclass[a4paper,12pt]{scrartcl}\n");
      fprintf(fi1,"\\usepackage[warn]{mathtext}\n");
      fprintf(fi1,"\\usepackage[T2A]{fontenc}\n");
      fprintf(fi1,"\\usepackage[utf8]{inputenc}\n");
      fprintf(fi1,"\\usepackage[english,russian]{babel}\n");
      fprintf(fi1,"\\usepackage{amsmath}\n");
      fprintf(fi1,"\\begin{document}\n");
    }

  fprintf(fi1, "$\n \\text{%s}: \\quad p_{\\rho}=%.3lf, \\mu = %.3lf \\\\ $\n", message, p_ro, mu);
  fprintf(fi1,"\\begin{tabular}{|p{0.6in}|p{1.2in}|p{1.2in}|p{1.2in}|} \\hline\n");
  fprintf(fi1,"$\\tau\\setminus h$ & $0.05000$ & $0.02500$& $0.01250$ \\\\ \\hline");
  it = 0;
  for(it_t = 0; it_t <= it_t_max; it_t++)
    {
      fprintf(fi1,"\n");
      fprintf(fi1,"$%.5lf$ & ",tauit[it_t]);
      for(it_sp = 0; it_sp <= it_sp_max; it_sp++)
        {
          fprintf(fi1, "$%.3le$ ",  data[it]);
          it++;
          if(it_sp < it_sp_max)
            fprintf(fi1,"&");
        }
      fprintf(fi1," \\\\ \\hline");
    }

    fprintf(fi1,"\n\\end{tabular}\\\\[20pt]\n");

  if(standalone)
    fprintf(fi1,"\\end{document}\n");

  fclose(fi1);
}

void make_tabletex (void)
{
  FILE *fout;

  fout = fopen (TABLETEX_SMOOTH, "a+");

  if (fout != NULL)
    {
      printhead (fout);
      fprintf (fout, "\\section* {Таблицы:}  \n");
      fprintf (fout, "\\input{gc.tex}  \n");
      fprintf (fout, "\\input{gl2.tex} \n");
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
