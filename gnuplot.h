
#include <stdio.h>
#include <math.h>

#define DIVISOR     5
#define OUTPUT      "results.txt"
#define OUTTEX      "theplot.tex"

char *plot_name (char *name, double tau, double h1, double h2, int j);

char *tex_name (char *name, double tau, double h1, double h2,  int t);

void print_plot (char* file_name, double *X, double *Y,
                 double *G, double *V1, double *V2,
                 int size, double tt);

int make_graph (char* texname, char* plotname, double h1, double h2,
                double tau, double t);

void printhead(FILE *fout);
void printtail(FILE *fout);
