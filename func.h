#define FIX_UNUSED(X) (void)(X)
#define FREE_ARRAY(X) if ((X)) free ((X)); (X) = NULL
#define CLOSE_FILE(X) if ((X)) fclose ((X)); (X) = NULL

#define MINIMAL_FOR_COMPARE   1.e-16
#define STAT_SOL_EPS          7 * 1.e-7

#define RHO_0   0.5
#define RHO_G   0.3
#define W       0.1

#define NO_SMOOTH      1
#define STAT_SOL_SRCH  0
#define EIG_FUNC_INIT  1
#define SQUARE         1

#define LIGHT_FUNCS   0
#define LIGHT_G       0
#define LIGHT_U1      0
#define LIGHT_U2      0


// Working area is a rectangle which is located in a first quadrant.
// Its sides are parallel to axes OXY and one vertex is located in (0, 0).
// Also it has a cut -- rectangle which sides are parallel to sides of the main one.
typedef struct
{
  double Segm_T;
  double Segm_X;    // x-coord of farthest vertex of main rect
  double Segm_Y;    // y-coord of farthest vertex of main rect
  double Segm_X_0;  // x-coord of farthest vertex of cut rect
  double Segm_Y_0;  // y-coord of farthest vertex of cut rect
  double p_ro;      // p_ro = p'_ro * ro
  double mu;
} P_dif;

// initialization of working area params
void param_dif (P_dif *p_d);

typedef struct
{
  int M_x;    // M_x   * h_x = Segm_X
  int M_y;    // M_y   * h_y = Segm_Y
  int M_x_0;  // M_x_0 * h_x = Segm_X_0
  int M_y_0;  // M_y_0 * h_y = Segm_Y_0
  int N;      // N     * tau = Segm_T
  int Dim;    // number of nodes
  double h_x;
  double h_y;
  double tau;
  double eta;
} P_she;

void param_she_step (P_she *p_s, P_dif *p_d, int it_t, int int_sp);

double Norm_c (double *A, int Dim, double *X, double *Y, double t, double (*f) (double tt, double x1, double x2));

double Norm_l2 (double *A, int Dim, double *X, double *Y, double t, double (*f) (double tt, double x1, double x2));

double ro (double t, double x, double y);
double gg (double t, double x, double y);
double p (double t, double x, double y, double p_ro);
double u1 (double t, double x, double y);
double u2 (double t, double x, double y);
double Func_g (double t, double x, double y);
double Func_v1 (double t, double x, double y, double p_ro, double mu);
double Func_v2 (double t, double x, double y, double p_ro, double mu);

void Setka (int *st, double *X, double *Y, int *M0L, int *M0R, P_she *p_s, P_dif *p_d);

int  Sxema (double *G, double *V1, double *V2,
            double *G_prev, double *V1_prev, double *V2_prev,
            const double *G_eigen, const double *V1_eigen, const double *V2_eigen,
            const double *G_stat, const double *V1_stat, const double *V2_stat,
            int *st, double *X, double *Y, int *M0L,
            int *M0R, P_she *p_s, P_dif *p_d);

double calc_sol_residual_norm (int Dim, const double *G, const double *V1,
                               const double *V2, const double *G_prev,
                               const double *V1_prev, const double *V2_prev);

void init_prev_with_curr (int Dim, double *G, double *V1, double *V2,
                          double *G_prev, double *V1_prev, double *V2_prev);

int is_equal (double x1, double x2);
