/* ########################
*/

#define n_grids 0
#define m_grids 0
#define equa 11

extern void init_diff_const(void);
double D[equa];
extern double reac(double*, int);
extern void eval(double *, double *);
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]);
extern double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l);
void calc_jac2(double S[equa], double Jac[equa][equa]);
