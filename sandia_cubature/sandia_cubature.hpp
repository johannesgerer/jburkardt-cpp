namespace webbur 
{
double c1_geg_monomial_integral ( double alpha, int expon );
double c1_jac_monomial_integral ( double alpha, double beta, int expon );
double c1_leg_monomial_integral ( int expon );
void cn_geg_01_1 ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_01_1_size ( int n, double alpha );
void cn_geg_02_xiu ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_02_xiu_size ( int n, double alpha );
void cn_geg_03_xiu ( int n, double alpha, int o, double x[], double w[] );
int cn_geg_03_xiu_size ( int n, double alpha );
double cn_geg_monomial_integral ( int n, double alpha, int expon[] );
void cn_jac_01_1 ( int n, double alpha, double beta, int o, double x[], 
  double w[] );
int cn_jac_01_1_size ( int n, double alpha, double beta );
double cn_jac_monomial_integral ( int n, double alpha, double beta, 
  int expon[] );
void cn_jac_02_xiu ( int n, double alpha, double beta, int o, double x[], 
  double w[] );
int cn_jac_02_xiu_size ( int n, double alpha, double beta );
void cn_leg_01_1 ( int n, int o, double x[], double w[] );
int cn_leg_01_1_size ( int n );
void cn_leg_02_xiu ( int n, int o, double x[], double w[] );
int cn_leg_02_xiu_size ( int n );
void cn_leg_03_1 ( int n, int o, double x[], double w[] );
int cn_leg_03_1_size ( int n );
void cn_leg_03_xiu ( int n, int o, double x[], double w[] );
int cn_leg_03_xiu_size ( int n );
void cn_leg_05_1 ( int n, int option, int o, double x[], double w[] );
int cn_leg_05_1_size ( int n );
void cn_leg_05_2 ( int n, int o, double x[], double w[] );
int cn_leg_05_2_size ( int n );
double cn_leg_monomial_integral ( int n, int expon[] );
void en_her_01_1 ( int n, int o, double x[], double w[] );
int en_her_01_1_size ( int n );
void en_her_02_xiu ( int n, int o, double x[], double w[] );
int en_her_02_xiu_size ( int n );
void en_her_03_1 ( int n, int o, double x[], double w[] );
int en_her_03_1_size ( int n );
void en_her_03_xiu ( int n, int o, double x[], double w[] );
int en_her_03_xiu_size ( int n );
void en_her_05_1 ( int n, int option, int o, double x[], double w[] );
int en_her_05_1_size ( int n );
void en_her_05_2 ( int n, int o, double x[], double w[] );
int en_her_05_2_size ( int n );
double en_her_monomial_integral ( int n, int alpha[] );
double ep1_glg_monomial_integral ( int expon, double alpha );
double ep1_lag_monomial_integral ( int expon );
void epn_glg_01_1 ( int n, double alpha, int o, double x[], double w[] );
int epn_glg_01_1_size ( int n, double alpha );
void epn_glg_02_xiu ( int n, double alpha, int o, double x[], double w[] );
int epn_glg_02_xiu_size ( int n, double alpha );
double epn_glg_monomial_integral ( int n, int expon[], double alpha );
void epn_lag_01_1 ( int n, int o, double x[], double w[] );
int epn_lag_01_1_size ( int n );
void epn_lag_02_xiu ( int n, int o, double x[], double w[] );
int epn_lag_02_xiu_size ( int n );
double epn_lag_monomial_integral ( int n, int expon[] );
void gw_02_xiu ( int n, int o, double gamma0, double delta0, double c1, 
  double volume_1d, double x[], double w[] );
int gw_02_xiu_size ( int n );
double *monomial_value ( int dim_num, int point_num, double x[], int expon[] );
}
