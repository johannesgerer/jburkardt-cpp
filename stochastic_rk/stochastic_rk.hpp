double rk1_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed );
double rk2_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed );
double rk3_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed );
double rk4_ti_step ( double x, double t, double h, double q, 
  double fi ( double x ), double gi ( double x ), int *seed );

double rk1_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed );
double rk2_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed );
double rk4_tv_step ( double x, double t, double h, double q, 
  double fv ( double t, double x ), double gv ( double t, double x ), 
  int *seed );

double r8_normal_01 ( int *seed );
double r8_uniform_01 ( int *seed );
void timestamp ( );
