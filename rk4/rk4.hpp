double rk4 ( double t0, double u0, double dt, double f ( double t, double u ) );
double *rk4vec ( double t0, int n, double u0[], double dt, 
  double *f ( double t, int n, double u[] ) );
void timestamp ( );
