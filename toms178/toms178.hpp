double best_nearby ( double delta[], double point[], double prevbest, 
  int nvars, double f ( double x[], int nvars ), int *funevals );
int hooke ( int nvars, double startpt[], double endpt[], double rho, double eps, 
  int itermax, double f ( double x[], int nvars ) );
double r8_abs ( double x );
void timestamp ( void );
