void ou_euler ( double theta, double mu, double sigma, double x0, double tmax, 
  int n, int &seed );
void ou_euler_maruyama ( double theta, double mu, double sigma, double x0, 
  double tmax, int n, int r, int &seed );
double r8_uniform_01 ( int &seed );
double *r8vec_linspace_new ( int n, double a, double b );
double *r8vec_normal_01_new ( int n, int &seed );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
