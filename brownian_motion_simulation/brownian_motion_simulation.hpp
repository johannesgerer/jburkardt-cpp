void brownian_displacement_display ( int k, int n, double d, double t, 
  double dsq[], string header );
double *brownian_displacement_simulation ( int k, int n, int m, double d, 
  double t, int &seed );
void brownian_motion_display ( int m, int n, double x[], string header );
double *brownian_motion_simulation ( int m, int n, double d, double t, 
  int &seed );
int *i4vec_uniform_new ( int n, int a, int b, int &seed );
double r8_normal_01 ( int &seed );
double r8_uniform_01 ( int &seed );
double *r8vec_normal_01_new ( int n, int &seed );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
