bool between ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
bool collinear ( double xa, double ya, double xb, double yb, double xc, 
  double yc );
bool diagonal ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool diagonalie ( int im1, int ip1, int n, int next[], double x[], double y[] );
bool in_cone ( int im1, int ip1, int n, int prev[], int next[], double x[], 
  double y[] );
bool intersect ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool intersect_prop ( double xa, double ya, double xb, double yb, double xc, 
  double yc, double xd, double yd );
bool l4_xor ( bool l1, bool l2 );
double *monomial_value ( int m, int n, int e[], double x[] );
double polygon_area ( int nv, double v[] );
double polygon_monomial_integral ( int nv, double v[], int e[] );
double *polygon_sample ( int nv, double v[], int n, int &seed );
int *polygon_triangulate ( int n, double x[], double y[] );
double r8_choose ( int n, int k );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int &seed );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
void timestamp ( );
double triangle_area ( double xa, double ya, double xb, double yb, double xc, 
  double yc );

