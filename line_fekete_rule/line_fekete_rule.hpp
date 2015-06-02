double *cheby_van1 ( int m, double a, double b, int n, double x[] );
double *legendre_van ( int m, double a, double b, int n, double x[] );
void line_fekete_chebyshev ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] );
void line_fekete_legendre ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] );
void line_fekete_monomial ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] );
double *line_monomial_moments ( double a, double b, int m );

