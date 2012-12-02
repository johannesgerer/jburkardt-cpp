double *correlation_besselj ( int n, double rho[], double rho0 );
double *correlation_besselk ( int n, double rho[], double rho0 );
double *correlation_brownian ( int m, int n, double s[], double t[], double rho0 );
void correlation_brownian_display ( );
double *correlation_circular ( int n, double rho[], double rho0 );
double *correlation_constant ( int n, double rho[], double rho0 );
double *correlation_cubic ( int n, double rho[], double rho0 );
double *correlation_damped_cosine ( int n, double rho[], double rho0 );
double *correlation_damped_sine ( int n, double rho[], double rho0 );
double *correlation_exponential ( int n, double rho[], double rho0 );
double *correlation_gaussian ( int n, double rho[], double rho0 );
double *correlation_hole ( int n, double rho[], double rho0 );
double *correlation_linear ( int n, double rho[], double rho0 );
double *correlation_matern ( int n, double rho[], double rho0 );
double *correlation_pentaspherical ( int n, double rho[], double rho0 );
void correlation_plot ( int n, double rho[], double c[], string header, 
  string title );
void correlation_plots ( int n, int n2, double rho[], double rho0[], double c[], 
  string header, string title );
double *correlation_power ( int n, double rho[], double rho0 );
double *correlation_rational_quadratic ( int n, double rho[], double rho0 );
double *correlation_spherical ( int n, double rho[], double rho0 );
double *correlation_to_covariance ( int n, double c[], double sigma[] );
double *correlation_white_noise ( int n, double rho[], double rho0 );
void covariance_to_correlation ( int n, double k[], double c[], double sigma[] );
int i4_abs ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
double *minij ( int m, int n );
void paths_plot ( int n, int n2, double rho[], double x[], string header, 
  string title );
double pythag ( double a, double b );
double r8_abs ( double x );
double r8_aint ( double x );
void r8_b0mp ( double x, double &ampl, double &theta );
double r8_besi1 ( double x );
double r8_besi1e ( double x );
double r8_besj0 ( double x );
double r8_besk ( double nu, double x );
double r8_besk1 ( double x );
double r8_besk1e ( double x );
double *r8_beskes ( double xnu, double x, int nin );
double *r8_besks ( double xnu, double x, int nin );
double r8_csevl ( double x, double a[], int n );
double r8_epsilon ( void );
void r8_gaml ( double &xmin, double &xmax );
double r8_gamma ( double x );
double r8_huge ( void );
int r8_inits ( double dos[], int nos, double eta );
void r8_knus ( double xnu, double x, double &bknu, double &bknu1, int &iswtch );
double r8_lgmc ( double x );
double r8_mach ( int i );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int &seed );
double *r8mat_cholesky_factor ( int n, double a[], int &flag );
double *r8mat_copy_new ( int m, int n, double a1[] );
double r8mat_is_symmetric ( int m, int n, double a[] );
double r8mat_max ( int m, int n, double a[] );
double r8mat_min ( int m, int n, double a[] );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
double *r8mat_normal_01_new ( int m, int n, int &seed );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_min ( int n, double r8vec[] );
double *r8vec_normal_01_new ( int n, int &seed );
void r8vec_print ( int n, double a[], string title );
double *r8vec_uniform_01_new ( int n, int &seed );
double *sample_paths_cholesky ( int n, int n2, double rhomax, double rho0, 
  double *correlation ( int n, double rho_vec[], double rho0 ), int &seed );
double *sample_paths_eigen ( int n, int n2, double rhomax, double rho0, 
  double *correlation ( int n, double rho_vec[], double rho0 ), int &seed );
double *sample_paths2_cholesky ( int n, int n2, double rhomin, double rhomax, 
  double rho0, double *correlation2 ( int m, int n, double s[], double t[], double rho0 ), 
  int &seed );
double *sample_paths2_eigen ( int n, int n2, double rhomin, double rhomax, 
  double rho0, double *correlation2 ( int m, int n, double s[], double t[], double rho0 ), 
  int &seed );
void timestamp ( void );
int tql2 ( int n, double d[], double e[], double z[] );
void tred2 ( int n, double a[], double d[], double e[], double z[] );

