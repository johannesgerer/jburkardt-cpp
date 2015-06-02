double circle_segment_angle_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] );
double circle_segment_angle_from_chord_angles ( double omega1, double omega2 );
double circle_segment_angle_from_height ( double r, double h );
double circle_segment_area_from_angle ( double r, double theta );
double circle_segment_area_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] );
double circle_segment_area_from_height ( double r, double h );
double circle_segment_area_from_sample ( double r, double c[2], double p1[2], 
  double p2[2], int n, int &seed );
double circle_segment_cdf ( double r, double h, double h2 );
double *circle_segment_centroid_from_chord ( double r, double c[2], 
  double p1[2], double p2[2] );
double *circle_segment_centroid_from_height ( double r, double h );
double *circle_segment_centroid_from_sample ( double r, double c[2], 
  double p1[2], double p2[2], int n, int &seed );
int circle_segment_contains_point ( double r, double c[2], double omega1, 
  double omega2, double xy[2] );
double circle_segment_height_from_angle ( double r, double angle );
double circle_segment_height_from_area ( double r, double area );
double circle_segment_height_from_chord ( double r, double c[2], double p1[2], 
  double p2[2] );
double circle_segment_rotation_from_chord ( double r, double c[], double p1[], 
  double p2[] );
void circle_segment_sample_from_chord ( double r, double c[2], double p1[2], 
  double p2[2], int n, int &seed, double x[], double y[] );
void circle_segment_sample_from_height ( double r, double h, int n, int &seed, 
  double x[], double y[] );
double circle_segment_width_from_height ( double r, double h );
void filename_inc ( string *filename );
void gauss ( int n, double alpha[], double beta[], double x[], double w[] );
int i4vec_sum ( int n, int a[] );
void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], 
  double d[], int &it_num, int &rot_num );
void r_jacobi ( int n, double a, double b, double alpha[], double beta[] );
double r8_abs ( double x );
double r8_acos ( double c );
double r8_asin ( double s );
double r8_atan ( double y, double x );
double r8_epsilon ( );
double r8_gamma ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int &seed );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
int s_len_trim ( string s );
void timestamp ( );

