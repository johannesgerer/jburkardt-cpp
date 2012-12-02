int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
void jed_to_weekday ( double jed, int *w, double *f );
double r8_abs ( double x );
double r8_mod ( double x, double y );
int r8_nint ( double x );
void timestamp ( );
string weekday_to_name_common ( int w );
void weekday_values ( int &n_data, int &y, int &m, int &d, int &w );
int y_common_to_astronomical ( int y );
string ymd_to_s_common ( int y, int m, int d );
int ymd_to_weekday_common ( int y, int m, int d );
int ymd_to_weekday_english ( int y, int m, int d );
int ymd_to_weekday_gregorian ( int y, int m, int d );
char ymdf_compare ( int y1, int m1, int d1, double f1, int y2, int m2, int d2, 
  double f2 );
double ymdf_to_jed_common ( int y, int m, int d, double f );
double ymdf_to_jed_english ( int y, int m, int d, double f );
double ymdf_to_jed_gregorian ( int y, int m, int d, double f );
double ymdf_to_jed_julian ( int y, int m, int d, double f );

