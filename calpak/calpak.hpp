double cws_to_jed_gps ( int c, int w, double s );
double datenum_to_jed ( double dn );

void day_borrow_alexandrian ( int &y, int &m, int &d );
void day_borrow_common ( int &y, int &m, int &d );
void day_borrow_eg_civil ( int &y, int &m, int &d );
void day_borrow_english ( int &y, int &m, int &d );
void day_borrow_gregorian ( int &y, int &m, int &d );
void day_borrow_hebrew ( int &y, int &m, int &d );
void day_borrow_islamic ( int &y, int &m, int &d );
void day_borrow_julian ( int &y, int &m, int &d );
void day_borrow_republican ( int &y, int &m, int &d );
void day_borrow_roman ( int &y, int &m, int &d );

void day_carry_alexandrian ( int &y, int &m, int &d );
void day_carry_common ( int &y, int &m, int &d );
void day_carry_eg_civil ( int &y, int &m, int &d );
void day_carry_english ( int &y, int &m, int &d );
void day_carry_gregorian ( int &y, int &m, int &d );
void day_carry_hebrew ( int &y, int &m, int &d );
void day_carry_islamic ( int &y, int &m, int &d );
void day_carry_julian ( int &y, int &m, int &d );
void day_carry_republican ( int &y, int &m, int &d );
void day_carry_roman ( int &y, int &m, int &d );

void day_list_common ( int y1, int m1, int d1, int y2, int m2, int d2 );

int days_before_month_common ( int y, int m );
int days_before_month_gregorian ( int y, int m );
int days_before_month_julian ( int y, int m );

void deflate_common ( int &y, int &m, int &d );
void deflate_english ( int &y, int &m, int &d );

char digit_to_ch ( int i );

void easter_ds ( int y, int &m, int &d );
void easter_egr ( int y, int &m, int &d );
void easter_egr2 ( int y, int &m, int &d );
void easter_julian ( int y, int &m, int &d );

double epoch_to_jed_akbar ( );
double epoch_to_jed_alexandrian ( );
double epoch_to_jed_armenian ( );
double epoch_to_jed_bahai ( );
double epoch_to_jed_bessel ( );
double epoch_to_jed_chinese ( );
double epoch_to_jed_common ( );
double epoch_to_jed_coptic ( );
double epoch_to_jed_deccan ( );
double epoch_to_jed_eg_civil ( );
double epoch_to_jed_eg_lunar ( );
double epoch_to_jed_english ( );
double epoch_to_jed_ethiopian ( );
double epoch_to_jed_gps ( );
double epoch_to_jed_greek ( );
double epoch_to_jed_gregorian ( );
double epoch_to_jed_hebrew ( );
double epoch_to_jed_hindu_lunar ( );
double epoch_to_jed_hindu_solar ( );
double epoch_to_jed_islamic_a ( );
double epoch_to_jed_islamic_b ( );
double epoch_to_jed_jed ( );
double epoch_to_jed_jelali ( );
double epoch_to_jed_julian ( );
double epoch_to_jed_khwarizmian ( );
double epoch_to_jed_macedonian ( );
double epoch_to_jed_matlab ( );
double epoch_to_jed_mayan_long ( );
double epoch_to_jed_mjd ( );
double epoch_to_jed_nyt ( );
double epoch_to_jed_persian ( );
double epoch_to_jed_persian_solar ( );
double epoch_to_jed_rd ( );
double epoch_to_jed_republican ( );
double epoch_to_jed_roman ( );
double epoch_to_jed_saka ( );
double epoch_to_jed_soor_san ( );
double epoch_to_jed_syrian ( );
double epoch_to_jed_unix ( );
double epoch_to_jed_y2k ( );
double epoch_to_jed_zoroastrian ( );

void frac_borrow_common ( int &y, int &m, int &d, double &f );
void frac_borrow_english ( int &y, int &m, int &d, double &f );
void frac_borrow_gregorian ( int &y, int &m, int &d, double &f );
void frac_borrow_hebrew ( int &y, int &m, int &d, double &f );
void frac_borrow_islamic ( int &y, int &m, int &d, double &f );
void frac_borrow_julian ( int &y, int &m, int &d, double &f );
void frac_borrow_republican ( int &y, int &m, int &d, double &f );
void frac_borrow_roman ( int &y, int &m, int &d, double &f );

void frac_carry_common ( int &y, int &m, int &d, double &f );
void frac_carry_english ( int &y, int &m, int &d, double &f );
void frac_carry_gregorian ( int &y, int &m, int &d, double &f );
void frac_carry_hebrew ( int &y, int &m, int &d, double &f );
void frac_carry_islamic ( int &y, int &m, int &d, double &f );
void frac_carry_julian ( int &y, int &m, int &d, double &f );
void frac_carry_republican ( int &y, int &m, int &d, double &f );
void frac_carry_roman ( int &y, int &m, int &d, double &f );

void frac_to_hms ( double f, int &h, int &m, int &s );

void hour_borrow_common ( int &y, int &m, int &d, int &h );
void hour_carry_common ( int &y, int &m, int &d, int &h );

int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
void i4_swap ( int &i, int &j );
char i4_to_a ( int i );
int i4_wrap ( int ival, int ilo, int ihi );

void inflate_common ( int &y, int &m, int &d );
void inflate_english ( int &y, int &m, int &d );

void j_borrow_common ( int &y, int &j );
void j_borrow_english ( int &y, int &j );
void j_borrow_gregorian ( int &y, int &j );
void j_borrow_hebrew ( int &y, int &j );
void j_borrow_islamic ( int &y, int &j );
void j_borrow_julian ( int &y, int &j );
void j_borrow_republican ( int &y, int &j );
void j_borrow_roman ( int &y, int &j );

void j_carry_common ( int &y, int &j );
void j_carry_english ( int &y, int &j );
void j_carry_gregorian ( int &y, int &j );
void j_carry_hebrew ( int &y, int &j );
void j_carry_islamic ( int &y, int &j );
void j_carry_julian ( int &y, int &j );
void j_carry_republican ( int &y, int &j );
void j_carry_roman ( int &y, int &j );

void jed_check ( double jed );
double jed_test ( int i );
void jed_to_cws_gps ( double jed, int &c, int &w, double &s );
double jed_to_datenum ( double jed );
void jed_to_mayan_long ( double jed, int &pictun, int &baktun, int &katun, 
  int &tun, int &uinal, int &kin, double &f );
void jed_to_mayan_round ( double jed, int &y, int &a, int &b, int &c, int &d, 
  double &f );
double jed_to_mjd ( double jed );
double jed_to_nearest_noon ( double jed1 );
double jed_to_next_noon ( double jed1 );
double jed_to_rd ( double jed );
double jed_to_ss_unix ( double jed );
void jed_to_weekday ( double jed, int &w, double &f );
int jed_to_year_hebrew ( double jed );
double jed_to_yearcount_bessel ( double jed );
double jed_to_yearcount_julian ( double jed );

void jed_to_yjf_common ( double jed, int &y, int &j, double &f );
void jed_to_yjf_english ( double jed, int &y, int &j, double &f );
void jed_to_yjf_gregorian ( double jed, int &y, int &j, double &f );
void jed_to_yjf_hebrew ( double jed, int &y, int &j, double &f );
void jed_to_yjf_islamic_a ( double jed, int &y, int &j, double &f );
void jed_to_yjf_islamic_b ( double jed, int &y, int &j, double &f );
void jed_to_yjf_julian ( double jed, int &y, int &j, double &f );
void jed_to_yjf_republican ( double jed, int &y, int &j, double &f );
void jed_to_yjf_roman ( double jed, int &y, int &j, double &f );

void jed_to_ymdf_alexandrian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_armenian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_bahai ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_common ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_coptic ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_eg_civil ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_eg_lunar ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_english ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_ethiopian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_gregorian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_gregorian2 ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_hebrew ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_hindu_solar ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_islamic_a ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_islamic_b ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_jelali ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_julian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_julian2 ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_julian3 ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_kwarizmian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_macedonian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_persian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_republican ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_roman ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_saka ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_soor_san ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_syrian ( double jed, int &y, int &m, int &d, double &f );
void jed_to_ymdf_zoroastrian ( double jed, int &y, int &m, int &d, double &f );

double mayan_long_to_jed ( int pictun, int baktun, int katun, int tun, 
  int uinal, int kin, double f );
double mayan_round_to_jed ( int y, int a, int b, int c, int d, double f );

void minute_borrow_common ( int &y, int &m, int &d, int &h, int &n );
void minute_carry_common ( int &y, int &m, int &d, int &h, int &n );

double mjd_to_jed ( double mjd );

void month_borrow_alexandrian ( int &y, int &m );
void month_borrow_bahai ( int &y, int &m );
void month_borrow_common ( int &y, int &m );
void month_borrow_eg_civil ( int &y, int &m );
void month_borrow_english ( int &y, int &m );
void month_borrow_gregorian ( int &y, int &m );
void month_borrow_hebrew ( int &y, int &m );
void month_borrow_islamic ( int &y, int &m );
void month_borrow_julian ( int &y, int &m );
void month_borrow_republican ( int &y, int &m );
void month_borrow_roman ( int &y, int &m );

void month_carry_alexandrian ( int &y, int &m );
void month_carry_bahai ( int &y, int &m );
void month_carry_common ( int &y, int &m );
void month_carry_eg_civil ( int &y, int &m );
void month_carry_english ( int &y, int &m );
void month_carry_gregorian ( int &y, int &m );
void month_carry_hebrew ( int &y, int &m );
void month_carry_islamic ( int &y, int &m );
void month_carry_julian ( int &y, int &m );
void month_carry_republican ( int &y, int &m );
void month_carry_roman ( int &y, int &m );

int month_length_alexandrian ( int y, int m );
int month_length_bahai ( int y, int m );
int month_length_common ( int y, int m );
int month_length_coptic ( int y, int m );
int month_length_eg_civil ( int y, int m );
int month_length_eg_lunar ( int y, int m );
int month_length_english ( int y, int m );
int month_length_ethiopian ( int y, int m );
int month_length_greek ( int y, int m );
int month_length_gregorian ( int y, int m );
int month_length_hebrew ( int y, int m );
double month_length_hindu_solar ( );
int month_length_iranian ( int y, int m );
int month_length_islamic ( int y, int m );
int month_length_julian ( int y, int m );
double month_length_lunar ( int y, int m );
int month_length_persian ( int y, int m );
int month_length_republican ( int y, int m );
int month_length_roman ( int y, int m );
double month_length_synodic ( );

string month_to_month_name_common ( int m );
string month_to_month_name_common3 ( int m );
int month_to_nones_roman ( int m );

void mothers_day ( int y, int &m, int &d );
double new_year_to_jed_hebrew ( int y );

double now_to_jed ( );
void now_to_yjf_common ( int &y, int &j, double &f );
void now_to_ymdf_common ( int &y, int &m, int &d, double &f );
void now_to_ymdhms_common ( int &y, int &m, int &d, int &h, int &n, int &s );

double nyt_to_jed ( int volume, int issue );
void nyt_to_ymd ( int volume, int issue, int &y, int &m, int &d );

double r8_abs ( double x );
double r8_mod ( double x, double y );
int r8_nint ( double x );
double r8_round ( double x );
void r8_swap ( double &x, double &y );
double r8_uniform_ab ( double b, double c, int &seed );

double rd_to_jed ( double rd );

void second_borrow_common ( int &y, int &m, int &d, int &h, int &n, int &s );
void second_carry_common ( int &y, int &m, int &d, int &h, int &n, int &s );

double ss_to_jed_unix ( double s );
void thanksgiving_canada ( int y, int &m, int &d );
void thanksgiving_us ( int y, int &m, int &d );
void timestamp ( );

double transition_to_jed_common ( );
double transition_to_jed_english ( );
double transition_to_jed_jed ( );
double transition_to_jed_mayan_long ( );

int weekday_check_common ( int w );
string weekday_to_name_bahai ( int w );
string weekday_to_name_common ( int w );
string weekday_to_name_common2 ( int w );
string weekday_to_name_common3 ( int w );
string weekday_to_name_french ( int w );
string weekday_to_name_german ( int w );
string weekday_to_name_hebrew ( int w );
string weekday_to_name_islamic ( int w );
string weekday_to_name_republican ( int w );
string weekday_to_name_roman ( int w );

int y_astronomical_to_common ( int y );

void y_check_alexandrian ( int y );
void y_check_bahai ( int y );
void y_check_common ( int y );
void y_check_eg_civil ( int y );
void y_check_english ( int y );
void y_check_greek ( int y );
void y_check_gregorian ( int y );
void y_check_hebrew ( int y );
void y_check_islamic ( int y );
void y_check_julian ( int y );
void y_check_republican ( int y );
void y_check_roman ( int y );

int y_common_to_astronomical ( int y );
int y_julian_to_roman ( int y );
int y_roman_to_julian ( int y );
string y_to_s_common ( int y );

bool year_is_embolismic_eg_lunar ( int y );
bool year_is_embolismic_greek ( int y );
bool year_is_embolismic_hebrew ( int y );

bool year_is_leap_alexandrian ( int y );
bool year_is_leap_bahai ( int y );
bool year_is_leap_common ( int y );
bool year_is_leap_coptic ( int y );
bool year_is_leap_eg_lunar ( int y );
bool year_is_leap_english ( int y );
bool year_is_leap_ethiopian ( int y );
bool year_is_leap_greek ( int y );
bool year_is_leap_gregorian ( int y );
bool year_is_leap_iranian ( int y );
bool year_is_leap_islamic ( int y );
bool year_is_leap_julian ( int y );
bool year_is_leap_persian ( int y );
bool year_is_leap_republican ( int y );
bool year_is_leap_roman ( int y );

int year_length_alexandrian ( int y );
int year_length_bahai ( int y );
int year_length_common ( int y );
int year_length_coptic ( int y );
int year_length_eg_civil ( int y );
int year_length_eg_lunar ( int y );
int year_length_english ( int y );
int year_length_ethiopian ( int y );
int year_length_greek ( int y );
int year_length_gregorian ( int y );
int year_length_hebrew ( int y );
double year_length_hindu_solar ( );
int year_length_islamic ( int y );
int year_length_julian ( int y );
double year_length_lunar ( int y );
int year_length_persian ( int y );
int year_length_republican ( int y );
int year_length_roman ( int y );
double year_length_solar ( int y );

int year_length_months_alexandrian ( int y );
int year_length_months_bahai ( int y );
int year_length_months_common ( int y );
int year_length_months_coptic ( int y );
int year_length_months_eg_civil ( int y );
int year_length_months_eg_lunar ( int y );
int year_length_months_english ( int y );
int year_length_months_ethiopian ( int y );
int year_length_months_greek ( int y );
int year_length_months_gregorian ( int y );
int year_length_months_hebrew ( int y );
int year_length_months_hindu_lunar ( int y );
int year_length_months_hindu_solar ( int y );
int year_length_months_islamic ( int y );
int year_length_months_julian ( int y );
int year_length_months_persian ( int y );
int year_length_months_republican ( int y );
int year_length_months_roman ( int y );

void year_to_dominical_common ( int y, int &n1, int &n2 );
void year_to_dominical_gregorian ( int y, int &n1, int &n2 );
void year_to_dominical_julian ( int y, int &n1, int &n2 );
int year_to_epact_gregorian ( int y );
int year_to_epact_julian ( int y );
int year_to_golden_number ( int y );
int year_to_indiction_common ( int y );
void year_to_scaliger_common ( int y, int &c1, int &c2, int &c3, int &r1, 
  int &r2, int &r3 );
int year_to_type_hebrew ( int y );

void yj_check_common ( int &y, int &j );
void yj_check_english ( int &y, int &j );
void yj_check_gregorian ( int &y, int &j );
void yj_check_hebrew ( int &y, int &j );
void yj_check_islamic ( int &y, int &j );
void yj_check_julian ( int &y, int &j );
void yj_check_republican ( int &y, int &j );
void yj_check_roman ( int &y, int &j );

void yjf_check_common ( int &y, int &j, double &f );
void yjf_check_english ( int &y, int &j, double &f );
char yjf_compare ( int y1, int j1, double f1, int y2, int j2, double f2 );
double yjf_dif_common ( int y1, int j1, double f1, int y2, int j2, double f2 );
void yjf_swap ( int &y1, int &j1, double &f1, int &y2, int &j2, double &f2 );

double yjf_to_jed_common ( int y, int j, double f );
double yjf_to_jed_english ( int y, int j, double f );
double yjf_to_jed_gregorian ( int y, int j, double f );
double yjf_to_jed_hebrew ( int y, int j, double f );
double yjf_to_jed_islamic_a ( int y, int j, double f );
double yjf_to_jed_islamic_b ( int y, int j, double f );
double yjf_to_jed_julian ( int y, int j, double f );
double yjf_to_jed_republican ( int y, int j, double f );
double yjf_to_jed_roman ( int y, int j, double f );

void yjf_to_ymdf_common ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_english ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_gregorian ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_hebrew ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_islamic ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_julian ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_republican ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );
void yjf_to_ymdf_roman ( int y1, int j1, double f1, int &y2, int &m2, int &d2, 
  double &f2 );

void yjf_uniform_common ( int y1, int j1, double f1, int y2, int j2, double f2, 
  int &seed, int &y, int &j, double &f );

void ym_check_alexandrian ( int &y, int &m );
void ym_check_bahai ( int &y, int &m );
void ym_check_common ( int &y, int &m );
void ym_check_eg_civil ( int &y, int &m );
void ym_check_english ( int &y, int &m );
void ym_check_gregorian ( int &y, int &m );
void ym_check_hebrew ( int &y, int &m );
void ym_check_islamic ( int &y, int &m );
void ym_check_julian ( int &y, int &m );
void ym_check_republican ( int &y, int &m );
void ym_check_roman ( int &y, int &m );

double ym_to_decimal ( int y, int m );

void ymd_check_alexandrian ( int &y, int &m, int &d );
void ymd_check_common ( int &y, int &m, int &d );
void ymd_check_eg_civil ( int &y, int &m, int &d );
void ymd_check_english ( int &y, int &m, int &d );
void ymd_check_gregorian ( int &y, int &m, int &d );
void ymd_check_hebrew ( int &y, int &m, int &d );
void ymd_check_islamic ( int &y, int &m, int &d );
void ymd_check_julian ( int &y, int &m, int &d );
void ymd_check_republican ( int &y, int &m, int &d );
void ymd_check_roman ( int &y, int &m, int &d );
char ymd_compare ( int y1, int m1, int d1, int y2, int m2, int d2 );
int ymd_dif_common ( int y1, int m1, int d1, int y2, int m2, int d2 );
void ymd_inc_ymd_common ( int y1, int m1, int d1, int yn, int mn, int dn, 
  int &y2, int &m2, int &d2 );
double ymd_to_decimal ( int y, int m, int d );
double ymd_to_jed_common ( int y, int m, int d );
double ymd_to_jed_gregorian ( int y, int m, int d );
double ymd_to_jed_julian ( int y, int m, int d );
string ymd_to_s_common ( int y, int m, int d );

void ymdf_check_common ( int &y, int &m, int &d, double &f );
void ymdf_check_julian ( int &y, int &m, int &d, double &f );

char ymdf_compare ( int y1, int m1, int d1, double f1, int y2, int m2, int d2, 
  double f2 );
void ymdf_next_common ( int y1, int m1, int d1, double f1, int &y2, int &m2, 
  int &d2, double &f2 );

double ymdf_to_jed_alexandrian ( int y, int m, int d, double f );
double ymdf_to_jed_armenian ( int y, int m, int d, double f );
double ymdf_to_jed_bahai ( int y, int m, int d, double f );
double ymdf_to_jed_common ( int y, int m, int d, double f );
double ymdf_to_jed_coptic ( int y, int m, int d, double f );
double ymdf_to_jed_eg_civil ( int y, int m, int d, double f );
double ymdf_to_jed_eg_lunar ( int y, int m, int d, double f );
double ymdf_to_jed_english ( int y, int m, int d, double f );
double ymdf_to_jed_ethiopian ( int y, int m, int d, double f );
double ymdf_to_jed_gregorian ( int y, int m, int d, double f );
double ymdf_to_jed_hebrew ( int y, int m, int d, double f );
double ymdf_to_jed_hindu_solar ( int y, int m, int d, double f );
double ymdf_to_jed_islamic_a ( int y, int m, int d, double f );
double ymdf_to_jed_islamic_a2 ( int y, int m, int d, double f );
double ymdf_to_jed_islamic_b ( int y, int m, int d, double f );
double ymdf_to_jed_jelali ( int y, int m, int d, double f );
double ymdf_to_jed_julian ( int y, int m, int d, double f );
double ymdf_to_jed_julian2 ( int y, int m, int d, double f );
double ymdf_to_jed_khwarizmian ( int y, int m, int d, double f );
double ymdf_to_jed_macedonian ( int y, int m, int d, double f );
double ymdf_to_jed_persian ( int y, int m, int d, double f );
double ymdf_to_jed_republican ( int y, int m, int d, double f );
double ymdf_to_jed_roman ( int y, int m, int d, double f );
double ymdf_to_jed_saka ( int y, int m, int d, double f );
double ymdf_to_jed_soor_san ( int y, int m, int d, double f );
double ymdf_to_jed_syrian ( int y, int m, int d, double f );
double ymdf_to_jed_zoroastrian ( int y, int m, int d, double f );

int ymdf_to_weekday_common ( int y, int m, int d, double f );
int ymdf_to_weekday_english ( int y, int m, int d, double f );
int ymdf_to_weekday_gregorian ( int y, int m, int d, double f );
int ymdf_to_weekday_julian ( int y, int m, int d, double f );

void ymdf_to_yjf_common ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );
void ymdf_to_yjf_english ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );
void ymdf_to_yjf_gregorian ( int y1, int m1, int d1, double f1, int &y2, 
  int &j2, double &f2 );
void ymdf_to_yjf_hebrew ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );
void ymdf_to_yjf_islamic ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );
void ymdf_to_yjf_julian ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );
void ymdf_to_yjf_republican ( int y1, int m1, int d1, double f1, int &y2, 
  int &j2, double &f2 );
void ymdf_to_yjf_roman ( int y1, int m1, int d1, double f1, int &y2, int &j2, 
  double &f2 );

void ymdf_uniform_common ( int y1, int m1, int d1, double f1, int y2, int m2, 
  int d2, double f2, int &seed, int &y, int &m, int &d, double &f );

double ymdhms_to_decimal ( int y, int m, int d, int h, int n, int s );
