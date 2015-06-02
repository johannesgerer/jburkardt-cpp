double angle_cdf ( double x, int n );
double angle_mean ( int a );
double angle_pdf ( double x, int a );

double anglit_cdf ( double x );
double anglit_cdf_inv ( double cdf );
double anglit_mean ( );
double anglit_pdf ( double x );
double anglit_sample ( int &seed );
double anglit_variance ( );

double arcsin_cdf ( double x, double a );
double arcsin_cdf_inv ( double cdf, double a );
bool arcsin_check ( double a );
double arcsin_mean ( double a );
double arcsin_pdf ( double x, double a );
double arcsin_sample ( double a, int &seed );
double arcsin_variance ( double a );

double benford_pdf ( double x );

double bernoulli_cdf ( int x, double a );
int bernoulli_cdf_inv ( double cdf, double a );
bool bernoulli_check ( double a );
double bernoulli_mean ( double a );
double bernoulli_pdf ( int x, double a );
int bernoulli_sample ( double a, int &seed );
double bernoulli_variance ( double a );

double bessel_i0 ( double arg );
void bessel_i0_values ( int &n_data, double &x, double &fx );
double bessel_i1 ( double arg );
void bessel_i1_values ( int &n_data, double &x, double &fx );
void bessel_ix_values ( int &n_data, double &nu, double &x, double &fx );

double beta ( double x, double y );

double beta_binomial_cdf ( int x, double a, double b, int c );
int beta_binomial_cdf_inv ( double cdf, double a, double b, int c );
bool beta_binomial_check ( double a, double b, int c );
double beta_binomial_mean ( double a, double b, int c );
double beta_binomial_pdf ( int x, double a, double b, int c );
int beta_binomial_sample ( double a, double b, int c, int &seed );
double beta_binomial_variance ( double a, double b, int c );

double beta_cdf ( double x, double a, double b );
double beta_cdf_inv ( double cdf, double a, double b );
double beta_cdf_inv_old ( double cdf, double a, double b );
void beta_cdf_values ( int &n_data, double &a, double &b, double &x, 
  double &fx );
bool beta_check ( double a, double b );
double beta_inc ( double a, double b, double x );
void beta_inc_values ( int &n_data, double &a, double &b, double &x, 
  double &fx );
double beta_mean ( double a, double b );
double beta_pdf ( double x, double a, double b );
double beta_sample ( double a, double b, int &seed );
double beta_variance ( double a, double b );

double binomial_cdf ( double x, int a, double b );
void binomial_cdf_values ( int &n_data, int &a, double &b, int &x, double &fx );
int binomial_cdf_inv ( double cdf, int a, double b );
bool binomial_check ( int a, double b );
int binomial_coef ( int n, int k );
double binomial_coef_log ( int n, int k );
double binomial_mean ( int a, double b );
double binomial_pdf ( int x, int a, double b );
int binomial_sample ( int a, double b, int &seed );
double binomial_variance ( int a, double b );

double birthday_cdf ( int n );
int birthday_cdf_inv ( double cdf );
double birthday_pdf ( int n );

double bradford_cdf ( double x, double a, double b, double c );
double bradford_cdf_inv ( double cdf, double a, double b, double c );
bool bradford_check ( double a, double b, double c );
double bradford_mean ( double a, double b, double c );
double bradford_pdf ( double x, double a, double b, double c );
double bradford_sample ( double a, double b, double c, int &seed );
double bradford_variance ( double a, double b, double c );

double buffon_laplace_pdf ( double a, double b, double l );
int buffon_laplace_simulate ( double a, double b, double l, int trial_num );

double buffon_pdf ( double a, double l );
int buffon_simulate ( double a, double l, int trial_num );

double burr_cdf ( double x, double a, double b, double c, double d );
double burr_cdf_inv ( double cdf, double a, double b, double c, double d );
bool burr_check ( double a, double b, double c, double d );
double burr_mean ( double a, double b, double c, double d );
double burr_pdf ( double x, double a, double b, double c, double d );
double burr_sample ( double a, double b, double c, double d, int &seed );
double burr_variance ( double a, double b, double c, double d );

double cardioid_cdf ( double x, double a, double b );
double cardioid_cdf_inv ( double cdf, double a, double b );
bool cardioid_check ( double a, double b );
double cardioid_mean ( double a, double b );
double cardioid_pdf ( double x, double a, double b );
double cardioid_sample ( double a, double b, int &seed );
double cardioid_variance ( double a, double b );

double cauchy_cdf ( double x, double a, double b );
double cauchy_cdf_inv ( double cdf, double a, double b );
void cauchy_cdf_values ( int &n_data, double &mu, double &sigma, double &x,
  double &fx );
bool cauchy_check ( double a, double b );
double cauchy_mean ( double a, double b );
double cauchy_pdf ( double x, double a, double b );
double cauchy_sample ( double a, double b, int &seed );
double cauchy_variance ( double a, double b );

double chi_cdf ( double x, double a, double b, double c );
double chi_cdf_inv ( double cdf, double a, double b, double c );
bool chi_check ( double a, double b, double c );
double chi_mean ( double a, double b, double c );
double chi_pdf ( double x, double a, double b, double c );
double chi_sample ( double a, double b, double c, int &seed );
double chi_variance ( double a, double b, double c );

double chi_square_cdf ( double x, double a );
double chi_square_cdf_inv ( double cdf, double a );
void chi_square_cdf_values ( int &n_data, int &a, double &x, double &fx );
bool chi_square_check ( double a );
double chi_square_mean ( double a );
double chi_square_pdf ( double x, double a );
double chi_square_sample ( double a, int &seed );
double chi_square_variance ( double a );

void chi_square_noncentral_cdf_values ( int &n_data, int &df, double &lambda, 
  double &x, double &cdf );
double chi_square_noncentral_check ( double a, double b );
double chi_square_noncentral_mean ( double a, double b );
double chi_square_noncentral_sample ( double a, double b, int &seed );
double chi_square_noncentral_variance ( double a, double b );

double *circle_sample ( double a, double b, double c, int &seed );

double *circular_normal_01_mean ( );
double circular_normal_01_pdf ( double x[2] );
double *circular_normal_01_sample ( int &seed );
double *circular_normal_01_variance ( );

double *circular_normal_mean ( double a[2], double b );
double circular_normal_pdf ( double x[2], double a[2], double b );
double *circular_normal_sample ( double a[2], double b, int &seed );
double *circular_normal_variance ( double a[2], double b );

int combinatorial ( int n, int k );

double cosine_cdf ( double x, double a, double b );
double cosine_cdf_inv ( double cdf, double a, double b );
bool cosine_check ( double a, double b );
double cosine_mean ( double a, double b );
double cosine_pdf ( double x, double a, double b );
double cosine_sample ( double a, double b, int &seed );
double cosine_variance ( double a, double b );

double coupon_complete_pdf ( int type_num, int box_num );
double coupon_mean ( int j, int n );
void coupon_simulate ( int n_type, int &seed, int coupon[], int *n_coupon );
double coupon_variance ( int j, int n );

double deranged_cdf ( int x, int a );
int deranged_cdf_inv ( double cdf, int a );
bool deranged_check ( int a );
int deranged_enum ( int n );
double deranged_mean ( int a );
double deranged_pdf ( int x, int a );
int deranged_sample ( int a, int &seed );
double deranged_variance ( int a );

double digamma ( double x );

double dipole_cdf ( double x, double a, double b );
double dipole_cdf_inv ( double cdf, double a, double b );
bool dipole_check ( double a, double b );
double dipole_pdf ( double x, double a, double b );
double dipole_sample ( double a, double b, int &seed );

bool dirichlet_check ( int n, double a[] );
double *dirichlet_mean ( int n, double a[] );
double *dirichlet_moment2 ( int n, double a[] );
double dirichlet_pdf ( double x[], int n, double a[] );
double *dirichlet_sample ( int n, double a[], int &seed );
double *dirichlet_variance ( int n, double a[] );

bool dirichlet_mix_check ( int comp_num, int elem_num, double a[], 
  double comp_weight[] );
double *dirichlet_mix_mean ( int comp_num, int elem_num, double a[], 
  double comp_weight[] );
double dirichlet_mix_pdf ( double x[], int comp_num, int elem_num, double a[], 
  double comp_weight[] );
double *dirichlet_mix_sample ( int comp_num, int elem_num, double a[], 
  double comp_weight[], int &seed, int *comp );

double dirichlet_multinomial_pdf ( int x[], int a, int b, double c[] );

double discrete_cdf ( int x, int a, double b[] );
int discrete_cdf_inv ( double cdf, int a, double b[] );
bool discrete_check ( int a, double b[] );
double discrete_mean ( int a, double b[] );
double discrete_pdf ( int x, int a, double b[] );
int discrete_sample ( int a, double b[], int &seed );
double discrete_variance ( int a, double b[] );

double e_constant ( );

double empirical_discrete_cdf ( double x, int a, double b[], double c[] );
double empirical_discrete_cdf_inv ( double cdf, int a, double b[], double c[] );
bool empirical_discrete_check ( int a, double b[], double c[] );
double empirical_discrete_mean ( int a, double b[], double c[] );
double empirical_discrete_pdf ( double x, int a, double b[], double c[] );
double empirical_discrete_sample ( int a, double b[], double c[], int &seed );
double empirical_discrete_variance ( int a, double b[], double c[] );

double english_sentence_length_cdf ( int x );
int english_sentence_length_cdf_inv ( double cdf );
double english_sentence_length_mean ( );
double english_sentence_length_pdf ( int x );
int english_sentence_length_sample ( int &seed );
double english_sentence_length_variance ( );

double english_word_length_cdf ( int x );
int english_word_length_cdf_inv ( double cdf );
double english_word_length_mean ( );
double english_word_length_pdf ( int x );
int english_word_length_sample ( int &seed );
double english_word_length_variance ( );

double erlang_cdf ( double x, double a, double b, int c );
double erlang_cdf_inv ( double cdf, double a, double b, int c );
bool erlang_check ( double a, double b, int c );
double erlang_mean ( double a, double b, int c );
double erlang_pdf ( double x, double a, double b, int c );
double erlang_sample ( double a, double b, int c, int &seed );
double erlang_variance ( double a, double b, int c );

double error_f ( double x );
double error_f_inverse ( double y );

double euler_constant ( );

double exponential_01_cdf ( double x );
double exponential_01_cdf_inv ( double cdf );
double exponential_01_mean ( );
double exponential_01_pdf ( double x );
double exponential_01_sample ( int &seed );
double exponential_01_variance ( );

double exponential_cdf ( double x, double a, double b );
double exponential_cdf_inv ( double cdf, double a, double b );
void exponential_cdf_values ( int &n_data, double &lambda, double &x, 
  double &fx );
bool exponential_check ( double a, double b );
double exponential_mean ( double a, double b );
double exponential_pdf ( double x, double a, double b );
double exponential_sample ( double a, double b, int &seed );
double exponential_variance ( double a, double b );

double extreme_values_cdf ( double x, double a, double b );
double extreme_values_cdf_inv ( double cdf, double a, double b );
void extreme_values_cdf_values ( int &n_data, double &alpha, double &beta, 
  double &x, double &fx );
bool extreme_values_check ( double a, double b );
double extreme_values_mean ( double a, double b );
double extreme_values_pdf ( double x, double a, double b );
double extreme_values_sample ( double a, double b, int &seed );
double extreme_values_variance ( double a, double b );

double f_cdf ( double x, int m, int n );
void f_cdf_values ( int &n_data, int &a, int &b, double &x, double &fx );
bool f_check ( int m, int n );
double f_mean ( int m, int n );
double f_pdf ( double x, int m, int n );
double f_sample ( int m, int n, int &seed );
double f_variance ( int m, int n );

void f_noncentral_cdf_values ( int &n_data, int &n1, int &n2, double &lambda, 
  double &x, double &fx );
bool f_noncentral_check ( double a, int m, int n );
double f_noncentral_mean ( double a, int m, int n );
double f_noncentral_variance ( double a, int m, int n );

double factorial_log ( int n );
double factorial_stirling ( int n );

double fisher_pdf ( double x[3], double kappa, double mu[3] );
double *fisher_sample ( double kappa, double xyz[], int n, int &seed );

double fisk_cdf ( double x, double a, double b, double c );
double fisk_cdf_inv ( double cdf, double a, double b, double c );
bool fisk_check ( double a, double b, double c );
double fisk_mean ( double a, double b, double c );
double fisk_pdf ( double x, double a, double b, double c );
double fisk_sample ( double a, double b, double c, int &seed );
double fisk_variance ( double a, double b, double c );

double folded_normal_cdf ( double x, double a, double b );
double folded_normal_cdf_inv ( double cdf, double a, double b );
bool folded_normal_check ( double a, double b );
double folded_normal_mean ( double a, double b );
double folded_normal_pdf ( double x, double a, double b );
double folded_normal_sample ( double a, double b, int &seed );
double folded_normal_variance ( double a, double b );

double frechet_cdf ( double x, double alpha );
double frechet_cdf_inv ( double cdf, double alpha );
double frechet_mean ( double alpha );
double frechet_pdf ( double x, double alpha );
double frechet_sample ( double alpha, int &seed );
double frechet_variance ( double alpha );

double gamma_cdf ( double x, double a, double b, double c );
void gamma_cdf_values ( int &n_data, double &mu, double &sigma, double &x,
  double &fx );
bool gamma_check ( double a, double b, double c );
double gamma_mean ( double a, double b, double c );
double gamma_pdf ( double x, double a, double b, double c );
double gamma_sample ( double a, double b, double c, int &seed );
double gamma_variance ( double a, double b, double c );

double gamma_inc ( double p, double x );
void gamma_inc_values ( int &n_data, double &a, double &x, double &fx );
double gamma_log ( double x );
double gamma_log_int ( int n );

double genlogistic_cdf ( double x, double a, double b, double c );
double genlogistic_cdf_inv ( double cdf, double a, double b, double c );
bool genlogistic_check ( double a, double b, double c );
double genlogistic_mean ( double a, double b, double c );
double genlogistic_pdf ( double x, double a, double b, double c );
double genlogistic_sample ( double a, double b, double c, int &seed );
double genlogistic_variance ( double a, double b, double c );

double geometric_cdf ( int x, double a );
int geometric_cdf_inv ( double cdf, double a );
void geometric_cdf_values ( int &n_data, int &x, double &p, double &cdf );
bool geometric_check ( double a );
double geometric_mean ( double a );
double geometric_pdf ( int x, double a );
int geometric_sample ( double a, int &seed );
double geometric_variance ( double a );

int get_seed ( );

double gompertz_cdf ( double x, double a, double b );
double gompertz_cdf_inv ( double cdf, double a, double b );
bool gompertz_check ( double a, double b );
double gompertz_pdf ( double x, double a, double b );
double gompertz_sample ( double a, double b, int &seed );

double gumbel_cdf ( double x );
double gumbel_cdf_inv ( double cdf );
double gumbel_mean ( );
double gumbel_pdf ( double x );
double gumbel_sample ( int &seed );
double gumbel_variance ( );

double half_normal_cdf ( double x, double a, double b );
double half_normal_cdf_inv ( double cdf, double a, double b );
bool half_normal_check ( double a, double b );
double half_normal_mean ( double a, double b );
double half_normal_pdf ( double x, double a, double b );
double half_normal_sample ( double a, double b, int &seed );
double half_normal_variance ( double a, double b );

double hypergeometric_cdf ( int x, int n, int m, int l );
void hypergeometric_cdf_values ( int &n_data, int &sam, int &suc, int &pop, 
  int &n, double &fx );
bool hypergeometric_check ( int n, int m, int l );
double hypergeometric_mean ( int n, int m, int l );
double hypergeometric_pdf ( int x, int n, int m, int l );
int hypergeometric_sample ( int n, int m, int l, int &seed );
double hypergeometric_variance ( int n, int m, int l );

double i4_factorial ( int n );
int i4_huge ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform_ab ( int b, int c, int &seed );

int *i4row_max ( int m, int n, int a[] );
double *i4row_mean ( int m, int n, int a[] );
int *i4row_min ( int m, int n, int a[] );
double *i4row_variance ( int m, int n, int a[] );

int i4vec_max ( int n, int x[] );
double i4vec_mean ( int n, int x[] );
int i4vec_min ( int n, int x[] );
void i4vec_print ( int n, int a[], char *title );
void i4vec_print ( int n, int a[], char *title );
int i4vec_sum ( int n, int a[] );
double i4vec_variance ( int n, int x[] );

double inverse_gaussian_cdf ( double x, double a, double b );
bool inverse_gaussian_check ( double a, double b );
double inverse_gaussian_mean ( double a, double b );
double inverse_gaussian_pdf ( double x, double a, double b );
double inverse_gaussian_sample ( double a, double b, int &seed );
double inverse_gaussian_variance ( double a, double b );

double laplace_cdf ( double x, double a, double b );
double laplace_cdf_inv ( double cdf, double a, double b );
void laplace_cdf_values ( int &n_data, double &mu, double &beta, double &x,
  double &fx );
bool laplace_check ( double a, double b );
double laplace_mean ( double a, double b );
double laplace_pdf ( double x, double a, double b );
double laplace_sample ( double a, double b, int &seed );
double laplace_variance ( double a, double b );

double lerch ( double a, double b, double c );

double levy_cdf ( double x, double a, double b );
double levy_cdf_inv ( double cdf, double a, double b );
double levy_pdf ( double x, double a, double b );
double levy_sample ( double a, double b, int &seed );

double log_normal_cdf ( double x, double a, double b );
double log_normal_cdf_inv ( double cdf, double a, double b );
void log_normal_cdf_values ( int &n_data, double &mu, double &sigma, 
  double &x, double &fx );
bool log_normal_check ( double a, double b );
double log_normal_mean ( double a, double b );
double log_normal_pdf ( double x, double a, double b );
double log_normal_sample ( double a, double b, int &seed );
double log_normal_variance ( double a, double b );

double log_series_cdf ( double x, double a );
int log_series_cdf_inv ( double cdf, double a );
void log_series_cdf_values ( int &n_data, double &t, int &n, double &fx );
bool log_series_check ( double a );
double log_series_mean ( double a );
double log_series_pdf ( int x, double a );
int log_series_sample ( double a, int &seed );
double log_series_variance ( double a );

double log_uniform_cdf ( double x, double a, double b );
double log_uniform_cdf_inv ( double cdf, double a, double b );
bool log_uniform_check ( double a, double b );
double log_uniform_mean ( double a, double b );
double log_uniform_pdf ( double x, double a, double b );
double log_uniform_sample ( double a, double b, int &seed );
double log_uniform_variance ( double a, double b );

double logistic_cdf ( double x, double a, double b );
double logistic_cdf_inv ( double cdf, double a, double b );
void logistic_cdf_values ( int &n_data, double &mu, double &beta, double &x,
  double &fx );
bool logistic_check ( double a, double b );
double logistic_mean ( double a, double b );
double logistic_pdf ( double x, double a, double b );
double logistic_sample ( double a, double b, int &seed );
double logistic_variance ( double a, double b );

double lorentz_cdf ( double x );
double lorentz_cdf_inv ( double cdf );
double lorentz_mean ( );
double lorentz_pdf ( double x );
double lorentz_sample ( int &seed );
double lorentz_variance ( );

double maxwell_cdf ( double x, double a );
double maxwell_cdf_inv ( double cdf, double a );
bool maxwell_check ( double a );
double maxwell_mean ( double a );
double maxwell_pdf ( double x, double a );
double maxwell_sample ( double a, int &seed );
double maxwell_variance ( double a );

bool multicoef_check ( int nfactor, int factor[] );
int multinomial_coef1 ( int nfactor, int factor[] );
int multinomial_coef2 ( int nfactor, int factor[] );

bool multinomial_check ( int a, int b, double c[] );
double *multinomial_covariance ( int a, int b, double c[] );
double *multinomial_mean ( int a, int b, double c[] );
double multinomial_pdf ( int x[], int a, int b, double c[] );
int *multinomial_sample ( int a, int b, double c[], int &seed );
double *multinomial_variance ( int a, int b, double c[] );

double *multivariate_normal_sample ( int n, double mean[], double covar_factor[],
  int &seed );

double nakagami_cdf ( double x, double a, double b, double c );
bool nakagami_check ( double a, double b, double c );
double nakagami_mean ( double a, double b, double c );
double nakagami_pdf ( double x, double a, double b, double c );
double nakagami_variance ( double a, double b, double c );

double negative_binomial_cdf ( int x, int a, double b );
int negative_binomial_cdf_inv ( double cdf, int a, double b );
void negative_binomial_cdf_values ( int &n_data, int &f, int &s, double &p, 
  double &cdf );
bool negative_binomial_check ( int a, double b );
double negative_binomial_mean ( int a, double b );
double negative_binomial_pdf ( int x, int a, double b );
int negative_binomial_sample ( int a, double b, int &seed );
double negative_binomial_variance ( int a, double b );

double normal_01_cdf ( double x );
double normal_01_cdf_inv ( double cdf );
void normal_01_cdf_values ( int &n_data, double &x, double &fx );
double normal_01_mean ( );
double normal_01_pdf ( double x );
double normal_01_sample ( int &seed );
double normal_01_variance ( );

double *normal_01_vector ( int n, int &seed );

double normal_cdf ( double x, double a, double b );
double normal_cdf_inv ( double cdf, double a, double b );
void normal_cdf_values ( int &n_data, double &mu, double &sigma, double &x, 
  double &fx );
bool normal_check ( double a, double b );
double normal_mean ( double a, double b );
double normal_pdf ( double x, double a, double b );
double normal_sample ( double a, double b, int &seed );
double normal_variance ( double a, double b );

double *normal_vector ( int n, double mean, double dev, int &seed );

double normal_truncated_ab_cdf ( double x, double mu, double s, double a, 
  double b );
double normal_truncated_ab_cdf_inv ( double cdf, double mu, double s, double a, 
  double b );
double normal_truncated_ab_mean ( double mu, double s, double a, double b );
double normal_truncated_ab_pdf ( double x, double mu, double s, double a, 
  double b );
double normal_truncated_ab_sample ( double mu, double s, double a, double b, 
  int &seed );
double normal_truncated_ab_variance ( double mu, double s, double a, double b );

double normal_truncated_a_cdf ( double x, double mu, double s, double a );
double normal_truncated_a_cdf_inv ( double cdf, double mu, double s, double a );
double normal_truncated_a_mean ( double mu, double s, double a );
double normal_truncated_a_pdf ( double x, double mu, double s, double a );
double normal_truncated_a_sample ( double mu, double s, double a, int &seed );
double normal_truncated_a_variance ( double mu, double s, double a );

double normal_truncated_b_cdf ( double x, double mu, double s, double b );
double normal_truncated_b_cdf_inv ( double cdf, double mu, double s, double b );
double normal_truncated_b_mean ( double mu, double s, double b );
double normal_truncated_b_pdf ( double x, double mu, double s, double b );
double normal_truncated_b_sample ( double mu, double s, double b, int &seed );
double normal_truncated_b_variance ( double mu, double s, double b );

void owen_values ( int &n_data, double &h, double &a, double &t );

double pareto_cdf ( double x, double a, double b );
double pareto_cdf_inv ( double cdf, double a, double b );
bool pareto_check ( double a, double b );
double pareto_mean ( double a, double b );
double pareto_pdf ( double x, double a, double b );
double pareto_sample ( double a, double b, int &seed );
double pareto_variance ( double a, double b );

bool pearson_05_check ( double a, double b, double c );
double pearson_05_mean ( double a, double b, double c );
double pearson_05_pdf ( double x, double a, double b, double c );
double pearson_05_sample ( double a, double b, double c, int &seed );

bool planck_check ( double a, double b );
double planck_mean ( double a, double b );
double planck_pdf ( double x, double a, double b );
double planck_sample ( double a, double b, int &seed );
double planck_variance ( double a, double b );

double point_distance_1d_pdf ( double x, int a, double b );
double point_distance_2d_pdf ( double x, int a, double b );
double point_distance_3d_pdf ( double x, int a, double b );

double poisson_cdf ( int k, double a );
int poisson_cdf_inv ( double cdf, double a );
void poisson_cdf_values ( int &n_data, double &a, int &x, double &fx );
bool poisson_check ( double a );
double poisson_kernel ( double r, int n, double c[], double x[], double y[] );
double poisson_mean ( double a );
double poisson_pdf ( int k, double a );
int poisson_sample ( double a, int &seed );
double poisson_variance ( double a );

double power_cdf ( double x, double a, double b );
double power_cdf_inv ( double cdf, double a, double b );
bool power_check ( double a, double b );
double power_mean ( double a, double b );
double power_pdf ( double x, double a, double b );
double power_sample ( double a, double b, int &seed );
double power_variance ( double a, double b );

void psi_values ( int &n_data, double &x, double &fx );

double quasigeometric_cdf ( int x, double a, double b );
int quasigeometric_cdf_inv ( double cdf, double a, double b );
bool quasigeometric_check ( double a, double b );
double quasigeometric_mean ( double a, double b );
double quasigeometric_pdf ( int x, double a, double b );
int quasigeometric_sample ( double a, double b, int &seed );
double quasigeometric_variance ( double a, double b );

float r4_abs ( float x );
int r4_nint ( float x );
float r4_uniform_01 ( int &seed );

double r8_abs ( double x );
int r8_ceiling ( double r );
double r8_csc ( double theta );
double r8_epsilon ( );
double r8_gamma ( double x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_modp ( double x, double y );
int r8_nint ( double x );
double r8_pi ( );
double r8_random ( double rlo, double rhi, int &seed );
double r8_sign ( double x );
double r8_uniform_01 ( int &seed );

void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );

double r8poly_value ( int n, double a[], double x );

double *r8row_max ( int m, int n, double a[] );
double *r8row_mean ( int m, int n, double a[] );
double *r8row_min ( int m, int n, double a[] );
double *r8row_variance ( int m, int n, double a[] );

double r8vec_circular_variance ( int n, double x[] );
double r8vec_diff_norm ( int n, double a[], double b[] );
double r8vec_dot ( int n, double a1[], double a2[] );
double r8vec_length ( int dim_num, double x[] );
double r8vec_max ( int n, double x[] );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double x[] );
void r8vec_print ( int n, double a[], char *title );
double *r8vec_random ( int n, double alo, double ahi, int &seed );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01 ( int n, int &seed );
void r8vec_unit_sum ( int n, double a[] );
double r8vec_variance ( int n, double x[] );

double rayleigh_cdf ( double x, double a );
double rayleigh_cdf_inv ( double cdf, double a );
void rayleigh_cdf_values ( int &n_data, double &sigma, double &x, double &fx );
bool rayleigh_check ( double a );
double rayleigh_mean ( double a );
double rayleigh_pdf ( double x, double a );
double rayleigh_sample ( double a, int &seed );
double rayleigh_variance ( double a );

double reciprocal_cdf ( double x, double a, double b );
double reciprocal_cdf_inv ( double cdf, double a, double b );
bool reciprocal_check ( double a, double b );
double reciprocal_mean ( double a, double b );
double reciprocal_pdf ( double x, double a, double b );
double reciprocal_sample ( double a, double b, int &seed );
double reciprocal_variance ( double a, double b );

int ribesl ( double kappa, double alpha, int nb, int ize, double b[] );

double runs_mean ( int m, int n );
double runs_pdf ( int m, int n, int r );
int runs_sample ( int m, int n, int &seed );
int *runs_simulate ( int m, int n, int &seed );
double runs_variance ( int m, int n );

int s_len_trim ( char *s );

double sech ( double x );

double sech_cdf ( double x, double a, double b );
double sech_cdf_inv ( double cdf, double a, double b );
bool sech_check ( double a, double b );
double sech_mean ( double a, double b );
double sech_pdf ( double x, double a, double b );
double sech_sample ( double a, double b, int &seed );
double sech_variance ( double a, double b );

double semicircular_cdf ( double x, double a, double b );
double semicircular_cdf_inv ( double cdf, double a, double b );
bool semicircular_check ( double a, double b );
double semicircular_mean ( double a, double b );
double semicircular_pdf ( double x, double a, double b );
double semicircular_sample ( double a, double b, int &seed );
double semicircular_variance ( double a, double b );

double sin_power_int ( double a, double b, int n );

double sphere_unit_area_nd ( int dim_num );

int stirling2_value ( int n, int m );

double student_cdf ( double x, double a, double b, double c );
void student_cdf_values ( int &n_data, double &c, double &x, double &fx );
bool student_check ( double a, double b, double c );
double student_mean ( double a, double b, double c );
double student_pdf ( double x, double a, double b, double c );
double student_sample ( double a, double b, double c, int &seed );
double student_variance ( double a, double b, double c );

double student_noncentral_cdf ( double x, int idf, double d );
void student_noncentral_cdf_values ( int &n_data, int &df, double &lambda, 
  double &x, double &fx );

double tfn ( double x, double fx );

void timestamp ( );

double triangle_cdf ( double x, double a, double b, double c );
double triangle_cdf_inv ( double cdf, double a, double b, double c );
bool triangle_check ( double a, double b, double c );
double triangle_mean ( double a, double b, double c );
double triangle_pdf ( double x, double a, double b, double c );
double triangle_sample ( double a, double b, double c, int &seed );
double triangle_variance ( double a, double b, double c );

double triangular_cdf ( double x, double a, double b );
double triangular_cdf_inv ( double cdf, double a, double b );
bool triangular_check ( double a, double b );
double triangular_mean ( double a, double b );
double triangular_pdf ( double x, double a, double b );
double triangular_sample ( double a, double b, int &seed );
double triangular_variance ( double a, double b );

double trigamma ( double x );

double uniform_01_cdf ( double x );
double uniform_01_cdf_inv ( double cdf );
double uniform_01_mean ( );
double uniform_01_pdf ( double x );
double uniform_01_sample ( int &seed );
double uniform_01_variance ( );

double *uniform_01_order_sample ( int n, int &seed );

double uniform_cdf ( double x, double a, double b );
double uniform_cdf_inv ( double cdf, double a, double b );
bool uniform_check ( double a, double b );
double uniform_mean ( double a, double b );
double uniform_pdf ( double x, double a, double b );
double uniform_sample ( double a, double b, int &seed );
double uniform_variance ( double a, double b );

double uniform_discrete_cdf ( int x, int a, int b );
int uniform_discrete_cdf_inv ( double cdf, int a, int b );
bool uniform_discrete_check ( int a, int b );
double uniform_discrete_mean ( int a, int b );
double uniform_discrete_pdf ( int x, int a, int b );
int uniform_discrete_sample ( int a, int b, int &seed );
double uniform_discrete_variance ( int a, int b );

double *uniform_nsphere_sample ( int n, int &seed );

double von_mises_cdf ( double x, double a, double b );
double von_mises_cdf_inv ( double cdf, double a, double b );
void von_mises_cdf_values ( int &n_data, double &a, double &b, double &x, 
  double &fx );
bool von_mises_check ( double a, double b );
double von_mises_mean ( double a, double b );
double von_mises_pdf ( double x, double a, double b );
double von_mises_sample ( double a, double b, int &seed );
double von_mises_circular_variance ( double a, double b );

double weibull_cdf ( double x, double a, double b, double c );
double weibull_cdf_inv ( double cdf, double a, double b, double c );
void weibull_cdf_values ( int &n_data, double &alpha, double &beta, 
  double &x, double &fx );
bool weibull_check ( double a, double b, double c );
double weibull_mean ( double a, double b, double c );
double weibull_pdf ( double x, double a, double b, double c );
double weibull_sample ( double a, double b, double c, int &seed );
double weibull_variance ( double a, double b, double c );

double weibull_discrete_cdf ( int x, double a, double b );
int weibull_discrete_cdf_inv ( double cdf, double a, double b );
bool weibull_discrete_check ( double a, double b );
double weibull_discrete_pdf ( int x, double a, double b );
int weibull_discrete_sample ( double a, double b, int &seed );

double zeta ( double p );

double zipf_cdf ( int x, double a );
bool zipf_check ( double a );
double zipf_mean ( double a );
double zipf_pdf ( int x, double a );
int zipf_sample ( double a, int &seed );
double zipf_variance ( double a );
