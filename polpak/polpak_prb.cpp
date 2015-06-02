# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polpak.hpp"

int main ( );

void agud_test ( );
void align_enum_test ( );
void bell_test ( );
void benford_test ( );
void bernoulli_number_test ( );
void bernoulli_number2_test ( );
void bernoulli_number3_test ( );
void bernoulli_poly_test ( );
void bernoulli_poly2_test ( );
void bernstein_poly_test ( );
void bpab_test ( );
void cardan_poly_test ( );
void cardan_poly_coef_test ( );
void cardinal_cos_test ( );
void cardinal_sin_test ( );
void catalan_test ( );
void catalan_row_next_test ( );
void charlier_test ( );
void cheby_t_poly_test ( );
void cheby_t_poly_coef_test ( );
void cheby_t_poly_zero_test ( );
void cheby_u_poly_test ( );
void cheby_u_poly_coef_test ( );
void cheby_u_poly_zero_test ( );
void chebyshev_discrete_test ( );
void collatz_count_test ( );
void collatz_count_max_test ( );
void comb_row_next_test ( );
void commul_test ( );
void complete_symmetric_poly_test ( );
void cos_power_int_test ( );
void euler_number_test ( );
void euler_number2_test ( );
void euler_poly_test ( );
void eulerian_test ( );
void f_hofstadter_test ( );
void fibonacci_direct_test ( );
void fibonacci_floor_test ( );
void fibonacci_recursive_test ( );
void g_hofstadter_test ( );
void gegenbauer_poly_test ( );
void gen_hermite_poly_test ( );
void gen_laguerre_poly_test ( );
void gud_test ( );
void hail_test ( );
void h_hofstadter_test ( );
void hermite_poly_phys_test ( );
void hermite_poly_phys_coef_test ( );
void i4_choose_test ( );
void i4_factor_test ( );
void i4_factorial_test ( );
void i4_factorial2_test ( );
void i4_is_triangular_test ( );
void i4_partition_distinct_count_test ( );
void i4_to_triangle_test ( );
void jacobi_poly_test ( );
void jacobi_symbol_test ( );
void krawtchouk_test ( );
void laguerre_associated_test ( );
void laguerre_poly_test ( );
void laguerre_poly_coef_test ( );
void legendre_poly_test ( );
void legendre_poly_coef_test ( );
void legendre_associated_test ( );
void legendre_associated_normalized_test ( );
void legendre_function_q_test ( );
void legendre_symbol_test ( );
void lerch_test ( );
void lgamma_test ( );
void lock_test ( );
void meixner_test ( );
void mertens_test ( );
void moebius_test ( );
void motzkin_test ( );
void normal_01_cdf_inverse_test ( );
void omega_test ( );
void pentagon_num_test ( );
void phi_test ( );
void plane_partition_num_test ( );
void poly_bernoulli_test ( );
void poly_coef_count_test ( );
void prime_test ( );
void pyramid_num_test ( );
void pyramid_square_num_test ( );
void r8_agm_test ( );
void r8_beta_test ( );
void r8_choose_test ( );
void r8_erf_test ( );
void r8_erf_inverse_test ( );
void r8_euler_constant_test ( );
void r8_factorial_test ( );
void r8_factorial_log_test ( );
void r8_hyper_2f1_test ( );
void r8_psi_test ( );
void r8poly_degree_test ( );
void r8poly_print_test ( );
void r8poly_value_horner_test ( );
void sigma_test ( );
void simplex_num_test ( );
void sin_power_int_test ( );
void slice_test ( );
void spherical_harmonic_test ( );
void stirling1_test ( );
void stirling2_test ( );
void tau_test ( );
void tetrahedron_num_test ( );
void triangle_num_test ( );
void triangle_to_i4_test ( );
void trinomial_test ( );
void v_hofstadter_test ( );
void vibonacci_test ( );
void zeckendorf_test ( );
void zernike_poly_test ( );
void zernike_poly_coef_test ( );
void zeta_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLPAK_PRB.
//
//  Discussion:
//
//    POLPAK_PRB tests the POLPAK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POLPAK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POLPAK library.\n";

  agud_test ( );
  align_enum_test ( );
  bell_test ( );
  benford_test ( );
  bernoulli_number_test ( );
  bernoulli_number2_test ( );
  bernoulli_number3_test ( );
  bernoulli_poly_test ( );
  bernoulli_poly2_test ( );
  bernstein_poly_test ( );
  bpab_test ( );
  cardan_poly_test ( );
  cardan_poly_coef_test ( );
  cardinal_cos_test ( );
  cardinal_sin_test ( );
  catalan_test ( );
  catalan_row_next_test ( );
  charlier_test ( );
  cheby_t_poly_test ( );
  cheby_t_poly_coef_test ( );
  cheby_t_poly_zero_test ( );
  cheby_u_poly_test ( );
  cheby_u_poly_coef_test ( );
  cheby_u_poly_zero_test ( );
  chebyshev_discrete_test ( );
  collatz_count_test ( );
  collatz_count_max_test ( );
  comb_row_next_test ( );
  commul_test ( );
  complete_symmetric_poly_test ( );
  cos_power_int_test ( );
  euler_number_test ( );
  euler_number2_test ( );
  euler_poly_test ( );
  eulerian_test ( );
  f_hofstadter_test ( );
  fibonacci_direct_test ( );
  fibonacci_floor_test ( );
  fibonacci_recursive_test ( );
  g_hofstadter_test ( );
  gegenbauer_poly_test ( );
  gen_hermite_poly_test ( );
  gen_laguerre_poly_test ( );
  gud_test ( );
  hail_test ( );
  h_hofstadter_test ( );
  hermite_poly_phys_test ( );
  hermite_poly_phys_coef_test ( );
  i4_choose_test ( );
  i4_factor_test ( );
  i4_factorial_test ( );
  i4_factorial2_test ( );
  i4_is_triangular_test ( );
  i4_partition_distinct_count_test ( );
  i4_to_triangle_test ( );
  jacobi_poly_test ( );
  jacobi_symbol_test ( );
  krawtchouk_test ( );
  laguerre_associated_test ( );
  laguerre_poly_test ( );
  laguerre_poly_coef_test ( );
  legendre_poly_test ( );
  legendre_poly_coef_test ( );
  legendre_associated_test ( );
  legendre_associated_normalized_test ( );
  legendre_function_q_test ( );
  legendre_symbol_test ( );
  lerch_test ( );
  lgamma_test ( );
  lock_test ( );
  meixner_test ( );
  mertens_test ( );
  moebius_test ( );
  motzkin_test ( );
  normal_01_cdf_inverse_test ( );
  omega_test ( );
  pentagon_num_test ( );
  phi_test ( );
  plane_partition_num_test ( );
  poly_bernoulli_test ( );
  poly_coef_count_test ( );
  prime_test ( );
  pyramid_num_test ( );
  pyramid_square_num_test ( );
  r8_agm_test ( );
  r8_beta_test ( );
  r8_choose_test ( );
  r8_erf_test ( );
  r8_erf_inverse_test ( );
  r8_euler_constant_test ( );
  r8_factorial_test ( );
  r8_factorial_log_test ( );
  r8_hyper_2f1_test ( );
  r8_psi_test ( );
  r8poly_degree_test ( );
  r8poly_print_test ( );
  r8poly_value_horner_test ( );
  sigma_test ( );
  simplex_num_test ( );
  sin_power_int_test ( );
  slice_test ( );
  spherical_harmonic_test ( );
  stirling1_test ( );
  stirling2_test ( );
  tau_test ( );
  tetrahedron_num_test ( );
  triangle_num_test ( );
  triangle_to_i4_test ( );
  trinomial_test ( );
  v_hofstadter_test ( );
  vibonacci_test ( );
  zeckendorf_test ( );
  zernike_poly_test ( );
  zernike_poly_coef_test ( );
  zeta_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLPAK_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void agud_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    AGUD_TEST tests AGUD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double g;
  int i;
  double x;
  double x2;

  cout << "\n";
  cout << "AGUD_TEST\n";
  cout << "  AGUD computes the inverse Gudermannian;\n";
  cout << "\n";
  cout << "         X     GUD(X)     AGUD(GUD(X))\n";
  cout << "\n";

  for ( i = 0; i <= 10; i++ )
  {
    x = 1.0 + ( ( double ) i ) / 5.0;
    g = gud ( x );
    x2 = agud ( g );

    cout << "  " << setw(10) << x
         << "  " << setw(10) << g
         << "  " << setw(10) << x2    << "\n";
  }

  return;
}
//****************************************************************************80

void align_enum_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ALIGN_ENUM_TEST tests ALIGN_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M_MAX 10
# define N_MAX 10

  int i;
  int j;

  cout << "\n";
  cout << "ALIGN_ENUM_TEST\n";
  cout << "  ALIGN_ENUM counts the number of possible\n";
  cout << "  alignments of two biological sequences.\n";

  cout << "\n";
  cout << "  Alignment enumeration table:\n";
  cout << "\n";

  cout << "      ";
  for ( j = 0; j <= 5; j++ )
  {
    cout << setw(8) << j << "  ";
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i <= M_MAX; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j <= 5; j++ )
    {
      cout << setw(8) << align_enum ( i, j ) << "  ";
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "      ";
  for ( j = 6; j <= N_MAX; j++ )
  {
    cout << setw(8) << j << "  ";
  }
  cout << "\n";
  cout << "\n";

  for ( i = 0; i <= M_MAX; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 6; j <= N_MAX; j++ )
    {
      cout << setw(8) << align_enum ( i, j ) << "  ";
    }
    cout << "\n";
  }
  return;
# undef M_MAX
# undef N_MAX
}
//****************************************************************************80

void bell_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_TEST tests BELL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "BELL_TEST\n";
  cout << "  BELL computes Bell numbers.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bell_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    bell ( n, c2 );

    cout                     << "  "
         << setw(4) << n     << "  "
         << setw(8) << c     << "  "
         << setw(8) << c2[n] << "\n";

    delete [] c2;

  }

  return;
}
//****************************************************************************80

void benford_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BENFORD_TEST tests BENFORD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "BENFORD_TEST\n";
  cout << "  BENFORD(I) is the Benford probability of the\n";
  cout << "  initial digit sequence I.\n";
  cout << "\n";
  cout << "     I  BENFORD(I)\n";
  cout << "\n";

  for ( i = 1; i <= 9; i++ )
  {
    cout                              << "  "
         << setw(4) << i              << "  "
         << setw(10) << benford ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void bernoulli_number_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER_TEST tests BERNOULLI_NUMBER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  cout << "\n";
  cout << "BERNOULLI_NUMBER_TEST\n";
  cout << "  BERNOULLI_NUMBER computes Bernoulli numbers;\n";

  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( n_data, n, c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number ( n, c1 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(10) << c0    << "  "
         << setw(10) << c1[n] << "\n";
  }

  return;
}
//****************************************************************************80

void bernoulli_number2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER2_TEST tests BERNOULLI_NUMBER2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1[31];
  int n;
  int n_data;

  cout << "\n";
  cout << "BERNOULLI_NUMBER2_TEST\n";
  cout << "  BERNOULLI_NUMBER2 computes Bernoulli numbers;\n";
  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( n_data, n, c0 );

    if ( n_data == 0 )
    {
      break;
    }

    bernoulli_number2 ( n, c1 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(10) << c0    << "  "
         << setw(10) << c1[n] << "\n";
  }

  return;
}
//****************************************************************************80

void bernoulli_number3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER3_TEST tests BERNOULLI_NUMBER3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double c0;
  double c1;
  int n;
  int n_data;

  cout << "\n";
  cout << "BERNOULLI_NUMBER3_TEST\n";
  cout << "  BERNOULLI_NUMBER3 computes Bernoulli numbers.\n";
  cout << "\n";
  cout << "   I      Exact     BERNOULLI_NUMBER3\n";
  cout << "\n";
  
  n_data = 0;

  for ( ; ; )
  {
    bernoulli_number_values ( n_data, n, c0 );

    if ( n_data == 0 )
    {
      break;
    }

    c1 = bernoulli_number3 ( n );

    cout                   << "  "
         << setw(4)  << n  << "  "
         << setw(14) << c0 << "  "
         << setw(14) << c1 << "\n";

  }
 
  return;
}
//****************************************************************************80

void bernoulli_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY_TEST tests BERNOULLI_POLY;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;

  cout << "\n";
  cout << "BERNOULLI_POLY_TEST\n";
  cout << "  BERNOULLI_POLY evaluates Bernoulli polynomials;\n";
  cout << "\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  I          BX\n";
  cout << "\n";

  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly ( i, x );

    cout                   << "  "
         << setw(6)  << i  << "  "
         << setw(10) << bx << "\n";
  }

  return;
}
//****************************************************************************80

void bernoulli_poly2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY2_TEST tests BERNOULLI_POLY2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double bx;
  int i;
  int n = 15;
  double x;

  x = 0.2;
 
  cout << "\n";
  cout << "BERNOULLI_POLY2_TEST\n";
  cout << "  BERNOULLI_POLY2 evaluates Bernoulli polynomials.\n";
  cout << "\n";
  cout << "  X = " << x << "\n";
  cout << "\n";
  cout << "  I          BX\n";
  cout << "\n";
 
  for ( i = 1; i <= n; i++ )
  {
    bx = bernoulli_poly2 ( i, x );

    cout                   << "  "
         << setw(2)  << i  << "  "
         << setw(16) << bx << "\n";
  }
 
  return;
}
//****************************************************************************80

void bernstein_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_TEST tests BERNSTEIN_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double bvec[11];
  int k;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "BERNSTEIN_POLY_TEST:\n";
  cout << "  BERNSTEIN_POLY evaluates the Bernstein polynomials.\n";
  cout << "\n";
  cout << "   N   K   X   Exact   B(N,K)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bernstein_poly_values ( n_data, n, k, x, b );

    if ( n_data == 0 )
    {
      break;
    }

    bernstein_poly ( n, x, bvec );

    cout << "  " << setw(4)  << n
         << "  " << setw(4)  << k
         << "  " << setw(7)  << x
         << "  " << setw(14) << b
         << "  " << setw(14) << bvec[k] << "\n";
  }

  return;
}
//****************************************************************************80

void bpab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BPAB_TEST tests BPAB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a;
  double b;
  double bern[N+1];
  int i;
  double x;

  cout << "\n";
  cout << "BPAB_TEST\n";
  cout << "  BPAB evaluates Bernstein polynomials.\n";
  cout << "\n";

  x = 0.3;
  a = 0.0;
  b = 1.0;

  bpab ( N, x, a, b, bern );

  cout << "  The Bernstein polynomials of degree " << N << "\n";
  cout << "  based on the interval from " << a << "\n";
  cout << "  to " << b << "\n";
  cout << "  evaluated at X = " << x << "\n";
  cout << "\n";

  for ( i = 0; i <= N; i++ )
  {
    cout << "  " << setw(4)  << i       
         << "  " << setw(14) << bern[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void cardan_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CARDAN_POLY_TEST tests CARDAN_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  double *cx2;
  int i;
  int n;
  double s;
  double x;

  cout << "\n";
  cout << "CARDAN_POLY_TEST\n";
  cout << "  CARDAN_POLY evaluates a Cardan polynomial directly.\n";
  cout << "\n";

  n = N_MAX;
  x = 0.25;
  s = 0.5;

  cout << "\n";
  cout << "  Compare CARDAN_POLY_COEF + R8POLY_VALUE_HORNER\n";
  cout << "  versus CARDAN_POLY alone.\n";
  cout << "\n";
  cout << "  Evaluate polynomials at X = " << x << "\n";
  cout << "  We use the parameter S = " << s << "\n";
  cout << "\n";
  cout << "  Order       Horner          Direct\n";
  cout << "\n";

  cx2 = cardan_poly ( n, x, s );

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );

    cx1 = r8poly_value_horner ( n, c, x );

    cout << "  " << setw(2)  << n
         << "  " << setw(14) << cx1
         << "  " << setw(14) << cx2[n] << "\n";
  }
  delete [] cx2;

  return;
# undef N_MAX
}
//****************************************************************************80

void cardan_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CARDAN_POLY_COEF_TEST tests CARDAN_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  double c[N_MAX+1];
  double cx1;
  double *cx2;
  int i;
  int n;
  double s;
  double x;

  s = 1.0;

  cout << "\n";
  cout << "CARDAN_POLY_COEF_TEST\n";
  cout << "  CARDAN_POLY_COEF returns the coefficients of a\n";
  cout << "  Cardan polynomial.\n";
  cout << "\n";
  cout << "  We use the parameter S = " << s << "\n";
  cout << "\n";
  cout << "  Table of polynomial coefficients:\n";
  cout << "\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    cardan_poly_coef ( n, s, c );
    cout << "  "
         << setw(2) << n << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(5) << c[i] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cardinal_cos_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CARDINAL_COS_TEST tests CARDINAL_COS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int m = 11;
  const double r8_pi = 3.141592653589793;
  double *t;

  cout << "\n";
  cout << "CARDINAL_COS_TEST\n";
  cout << "  CARDINAL_COS evaluates cardinal cosine functions.\n";
  cout << "  Ci(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).\n";
  cout << "  A simple check of all pairs should form the identity matrix.\n";

  cout << "\n";
  cout << "  The CARDINAL_COS test matrix:\n";
  cout << "\n";

  t = r8vec_linspace_new ( m + 2, 0.0, r8_pi );

  for ( j = 0; j <= m + 1; j++ )
  {
    c = cardinal_cos ( j, m, m + 2, t );
    for ( i = 0; i <= m + 1; i++ )
    {
      cout << "  " << setw(4) << c[i];
    }
    cout << "\n";
    delete [] c;
  }

  delete [] t;

  return;
}
//****************************************************************************80

void cardinal_sin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CARDINAL_SIN_TEST tests CARDINAL_SIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 May 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m = 11;
  const double r8_pi = 3.141592653589793;
  double *s;
  double *t;

  cout << "\n";
  cout << "CARDINAL_SIN_TEST\n";
  cout << "  CARDINAL_SIN evaluates cardinal sine functions.\n";
  cout << "  Si(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).\n";
  cout << "  A simple check of all pairs should form the identity matrix.\n";

  t = r8vec_linspace_new ( m + 2, 0.0, r8_pi );

  cout << "\n";
  cout << "  The CARDINAL_SIN test matrix:\n";
  cout << "\n";
  for ( j = 0; j <= m + 1; j++ )
  {
    s = cardinal_sin ( j, m, m + 2, t );
    for ( i = 0; i <= m + 1; i++ )
    {
      cout << "  " << setw(4) << s[i];
    }
    cout << "\n";
    delete [] s;
  }

  delete [] t;

  return;
}
//****************************************************************************80

void catalan_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_TEST tests CATALAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "CATALAN_TEST\n";
  cout << "  CATALAN computes Catalan numbers.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = new int[n+1];

    catalan ( n, c2 );

    cout << "  " << setw(4) << n
         << "  " << setw(8) << c
         << "  " << setw(8) << c2[n] << "\n";

    delete [] c2;
  }

  return;
}
//****************************************************************************80

void catalan_row_next_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_ROW_NEXT_TEST tests CATALAN_ROW_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;
  bool next;

  cout << "\n";
  cout << "CATALAN_ROW_NEXT_TEST\n";
  cout << "  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.\n";
  cout << "\n";
  cout << "  First, compute row 7:\n";
  cout << "\n";

  next = false;
  n = 7;
  catalan_row_next ( next, n, c );

  cout << setw(4) << n << "  ";
  for ( i = 0; i <= n; i++ )
  {
    cout << setw(8) << c[i] << "  ";
  }
  cout << "\n";

  cout << "\n";
  cout << "  Now compute rows consecutively, one at a time:\n";
  cout << "\n";

  next = false;

  for ( n = 0; n <= N_MAX; n++ )
  {
    catalan_row_next ( next, n, c );
    next = true;

    cout << setw(4) << i << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(8) << c[i] << "  ";
    }
    cout << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void charlier_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHARLIER_TEST tests CHARLIER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5
# define N 5

  double a;
  double a_test[TEST_NUM] = { 0.25, 0.5, 1.0, 2.0, 10.0 };
  int i;
  int j;
  int n;
  int test;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "CHARLIER_TEST:\n";
  cout << "  CHARLIER evaluates Charlier polynomials.\n";
  cout << "\n";
  cout << "       N      A         X        P(N,A,X)\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    a = a_test[test];

    cout << "\n";

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      charlier ( n, a, x, value );

      cout << "\n";
      for ( i = 0; i <= 5; i++ )
      {

        cout << "  " << setw(6)  << i     
             << "  " << setw(8)  << a
             << "  " << setw(8)  << x
             << "  " << setw(14) << value[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void cheby_t_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_TEST tests CHEBY_T_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "CHEBY_T_POLY_TEST:\n";
  cout << "  CHEBY_T_POLY evaluates the Chebyshev T polynomial.\n";
  cout << "\n";
  cout << "     N      X        Exact F       T(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_t_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_t_poly ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_t_poly_zero_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_ZERO_TEST tests CHEBY_T_POLY_ZERO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  cout << "\n";
  cout << "CHEBY_T_POLY_ZERO_TEST:\n";
  cout << "  CHEBY_T_POLY_ZERO returns zeroes of T(N,X).\n";
  cout << "\n";
  cout << "       N      X        T(N,X)\n";
  cout << "\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_t_poly_zero ( n );
    fx = cheby_t_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[i+n*n] << "\n";
    }
    cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_t_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_COEF_TEST tests CHEBY_T_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int j;
  int n = 5;

  cout << "\n";
  cout << "CHEBY_T_POLY_COEF_TEST\n";
  cout << "  CHEBY_T_POLY_COEF determines the  polynomial coefficients\n";
  cout << "  of the Chebyshev polynomial T(n,x).\n";

  c = cheby_t_poly_coef ( n );
 
  for ( i = 0; i <= n; i++ )
  {
    cout << "\n";
    cout << "  T(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          cout << setw(14) << c[i+j*(n+1)] << "\n";;
        }
        else if ( j == 1 )
        {
          cout << setw(14) << c[i+j*(n+1)] << " * x\n";
        }
        else
        {
          cout << setw(14) << c[i+j*(n+1)] << " * x^" << j << "\n";
        }
      }
    }
  }
 
  delete [] c;

  return;
}
//****************************************************************************80

void cheby_u_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_TEST tests CHEBY_U_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2015
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double *fx2;
  int n;
  int n_data;
  double x;
  double x_vec[1];

  cout << "\n";
  cout << "CHEBY_U_POLY_TEST:\n";
  cout << "  CHEBY_U_POLY evaluates the Chebyshev U polynomial.\n";
  cout << "\n";
  cout << "     N      X        Exact F       U(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cheby_u_poly_values ( n_data, n, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2 = cheby_u_poly ( 1, n, x_vec );

    cout << "  " << setw(8)  << n
         << "  " << setw(8)  << x
         << "  " << setw(14) << fx
         << "  " << setw(14) << fx2[n] << "\n";

    delete [] fx2;

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_u_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_COEF_TEST tests CHEBY_U_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2012
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "CHEBY_U_POLY_COEF_TEST\n";
  cout << "  CHEBY_U_POLY_COEF determines the polynomial coefficients\n";
  cout << "  of the Chebyshev polynomial U(n,x).\n";

  cheby_u_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  U(" << i << ",x)\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(N+1)] != 0.0 )
      {
        if ( j == 0 )
        {
          cout << setw(14) << c[i+j*(N+1)] << "\n";
        }
        else if ( j == 1 )
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x\n";
        }
        else
        {
          cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
        }
      }
    }
  }
 
  return;
# undef N
}
//****************************************************************************80

void cheby_u_poly_zero_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_ZERO_TEST tests CHEBY_U_POLY_ZERO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2015
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 4

  double *fx;
  int i;
  int n;
  double *z;

  cout << "\n";
  cout << "CHEBY_U_POLY_ZERO_TEST:\n";
  cout << "  CHEBY_U_POLY_ZERO returns zeroes of U(N,X).\n";
  cout << "\n";
  cout << "       N      X        U(N,X)\n";
  cout << "\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    z = cheby_u_poly_zero ( n );
    fx = cheby_u_poly ( n, n, z );
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << n
           << "  " << setw(8) << z[i]
           << "  " << setw(14) << fx[i+n*n] << "\n";
    }
    cout << "\n";
    delete [] fx;
    delete [] z;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void chebyshev_discrete_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_DISCRETE_TEST tests CHEBYSHEV_DISCRETE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 5
# define N 5

  int i;
  int j;
  int m;
  int n;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "CHEBYSHEV_DISCRETE_TEST:\n";
  cout << "  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.\n";
  cout << "\n";
  cout << "       N      M         X        T(N,M,X)\n";

  m = 5;
  n = N;

  for ( j = 0; j <= 5; j++ )
  {
    x = ( double ) ( j ) / 2.0;

    chebyshev_discrete ( n, m, x, value );

    cout << "\n";
    for ( i = 0; i <= 5; i++ )
    {
      cout << "  " << setw(6)  << i     
           << "  " << setw(8)  << m
           << "  " << setw(8)  << x
           << "  " << setw(14) << value[i] << "\n";
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void collatz_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_TEST tests COLLATZ_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int count;
  int count2;
  int n;
  int n_data;

  cout << "\n";
  cout << "COLLATZ_COUNT_TEST:\n";
  cout << "  COLLATZ_COUNT(N) counts the length of the\n";
  cout << "  Collatz sequence beginning with N.\n";
  cout << "\n";
  cout << "       N       COUNT(N)     COUNT(N)\n";
  cout << "              (computed)    (table)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    collatz_count_values ( &n_data, &n, &count );

    if ( n_data == 0 )
    {
      break;
    }

    count2 = collatz_count ( n );

    cout << "  " << setw(8) << n
         << "  " << setw(8) << count
         << "  " << setw(8) << count2 << "\n";
  }

  return;
}
//****************************************************************************80

void collatz_count_max_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_MAX_TEST tests COLLATZ_COUNT_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i_max;
  int j_max;
  int n;

  cout << "\n";
  cout << "COLLATZ_COUNT_MAX_TEST:\n";
  cout << "  COLLATZ_COUNT_MAX(N) returns the length of the\n";
  cout << "  longest Collatz sequence from 1 to N.\n";
  cout << "\n";
  cout << "         N     I_MAX     J_MAX\n";
  cout << "\n";

  n = 10;

  while ( n <= 100000 )
  {
    collatz_count_max ( n, &i_max, &j_max );

    cout << "  " << setw(8) << n
         << "  " << setw(8) << i_max
         << "  " << setw(8) << j_max << "\n";

    n = n * 10;
  }

  return;
}
//****************************************************************************80

void comb_row_next_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMB_ROW_NEXT_TEST tests COMB_ROW_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 December 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 10

  int c[N_MAX+1];
  int i;
  int n;

  cout << "\n";
  cout << "COMB_ROW_NEXT_TEST\n";
  cout << "  COMB_ROW_NEXT computes the next row of Pascal's triangle.\n";
  cout << "\n";

  for ( n = 0; n <= N_MAX; n++ )
  {
    comb_row_next ( n, c );
    cout << "  " << setw(2) << n << "  ";
    for ( i = 0; i <= n; i++ )
    {
      cout << setw(5) << c[i];
    }
    cout << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void commul_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMMUL_TEST tests COMMUL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int factor[4];
  int i;
  int ncomb;
  int nfactor;

  cout << "\n";
  cout << "COMMUL_TEST\n";
  cout << "  COMMUL computes a multinomial coefficient.\n";
  cout << "\n";

  n = 8;
  nfactor = 2;
  factor[0] = 6;
  factor[1] = 2;
  ncomb = commul ( n, nfactor, factor );
  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of factors = " << factor << "\n";
  for ( i = 0; i < nfactor; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(8) << factor[i] << "\n";
  }
  cout << "  Value of coefficient = " << ncomb << "\n";

  n = 8;
  nfactor = 3;
  factor[0] = 2;
  factor[1] = 2;
  factor[2] = 4;
  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of factors = " << factor << "\n";
  for ( i = 0; i < nfactor; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(8) << factor[i] << "\n";
  }
  cout << "  Value of coefficient = " << ncomb << "\n";

  n = 13;
  nfactor = 4;
  factor[0] = 5;
  factor[1] = 3;
  factor[2] = 3;
  factor[3] = 2;
  ncomb = commul ( n, nfactor, factor );
  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  Number of factors = " << factor << "\n";
  for ( i = 0; i < nfactor; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(8) << factor[i] << "\n";
  }
  cout << "  Value of coefficient = " << ncomb << "\n";

  return;
}
//****************************************************************************80

void complete_symmetric_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COMPLETE_SYMMETRIC_POLY_TEST tests COMPLETE_SYMMETRIC_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int n = 5;
  int nn;
  int r;
  int rr;
  double value;
  double x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  cout << "\n";
  cout << "COMPLETE_SYMMETRIC_POLY_TEST\n";
  cout << "  COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric.\n";
  cout << "  polynomial in a given set of variables X.\n";
 
  r8vec_print ( n, x, "  Variable vector X:" );

  cout << "\n";
  cout << "   N\\R     0       1       2       3       4       5\n";
  cout << "\n";

  for ( nn = 0; nn <= n; nn++ )
  {
    cout << "  " << setw(2) <<  nn;
    for ( rr = 0; rr <= 5; rr++ )
    {
      value = complete_symmetric_poly ( nn, rr, x );
      cout << "  " << setw(6) << value;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void cos_power_int_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    COS_POWER_INT_TEST tests COS_POWER_INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  cout << "\n";
  cout << "COS_POWER_INT_TEST:\n";
  cout << "  COS_POWER_INT computes the integral of the N-th power\n";
  cout << "  of the cosine function.\n";
  cout << "\n";
  cout << "         A         B       N        Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    cos_power_int_values ( n_data, a, b, n, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = cos_power_int ( a, b, n );

    cout                    << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(6)  << n   << "  "
         << setw(12) << fx  << "  "
         << setw(12) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void euler_number_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER_TEST tests EULER_NUMBER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  int c2[13];
  int n;
  int n_data;

  cout << "\n";
  cout << "EULER_NUMBER_TEST\n";
  cout << "  EULER_NUMBER computes Euler numbers.\n";
  cout << "\n";
  cout << "  N  exact   EULER_NUMBER\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    euler_number ( n, c2 );

    cout                      << "  "
         << setw(4)  << n     << "  "
         << setw(12) << c1    << "  "
         << setw(12) << c2[n] << "\n";

  }
 
  return;
}
//****************************************************************************80

void euler_number2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER2_TEST tests EULER_NUMBER2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c1;
  double c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "EULER_NUMBER2_TEST\n";
  cout << "  EULER_NUMBER2 computes Euler numbers.\n";
  cout << "\n";
  cout << "  N  exact   EULER_NUMBER2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    euler_number_values ( &n_data, &n, &c1 );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = euler_number2 ( n );

    cout                   << "  "
         << setw(4)  << n  << "  "
         << setw(12) << c1 << "  "
         << setw(14) << c2 << "\n";

  }
 
  return;
}
//****************************************************************************80

void euler_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_POLY_TEST tests EULER_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  int i;
  int n = 15;
  double x;

  x = 0.5;
 
  cout << "\n";
  cout << "EULER_POLY_TEST\n";
  cout << "  EULER_POLY evaluates Euler polynomials.\n";
  cout << "\n";
  cout << "  N         X              F(X)\n";
  cout << "\n";
   
  for ( i = 0; i <= n; i++ )
  {
    f = euler_poly ( i, x );

    cout                  << "  "
         << setw(2)  << i << "  "
         << setw(14) << x << "  "
         << setw(14) << f << "\n";
  }
 
  return;
}
//****************************************************************************80

void eulerian_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EULERIAN_TEST tests EULERIAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 7

  int e[N*N];
  int i;
  int j;

  cout << "\n";
  cout << "EULERIAN_TEST\n";
  cout << "  EULERIAN evaluates Eulerian numbers.\n";
  cout << "\n";
 
  eulerian ( N, e );

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      cout << setw(6) << e[i+j*N] << "  ";
    }
    cout << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void f_hofstadter_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    F_HOFSTADTER_TEST tests F_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;

  cout << "\n";
  cout << "F_HOFSTADTER_TEST\n";
  cout << "  F_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  F function.\n";
  cout << "\n";
  cout << "     N   F(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    f = f_hofstadter ( i );

    cout                 << "  "
         << setw(6) << i << "  "
         << setw(6) << f << "\n";
  }

  return;
}
//****************************************************************************80

void fibonacci_direct_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_DIRECT_TEST tests FIBONACCI_DIRECT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;
  int n = 20;

  cout << "\n";
  cout << "FIBONACCI_DIRECT_TEST\n";
  cout << "  FIBONACCI_DIRECT evalutes a Fibonacci number directly.\n";
  cout << "\n";
  
  for ( i = 1; i <= n; i++ )
  {
    f = fibonacci_direct ( i );

    cout                  << "  "
         << setw(6)  << i << "  "
         << setw(10) << f << "\n";
  }
 
  return;
}
//****************************************************************************80

void fibonacci_floor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_FLOOR_TEST tests FIBONACCI_FLOOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int f;
  int i;
  int n;

  cout << "\n";
  cout << "FIBONACCI_FLOOR_TEST\n";
  cout << "  FIBONACCI_FLOOR computes the largest Fibonacci number\n";
  cout << "  less than or equal to a given positive integer.\n";
  cout << "\n";
  cout << "     N  Fibonacci  Index\n";
  cout << "\n";

  for ( n = 1; n <= 20; n++ )
  {
    fibonacci_floor ( n, &f, &i );

    cout                 << "  "
         << setw(6) << n << "  "
         << setw(6) << f << "  "
         << setw(6) << i << "\n";
  }
 
  return;
}
//****************************************************************************80

void fibonacci_recursive_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_RECURSIVE_TEST tests FIBONACCI_RECURSIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20

  int f[N];
  int i;

  cout << "\n";
  cout << "FIBONACCI_RECURSIVE_TEST\n";
  cout << "  FIBONACCI_RECURSIVE computes the Fibonacci sequence.\n";
  cout << "\n";
 
  fibonacci_recursive ( N, f );
 
  for ( i = 1; i <= N; i++ )
  {
    cout                       << "  "
         << setw(6)  << i      << "  "
         << setw(10) << f[i-1] << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void g_hofstadter_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    G_HOFSTADTER_TEST tests G_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "G_HOFSTADTER_TEST\n";
  cout << "  G_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  G function.\n";
  cout << "\n";
  cout << "     N   G(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << g_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void gegenbauer_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLY_TEST tests GEGENBAUER_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "GEGENBAUER_POLY_TEST\n";
  cout << "  GEGENBAUER_POLY evaluates the Gegenbauer polynomials.\n";
  cout << "\n";
  cout << "        N       A       X       GPV      GEGENBAUER\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {

    gegenbauer_poly_values ( &n_data, &n, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = new double[n+1];

    gegenbauer_poly ( n, a, x, c );
    fx2 = c[n];

    cout                    << "  "
         << setw(6)  << n   << "  "
         << setw(10) << a   << "  "
         << setw(10) << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";

    delete [] c;
  }
 
  return;
}
//****************************************************************************80

void gen_hermite_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_POLY_TEST tests GEN_HERMITE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2015
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  double mu;
  double mu_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  cout << "\n";
  cout << "GEN_HERMITE_POLY_TEST\n";
  cout << "  GEN_HERMITE_POLY evaluates the generalized Hermite\n";
  cout << "  polynomial.\n";

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    mu = mu_test[i];

    cout << "\n";
    cout << "  Table of H(N,MU)(X) for\n";
    cout << "\n";
    cout << "    N(max) = " << N << "\n";
    cout << "    MU =     " << mu << "\n";
    cout << "    X =      " << x << "\n";
    cout << "\n";
  
    gen_hermite_poly ( N, x, mu, c );
 
    for ( j = 0; j <= N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(14) << c[j] << "\n";
    }
  }
 
  return;
# undef N
# undef N_TEST
}
//****************************************************************************80

void gen_laguerre_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_POLY_TEST tests GEN_LAGUERRE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define N_TEST 6

  double alpha;
  double alpha_test[N_TEST] = { 0.0, 0.0, 0.1, 0.1, 0.5, 1.0 };
  double c[N+1];
  int i;
  int j;
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  cout << "\n";
  cout << "GEN_LAGUERRE_POLY_TEST\n";
  cout << "  GEN_LAGUERRE_POLY evaluates the generalized Laguerre\n";
  cout << "  functions.\n";

  for ( i = 0; i < N_TEST; i++ )
  {

    x = x_test[i];
    alpha = alpha_test[i];

    cout << "\n";
    cout << "  Table of L(N,ALPHA)(X) for\n";
    cout << "\n";
    cout << "    N(max) = " << N << "\n";
    cout << "    ALPHA =  " << alpha << "\n";
    cout << "    X =      " << x << "\n";
    cout << "\n";
  
    gen_laguerre_poly ( N, alpha, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(14) << c[j] << "\n";
    }
  }
 
  return;
# undef N
# undef N_TEST
}
//****************************************************************************80

void gud_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GUD_TEST tests GUD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "GUD_TEST:\n";
  cout << "  GUD evaluates the Gudermannian function.\n";
  cout << "\n";
  cout << "     X      Exact F       GUD(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gud_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gud ( x );

    cout                    << "  "
         << setw(10) << x   << "  "
         << setw(10) << fx  << "  "
         << setw(10) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void hail_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HAIL_TEST tests HAIL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "HAIL_TEST\n";
  cout << "  HAIL(I) computes the length of the hail sequence\n";
  cout << "  for I, also known as the 3*N+1 sequence.\n";
  cout << "\n";
  cout << "  I,  HAIL(I)\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    cout                          << "  "
         << setw(4) << i          << "  "
         << setw(6) << hail ( i ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void h_hofstadter_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    H_HOFSTADTER_TEST tests H_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;

  cout << "\n";
  cout << "H_HOFSTADTER_TEST\n";
  cout << "  H_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  H function.\n";

  cout << "\n";
  cout << "     N   H(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << h_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void hermite_poly_phys_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PHYS_TEST tests HERMITE_POLY_PHYS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "HERMITE_POLY_PHYS_TEST:\n";
  cout << "  HERMITE_POLY_PHYS evaluates the physicist's Hermite polynomial.\n";
  cout << "\n";
  cout << "     N      X        Exact F       H(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hermite_poly_phys_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    hermite_poly_phys ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hermite_poly_phys_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PHYS_COEF_TEST tests HERMITE_POLY_PHYS_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "HERMITE_POLY_PHYS_COEF_TEST\n";
  cout << "  HERMITE_POLY_PHYS_COEF: physicist's Hermite polynomial coefficients.\n";

  hermite_poly_phys_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  H(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }

  return;
# undef N
}
//****************************************************************************80

void i4_choose_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE_TEST tests I4_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int cnk;
  int k;
  int n;

  cout << "\n";
  cout << "I4_CHOOSE_TEST\n";
  cout << "  I4_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "   N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i4_choose ( n, k );

      cout                   << "  "
           << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void i4_factor_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTOR_TEST tests I4_FACTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int maxfactor = 10;
  int n;
  int n_test[3] = { 60, 664048, 8466763 };
  int nfactor;
  int nleft;
  int factor[10];
  int power[10];

  cout << "\n";
  cout << "I4_FACTOR_TEST:\n";
  cout << "  I4_FACTOR tries to factor an I4\n";

  for ( i = 0; i < 3; i++ )
  {
    n = n_test[i];
    i4_factor ( n, maxfactor, nfactor, factor, power, nleft );
    cout << "\n";
    cout << "  Factors of N = " << n << "\n";
    for ( j = 0; j < nfactor; j++ )
    {
      cout << "    " << factor[j] << "^" <<  power[j] << "\n";
    }
    if ( nleft != 1 )
    {
      cout << "  Unresolved factor NLEFT = " << nleft << "\n";
    }
  }

  return;
}
//****************************************************************************80

void i4_factorial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL_TEST tests I4_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int fn2;
  int n;
  int n_data;

  cout << "\n";
  cout << "I4_FACTORIAL_TEST:\n";
  cout << "  I4_FACTORIAL evaluates the factorial function.\n";
  cout << "\n";
  cout << "     X       Exact F       I4_FACTORIAL(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial ( n );

    cout                    << "  "
         << setw(4)  << n   << "  "
         << setw(12) << fn  << "  "
         << setw(12) << fn2 << "\n";

  }

  return;
}
//****************************************************************************80

void i4_factorial2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL2_TEST tests I4_FACTORIAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int fn;
  int fn2;
  int n;
  int n_data;

  cout << "\n";
  cout << "I4_FACTORIAL2_TEST:\n";
  cout << "  I4_FACTORIAL2 evaluates the double factorial function.\n";
  cout << "\n";
  cout << "   N   Exact  I4_FACTORIAL2(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    i4_factorial2_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    fn2 = i4_factorial2 ( n );

    cout                   << "  "
         << setw(4) << n   << "  "
         << setw(8) << fn  << "  "
         << setw(8) << fn2 << "\n";
  }

  return;
}
//****************************************************************************80

void i4_is_triangular_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_IS_TRIANGULAR_TEST tests I4_IS_TRIANGULAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i2;
  int j;
  int k;
  int k2;
  bool l;

  cout << "\n";
  cout << "I4_IS_TRIANGULAR_TEST\n";
  cout << "  I4_IS_TRIANGULAR returns 0 or 1 depending on\n";
  cout << "  whether I is triangular.\n";
  cout << "\n";
  cout << "   I  =>   0/1\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    l = i4_is_triangular ( i );

    cout << "  " << setw(4) << i  
         << "  " << setw(1) << l << "\n";
  }
 
  return;
}
//****************************************************************************80

void i4_partition_distinct_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_PARTITION_DISTINCT_COUNT_TEST tests I4_PARTITION_DISTINCT_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int c2;
  int n;
  int n_data;
  int n_max = 20;

  cout << "\n";
  cout << "I4_PARTITION_DISTINCT_COUNT_TEST:\n";
  cout << "  For the number of partitions of an integer\n";
  cout << "  into distinct parts,\n";
  cout << "  I4_PARTITION_DISTINCT_COUNT computes any value.\n";
  cout << "\n";
  cout << "     N       Exact F    Q(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    partition_distinct_count_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    if ( n_max < n )
    {
      continue;
    }

    c2 = i4_partition_distinct_count ( n );

    cout                   << "  "
         << setw(10) << n  << "  "
         << setw(10) << c  << "  "
         << setw(10) << c2 << "\n";

  }

  return;
}
//****************************************************************************80

void i4_to_triangle_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_TRIANGLE_TEST tests I4_TO_TRIANGLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "I4_TO_TRIANGLE_TEST\n";
  cout << "  I4_TO_TRIANGLE converts a linear index to a\n";
  cout << "  triangular one.\n";
  cout << "\n";
  cout << "     K  => I     J\n";
  cout << "\n";

  for ( k = 0; k <= 20; k++ )
  {
    i4_to_triangle ( k, &i, &j );

    cout << "  " << setw(4) << k  
         << "  " << setw(4) << i
         << "  " << setw(4) << j << "\n";
  }
 
  return;
}
//****************************************************************************80

void jacobi_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY_TEST tests JACOBI_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *c;
  double fx;
  double fx2;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "JACOBI_POLY_TEST\n";
  cout << "  JACOBI_POLY evaluates the Jacobi polynomials.\n";
  cout << "  the Jacobi polynomials.\n";
  cout << "\n";
  cout << "        N       A       B      X       JPV      JACOBI\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {

    jacobi_poly_values ( n_data, n, a, b, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    c = jacobi_poly ( n, a, b, x );
    fx2 = c[n];

    cout                    << "  "
         << setw(6)  << n   << "  "
         << setw(6)  << a   << "  "
         << setw(6)  << b   << "  "
         << setw(10) << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";

    delete [] c;
  }
 
  return;
}
//****************************************************************************80

void jacobi_symbol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SYMBOL_TEST tests JACOBI_SYMBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 4

  int i;
  int p;
  int ptest[N_TEST] = { 3, 9, 10, 12 };
  int q;

  cout << "\n";
  cout << "JACOBI_SYMBOL_TEST\n";
  cout << "  JACOBI_SYMBOL computes the Jacobi symbol\n";
  cout << "  (Q/P), which records if Q is a quadratic\n";
  cout << "  residue modulo the number P.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    cout << "\n";
    cout << "Jacobi Symbols for P = " << p << "\n";
    cout << "\n";
    for ( q = 0; q <= p; q++ )
    {
      cout << "  " << setw(8) << p
           << "  " << setw(8) << q
           << "  " << setw(8) << jacobi_symbol ( q, p ) << "\n";
    }
  }

  return;
# undef N_TEST
}
//****************************************************************************80

void krawtchouk_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KRAWTCHOUK_TEST tests KRAWTCHOUK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 2
# define N 5

  int i;
  int j;
  int m;
  int n;
  double p;
  double p_test[TEST_NUM] = { 0.25, 0.5 };
  int test;
  double x;
  double value[N+1];

  cout << "\n";
  cout << "KRAWTCHOUK_TEST:\n";
  cout << "  KRAWTCHOUK evaluates Krawtchouk polynomials.\n";
  cout << "\n";
  cout << "        N         P         X          M      K(N,P,X,M)\n";
  cout << "\n";

  m = 5;
  n = N;

  for ( test = 0; test < TEST_NUM; test++ )
  {
    p = p_test[test];

    cout << "\n";

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      krawtchouk ( n, p, x, m, value );

      cout << "\n";
      for ( i = 0; i <= 5; i++ )
      {

        cout << "  " << setw(8)  << i     
             << "  " << setw(8)  << p
             << "  " << setw(8)  << x
             << "  " << setw(8)  << m
             << "  " << setw(14) << value[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void laguerre_associated_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_ASSOCIATED_TEST tests LAGUERRE_ASSOCIATED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
# define N_TEST 6

  double c[N+1];
  int i;
  int j;
  int m;
  int m_test[N_TEST] = { 0, 0, 1, 2, 3, 1 };
  double x;
  double x_test[N_TEST] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.5 };

  cout << "\n";
  cout << "LAGUERRE_ASSOCIATED_TEST\n";
  cout << "  LAGUERRE_ASSOCIATED evaluates the associated Laguerre\n";
  cout << "  polynomials.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    m = m_test[i];
    x = x_test[i];

    cout << "\n";
    cout << "  Table of L(N,M)(X) for\n";
    cout << "\n";
    cout << "  N(max) = " << N << "\n";
    cout << "  M      = " << m << "\n";
    cout << "  X =      " << x << "\n";
    cout << "\n";
 
    laguerre_associated ( N, m, x, c );
 
    for ( j = 0; j <= N; j++ )
    {
      cout                     << "  "
           << setw(6)  << j    << "  "
           << setw(14) << c[j] << "\n";
    }
 
  }

  return;
# undef N
# undef N_TEST
}
//****************************************************************************80

void laguerre_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLY_TEST tests LAGUERRE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "LAGUERRE_POLY_COEF:\n";
  cout << "  LAGUERRE_POLY evaluates the Laguerre polynomial.\n";
  cout << "\n";
  cout << "     N      X        Exact F       L(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    laguerre_polynomial_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    laguerre_poly ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void laguerre_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests LAGUERRE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  double fact;
  int i;
  int j;

  cout << "\n";
  cout << "LAGUERRE_POLY_COEF_TEST\n";
  cout << "  LAGUERRE_POLY_COEF determines Laguerre \n";
  cout << "  polynomial coefficients.\n";

  laguerre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  L(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
 
  for ( i = 0; i <= N; i++ )
  {
    fact = r8_factorial ( i );
    cout << "\n";
    cout << "  Factorially scaled L(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << fact *c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
  return;
# undef N
}
//****************************************************************************80

void legendre_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_TEST tests LEGENDRE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fp2[N_MAX+1];
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "LEGENDRE_POLY_TEST:\n";
  cout << "  LEGENDRE_POLY evaluates the Legendre PN function.\n";
  cout << "\n";
  cout << "     N      X        Exact F       P(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_poly_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_poly ( n, x, fx2, fp2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_COEF_TEST tests LEGENDRE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double c[(N+1)*(N+1)];
  int i;
  int j;

  cout << "\n";
  cout << "LEGENDRE_POLY_COEF_TEST\n";
  cout << "  LEGENDRE_POLY_COEF determines the Legendre P \n";
  cout << "  polynomial coefficients.\n";

  legendre_poly_coef ( N, c );
 
  for ( i = 0; i <= N; i++ )
  {
    cout << "\n";
    cout << "  P(" << i << ")\n";
    cout << "\n";
    for ( j = i; 0 <= j; j-- )
    {
      if ( j == 0 )
      {
        cout << setw(14) << c[i+j*(N+1)] << "\n";;
      }
      else if ( j == 1 )
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x\n";
      }
      else
      {
        cout << setw(14) << c[i+j*(N+1)] << " * x^" << j << "\n";
      }
    }
  }
 
  return;
# undef N
}
//****************************************************************************80

void legendre_associated_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_TEST tests LEGENDRE_ASSOCIATED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "LEGENDRE_ASSOCIATED_TEST:\n";
  cout << "  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.\n";
  cout << "\n";
  cout << "      N       M    X     Exact F     PNM(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_values ( &n_data, &n, &m, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated ( n, m, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << m      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated_normalized_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_TEST tests LEGENDRE_ASSOCIATED_NORMALIZED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 February 2015
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double fx2[N_MAX+1];
  double fx;
  int m;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "LEGENDRE_ASSOCIATED_NORMALIZED_TEST:\n";
  cout << "  LEGENDRE_ASSOCIATED_NORMALIZED evaluates \n";
  cout << "  normalized associated Legendre functions.\n";
  cout << "\n";
  cout << "      N       M    X     Exact F     PNM(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_associated_normalized ( n, m, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << m      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_function_q_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_FUNCTION_Q_TEST tests LEGENDRE_FUNCTION_Q.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 12

  double fx;
  double fx2[N_MAX+1];
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "LEGENDRE_FUNCTION_Q_TEST:\n";
  cout << "  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.\n";
  cout << "\n";
  cout << "     N      X        Exact F       Q(N)(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    legendre_function_q_values ( &n_data, &n, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    legendre_function_q ( n, x, fx2 );

    cout                       << "  "
         << setw(8)  << n      << "  "
         << setw(8)  << x      << "  "
         << setw(14) << fx     << "  "
         << setw(14) << fx2[n] << "\n";

  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_symbol_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SYMBOL_TEST tests LEGENDRE_SYMBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TEST 4

  int i;
  int l;
  int p;
  int ptest[N_TEST] = { 7, 11, 13, 17 };
  int q;

  cout << "\n";
  cout << "LEGENDRE_SYMBOL_TEST\n";
  cout << "  LEGENDRE_SYMBOL computes the Legendre\n";
  cout << "  symbol (Q/P) which records whether Q is \n";
  cout << "  a quadratic residue modulo the prime P.\n";

  for ( i = 0; i < N_TEST; i++ )
  {
    p = ptest[i];
    cout << "\n";
    cout << "  Legendre Symbols for P = " << p << "\n";
    cout << "\n";
    for ( q = 0; q <= p; q++ )
    {
      cout                                        << "  "
           << setw(8) << p                        << "  "
           << setw(8) << q                        << "  "
           << setw(8) << legendre_symbol ( q, p ) << "\n";
    }
  }

  return;
# undef N_TEST
}
//****************************************************************************80

void lerch_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH_TEST tests LERCH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int s;
  double z;

  cout << "\n";
  cout << "LERCH_TEST:\n";
  cout << "  LERCH evaluates the Lerch function.\n";
  cout << "\n";
  cout << "       Z       S       A         Lerch           Lerch\n";
  cout << "                             Tabulated        Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    lerch_values ( &n_data, &z, &s, &a, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lerch ( z, s, a );

    cout                    << "  "
         << setw(8)  << z   << "  "
         << setw(4)  << s   << "  "
         << setw(8)  << a   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void lgamma_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LGAMMA_TEST tests LGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "LGAMMA_TEST:\n";
  cout << "  LGAMMA is a C math library function which evaluates\n";
  cout << "  the logarithm of the Gamma function.\n";
  cout << "\n";
  cout << "     X       Exact F       LGAMMA(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_log_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = lgamma ( x );

    cout                    << "  "
         << setw(8)  << x   << "  "
         << setw(10) << fx  << "  "
         << setw(10) << fx2 << "\n";

  }

  return;
}
//****************************************************************************80

void lock_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOCK_TEST tests LOCK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N+1];
  int i;

  cout << "\n";
  cout << "LOCK_TEST\n";
  cout << "  LOCK counts the combinations on a button lock.\n";
  cout << "\n";
  cout << "  I,  LOCK(I)\n";
  cout << "\n";

  lock ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    cout                     << "  "
         << setw(4)  << i    << "  "
         << setw(10) << a[i] << "\n";
  }
 
  return;
# undef N
}
//****************************************************************************80

void meixner_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MEIXNER_TEST tests MEIXNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2009
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define TEST_NUM 3

  double beta;
  double beta_test[TEST_NUM] = { 0.5, 1.0, 2.0 };
  double c;
  double c_test[TEST_NUM] = { 0.125, 0.25, 0.5 };
  int i;
  int j;
  int n;
  int test;
  double v[N+1];
  double x;

  cout << "\n";
  cout << "MEIXNER_TEST:\n";
  cout << "  MEIXNER evaluates Meixner polynomials.\n";
  cout << "\n";
  cout << "       N      BETA         C         X        M(N,BETA,C,X)\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    n = N;
    beta = beta_test[test];
    c = c_test[test];

    for ( j = 0; j <= 5; j++ )
    {
      x = ( double ) ( j ) / 2.0;

      meixner ( n, beta, c, x, v );

      cout << "\n";

      for ( i = 0; i <= n; i++ )
      {
        cout << "  " << setw(8) << i
             << "  " << setw(8) << beta
             << "  " << setw(8) << c
             << "  " << setw(8) << x
             << "  " << setw(14) << v[i] << "\n";
      }
    }
  }

  return;
# undef N
# undef TEST_NUM
}
//****************************************************************************80

void mertens_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MERTENS_TEST tests MERTENS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "MERTENS_TEST\n";
  cout << "  MERTENS computes the Mertens function.\n";
  cout << "\n";
  cout << "      N   Exact   MERTENS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     mertens_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << mertens ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void moebius_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_TEST tests MOEBIUS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "MOEBIUS_TEST\n";
  cout << "  MOEBIUS computes the Moebius function.\n";
  cout << "\n";
  cout << "      N   Exact   MOEBIUS(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
     moebius_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                              << "  "
         << setw(8)  << n             << "  "
         << setw(10) << c             << "  "
         << setw(10) << moebius ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void motzkin_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MOTZKIN_TEST tests MOTZKIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int a[N+1];
  int i;

  cout << "\n";
  cout << "MOTZKIN_TEST\n";
  cout << "  MOTZKIN computes the Motzkin numbers A(0:N).\n";
  cout << "  A(N) counts the paths from (0,0) to (N,0).\n";
  cout << "\n";
  cout << "  I,  A(I)\n";
  cout << "\n";

  motzkin ( N, a );

  for ( i = 0; i <= N; i++ )
  {
    cout                     << "  "
         << setw(4)  << i    << "  "
         << setw(10) << a[i] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void normal_01_cdf_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_INVERSE_TEST tests NORMAL_01_CDF_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x;
  double x2;

  cout << "\n";
  cout << "NORMAL_01_CDF_INVERSE_TEST:\n";
  cout << "  NORMAL_01_CDF_INVERSE inverts the normal 01 CDF.\n";
  cout << "\n";
  cout << "    FX      X    NORMAL_01_CDF_INVERSE(FX)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = normal_01_cdf_inverse ( fx );

    cout                    << "  "
         << setw(8)  << fx   << "  "
         << setw(14) << x  << "  "
         << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void omega_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    OMEGA_TEST tests OMEGA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "OMEGA_TEST\n";
  cout << "  OMEGA computes the OMEGA function.\n";
  cout << "\n";
  cout << "          N   Exact   OMEGA(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    omega_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                            << "  "
         << setw(12) << n           << "  "
         << setw(10) << c           << "  "
         << setw(10) << omega ( n ) << "\n";

  }
 
  return;
}
//****************************************************************************80

void pentagon_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PENTAGON_NUM_TEST tests PENTAGON_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "PENTAGON_NUM_TEST\n";
  cout << "  PENTAGON_NUM computes the pentagonal numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout << "  " << setw(4) << n
         << "  " << setw(6) << pentagon_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void phi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PHI_TEST tests PHI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "PHI_TEST\n";
  cout << "  PHI computes the PHI function.\n";
  cout << "\n";
  cout << "  N   Exact   PHI(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    phi_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                          << "  "
         << setw(4)  << n         << "  "
         << setw(10) << c         << "  "
         << setw(10) << phi ( n ) << "\n";

  }
 
  return;
}
//****************************************************************************80

void plane_partition_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PLANE_PARTITION_NUM_TEST tests PLANE_PARTITION_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "PLANE_PARTITION_NUM_TEST\n";
  cout << "  PLANE_PARTITION_NUM computes the number of plane\n";
  cout << "  partitions of an integer.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout << "  " << setw(4) << n
         << "  " << setw(6) << plane_partition_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void poly_bernoulli_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_BERNOULLI_TEST tests POLY_BERNOULLI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int b;
  int k;
  int n;

  cout << "\n";
  cout << "POLY_BERNOULLI_TEST\n";
  cout << "  POLY_BERNOULLI computes the poly-Bernoulli numbers\n";
  cout << "  of negative index, B_n^(-k)\n";
  cout << "\n";
  cout << "   N   K    B_N^(-K)\n";
  cout << "\n";

  for ( k = 0; k <= 6; k++ )
  {
    cout << "\n";
    for ( n = 0; n <= 6; n++ )
    {
      b = poly_bernoulli ( n, k );

      cout << "  " << setw(2)  << n
           << "  " << setw(2)  << k
           << "  " << setw(12) << b << "\n";
    }
  }

  return;
}
//****************************************************************************80

void poly_coef_count_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POLY_COEF_COUNT_TEST tests POLY_COEF_COUNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int degree;
  int dim;
  int n;

  cout << "\n";
  cout << "POLY_COEF_COUNT_TEST\n";
  cout << "  POLY_COEF_COUNT counts the number of coefficients\n";
  cout << "  in a polynomial of degree DEGREE and dimension DIM.\n";
  cout << "\n";
  cout << " Dimension    Degree     Count\n";

  for ( dim = 1; dim <= 10; dim = dim + 3 )
  {
    cout << "\n";
    for ( degree = 0; degree <= 5; degree++ )
    {
      cout << "  " << setw(8) << dim
           << "  " << setw(8) << degree
           << "  " << setw(8) << poly_coef_count ( dim, degree ) << "\n";
    }
  }

  return;
}
//****************************************************************************80

void prime_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_TEST tests PRIME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int prime_max;

  cout << "\n";
  cout << "PRIME_TEST\n";
  cout << "  PRIME returns primes from a table.\n";

  n = -1;
  prime_max = prime ( n );
  cout << "\n";
  cout << "  Number of primes stored is " << prime_max << "\n";
  cout << "\n";
  cout << "     I    Prime(I)\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  "
         << setw(4) << i << "  "
         << setw(6) << prime ( i ) << "\n";
  }
  cout << "\n";
  for ( i = prime_max - 10; i <= prime_max; i++ )
  {
    cout << "  "
         << setw(4) << i << "  "
         << setw(6) << prime ( i ) << "\n";
  }
  
  return;
}
//****************************************************************************80

void pyramid_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_NUM_TEST tests PYRAMID_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "PYRAMID_NUM_TEST\n";
  cout << "  PYRAMID_NUM computes the pyramidal numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                 << "  "
         << setw(4) << n                 << "  "
         << setw(6) << pyramid_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void pyramid_square_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_SQUARE_NUM_TEST tests PYRAMID_SQUARE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "PYRAMID_SQUARE_NUM_TEST\n";
  cout << "  PYRAMID_SQUARE_NUM computes the pyramidal square numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout << "  "
         << setw(6) << n << "  "
         << setw(6) << pyramid_square_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void r8_agm_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AGM_TEST tests R8_AGM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;

  cout << "\n";
  cout << "R8_AGM_TEST\n";
  cout << "  R8_AGM computes the arithmetic geometric mean.\n";
  cout << "\n";
  cout << "           A           B         "
       << "   AGM                       AGM               Diff\n";
  cout << "                             "
       << "      (Tabulated)             R8_AGM(A,B)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    agm_values ( n_data, a, b, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_agm ( a, b );

    cout << "  " << setprecision(6)  << setw(10) << a   
         << "  " << setprecision(6)  << setw(10) << b   
         << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << fx2 
         << "  " << setprecision(6)  << setw(10) << fabs ( fx - fx2 ) << "\n";
  }

  return;
}
//****************************************************************************80

void r8_beta_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_TEST tests R8_BETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fxy;
  double fxy2;
  int n_data;
  double x;
  double y;

  cout << "\n";
  cout << "R8_BETA_TEST:\n";
  cout << "  R8_BETA evaluates the Beta function.\n";
  cout << "\n";
  cout << "     X      Y        Exact F       R8_BETA(X,Y)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_values ( &n_data, &x, &y, &fxy );

    if ( n_data == 0 )
    {
      break;
    }

    fxy2 = r8_beta ( x, y );

    cout                     << "  "
         << setw(10) << x    << "  "
         << setw(10) << y    << "  "
         << setw(10) << fxy  << "  "
         << setw(10) << fxy2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_choose_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE_TEST tests R8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cnk;
  int k;
  int n;

  cout << "\n";
  cout << "R8_CHOOSE_TEST\n";
  cout << "  R8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "   N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );

      cout                   << "  "
           << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void r8_erf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERF_TEST tests R8_ERF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "R8_ERF_TEST:\n";
  cout << "  R8_ERF evaluates the error function.\n";
  cout << "\n";
  cout << "     X      Exact F     R8_ERF(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( n_data, x, fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_erf ( x );

    cout                    << "  "
         << setw(8)  << x   << "  "
         << setw(14) << fx  << "  "
         << setw(14) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_erf_inverse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERF_INVERSE_TEST tests R8_ERF_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  int n_data;
  double x1;
  double x2;

  cout << "\n";
  cout << "R8_ERF_INVERSE_TEST\n";
  cout << "  R8_ERF_INVERSE inverts the error function.\n";
  cout << "\n";
  cout << "    FX           X1           X2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    erf_values ( n_data, x1, fx );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = r8_erf_inverse ( fx );

    cout                   << "  "
         << setw(8)  << fx << "  "
         << setw(14) << x1 << "  "
         << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_euler_constant_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EULER_CONSTANT_TEST tests R8_EULER_CONSTANT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double g;
  double g_approx;
  int i;
  int n;
  double n_r8;
  int test;

  g = r8_euler_constant ( );

  cout << "\n";
  cout << "R8_EULER_CONSTANT_TEST:\n";
  cout << "  R8_EULER_CONSTANT returns the Euler-Mascheroni constant\n";
  cout << "  sometimes denoted by 'gamma'.\n";
  cout << "\n";
  cout << "  gamma = limit ( N -> oo ) ( sum ( 1 <= I <= N ) 1 / I ) - log ( N )\n";
  cout << "\n";
  cout << "  Numerically, g = " << g << "\n";
  cout << "\n";
  cout << "         N      Partial Sum    |gamma - partial sum|\n";
  cout << "\n";

  n = 1;
  for ( test = 0; test <= 20; test++ )
  {
    n_r8 = ( double ) ( n );
    g_approx = - log ( n_r8 );
    for ( i = 1; i <= n; i++ )    
    {
      g_approx = g_approx + 1.0 / ( double ) ( i );
    }
    cout << "  " << setw(8) << n
         << "  " << setw(14) << g_approx
         << "  " << setw(14) <<  fabs ( g_approx - g ) << "\n";
    n = n * 2;
  }

  return;
}
//****************************************************************************80

void r8_factorial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_TEST tests R8_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n_data;
  int n;

  cout << "\n";
  cout << "R8_FACTORIAL_TEST:\n";
  cout << "  R8_FACTORIAL evaluates the factorial function.\n";
  cout << "\n";
  cout << "     N       Exact F       R8_FACTORIAL(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    cout                                  << "  "
         << setw(4)  << n                 << "  "
         << setw(14) << fn                << "  "
         << setw(14) << r8_factorial ( n ) << "\n";
  }

  return;
}
//****************************************************************************80

void r8_factorial_log_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_LOG_TEST tests R8_FACTORIAL_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fn;
  int n_data;
  int n;

  cout << "\n";
  cout << "R8_FACTORIAL_LOG_TEST:\n";
  cout << "  R8_FACTORIAL_LOG evaluates the logarithm of the\n";
  cout << "  factorial function.\n";
  cout << "\n";
  cout << "     N	   Exact F	 R8_FACTORIAL_LOG(N)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial_log_values ( &n_data, &n, &fn );

    if ( n_data == 0 )
    {
      break;
    }

    cout                                      << "  "
         << setw(5)  << n                     << "  "
         << setw(14) << fn                    << "  "
         << setw(14) << r8_factorial_log ( n ) << "\n";

  }

  return;
}
//****************************************************************************80

void r8_hyper_2f1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HYPER_2F1_TEST tests R8_HYPER_2F1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << " R8_HYPER_2F1_TEST:\n";
  cout << "   R8_HYPER_2F1 evaluates the hypergeometric function 2F1.\n";
  cout << "\n";
  cout << "      A       B       C       X      ";
  cout << " 2F1                       2F1                     DIFF\n";
  cout << "                                     ";
  cout << "(tabulated)               (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hyper_2f1_values ( &n_data, &a, &b, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_hyper_2f1 ( a, b, c, x );

    cout << "  " << setw(6)  << setprecision(2)  << a
         << "  " << setw(6)  << setprecision(2)  << b  
         << "  " << setw(6)  << setprecision(2)  << c  
         << "  " << setw(6)  << setprecision(2)  << x  
         << "  " << setw(24) << setprecision(16) << fx
         << "  " << setw(24) << setprecision(16) << fx2
         << "  " << setw(10) << setprecision(4)  << fabs ( fx - fx2 ) << "\n";
  }
  return;
}
//****************************************************************************80

void r8_psi_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PSI_TEST tests R8_PSI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "R8_PSI_TEST:\n";
  cout << "  R8_PSI evaluates the Psi function.\n";
  cout << "\n";
  cout << "         X                  Psi(X)           " 
       << "         Psi(X)          DIFF\n";
  cout << "                         (Tabulated)         " 
       << "       (R8_PSI)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = r8_psi ( x );

    cout << "  " << setprecision(2) << setw(8)  << x   
         << "  " << setprecision(16) << setw(24) << fx  
         << "  " << setprecision(16) << setw(24) << fx2 
         << "  " << setprecision(4) << setw(10) << fabs ( fx - fx2 ) << "\n";

  }

  return;
}
//****************************************************************************80

void r8poly_degree_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_DEGREE_TEST tests R8POLY_DEGREE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c1[4] = { 1.0, 2.0, 3.0, 4.0 }; 
  double c2[4] = { 1.0, 2.0, 3.0, 0.0 };
  double c3[4] = { 1.0, 2.0, 0.0, 4.0 };
  double c4[4] = { 1.0, 0.0, 0.0, 0.0 };
  double c5[4] = { 0.0, 0.0, 0.0, 0.0 };
  int d;
  int m;
 
  cout << "\n";
  cout << "R8POLY_DEGREE_TEST\n";
  cout << "  R8POLY_DEGREE determines the degree of an R8POLY.\n";

  m = 3;

  r8poly_print ( m, c1, "  The R8POLY:" );
  d = r8poly_degree ( m, c1 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c2, "  The R8POLY:" );
  d = r8poly_degree ( m, c2 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c3, "  The R8POLY:" );
  d = r8poly_degree ( m, c3 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c4, "  The R8POLY:" );
  d = r8poly_degree ( m, c4 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  r8poly_print ( m, c5, "  The R8POLY:" );
  d = r8poly_degree ( m, c5 );
  cout << "  Dimensioned degree = " << m << ",  Actual degree = " << d << "\n";

  return;
}
//****************************************************************************80

void r8poly_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT_TEST tests R8POLY_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  cout << "\n";
  cout << "R8POLY_PRINT_TEST\n";
  cout << "  R8POLY_PRINT prints an R8POLY.\n";

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
//****************************************************************************80

void r8poly_value_horner_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  cout << "\n";
  cout << "R8POLY_VALUE_HORNER_TEST\n";
  cout << "  R8POLY_VALUE_HORNER evaluates a polynomial at\n";
  cout << "  one point, using Horner's method.\n";

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  cout << "\n";
  cout << "   I    X    P(X)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    cout << "  " << setw(2) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(14) << p << "\n";
  }

  delete [] x;

  return;
}
//****************************************************************************80

void sigma_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA_TEST tests SIGMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "SIGMA_TEST\n";
  cout << "  SIGMA computes the SIGMA function.\n";
  cout << "\n";
  cout << "  N   Exact   SIGMA(N)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    sigma_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                            << "  "
         << setw(4)  << n           << "  "
         << setw(10) << c           << "  "
         << setw(10) << sigma ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void simplex_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_NUM_TEST tests SIMPLEX_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int value;

  cout << "\n";
  cout << "SIMPLEX_NUM_TEST\n";
  cout << "  SIMPLEX_NUM computes the N-th simplex number\n";
  cout << "  in M dimensions.\n";
  cout << "\n";
  cout << "      M: 0     1     2     3     4     5\n";
  cout << "   N\n";
 
  for ( n = 0; n <= 10; n++ )
  {
    cout << "  " << setw(2) << n;
    for ( m = 0; m <= 5; m++ )
    {
      value = simplex_num ( m, n );
      cout << "  " << setw(4) << value;
    }
    cout << "\n";
  } 
  return;
}
//****************************************************************************80

void sin_power_int_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT_TEST tests SIN_POWER_INT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n;
  int n_data;

  cout << "\n";
  cout << "SIN_POWER_INT_TEST:\n";
  cout << "  SIN_POWER_INT computes the integral of the N-th power\n";
  cout << "  of the sine function.\n";
  cout << "\n";
  cout << "         A         B       N        Exact    Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
   sin_power_int_values ( &n_data, &a, &b, &n, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = sin_power_int ( a, b, n );

    cout                    << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(6)  << n   << "  "
         << setw(12) << fx  << "  "
         << setw(12) << fx2 << "\n";
  }
  return;
}
//****************************************************************************80

void slice_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SLICE_TEST tests SLICE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2011
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 5
# define SLICE_MAX 8

  int dim_max = DIM_MAX;
  int dim_num;
  int p[DIM_MAX*SLICE_MAX];
  int piece_num;
  int slice_max = SLICE_MAX;
  int slice_num;

  cout << "\n";
  cout << "SLICE_TEST:\n";
  cout << "  SLICE determines the maximum number of pieces created\n";
  cout << "  by SLICE_NUM slices in a DIM_NUM space.\n";

  for ( dim_num = 1; dim_num <= dim_max; dim_num++ )
  {
    for ( slice_num = 1; slice_num <= slice_max; slice_num++ )
    {
      piece_num = slice ( dim_num, slice_num );
      p[dim_num-1+(slice_num-1)*dim_max] = piece_num;
    }
  }

  i4mat_print ( dim_max, slice_max, p, "  Slice Array:" );

  return;
# undef DIM_MAX
# undef SLICE_MAX
}
//****************************************************************************80

void spherical_harmonic_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERICAL_HARMONIC_TEST tests SPHERICAL_HARMONIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 20

  double c[N_MAX+1];
  int l;
  int m;
  int n_data;
  double phi;
  double s[N_MAX+1];
  double theta;
  double yi;
  double yi2;
  double yr;
  double yr2;

  cout << "\n";
  cout << "SPHERICAL_HARMONIC_TEST:\n";
  cout << "  SPHERICAL_HARMONIC evaluates spherical harmonic functions.\n";
  cout << "\n";
  cout << 
    "         N         M    THETA      PHI            YR            YI\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    spherical_harmonic_values ( &n_data, &l, &m, &theta, &phi, &yr, &yi );

    if ( n_data == 0 )
    {
      break;
    }

    spherical_harmonic ( l, m, theta, phi, c, s );

    yr2 = c[l];
    yi2 = s[l];

    cout                      << "  "
         << setw(8)  << l     << "  "
         << setw(8)  << m     << "  "
         << setw(8)  << theta << "  "
         << setw(8)  << phi   << "  "
         << setw(14) << yr    << "  "
         << setw(14) << yi    << "\n";

    cout                      << "  "
         << "        "        << "  "
         << "        "        << "  "
         << "        "        << "  "
         << "        "        << "  "
         << setw(14) << yr2   << "  "
         << setw(14) << yi2   << "\n";
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void stirling1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING1_TEST tests STIRLING1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m = 8;
  int n = 8;
  int *s1;

  cout << "\n";
  cout << "STIRLING1_TEST\n";
  cout << "  STIRLING1: Stirling numbers of first kind.\n";
  cout << "  Get rows 1 through " << m << "\n";
  cout << "\n";
 
  s1 = stirling1 ( m, n );
 
  for ( i = 0; i < m; i++ )
  {
    cout << setw(6) << i+1 << "  ";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(6) << s1[i+j*m] << "  ";
    }
    cout << "\n";
  }

  delete [] s1;
 
  return;
}
//****************************************************************************80

void stirling2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING2_TEST tests STIRLING2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 June 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define N 8

  int i;
  int j;
  int *s2;

  cout << "\n";
  cout << "STIRLING2_TEST\n";
  cout << "  STIRLING2: Stirling numbers of second kind.\n";
  cout << "  Get rows 1 through " << M << "\n";
  cout << "\n";
 
  s2 = stirling2 ( M, N );
 
  for ( i = 0; i < M; i++ )
  {
    cout << setw(6) << i+1 << "  ";
    for ( j = 0; j < N; j++ )
    {
      cout << setw(6) << s2[i+j*M] << "  ";
    }
    cout << "\n";
  }
 
  delete [] s2;

  return;
# undef M
# undef N
}
//****************************************************************************80

void tau_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TAU_TEST tests TAU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int n;
  int n_data;

  cout << "\n";
  cout << "TAU_TEST\n";
  cout << "  TAU computes the Tau function.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";
 
  n_data = 0;

  for ( ; ; )
  {
    tau_values ( &n_data, &n, &c );

    if ( n_data == 0 )
    {
      break;
    }

    cout                          << "  "
         << setw(4)  << n         << "  "
         << setw(10) << c         << "  "
         << setw(10) << tau ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void tetrahedron_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_NUM_TEST tests TETRAHEDRON_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TETRAHEDRON_NUM_TEST\n";
  cout << "  TETRAHEDRON_NUM computes the tetrahedron numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                     << "  "
         << setw(4) << n                     << "  "
         << setw(6) << tetrahedron_num ( n ) << "\n";
  }
 
  return;
}
//****************************************************************************80

void triangle_num_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NUM_TEST tests TRIANGLE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;

  cout << "\n";
  cout << "TRIANGLE_NUM_TEST\n";
  cout << "  TRIANGLE_NUM computes the triangular numbers.\n";
  cout << "\n";
 
  for ( n = 1; n <= 10; n++ )
  {
    cout                                  << "  "
         << setw(4) << n                  << "  "
         << setw(6) << triangle_num ( n ) << "\n";;
  }
 
  return;
}
//****************************************************************************80

void triangle_to_i4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_TO_I4_TEST tests TRIANGLE_TO_I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TRIANGLE_TO_I4_TEST\n";
  cout << "  TRIANGLE_TO_I4 converts a triangular index to a\n";
  cout << "  linear one.\n";
  cout << "\n";
  cout << "     I     J ==>   K\n";
  cout << "\n";

  for ( i = 0; i <= 4; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      k = triangle_to_i4 ( i, j );

      cout << "  " << setw(4) << i
           << "  " << setw(4) << j
           << "    " << setw(4) << k << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void trinomial_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRINOMIAL_TEST tests TRINOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int t;

  cout << "\n";
  cout << "TRINOMIAL_TEST\n";
  cout << "  TRINOMIAL evaluates the trinomial coefficient:\n";
  cout << "\n";
  cout << "  T(I,J,K) = (I+J+K)! / I! / J! / K!\n";
  cout << "\n";
  cout << "     I     J     K    T(I,J,K)\n";
  cout << "\n";
 
  for ( k = 0; k <= 4; k++ )
  {
    for ( j = 0; j <= 4; j++ )
    {
      for ( i = 0; i <= 4; i++ )
      {
        t = trinomial ( i, j, k );
        cout << "  " << setw(4) << i
             << "  " << setw(4) << j
             << "  " << setw(4) << k
             << "  " << setw(8) << t << "\n";
      }
    }
  }
 
  return;
}
//****************************************************************************80

void v_hofstadter_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    V_HOFSTADTER_TEST tests V_HOFSTADTER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int v;

  cout << "\n";
  cout << "V_HOFSTADTER_TEST\n";
  cout << "  V_HOFSTADTER evaluates Hofstadter's recursive\n";
  cout << "  V function.\n";
  cout << "\n";
  cout << "     N   V(N)\n";
  cout << "\n";

  for ( i = 0; i <= 30; i++ )
  {
    cout                                  << "  "
         << setw(6) << i                  << "  "
         << setw(6) << v_hofstadter ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void vibonacci_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    VIBONACCI_TEST tests VIBONACCI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20
  int i;
  int j;
  int seed;
  int v1[N];
  int v2[N];
  int v3[N];

  cout << "\n";
  cout << "VIBONACCI_TEST\n";
  cout << "  VIBONACCI computes a Vibonacci sequence.\n";
  cout << "\n";
  cout << "  We compute the series 3 times.\n";
  cout << "\n";
  cout << "     I      V1      V2      V3\n";
  cout << "\n";

  seed = 123456789;

  vibonacci ( N, seed, v1 );
  vibonacci ( N, seed, v2 );
  vibonacci ( N, seed, v3 );

  for ( i = 0; i < N; i++ )
  {
    cout                     << "  "
         << setw(6) << i     << "  "
         << setw(6) << v1[i] << "  "
         << setw(6) << v2[i] << "  "
         << setw(6) << v3[i] << "\n";
  } 

  return;
# undef N
}
//****************************************************************************80

void zeckendorf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ZECKENDORF_TEST tests ZECKENDORF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define M_MAX 20

  int i;
  int i_list[M_MAX];
  int j;
  int f_list[M_MAX];
  int f_sum;
  int m;
  int n;

  cout << "\n";
  cout << "ZECKENDORF_TEST\n";
  cout << "  ZECKENDORF computes the Zeckendorf decomposition of\n";
  cout << "  an integer N into nonconsecutive Fibonacci numbers.\n";
  cout << "\n";
  cout << "   N Sum M Parts\n";
  cout << "\n";

  for ( n = 1; n <= 100; n++ )
  {
    zeckendorf ( n, M_MAX, &m, i_list, f_list );

    cout << setw(4) << n << "  ";
    for ( j = 0; j < m; j++ )
    {
      cout << setw(4) << f_list[j] << "  ";
    }
    cout << "\n";

  }

  return;
# undef M_MAX
}
//****************************************************************************80

void zernike_poly_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ZERNIKE_POLY_TEST tests ZERNIKE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int i;
  int m;
  int n;
  double rho;
  double z1;
  double z2;

  cout << "\n";
  cout << "ZERNIKE_POLY_TEST\n";
  cout << "  ZERNIKE_POLY evaluates a Zernike polynomial directly.\n";
  cout << "\n";
  cout << "  Table of polynomial coefficients:\n";
  cout << "\n";
  cout << "   N   M\n";
  cout << "\n";

  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      cout << "  " << setw(2) << n
           << "  " << setw(2) << m;
      for ( i = 0; i <= n; i++ )
      {
        cout << "  " << setw(7) << c[i];
      }
      cout << "\n";
      delete [] c;
    }
  }

  rho = 0.987654321;

  cout << "\n";
  cout << "  Z1: Compute polynomial coefficients,\n";
  cout << "  then evaluate by Horner's method;\n";
  cout << "  Z2: Evaluate directly by recursion.\n";
  cout << "\n";
  cout << "   N   M       Z1              Z2\n";
  cout << "\n";

  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( m = 0; m <= n; m++ )
    {
      c = zernike_poly_coef ( m, n );
      z1 = r8poly_value_horner ( n, c, rho );

      z2 = zernike_poly ( m, n, rho );
      cout << "  " << setw(2)  << n
           << "  " << setw(2)  << m
           << "  " << setw(16) << z1
           << "  " << setw(16) << z2 << "\n";

      delete [] c;
    }
  }

  return;
}
//****************************************************************************80

void zernike_poly_coef_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ZERNIKE_POLY_COEF_TEST tests ZERNIKE_POLY_COEF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  int m;
  int n;

  cout << "\n";
  cout << "ZERNIKE_POLY_COEF_TEST\n";
  cout << "  ZERNIKE_POLY_COEF determines the Zernike\n";
  cout << "  polynomial coefficients.\n";

  n = 5;

  for ( m = 0; m <= n; m++ )
  {
    c = zernike_poly_coef ( m, n );
    r8poly_print ( n, c, "  Zernike polynomial" );
    delete [] c;
  }

  return;
}
//****************************************************************************80

void zeta_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA_TEST tests ZETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2005
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;
  double n_real;
  double z1;
  double z2;

  cout << "\n";
  cout << "ZETA_TEST\n";
  cout << "  ZETA computes the Zeta function.\n";
  cout << "\n";
  cout << "       N            exact Zeta         computed Zeta\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    zeta_values ( &n_data, &n, &z1 );

    if ( n_data == 0 )
    {
      break;
    }

    n_real = ( double ) n;

    z2 = zeta ( n_real );

    cout                   << "  "
         << setw(6)  << n  << "  "
         << setw(20) << z1 << "  "
         << setw(20) << z2 << "\n";

  }

  return;
}
