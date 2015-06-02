# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sparse_grid_hw.hpp"

int main ( );
void ccl_test ( );
void ccl_sparse_test ( );
void ccs_test ( );
void ccs_sparse_test ( );
void cce_test ( );
void cce_sparse_test ( );
void get_seq_test ( );
void gqn_test ( );
void gqn_sparse_test ( );
void gqn2_sparse_test ( );
void gqu_test ( );
void gqu_sparse_test ( );
void kpn_test ( );
void kpn_sparse_test ( );
void kpu_test ( );
void kpu_sparse_test ( );
void nwspgr_size_test ( );
void nwspgr_time_test ( );
void nwspgr_test ( );
void order_report ( );
void symmetric_sparse_size_test ( );
void tensor_product_test ( );
void tensor_product_cell_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_HW_PRB.
//
//  Discussion:
//
//    SPARSE_GRID_HW_PRB tests the SPARSE_GRID_HW library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_HW_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPARSE_GRID_HW library.\n";

  ccl_test ( );
  ccl_sparse_test ( );

  ccs_test ( );  
  ccs_sparse_test ( );

  cce_test ( );
  cce_sparse_test ( );

  get_seq_test ( );

  gqn_test ( );
  gqn_sparse_test ( );
  gqn2_sparse_test ( );

  gqu_test ( );
  gqu_sparse_test ( );

  kpn_test ( );
  kpn_sparse_test ( );

  kpu_test ( );
  kpu_sparse_test ( );

  nwspgr_size_test ( );
  nwspgr_time_test ( );
  nwspgr_test ( );

  order_report ( );

  symmetric_sparse_size_test ( );

  tensor_product_test ( );
  tensor_product_cell_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPARSE_GRID_HW_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void ccl_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCL_TEST uses CCL_ORDER + CC for 1D quadrature over [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "CCL_TEST:\n";
  cout << "  CCL_ORDER + CC\n";
  cout << "  Clenshaw Curtis Linear (CCL) quadrature over [0,1]:\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = ccl_order ( l );

    x = new double[n];
    w = new double[n];

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void ccl_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCL_SPARSE_TEST uses CCL_ORDER + CC for a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  cout << "\n";
  cout << "CCL_SPARSE_TEST:\n";
  cout << "  CCL_ORDER + CC\n";
  cout << "  Sparse Clenshaw Curtis Linear quadrature over [0,1].\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( ccl_order, d, k );

    x = new double[d*n];
    w = new double[n];

    nwspgr ( cc, ccl_order, d, k, n, n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }

  return;
}
//****************************************************************************80

void ccs_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCS_TEST uses CCS_ORDER + CC for 1D quadrature over [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "CCS_TEST:\n";
  cout << "  CCS_ORDER + CC\n";
  cout << "  Clenshaw Curtis Slow quadrature over [0,1]:\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = ccs_order ( l );

    x = new double[n];
    w = new double[n];

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void ccs_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCS_SPARSE_TEST uses CCS_ORDER + CC for a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  cout << "\n";
  cout << "CCS_SPARSE_TEST:\n";
  cout << "  CCS_ORDER + CC\n";
  cout << "  Sparse Clenshaw Curtis Slow quadrature over [0,1].\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( ccs_order, d, k );

    x = new double[d*n];
    w = new double[n];

    nwspgr ( cc, ccs_order, d, k, n, n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }

  return;
}
//****************************************************************************80

void cce_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCE_TEST uses CCE_ORDER + CC for 1D quadrature over [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "CCE_TEST:\n";
  cout << "  CCE_ORDER + CC\n";
  cout << "  Clenshaw Curtis Exponential quadrature over [0,1]:\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = cce_order ( l );

    x = new double[n];
    w = new double[n];

    cc ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void cce_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    CCE_SPARSE_TEST uses CCE_ORDER + CC for a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  cout << "\n";
  cout << "CCE_SPARSE_TEST:\n";
  cout << "  CCE_ORDER + CC\n";
  cout << "  Sparse Clenshaw Curtis Exponential quadrature over [0,1].\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( cce_order, d, k );

    x = new double[d*n];
    w = new double[n];

    nwspgr ( cc, cce_order, d, k, n, n2, x, w );

    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }

  return;
}
//****************************************************************************80

void get_seq_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEQ_TEST tests GET_SEQ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;

  int *fs;
  int norm;
  int seq_num;

  cout << "\n";
  cout << "GET_SEQ_TEST\n";
  cout << "  GET_SEQ returns all D-dimensional vectors that sum to NORM.\n";

  d = 3;
  norm = 6;

  cout << "\n";
  cout << "  D = " << d << "\n";
  cout << "  NORM = " << norm << "\n";

  seq_num = num_seq ( norm - d, d );

  fs = get_seq ( d, norm, seq_num );

  i4mat_print ( seq_num, d, fs, "  The compositions" );

  delete [] fs;

  return;
}
//****************************************************************************80

void gqn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GQN_TEST uses the GQN function for 1D quadrature over (-oo,+oo).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "\n";
  cout << "GQN_TEST:\n";
  cout << "  Gauss-Hermite quadrature over (-oo,+oo):\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fn_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = l;
    x = new double[n];
    w = new double[n];

    gqn ( n, x, w );

    fx = fn_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;

  }
  return;
}
//****************************************************************************80

void gqn_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GQN_SPARSE_TEST uses the GQN function to build a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fn_integral ( d );

  cout << "\n";
  cout << "GQN_SPARSE_TEST:\n";
  cout << "  GQN sparse grid:\n";
  cout << "  Sparse Gaussian quadrature with Hermite weight over (-oo,+oo).\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( gqn_order, d, k );
    x = new double[d*n];
    w = new double[n];
    nwspgr ( gqn, gqn_order, d, k, n, n2, x, w );
    fx = fn_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );

    error_sg = r8_abs ( ( estimate - trueval ) / trueval );

    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_normal_01_new ( d, n2, seed );
      fx = fn_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }
  return;
}
//****************************************************************************80

void gqn2_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GQN2_SPARSE_TEST uses the GQN and GQN2_ORDER functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  int j;
  int k;
  int maxk;
  int n;
  int n2;
  double *w;
  double *x;

  d = 2;
  maxk = 4;

  cout << "\n";
  cout << "GQN2_SPARSE_TEST:\n";
  cout << "  GQN sparse grid:\n";
  cout << "  Gauss-Hermite sparse grids over (-oo,+oo).\n";
  cout << "  Use GQN2_ORDER, the growth rule N = 2 * L - 1.\n";

  for ( k = 2; k <= maxk; k++ )
  {
    cout << "\n";
    cout << "     J      W                X               Y\n";
    cout << "\n";

    n = nwspgr_size ( gqn2_order, d, k );

    x = new double[d*n];
    w = new double[n];
    nwspgr ( gqn, gqn2_order, d, k, n, n2, x, w );

    for ( j = 0; j < n2; j++ )
    {
      cout << "  " << setw(4) << j
           << "  " << setw(14) << w[j]
           << "  " << setw(14) << x[0+j*d]
           << "  " << setw(14) << x[1+j*d] << "\n";
    }

    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void gqu_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GQU_TEST uses the GQU function for 1D quadrature over [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "GQU_TEST:\n";
  cout << "  Gauss-Legendre quadrature over [0,1]:\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = l;
    x = new double[n];
    w = new double[n];

    gqu ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }

  return;
}
//****************************************************************************80

void gqu_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GQU_SPARSE_TEST uses the GQU function to build a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  cout << "\n";
  cout << "GQU_SPARSE_TEST:\n";
  cout << "  GQU sparse grid:\n";
  cout << "  Sparse Gaussian unweighted quadrature over [0,1].\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( gqu_order, d, k );
    x = new double[d*n];
    w = new double[n];
    nwspgr ( gqu, gqu_order, d, k, n, n2, x, w );
    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }

  return;
}
//****************************************************************************80

void kpn_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KPN_TEST uses the KPN function for 1D quadrature over (-oo,+oo).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "KPN_TEST:\n";
  cout << "  Kronrod-Patterson-Hermite quadrature over (-oo,+oo):\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fn_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = kpn_order ( l );
    x = new double[n];
    w = new double[n];
    kpn ( n, x, w );

    fx = fn_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void kpn_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KPN_SPARSE_TEST uses the KPN function to build a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fn_integral ( d );

  cout << "\n";
  cout << "KPN_SPARSE_TEST:\n";
  cout << "  KPN sparse grid:\n";
  cout << "  Sparse Kronrod-Patterson quadrature with Hermite weight over (-oo,+oo).\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( kpn_order, d, k );
    x = new double[d*n];
    w = new double[n];
    nwspgr ( kpn, kpn_order, d, k, n, n2, x, w );
    fx = fn_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );

    error_sg = r8_abs ( ( estimate - trueval ) / trueval );

    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_normal_01_new ( d, n2, seed );
      fx = fn_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }
  return;
}
//****************************************************************************80

void kpu_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KPU_TEST uses the KPU function for 1D quadrature over [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  double e;
  double exact;
  double *fx;
  int i;
  int l;
  int n;
  double q;
  double *w;
  double *x;

  cout << "\n";
  cout << "KPU_TEST:\n";
  cout << "  Kronrod-Patterson quadrature over [0,1]:\n";
  cout << "\n";
  cout << "   Level   Nodes    Estimate  Error\n";
  cout << "\n";

  d = 1;
  exact = fu_integral ( d );

  for ( l = 1; l <= 5; l++ )
  {
    n = kpu_order ( l );
    x = new double[n];
    w = new double[n];
    kpu ( n, x, w );

    fx = fu_value ( d, n, x );

    q = r8vec_dot_product ( n, w, fx );

    e = r8_abs ( q - exact ) / exact;

    cout << "  " << setw(2) << l
         << "  " << setw(6) << n
         << "  " << setw(14) << q
         << "  " << setw(14) << e << "\n";

    delete [] fx;
    delete [] w;
    delete [] x;
  }
  return;
}
//****************************************************************************80

void kpu_sparse_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    KPU_SPARSE_TEST uses the KPU function to build a sparse grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int d;
  double error_mc;
  double error_sg;
  double estimate;
  double *fx;
  int k;
  int maxk;
  int n;
  int n2;
  int r;
  double *s;
  int s_num;
  int seed;
  double trueval;
  double *w;
  double *x;

  d = 10;
  maxk = 7;

  trueval = fu_integral ( d );

  cout << "\n";
  cout << "KPU_SPARSE_TEST:\n";
  cout << "  KPU sparse grid:\n";
  cout << "  Sparse Kronrod-Patterson unweighted quadrature over [0,1].\n";
  cout << "\n";
  cout << "   D  Level   Nodes    SG error    MC error\n";
  cout << "\n";

  for ( k = 2; k <= maxk; k++ )
  {
//
//  Compute sparse grid estimate.
//
    n = nwspgr_size ( kpu_order, d, k );
    x = new double[d*n];
    w = new double[n];
    nwspgr ( kpu, kpu_order, d, k, n, n2, x, w );
    fx = fu_value ( d, n2, x );
    estimate = r8vec_dot_product ( n2, w, fx );
    error_sg = r8_abs ( ( estimate - trueval ) / trueval );
    delete [] fx;
    delete [] w;
    delete [] x;
//
//  Compute 1000 Monte Carlo estimates with same number of points, and average.
//
    s_num = 1000;
    s = new double[s_num];
    seed = 123456789;
    for ( r = 0; r < 1000; r++ )
    {
      x = r8mat_uniform_01_new ( d, n2, seed );
      fx = fu_value ( d, n2, x );
      s[r] = r8vec_sum ( n2, fx ) / ( double ) ( n2 );
      delete [] fx;
      delete [] x;
    }

    error_mc = 0.0;
    for ( r = 0; r < s_num; r++ )
    {
      error_mc = error_mc + pow ( s[r] - trueval, 2 );
    }
    error_mc = sqrt ( error_mc / ( double ) ( s_num ) ) / trueval;

    cout << "  " << setw(2) << d
         << "  " << setw(5) << k
         << "  " << setw(6) << n2
         << "  " << setw(10) << error_sg
         << "  " << setw(10) << error_mc << "\n";

    delete [] s;
  }

  return;
}
//****************************************************************************80

void nwspgr_size_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NWSPGR_SIZE_TEST tests NWSPGR_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 January 2013
//
//  Author:
//
//    John Burkardt.
//
{
  int dim;
  int k;
  int r_size;

  cout << "\n";
  cout << "NWSPGR_SIZE_TEST:\n";
  cout << "  NWSPGR_SIZE returns the size of a sparse grid, based on either:\n";
  cout << "  one of the built-in 1D rules, or a family of 1D rules\n";
  cout << "  supplied by the user.\n";

  dim = 2;
  k = 3;
  cout << "\n";
  cout << "  Kronrod-Patterson, [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( kpu_order, dim, k );
  cout << "  Full          " << r_size << "\n";

  dim = 2;
  k = 3;
  cout << "\n";
  cout << "  Kronrod-Patterson, (-oo,+oo), Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( kpn_order, dim, k );
  cout << "  Full          " << r_size << "\n";

  dim = 2;
  k = 3;
  cout << "\n";
  cout << "  Gauss-Legendre, [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( gqu_order, dim, k );
  cout << "  Full          " << r_size << "\n";

  dim = 2;
  k = 3;
  cout << "\n";
  cout << "  Gauss Hermite, (-oo,+oo), [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( gqn_order, dim, k );
  cout << "  Full          " << r_size << "\n";

  dim = 2;
  k = 3;
  cout << "\n";
  cout << "  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim " << dim << ", Level " << k << ", Unsymmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( cce_order, dim, k );
  cout << "  Full          " << r_size << "\n";
//
//  Do a table.
//
  cout << "\n";
  cout << "  Dimension / Level table for Clenshaw Curtis Exponential\n";
  cout << "\n";
  cout << " Dim: ";
  for ( dim = 1; dim <= 10; dim++ )
  {
    cout << "  " << setw(6) << dim;
  }
  cout << "\n";
  cout << "Level\n";
  for ( k = 1; k <= 5; k++ )
  {
    cout << "  " << setw(2) << k << "  ";
    for ( dim = 1; dim <= 10; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      cout << "  " << setw(6) << r_size;
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void nwspgr_time_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NWSPGR_TIME_TEST times NWSPGR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2013
//
//  Author:
//
//    John Burkardt.
//
{
  int dim;
  int k;
  double *nodes;
  int r_size;
  int s_size;
  double t1;
  double t2;
  double *weights;

  cout << "\n";
  cout << "  This function measures the time in seconds required by NWSPGR\n";
  cout << "  to compute a sparse grid, based on either:\n";
  cout << "  one of the built-in 1D rules, or a family of 1D rules\n";
  cout << "  supplied by the user.\n";

  dim = 20;
  k = 5;
  cout << "\n";
  cout << "  Kronrod-Patterson, [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( kpu_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  t1 = cpu_time ( );
  nwspgr ( kpu, kpu_order, dim, k, r_size, s_size, nodes, weights );
  t2 = cpu_time ( );
  delete [] nodes;
  delete [] weights;
  cout << "  Full          " << t2 - t1 << "\n";

  dim = 20;
  k = 5;
  cout << "\n";
  cout << "  Kronrod-Patterson, (-oo,+oo), Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( kpn_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  t1 = cpu_time ( );
  nwspgr ( kpn, kpn_order, dim, k, r_size, s_size, nodes, weights );
  t2 = cpu_time ( );
  delete [] nodes;
  delete [] weights;
  cout << "  Full          " << t2 - t1 << "\n";

  dim = 20;
  k = 5;
  cout << "\n";
  cout << "  Gauss-Legendre, [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( gqu_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  t1 = cpu_time ( );
  nwspgr ( gqu, gqu_order, dim, k, r_size, s_size, nodes, weights );
  t2 = cpu_time ( );
  delete [] nodes;
  delete [] weights;
  cout << "  Full          " << t2 - t1 << "\n";

  dim = 20;
  k = 5;
  cout << "\n";
  cout << "  Gauss Hermite, (-oo,+oo), [0,1], Dim " << dim << ", Level " << k << ", Symmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( gqn_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  t1 = cpu_time ( );
  nwspgr ( gqn, gqn_order, dim, k, r_size, s_size, nodes, weights );
  t2 = cpu_time ( );
  delete [] nodes;
  delete [] weights;
  cout << "  Full          " << t2 - t1 << "\n";

  dim = 20;
  k = 5;
  cout << "\n";
  cout << "  Clenshaw Curtis Exponential, [-1,+1], [0,1], Dim " << dim << ", Level " << k << ", Unsymmetric\n";
  cout << "\n";
  r_size = nwspgr_size ( cce_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  t1 = cpu_time ( );
  nwspgr ( cc, cce_order, dim, k, r_size, s_size, nodes, weights );
  t2 = cpu_time ( );
  delete [] nodes;
  delete [] weights;
  cout << "  Full          " << t2 - t1 << "\n";
/*
  Do a table.
*/
  cout << "\n";
  cout << "  Dimension / Level table for Clenshaw Curtis Exponential\n";
  cout << "\n";
  cout << " Dim: ";
  for ( dim = 1; dim <= 10; dim++ )
  {
    cout << "  " << setw(6) << dim;
  }
  cout << "\n";
  cout << "Level\n";
  for ( k = 1; k <= 5; k++ )
  {
    cout << "  " << setw(2) << k << "  ";
    for ( dim = 1; dim <= 10; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      nodes = new double[dim*r_size];
      weights = new double[r_size];
      t1 = cpu_time ( );
      nwspgr ( cc, cce_order, dim, k, r_size, s_size, nodes, weights );
      t2 = cpu_time ( );
      delete [] nodes;
      delete [] weights;
      cout << "  " << setw(10) << t2 - t1;
    }
    cout << "\n";
  }
  cout << "\n";
  cout << " Dim: ";
  for ( dim = 11; dim <= 20; dim++ )
  {
    cout << "  " << setw(6) << dim;
  }
  cout << "\n";
  cout << "Level\n";
  for ( k = 1; k <= 5; k++ )
  {
    cout << "  " << setw(2) << k << "  ";
    for ( dim = 11; dim <= 20; dim++ )
    {
      r_size = nwspgr_size ( cce_order, dim, k );
      nodes = new double[dim*r_size];
      weights = new double[r_size];
      t1 = cpu_time ( );
      nwspgr ( cc, cce_order, dim, k, r_size, s_size, nodes, weights );
      t2 = cpu_time ( );
      delete [] nodes;
      delete [] weights;
      cout << "  " << setw(10) << t2 - t1;
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void nwspgr_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NWSPGR_TEST tests NWSPGR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt.
//
{
  int dim;
  int k;
  double *nodes;
  int r_size;
  int s_size;
  double *weights;

  cout << "\n";
  cout << "NWSPGR_TEST:\n";
  cout << "  NWSPGR generates a sparse grid, based on either:\n";
  cout << "  one of the built-in 1D rules, or a family of 1D rules\n";
  cout << "  supplied by the user.\n";

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( kpu_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  nwspgr ( kpu, kpu_order, dim, k, r_size, s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Kronrod-Patterson, [0,1], Dim 2, Level 3" );
  delete [] nodes;
  delete [] weights;

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( kpn_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  nwspgr ( kpn, kpn_order, dim, k, r_size, s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Kronrod-Patterson, (-oo,+oo), Dim 2, Level 3" );
  delete [] nodes;
  delete [] weights;

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( gqu_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  nwspgr ( gqu, gqu_order, dim, k, r_size, s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, "  Gauss-Legendre, [0,1], Dim 2, Level 3" );
  delete [] nodes;
  delete [] weights;

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( gqn_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  nwspgr ( gqn, gqn_order, dim, k, r_size, s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, 
    "  Gauss Hermite, (-oo,+oo), Dim 2, Level 3" );
  delete [] nodes;
  delete [] weights;

  dim = 2;
  k = 3;
  r_size = nwspgr_size ( cce_order, dim, k );
  nodes = new double[dim*r_size];
  weights = new double[r_size];
  nwspgr ( cc, cce_order, dim, k, r_size, s_size, nodes, weights );
  quad_rule_print ( dim, s_size, nodes, weights, 
    "  Clenshaw Curtis Exponential, [-1,+1], Dim 2, Level 3" );
  delete [] nodes;
  delete [] weights;

  return;
}
//****************************************************************************80

void order_report ( )

//****************************************************************************80
//
//  Purpose:
//
//    ORDER_REPORT reports on the order of each family of rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int ap;
  int k;
  int kpn_order[5] = {
    1, 3, 9, 19, 35 };
  int l;
  int o;
  int rp;

  cout << "\n";
  cout << "ORDER_REPORT\n";
  cout << "  For each family of rules, report:\n";
  cout << "  L,  the level index,\n";
  cout << "  RP, the required polynomial precision,\n";
  cout << "  AP, the actual polynomial precision,\n";
  cout << "  O,  the rule order (number of points).\n";

  cout << "\n";
  cout << "  GQN family\n";
  cout << "  Gauss quadrature, exponential weight, (-oo,+oo)\n";
  cout << "\n";
  cout << "   L  RP  AP   O\n";
  cout << "\n";

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = l;
    ap = 2 * o - 1;
    cout << "  " << setw(2) << l
         << "  " << setw(2) << rp
         << "  " << setw(2) << ap
         << "  " << setw(2) << o << "\n";
  }

  cout << "\n";
  cout << "  GQU family\n";
  cout << "  Gauss quadrature, unit weight, [0,1]\n";
  cout << "\n";
  cout << "   L  RP  AP   O\n";
  cout << "\n";

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = l;
    ap = 2 * o - 1;
    cout << "  " << setw(2) << l
         << "  " << setw(2) << rp
         << "  " << setw(2) << ap
         << "  " << setw(2) << o << "\n";
  }

  cout << "\n";
  cout << "  KPN family\n";
  cout << "  Gauss-Kronrod-Patterson quadrature, exponential weight, (-oo,+oo)\n";
  cout << "\n";
  cout << "   L  RP  AP   O\n";
  cout << "\n";

  k = 1;
  o = 1;
  ap = 1;

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;

    while ( ap < rp )
    {
      if ( k == 5 )
      {
        cout << "\n";
        cout << "  No higher order rule is available!\n";
        break;
      }
//
//  Can we use a simple rule?
//
      if ( rp < kpn_order[k] )
      {
        o = rp;
        ap = rp;
      }
//
//  Otherwise, move to next higher rule.
//
      else
      {
        k = k + 1;
        ap = 2 * kpn_order[k-1] - kpn_order[k-2];
        o = kpn_order[k-1];
      }
    }
    cout << "  " << setw(2) << l
         << "  " << setw(2) << rp
         << "  " << setw(2) << ap
         << "  " << setw(2) << o << "\n";
  }

  cout << "\n";
  cout << "  KPU family\n";
  cout << "  Gauss-Kronrod-Patterson quadrature, unit weight, [0,1]\n";
  cout << "\n";
  cout << "   L  RP  AP   O\n";
  cout << "\n";

  for ( l = 1; l <= 25; l++ )
  {
    rp = 2 * l - 1;
    o = 1;
    ap = 1;
    while ( ap < rp )
    {
      o = 2 * ( o + 1 ) - 1;
      ap = ( 3 * o + 1 ) / 2;
    }
    cout << "  " << setw(2) << l
         << "  " << setw(2) << rp
         << "  " << setw(2) << ap
         << "  " << setw(2) << o << "\n";
  }

  return;
}
//****************************************************************************80

void symmetric_sparse_size_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    SYMMETRIC_SPARSE_SIZE_TEST tests SYMMETRIC_SPARSE_SIZE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt.
//
//  Local parameters:
//
//    Local, int D, the spatial dimension.
//
//    Local, int MAXK, the maximum level to check.
//
{
  int test_num = 3;

  int dim;
  int dim_test[3] = { 5, 5, 3 };
  double nodes1[6*5] = {
   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
   0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  double nodes2[21*5] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
    0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
  double nodes3[23*3] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 
    1.0, 1.0, 1.0, 1.0, 1.0, 1.73205, 1.73205, 1.73205, 2.33441, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.741964, 1.0, 1.0, 1.0, 1.73205, 1.73205, 2.33441, 
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.73205, 0.0, 0.0, 1.0, 0.0, 
    0.0, 0.741964, 1.0, 1.73205, 2.33441, 0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0, 1.73205, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
  int r;
  int r_test[3] = { 6, 21, 23 };
  int r2;
  int test;
  double x0;

  cout << "\n";
  cout << "SYMMETRIC_SPARSE_SIZE_TEST\n";
  cout << "  Given a symmetric sparse grid rule represented only by\n";
  cout << "  the points with positive values, determine the total number\n";
  cout << "  of points in the grid.\n";
  cout << "\n";
  cout << "  For dimension DIM, we report\n";
  cout << "  R, the number of points in the positive orthant, and\n";
  cout << "  R2, the total number of points.\n";
  cout << "\n";
  cout << "       DIM         R        R2\n";
  cout << "\n";

  x0 = 0.0;

  for ( test = 0; test < test_num; test++ )
  {
    r = r_test[test];
    dim = dim_test[test];
    if ( test == 0 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes1, x0 );
    }
    else if ( test == 1 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes2, x0 );
    }
    else if ( test == 2 )
    {
      r2 = symmetric_sparse_size ( r, dim, nodes3, x0 );
    }
    cout << "  " << setw(8) << dim
         << "  " << setw(8) << r
         << "  " << setw(8) << r2 << "\n";
  }

  return;
}
//****************************************************************************80

void tensor_product_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TENSOR_PRODUCT_TEST tests TENSOR_PRODUCT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int i1;
  int i2;
  int j;
  int n;
  int n1d;
  int order1 = 2;
  int order2 = 3;
  int order3 = 2;
  int *order1d;
  double *w1d;
  double *wnd;
  double w1_1d[2] = { 1.0, 1.0 };
  double w2_1d[3] = { 0.25, 0.50, 0.25 };
  double w3_1d[2] = { 2.50, 2.50 };
  double x1_1d[2] = { -1.0, +1.0 };
  double x2_1d[3] = { 2.0, 2.5, 3.0 };
  double x3_1d[2] = { 10.0, 15.0 };
  double *x1d;
  double *xnd;

  cout << "\n";
  cout << "TENSOR_PRODUCT_TEST:\n";
  cout << "  Given a sequence of 1D quadrature rules, construct the\n";
  cout << "  tensor product rule.\n";
//
//  1D rule.
//
  d = 1;
  order1d = new int[d];

  order1d[0] = order1;

  n1d = i4vec_sum ( d, order1d );
  x1d = new double[n1d];
  w1d = new double[n1d];

  n = i4vec_product ( d, order1d );
  xnd = new double[d*n];
  wnd = new double[n];

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );
  
  quad_rule_print ( d, n, xnd, wnd, "  A 1D rule over [-1,+1]:" );

  delete [] order1d;
  delete [] w1d;
  delete [] wnd;
  delete [] x1d;
  delete [] xnd;
//
//  2D rule.
//
  d = 2;
  order1d = new int[d];

  order1d[0] = order1;
  order1d[1] = order2;

  n1d = i4vec_sum ( d, order1d );
  x1d = new double[n1d];
  w1d = new double[n1d];

  n = i4vec_product ( d, order1d );
  xnd = new double[d*n];
  wnd = new double[n];

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order2; i++ )
  {
    x1d[j] = x2_1d[i];
    w1d[j] = w2_1d[i];
    j = j + 1;
  }

  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );
  
  quad_rule_print ( d, n, xnd, wnd, "  A 2D rule over [-1,+1] x [2.0,3.0]:" );

  delete [] order1d;
  delete [] w1d;
  delete [] wnd;
  delete [] x1d;
  delete [] xnd;
//
//  3D rule.
//
  d = 3;
  order1d = new int[d];

  order1d[0] = order1;
  order1d[1] = order2;
  order1d[2] = order3;

  n1d = i4vec_sum ( d, order1d );
  x1d = new double[n1d];
  w1d = new double[n1d];

  n = i4vec_product ( d, order1d );
  xnd = new double[d*n];
  wnd = new double[n];

  j = 0;
  for ( i = 0; i < order1; i++ )
  {
    x1d[j] = x1_1d[i];
    w1d[j] = w1_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order2; i++ )
  {
    x1d[j] = x2_1d[i];
    w1d[j] = w2_1d[i];
    j = j + 1;
  }
  for ( i = 0; i < order3; i++ )
  {
    x1d[j] = x3_1d[i];
    w1d[j] = w3_1d[i];
    j = j + 1;
  }

  tensor_product ( d, order1d, n1d, x1d, w1d, n, xnd, wnd );

  quad_rule_print ( d, n, xnd, wnd, 
    "  A 3D rule over [-1,+1] x [2.0,3.0] x [10.0,15.0]:" );

  delete [] order1d;
  delete [] w1d;
  delete [] wnd;
  delete [] x1d;
  delete [] xnd;

  return;
}
//****************************************************************************80

void tensor_product_cell_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TENSOR_PRODUCT_CELL_TEST tests TENSOR_PRODUCT_CELL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2013
//
//  Author:
//
//   John Burkardt
//
{
  int d;
  int i1;
  int i2;
  int n1d;
  int nc;
  int np;
  int nr[3] = { 2, 3, 2 };
  int order1 = 2;
  int order2 = 3;
  int order3 = 2;
  int *order1d;
  int *roff;
  double *w1d;
  double *wc;
  double *wp;
  double w1_1d[2] = { 1.0, 1.0 };
  double w2_1d[3] = { 0.25, 0.50, 0.25 };
  double w3_1d[2] = { 2.50, 2.50 };
  double x1_1d[2] = { -1.0, +1.0 };
  double x2_1d[3] = { 2.0, 2.5, 3.0 };
  double x3_1d[2] = { 10.0, 15.0 };
  double *x1d;
  double *xc;
  double *xp;

  cout << "\n";
  cout << "TENSOR_PRODUCT_TEST_CELL:\n";
  cout << "  Given a set of 1D quadrature rules stored in a cell array,\n";
  cout << "  construct the tensor product rule.\n";
//
//  We can construct ROFF once and for all.
//
  roff = r8cvv_offset ( 3, nr );
//
//  1D rule.
//
  d = 1;
  nc = i4vec_sum ( d, nr );
  xc = new double[nc];
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  wc = new double[nc];
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  np = i4vec_product ( d, nr );
  xp = new double[d * np];
  wp = new double[np];
  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );
  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  delete [] wc;
  delete [] wp;
  delete [] xc;
  delete [] xp;
//
//  2D rule.
//
  d = 2;
  nc = i4vec_sum ( d, nr );
  xc = new double[nc];
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  r8cvv_rset ( nc, xc, d, roff, 1, x2_1d );
  wc = new double[nc];
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  r8cvv_rset ( nc, wc, d, roff, 1, w2_1d );
  np = i4vec_product ( d, nr );
  xp = new double[d * np];
  wp = new double[np];

  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );

  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  delete [] wc;
  delete [] wp;
  delete [] xc;
  delete [] xp;
//
//  3D rule.
//
  d = 3;
  nc = i4vec_sum ( d, nr );
  xc = new double[nc];
  r8cvv_rset ( nc, xc, d, roff, 0, x1_1d );
  r8cvv_rset ( nc, xc, d, roff, 1, x2_1d );
  r8cvv_rset ( nc, xc, d, roff, 2, x3_1d );
  wc = new double[nc];
  r8cvv_rset ( nc, wc, d, roff, 0, w1_1d );
  r8cvv_rset ( nc, wc, d, roff, 1, w2_1d );
  r8cvv_rset ( nc, wc, d, roff, 2, w3_1d );
  np = i4vec_product ( d, nr );
  xp = new double[d * np];
  wp = new double[np];

  tensor_product_cell ( nc, xc, wc, d, nr, roff, np, xp, wp );

  quad_rule_print ( d, np, xp, wp, "  A 1D rule over [-1,+1]:" );

  delete [] roff;
  delete [] wc;
  delete [] wp;
  delete [] xc;
  delete [] xp;

  return;
}
