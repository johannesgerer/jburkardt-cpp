# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

# include "sandia_rules.hpp"
# include "sandia_sgmgg.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void sandia_sgmgg_coef_naive_test ( int dim_num, int point_num,
  int sparse_index[] );
void sandia_sgmgg_neighbor_naive_test ( int dim_num, int point_num,
  int sparse_index[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_PRB tests generalized sparse grid routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  webbur::timestamp ( );
  std::cout << "\n";
  std::cout << "SANDIA_SGMGG_PRB:\n";
  std::cout << "  C++ version.\n";
  std::cout << "  Test the SANDIA_SGMGG library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  std::cout << "\n";
  std::cout << "SANDIA_SGMGG_PRB:\n";
  std::cout << "  Normal end of execution.\n";
  std::cout << "\n";
  webbur::timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 demonstrates naive coefficient calculation.
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
  int dim_num;
  int point_num;
  int sparse_index1[2*7] = {
    0, 2,
    0, 3,
    1, 1,
    1, 2,
    2, 0,
    2, 1,
    3, 0 };
  int sparse_index2[3*19] = {
    0, 1, 0,
    0, 2, 0,
    0, 3, 0,
    1, 0, 0,
    1, 1, 0,
    1, 2, 0,
    2, 0, 0,
    2, 1, 0,
    3, 0, 0,
    0, 0, 1,
    0, 1, 1,
    0, 2, 1,
    1, 0, 1,
    1, 1, 1,
    2, 0, 1,
    0, 0, 2,
    0, 1, 2,
    1, 0, 2,
    0, 0, 3 };
  int sparse_index3[2*8] = {
    0, 2,
    1, 1,
    1, 2,
    2, 1,
    3, 0,
    3, 1,
    4, 0,
    5, 0 };
  int sparse_index4[2*8] = {
    0, 0,
    0, 1,
    0, 2,
    0, 3,
    1, 0,
    1, 1,
    2, 0,
    3, 0 };

  std::cout << "\n";
  std::cout << "TEST01:\n";
  std::cout << "  Demonstrate naive coefficient calculations.\n";
//
//  Isotropic grid in 2D.
//
  std::cout << "\n";
  std::cout << "  1) Isotropic grid in 2D:\n";
  dim_num = 2;
  point_num = 7;

  sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index1 );
//
//  Isotropic grid in 3D.
//
  std::cout << "\n";
  std::cout << "  2) Isotropic grid in 3D:\n";
  dim_num = 3;
  point_num = 19;

  sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index2 );
//
//  Anisotropic grid in 2D.
//
  std::cout << "\n";
  std::cout << "  3) Ansotropic grid in 2D:\n";
  dim_num = 2;
  point_num = 8;

  sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index3 );
//
//  Generalized grid in 2D.
//
  std::cout << "\n";
  std::cout << "  4) Generalized grid in 2D:\n";
  dim_num = 2;
  point_num = 8;

  sandia_sgmgg_coef_naive_test ( dim_num, point_num, sparse_index4 );

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates naive neighbor calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int point_num;
  int *sparse_index;
  int sparse_index1[2*7] = {
    0, 2,
    0, 3,
    1, 1,
    1, 2,
    2, 0,
    2, 1,
    3, 0 };
  int sparse_index2[3*19] = {
    0, 1, 0,
    0, 2, 0,
    0, 3, 0,
    1, 0, 0,
    1, 1, 0,
    1, 2, 0,
    2, 0, 0,
    2, 1, 0,
    3, 0, 0,
    0, 0, 1,
    0, 1, 1,
    0, 2, 1,
    1, 0, 1,
    1, 1, 1,
    2, 0, 1,
    0, 0, 2,
    0, 1, 2,
    1, 0, 2,
    0, 0, 3 };
  int sparse_index3[2*8] = {
    0, 2,
    1, 1,
    1, 2,
    2, 1,
    3, 0,
    3, 1,
    4, 0,
    5, 0 };
  int sparse_index4[2*8] = {
    0, 0,
    0, 1,
    0, 2,
    0, 3,
    1, 0,
    1, 1,
    2, 0,
    3, 0 };

  std::cout << "\n";
  std::cout << "TEST02:\n";
  std::cout << "  Demonstrate naive neighbor calculations.\n";
//
//  Isotropic grid in 2D.
//
  std::cout << "\n";
  std::cout << "  1) Isotropic grid in 2D:\n";
  dim_num = 2;
  point_num = 7;

  sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index1  );
//
//  Isotropic grid in 3D.
//
  std::cout << "\n";
  std::cout << "  2) Isotropic grid in 3D:\n";
  dim_num = 3;
  point_num = 19;

  sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index2  );
//
//  Anisotropic grid in 2D.
//
  std::cout << "\n";
  std::cout << "  3) Ansotropic grid in 2D:\n";
  dim_num = 2;
  point_num = 8;

  sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index3  );
//
//  Generalized grid in 2D.
//
  std::cout << "\n";
  std::cout << "  4) Generalized grid in 2D:\n";
  dim_num = 2;
  point_num = 8;

  sandia_sgmgg_neighbor_naive_test ( dim_num, point_num, sparse_index4  );

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 sets up the GG data structure.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
{
# define GG_NA_01 4
# define GG_ND_01 2
# define GG_NI_01 10
# define GG_NO_01 6

  int *gg_a;
  int gg_a_01[GG_NA_01] = { 3, 6, 8, 9 };
  int *gg_b;
  int gg_b_01[GG_ND_01*GG_NI_01] = {
    -1, -1,
    -1,  0,
    -1,  1,
    -1,  2,
     0, -1,
     1,  4,
     2,  5,
     4, -1,
     5,  7,
     7, -1 };
  int *gg_f;
  int gg_f_01[GG_ND_01*GG_NI_01] = {
     4,  1,
     5,  2,
     6,  3,
    -1, -1,
     7,  5,
     8,  6,
    -1, -1,
     9,  8,
    -1, -1,
    -1, -1 };
  double *gg_g;
  double gg_g_01[GG_NI_01] = {
    0.1, 1.1, 2.2, 3.0, 1.0,
    2.1, 3.2, 2.0, 3.3, 3.1 };
  int *gg_i;
  int gg_i_01[GG_ND_01 * GG_NI_01] = {
    0, 0,
    0, 1,
    0, 2,
    0, 3,
    1, 0,
    1, 1,
    1, 2,
    2, 0,
    2, 1,
    3, 0 };
  int gg_ma = 100;
  int gg_mi = 100;
  int gg_mo = 100;
  int gg_na;
  int gg_nd;
  int gg_ni;
  int gg_no;
  int *gg_o;
  int gg_o_01[GG_NO_01] = { 0, 1, 2, 4, 5, 7 };
  int i;
  int max_index;

  std::cout << "\n";
  std::cout << "TEST03:\n";
  std::cout << "  Set up examples of the GG (Gerstner-Griebel)\n";
  std::cout << "  data structure for generalized sparse grids.\n";
//
//  Isotropic grid in 2D.
//
//  3 | 4
//  2 | 3  7
//  1 | 2  6  9
//  0 | 1  5  8 10
//    +-----------
//      0  1  2  3
//
  std::cout << "\n";
  std::cout << "  1) Isotropic grid in 2D\n";
//
//  Get sizes for this problem.
//
  gg_na = GG_NA_01;
  gg_nd = GG_ND_01;
  gg_ni = GG_NI_01;
  gg_no = GG_NO_01;
//
//  Create work arrays.
//  Because these will grow and shrink during the computation, we must give
//  them an initial size big enough that we don't need to resize.
//
  gg_a = new int    [         gg_ma ];
  gg_b = new int    [ gg_nd * gg_mi ];
  gg_f = new int    [ gg_nd * gg_mi ];
  gg_g = new double [         gg_mi ];
  gg_i = new int    [ gg_nd * gg_mi ];
  gg_o = new int    [         gg_mo ];
//
//  Copy the initial problem data vectors.
//
  webbur::i4vec_copy ( gg_na,        gg_a_01, gg_a );
  webbur::i4mat_copy ( gg_nd, gg_ni, gg_b_01, gg_b );
  webbur::i4mat_copy ( gg_nd, gg_ni, gg_f_01, gg_f );
  webbur::r8vec_copy ( gg_ni,        gg_g_01, gg_g );
  webbur::i4mat_copy ( gg_nd, gg_ni, gg_i_01, gg_i );
  webbur::i4vec_copy ( gg_no,        gg_o_01, gg_o );
//
//  Implicit decreasing heap sort on GG_G and GG_A.
//
  std::cout << "\n";
  std::cout << "  Before Heap:\n";
  std::cout << "     I     A      G\n";
  std::cout << "\n";

  for ( i = 0; i < gg_na; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(4) << gg_a[i]
              << "  " << std::setw(14) << gg_g[gg_a[i]] << "\n";
  }

  webbur::r8vec_indexed_heap_d ( gg_na, gg_g, gg_a );

  std::cout << "\n";
  std::cout << "  After Heap:\n";
  std::cout << "     I     A      G\n";
  std::cout << "\n";

  for ( i = 0; i < gg_na; i++ )
  {
    std::cout << "  " << std::setw(4) << i
              << "  " << std::setw(4) << gg_a[i]
              << "  " << std::setw(14) << gg_g[gg_a[i]] << "\n";
  }
//
//  One step of the algorithm is to identify MAX_INDEX, the active index with
//  maximum G value, and move it to the OLD array.
//
  max_index = webbur::r8vec_indexed_heap_d_extract ( &gg_na, gg_g, gg_a );
  std::cout << "\n";
  std::cout << "  Transferring index " << max_index
            << " from active to old set.\n";
//
//  GG_NA was decremented.
//  The heap structure of GG_A was restored.
//  To complete the update, we must
//    append MAX_INDEX to GG_O.
//    update the active set with the forward neighbors of MAX_INDEX.
//
  gg_o[gg_no] = max_index;
  gg_no = gg_no + 1;
//
//  Compute forward neighbors of MAX_INDEX.
//

//
//  DISPLAY CURRENT DATA.
//
  delete [] gg_a;
  delete [] gg_b;
  delete [] gg_f;
  delete [] gg_g;
  delete [] gg_i;
  delete [] gg_o;

  return;
# undef GG_NA_01
# undef GG_ND_01
# undef GG_NI_01
# undef GG_NO_01
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 simulates a set of incremental coefficient calculations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 January 2011
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int point_num;
  int point_num2;
  int sparse_index[2*8] = {
    0, 0,
    0, 1,
    0, 2,
    0, 3,
    1, 0,
    1, 1,
    2, 0,
    3, 0 };

  std::cout << "\n";
  std::cout << "TEST04:\n";
  std::cout << "  Simulate incremental coefficient calculations.\n";
//
//  Generalized grid in 2D.
//
  std::cout << "\n";
  std::cout << "  Generalized grid in 2D:\n";
  dim_num = 2;
  point_num = 8;

  for ( point_num2 = 1; point_num2 <= point_num; point_num2++ )
  {
    sandia_sgmgg_coef_naive_test ( dim_num, point_num2, sparse_index );
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 demonstrates SANDIA_SGMGG_COEF_INC2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N1 5

  int c1[N1] = {
    +1,
    +1,
    +1,
    -1,
    -1
  };
  int *c3;
  int i;
  int j;
  int m = M;
  int n1 = N1;
  int s1[M*N1] = {
    0, 2,
    1, 1,
    2, 0,
    0, 1,
    1, 0 };
  int s2[M] = {
    1, 2 };

  std::cout << "\n";
  std::cout << "TEST05:\n";
  std::cout << "  Predict new coefficients given candidate index.\n";

  std::cout << "\n";
  std::cout << "  Spatial dimension M = " << m << "\n";
  std::cout << "  Number of items in active set N1 = " << n1 << "\n";
  std::cout << "\n";
  std::cout << "  Index  Coef   Indices\n";
  std::cout << "\n";
  for ( j = 0; j < n1; j++ )
  {
    std::cout << "    " << std::setw(2) << j << ":"
              << "  " << std::setw(4) << c1[j] << "  ";
    for ( i = 0; i < m; i++ )
    {
      std::cout << "  " << std::setw(2) << s1[i+j*m];
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  Candidate:   ";
  for ( i = 0; i < m; i++ )
  {
    std::cout << "  " << std::setw(2) << s2[i];
  }
  std::cout << "\n";
//
//  Generalized grid in 2D.
//
  c3 = new int[n1+1];

  webbur::sandia_sgmgg_coef_inc2 ( m, n1, s1, c1, s2, c3 );

  std::cout << "\n";
  std::cout << "  Index  Coef  Coef\n";
  std::cout << "          Old   New\n";
  std::cout << "\n";
  for ( i = 0; i < n1; i++ )
  {
    std::cout << "    " << std::setw(2) << i << ":"
              << "  " << std::setw(4) << c1[i]
              << "  " << std::setw(4) << c3[i] << "\n";
  }

  std::cout << "    " << std::setw(2) << n1 << ":"
            << "  " << "    "
            << "  " << std::setw(4) << c3[n1] << "\n";

  std::cout << "    --   ----  ----\n";
  std::cout << "   Sum:"
            << "  " << std::setw(4) << webbur::i4vec_sum ( n1,     c1 )
            << "  " << std::setw(4) << webbur::i4vec_sum ( n1 + 1, c3 ) << "\n";
  delete [] c3;

  return;
# undef M
# undef N1
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests a special case for SANDIA_SGMGG_COEF_INC2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2011
//
//  Author:
//
//    John Burkardt
//
{
# define M 2
# define N1 5

  int c1[N1] = {
    +1,
    +1,
    +1,
    -1,
    -1
  };
  int *c3;
  int i;
  int j;
  int m = M;
  int n1 = N1;
  int s1[M*N1] = {
    2, 0,
    1, 1,
    0, 2,
    1, 0,
    0, 1 };
  int s2[M] = {
    3, 0 };

  std::cout << "\n";
  std::cout << "TEST06:\n";
  std::cout << "  Predict new coefficients given candidate index.\n";

  std::cout << "\n";
  std::cout << "  Spatial dimension M = " << m << "\n";
  std::cout << "  Number of items in active set N1 = " << n1 << "\n";
  std::cout << "\n";
  std::cout << "  Index  Coef   Indices\n";
  std::cout << "\n";
  for ( j = 0; j < n1; j++ )
  {
    std::cout << "    " << std::setw(2) << j << ":"
              << "  " << std::setw(4) << c1[j] << "  ";
    for ( i = 0; i < m; i++ )
    {
      std::cout << "  " << std::setw(2) << s1[i+j*m];
    }
    std::cout << "\n";
  }
  std::cout << "\n";
  std::cout << "  Candidate:   ";
  for ( i = 0; i < m; i++ )
  {
    std::cout << "  " << std::setw(2) << s2[i];
  }
  std::cout << "\n";
//
//  Generalized grid in 2D.
//
  c3 = new int[n1+1];

  webbur::sandia_sgmgg_coef_inc2 ( m, n1, s1, c1, s2, c3 );

  std::cout << "\n";
  std::cout << "  Index  Coef  Coef\n";
  std::cout << "          Old   New\n";
  std::cout << "\n";
  for ( i = 0; i < n1; i++ )
  {
    std::cout << "    " << std::setw(2) << i << ":"
              << "  " << std::setw(4) << c1[i]
              << "  " << std::setw(4) << c3[i] << "\n";
  }

  std::cout << "    " << std::setw(2) << n1 << ":"
            << "  " << "    "
            << "  " << std::setw(4) << c3[n1] << "\n";

  std::cout << "    --   ----  ----\n";
  std::cout << "   Sum:"
            << "  " << std::setw(4) << webbur::i4vec_sum ( n1,     c1 )
            << "  " << std::setw(4) << webbur::i4vec_sum ( n1 + 1, c3 ) << "\n";
  delete [] c3;

  return;
# undef M
# undef N1
}
//****************************************************************************80

void sandia_sgmgg_coef_naive_test ( int dim_num, int point_num,
  int sparse_index[]  )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_COEF_NAIVE_TEST demonstrates "naive" coefficient computations.
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
  int *coef;
  int coef_sum;

  webbur::i4mat_transpose_print ( dim_num, point_num, sparse_index,
    "  SPARSE_INDEX matrix:" );

  coef = new int[point_num];

  webbur::sandia_sgmgg_coef_naive ( dim_num, point_num, sparse_index, coef );

  webbur::i4vec_print ( point_num, coef, "  COEF vector:" );

  coef_sum = webbur::i4vec_sum ( point_num, coef );

  std::cout << "       ---         -\n";
  std::cout << "       Sum:        " << coef_sum << "\n";

  delete [] coef;

  return;
}
//****************************************************************************80

void sandia_sgmgg_neighbor_naive_test ( int dim_num, int point_num,
  int sparse_index[]  )

//****************************************************************************80
//
//  Purpose:
//
//    SANDIA_SGMGG_NEIGHBOR_NAIVE_TEST demonstrates "naive" neighbor computations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *neighbor;

  webbur::i4mat_transpose_print ( dim_num, point_num, sparse_index,
    "  SPARSE_INDEX matrix:" );

  neighbor = new int[dim_num*point_num];

  webbur::sandia_sgmgg_neighbor_naive ( dim_num, point_num, sparse_index,
    neighbor );

  webbur::i4mat_transpose_print ( dim_num, point_num, neighbor,
    "  NEIGHBOR matrix:" );

  delete [] neighbor;

  return;
}
