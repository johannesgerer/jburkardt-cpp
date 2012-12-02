# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void basis_mn_tet4 ( double t[3*4], int n, double p[], double phi[] );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
double *fem3d_transfer ( int sample_node_num, int sample_element_order, 
  int sample_element_num, int sample_value_dim, int sample_value_num, 
  double sample_node_xyz[], int sample_element_node[], 
  int sample_element_neighbor[], double sample_value[], int fem_node_num, 
  int fem_element_order, int fem_element_num, int fem_value_dim, 
  int fem_value_num, double fem_node_xyz[], int fem_element_node[] );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
int i4mat_min ( int m, int n, int a[] );
double *projection ( int fem_node_num, double fem_node_xyz[], 
  int fem_element_order, int fem_element_num, int fem_element_node[], 
  int fem_element_neighbor[], int fem_value_dim, double fem_value[], 
  int sample_node_num, double sample_node_xyz[] );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_abs ( double x );
double r8_min ( double x, double y );
double *r8ge_fss_new ( int n, double a[], int nb, double b[] );
double *r8mat_data_read ( string input_filename, int m, int n );
double r8mat_det_4d ( double a[] );
void r8mat_header_read ( string input_filename, int *m, int *n );
int r8mat_solve ( int n, int rhs_num, double a[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8mat_zero_new ( int m, int n );
bool r8vec_is_nonnegative ( int n, double x[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
int *tet_mesh_neighbor_tets ( int tetra_order, int tetra_num, 
  int tetra_node[] );
int tet_mesh_search_delaunay ( int node_num, double node_xyz[], int tet_order, 
  int tet_num, int tet_node[], int tet_neighbor[], double p[], int *face, 
  int *step_num );
int tet_mesh_search_naive ( int node_num, double node_xyz[],
  int tet_order, int tet_num, int tet_node[], double p[], int *step_num );
double *tetrahedron_barycentric ( double tetra[3*4], double p[3] );
double tetrahedron_volume ( double tetra[3*4] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM3D_PROJECT.
//
//  Discussion:
//
//    FEM3D_PROJECT reads files defining a sampling of a (scalar or vector)
//    function of 3 arguments, and a list of nodes and tetrahedral elements 
//    to use for a finite element representation of the data.
//
//    It computes a set of finite element coefficients to be associated with
//    the given finite element mesh, and writes that information to a file
//    so that an FEM representation is formed by the node, element and value 
//    files.
//
//  Usage:
//
//    fem3d_project sample_prefix fem_prefix
//
//    where 'sample_prefix' is the common prefix for the SAMPLE files:
//
//    * sample_prefix_nodes.txt,     the node coordinates where samples were taken,
//    * sample_prefix_elements.txt,  the 4 nodes that make up each element;
//    * sample_prefix_values.txt,    the sample values.
//
//    and 'fem_prefix' is the common prefix for the FEM files:
//
//    * fem_prefix_nodes.txt,    the node coordinates.
//    * fem_prefix_elements.txt, the 4 nodes that make up each element;
//    * fem_prefix_values.txt,   the values defined at each node (output).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int element_min;
  string fem_element_filename;
  int *fem_element_node;
  int fem_element_num;
  int fem_element_order;
  int fem_node_dim;
  string fem_node_filename;
  int fem_node_num;
  double *fem_node_xyz;
  string fem_prefix;
  double *fem_value;
  int fem_value_dim;
  string fem_value_filename;
  int fem_value_num;
  int i;
  int j;
  string sample_element_filename;
  int *sample_element_neighbor;
  int *sample_element_node;
  int sample_element_num;
  int sample_element_order;
  string sample_prefix;
  int sample_node_dim;
  string sample_node_filename;
  int sample_node_num;
  double *sample_node_xyz;
  int sample_value_dim;
  int sample_value_num;
  double *sample_value;
  string sample_value_filename;

  timestamp ( );

  cout << "\n";
  cout << "FEM3D_PROJECT\n";
  cout << "  C++ version.\n";
  cout << "\n";
  cout << "  Read files defining a sampling of a function of 3 arguments.\n";
  cout << "  Read files defining a finite element mesh.\n";
  cout << "  Project the sample data onto the mesh, and\n";
  cout << "  write a file of FEM coefficient values.\n";
//
//  Get the number of command line arguments.
//
  if ( 1 < argc )
  {
    sample_prefix = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "Enter the sample file prefix:\n";
    cin >> sample_prefix;
  }

  if ( 2 < argc )
  {
    fem_prefix = argv[2];
  }
  else
  {
    cout << "\n";
    cout << "Enter the FEM file prefix:\n";
    cin >> fem_prefix;
  }
//
//  Create the filenames.
//
  sample_node_filename = sample_prefix + "_nodes.txt";
  sample_element_filename = sample_prefix + "_elements.txt";
  sample_value_filename = sample_prefix + "_values.txt";

  fem_node_filename = fem_prefix + "_nodes.txt";
  fem_element_filename = fem_prefix + "_elements.txt";
  fem_value_filename = fem_prefix + "_values.txt";
//
//  Read the SAMPLE NODE, ELEMENT and VALUE data.
//
  r8mat_header_read ( sample_node_filename, &sample_node_dim, 
    &sample_node_num );

  cout << "\n";
  cout << "  Sample node spatial dimension is " << sample_node_dim << "\n";
  cout << "  Sample node number is            " << sample_node_num << "\n";

  if ( sample_node_dim != 3 )
  {
    cout << "\n";
    cout << "FEM3D_PROJECT - Fatal error!\n";
    cout << "  Spatial dimension of the sample nodes is not 3.\n";
    exit ( 1 );
  }

  sample_node_xyz = r8mat_data_read ( sample_node_filename, sample_node_dim, 
    sample_node_num );

  i4mat_header_read ( sample_element_filename, &sample_element_order, 
    &sample_element_num );

  cout << "\n";
  cout << "  Sample element order is  " << sample_element_order << "\n";
  cout << "  Sample element number is " << sample_element_num << "\n";

  if ( sample_element_order != 4 )
  {
    cout << "\n";
    cout << "FEM3D_PROJECT - Fatal error!\n";
    cout << "  The sample element order must be 4.\n";
    exit ( 1 );
  }

  sample_element_node = new int[sample_element_order*sample_element_num];

  sample_element_node = i4mat_data_read ( sample_element_filename, 
    sample_element_order, sample_element_num );

  element_min = i4mat_min ( sample_element_order, sample_element_num, 
    sample_element_node );

  if ( element_min == 1 )
  {
    cout << "\n";
    cout << "  Converting 1-based sample element array to 0 base.\n";
    for ( j = 0; j < sample_element_num; j++ )
    {
      for ( i = 0; i < sample_element_order; i++ )
      {
        sample_element_node[i+j*sample_element_order] = 
        sample_element_node[i+j*sample_element_order] - 1;
      }
    }
  }

  r8mat_header_read ( sample_value_filename, &sample_value_dim,
    &sample_value_num );

  cout << "\n";
  cout << "  The sample value dimension is    " << sample_value_dim << "\n";
  cout << "  The sample value number is        " << sample_value_num << "\n";

  if ( sample_value_num != sample_node_num )
  {
    cout << "\n";
    cout << "FEM3D_PROJECT - Fatal error!\n";
    cout << "  Number of sample values and nodes differ.\n";
    exit ( 1 );
  }
  sample_value = r8mat_data_read ( sample_value_filename, sample_value_dim, 
    sample_value_num );
//
//  Create the sample element neighbor array.
//
  sample_element_neighbor = tet_mesh_neighbor_tets ( sample_element_order,
    sample_element_num, sample_element_node );

  cout << "\n";
  cout << "  The element neighbor array has been computed.\n";
//
//  Read the FEM NODE and ELEMENT data.
//
  r8mat_header_read ( fem_node_filename, &fem_node_dim, &fem_node_num );

  cout << "\n";
  cout << "  The FEM node dimension is        " << fem_node_dim << "\n";
  cout << "  The FEM node number is           " << fem_node_num << "\n";

  if ( fem_node_dim != 3 )
  {
    cout << "\n";
    cout << "FEM3D_PROJECT - Fatal error!\n";
    cout << "  Spatial dimension of the nodes is not 3.\n";
    exit ( 1 );
  }

  fem_node_xyz = r8mat_data_read ( fem_node_filename, fem_node_dim, fem_node_num );

  i4mat_header_read ( fem_element_filename, &fem_element_order, &fem_element_num );

  cout << "  The FEM element order is         " << fem_element_order << "\n";
  cout << "  The FEM element number is        " << fem_element_num << "\n";

  if ( fem_element_order != 4 )
  {
    cout << "\n";
    cout << "FEM3D_PROJECT - Fatal error!\n";
    cout << "  The FEM element order is not 4.\n";
    exit ( 1 );
  }

  fem_element_node = i4mat_data_read ( fem_element_filename, fem_element_order, 
    fem_element_num );

  element_min = i4mat_min ( fem_element_order, fem_element_num, 
    fem_element_node );

  if ( element_min == 1 )
  {
    cout << "\n";
    cout << "  Converting 1-based FEM element array to 0 base.\n";
    for ( j = 0; j < fem_element_num; j++ )
    {
      for ( i = 0; i < fem_element_order; i++ )
      {
        fem_element_node[i+j*fem_element_order] = 
        fem_element_node[i+j*fem_element_order] - 1;
      }
    }
  }
//
//  Compute the FEM values.
//
  fem_value_dim = sample_value_dim;
  fem_value_num = fem_node_num;

  fem_value = fem3d_transfer ( sample_node_num, sample_element_order, 
    sample_element_num, sample_value_dim, sample_value_num, 
    sample_node_xyz, sample_element_node, sample_element_neighbor, sample_value, 
    fem_node_num, fem_element_order, 
    fem_element_num, fem_value_dim, fem_value_num, 
    fem_node_xyz, fem_element_node );
//
//  Write the FEM values.
//
  r8mat_write ( fem_value_filename, fem_value_dim, fem_value_num, 
    fem_value );

  cout << "\n";
  cout << "  FEM value data written to \"" << fem_value_filename << "\"\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM3D_PROJECT\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  delete [] fem_element_node;
  delete [] fem_node_xyz;
  delete [] fem_value;

  delete [] sample_element_neighbor;
  delete [] sample_element_node;
  delete [] sample_node_xyz;
  delete [] sample_value;

  return 0;
}
//****************************************************************************80

void basis_mn_tet4 ( double t[3*4], int n, double p[], double phi[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_TET4: all bases at N points for a TET4 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a tetrahedron.
//
//    It works directly with these coordinates, and does not refer to a
//    reference element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54.
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], the points where the basis functions
//    are to be evaluated.
//
//    Output, double PHI[4*N], the value of the basis functions
//    at the evaluation points.
//
{
  int j;
  double volume;
//
//           | x1 x2 x3 x4 |
//  Volume = | y1 y2 y3 y4 |
//           | z1 z2 z3 z4 |
//           |  1  1  1  1 |
//
  volume =                             
      t[0+0*3] * (                       
        t[1+1*3] * ( t[2+2*3] - t[2+3*3] )   
      - t[1+2*3] * ( t[2+1*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+1*3] - t[2+2*3] ) ) 
    - t[0+1*3] * (                       
        t[1+0*3] * ( t[2+2*3] - t[2+3*3] )   
      - t[1+2*3] * ( t[2+0*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+0*3] - t[2+2*3] ) ) 
    + t[0+2*3] * (                       
        t[1+0*3] * ( t[2+1*3] - t[2+3*3] )   
      - t[1+1*3] * ( t[2+0*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+0*3] - t[2+1*3] ) ) 
    - t[0+3*3] * (                       
        t[1+0*3] * ( t[2+1*3] - t[2+2*3] )   
      - t[1+1*3] * ( t[2+0*3] - t[2+2*3] )   
      + t[1+2*3] * ( t[2+0*3] - t[2+1*3] ) );

  if ( volume == 0.0 )
  {
    cerr << "\n";
    cerr << "BASIS_MN_TET4 - Fatal error!\n";
    cerr << "  Element has zero volume.\n";
    exit ( 1 );
  }
//
//             | xp x2 x3 x4 |
//  Phi(1,P) = | yp y2 y3 y4 | / volume
//             | zp z2 z3 z4 |
//             |  1  1  1  1 |
//
  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] = (                           
        p[0+j*3] * (                         
          t[1+1*3] * ( t[2+2*3] - t[2+3*3] )   
        - t[1+2*3] * ( t[2+1*3] - t[2+3*3] )   
        + t[1+3*3] * ( t[2+1*3] - t[2+2*3] ) ) 
      - t[0+1*3] * (                           
          p[1+j*3] * ( t[2+2*3] - t[2+3*3] )   
        - t[1+2*3] * ( p[2+j*3] - t[2+3*3] )   
        + t[1+3*3] * ( p[2+j*3] - t[2+2*3] ) ) 
      + t[0+2*3] * (                           
          p[1+j*3] * ( t[2+1*3] - t[2+3*3] )   
        - t[1+1*3] * ( p[2+j*3] - t[2+3*3] )   
        + t[1+3*3] * ( p[2+j*3] - t[2+1*3] ) ) 
      - t[0+3*3] * (                           
          p[1+j*3] * ( t[2+1*3] - t[2+2*3] )   
        - t[1+1*3] * ( p[2+j*3] - t[2+2*3] )   
        + t[1+2*3] * ( p[2+j*3] - t[2+1*3] ) ) ) / volume;
//
//             | x1 xp x3 x4 |
//  Phi(2,P) = | y1 yp y3 y4 | / volume
//             | z1 zp z3 z4 |
//             |  1  1  1  1 |
//
    phi[1+j*4] = (                             
        t[0+0*3] * (                             
          p[1+j*3] * ( t[2+2*3] - t[2+3*3] )     
        - t[1+2*3] * ( p[2+j*3] - t[2+3*3] )     
        + t[1+3*3] * ( p[2+j*3] - t[2+2*3] ) )   
      - p[0+j*3]   * (                         
          t[1+0*3] * ( t[2+2*3] - t[2+3*3] )     
        - t[1+2*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - t[2+2*3] ) )   
      + t[0+2*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - p[2+j*3] ) ) 
      - t[0+3*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+2*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+2*3] )     
        + t[1+2*3] * ( t[2+0*3] - p[2+j*3] ) ) ) / volume;
//
//             | x1 x2 xp x4 |
//  Phi(3,P) = | y1 y2 yp y4 | / volume
//             | z1 z2 zp z4 |
//             |  1  1  1  1 |
//
    phi[2+j*4] = (                              
        t[0+0*3] * (                             
          t[1+1*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+1*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+1*3] - p[2+j*3] ) ) 
      - t[0+1*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - p[2+j*3] ) ) 
      + p[0+j*3] * (                           
          t[1+0*3] * ( t[2+1*3] - t[2+3*3] )     
        - t[1+1*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - t[2+1*3] ) )   
      - t[0+3*3] * (                             
          t[1+0*3] * ( t[2+1*3] - p[2+j*3] )   
        - t[1+1*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+1*3] ) ) ) / volume;
//
//             | x1 x2 x3 xp |
//  Phi(4,P) = | y1 y2 y3 yp | / volume
//             | z1 z2 z3 zp |
//             |  1  1  1  1 |
//
    phi[3+j*4] = (                             
        t[0+0*3] * (                             
          t[1+1*3] * ( t[2+2*3] - p[2+j*3] )   
        - t[1+2*3] * ( t[2+1*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+1*3] - t[2+2*3] ) )   
      - t[0+1*3] * (                             
          t[1+0*3] * ( t[2+2*3] - p[2+j*3] )   
        - t[1+2*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+2*3] ) )   
      + t[0+2*3] * (                             
          t[1+0*3] * ( t[2+1*3] - p[2+j*3] )   
        - t[1+1*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+1*3] ) )   
      - p[0+j*3] * (                           
          t[1+0*3] * ( t[2+1*3] - t[2+2*3] )     
        - t[1+1*3] * ( t[2+0*3] - t[2+2*3] )     
        + t[1+2*3] * ( t[2+0*3] - t[2+1*3] ) ) ) / volume;
  }
  return;
}
//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  return ( ch1 == ch2 );
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

double *fem3d_transfer ( int sample_node_num, int sample_element_order, 
  int sample_element_num, int sample_value_dim, int sample_value_num, 
  double sample_node_xyz[], int sample_element_node[], 
  int sample_element_neighbor[], double sample_value[], int fem_node_num, 
  int fem_element_order, int fem_element_num, int fem_value_dim, 
  int fem_value_num, double fem_node_xyz[], int fem_element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEM3D_TRANSFER "transfers" from one finite element mesh to another.
//
//  BAD THINGS:
//
//    1) the linear system A*X=B is defined with A being a full storage matrix.
//    2) the quadrature rule used is low order.
//    3) the elements are assumed to be linear.
//
//  Discussion:
//
//    We are also given a set of "sample" finite element function defined
//    by SAMPLE_NODE_XYZ, SAMPLE_ELEMENT, and SAMPLE_VALUE.
//
//    We are given a second finite element mesh, FEM_NODE_XYZ and 
//    FEM_ELEMENT_NODE.
//
//    Our aim is to "project" the sample data values into the finite element 
//    space, that is, to come up with a finite element function FEM_VALUE which
//    well approximates the sample data.
//
//    Now let W(x,y,z) represent a function interpolating the sample data, and
//    let Vijk(x,y,z) represent the finite element basis function associated with
//    node IJK.
//
//    Then we seek the coefficient vector U corresponding to a finite element
//    function U(x,y,z) of the form:
//
//      U(x,y,z) = sum ( 1 <= IJK <= N ) Uijk * Vijk(x,y,z)
//
//    To determine the coefficent vector entries U, we form a set of
//    projection equations.  For node IJK at grid point (I,J,K), the associated
//    basis function Vk(x,y,z) is used to pose the equation:
//
//      Integral U(x,y,z) Vijk(x,y,z) dx dy dz 
//        = Integral W(x,y,z) Vijk(x,y,z) dx dy dz
//
//    The left hand side is the usual stiffness matrix times the desired
//    coefficient vector U.  To complete the system, we simply need to
//    determine the right hand side, that is, the integral of the data function
//    W against the basis function Vk.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SAMPLE_NODE_NUM, the number of nodes.
//
//    Input, int SAMPLE_ELEMENT_ORDER, the element order.
//
//    Input, int SAMPLE_ELEMENT_NUM, the number of elements.
//
//    Input, int SAMPLE_VALUE_DIM, the value dimension.
//
//    Input, int SAMPLE_VALUE_NUM, the number of values.
//
//    Input, double SAMPLE_NODE_XYZ[3*SAMPLE_NODE_NUM], the nodes.
//
//    Input, int SAMPLE_ELEMENT_NODE[SAMPLE_ELEMENT_ORDER*SAMPLE_ELEMENT_NUM],
//    the nodes that make up each element.
//
//    Input, int SAMPLE_ELEMENT_NEIGHBOR[3*SAMPLE_ELEMENT_NUM],
//    the neighbor triangles.
//
//    Input, double SAMPLE_VALUE[SAMPLE_VALUE_DIM*SAMPLE_NODE_NUM],
//    the values.
//
//    Input, int FEM_NODE_NUM, the number of nodes.
//
//    Input, int FEM_ELEMENT_ORDER, the element order.
//
//    Input, int FEM_ELEMENT_NUM, the number of elements.
//
//    Input, int FEM_VALUE_DIM, the value dimension.
//
//    Input, int FEM_VALUE_NUM, the number of values.
//
//    Input, double FEM_NODE_XYZ[3*FEM_NODE_NUM], the nodes.  
//
//    Input, int FEM_ELEMENT_NODE[FEM_ELEMENT_ORDER*FEM_ELEMENT_NUM],
//    the nodes that make up each element.
//
//    Output, double FEM3D_TRANSFER[FEM_VALUE_DIM*FEM_VALUE_NUM], 
//    the values.
//
{
  double *a;
  double *b;
  int element;
  double *fem_value;
  int i;
  int j;
  int j2;
  int k;
  int ni;
  int nj;
  double *phi;
  int project_node_num = 1;
  double project_node_xyz[3*1];
  double *project_value;
  int quad;
  int quad_num = 4;
  double *ref_quad;
  double *ref_weight;
  double *tet_quad;
  double *tet_xyz;
  double volume;
  double *x;
//
//  Assemble the coefficient matrix A and the right-hand side B.
//
  b = r8mat_zero_new ( fem_node_num, fem_value_dim );
  a = r8mat_zero_new ( fem_node_num, fem_node_num );

  phi = new double[4];
  ref_weight = new double[quad_num];
  ref_quad = new double[4*quad_num];
  tet_quad = new double[3*quad_num];
  tet_xyz = new double[3*4];

  for ( element = 0; element < fem_element_num; element++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        j2 = fem_element_node[j+element*4];
        tet_xyz[i+j*3] = fem_node_xyz[i+j2*3];
      }
    }
 
    volume = tetrahedron_volume ( tet_xyz );

    for ( j = 0; j < quad_num; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        tet_quad[i+j*3] = 0.0;
        for ( k = 0; k < 4; k++ )
        {
          tet_quad[i+j*3] = tet_quad[i+j*3] + tet_xyz[i+k*3] * ref_quad[k+j*4];
        }
      }
    }
//
//  Consider each quadrature point.
//  Here, we use the midside nodes as quadrature points.
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        project_node_xyz[i+0*3] = tet_quad[i+quad*3];
      }
      basis_mn_tet4 ( tet_xyz, 1, project_node_xyz, phi );

      for ( i = 0; i < 4; i++ )
      {
        ni = fem_element_node[i+element*fem_element_order];
//
//  The projection takes place here.  The finite element code needs the value
//  of the sample function at the point (XQ,YQ).  The call to PROJECTION
//  locates (XQ,YQ) in the triangulated mesh of sample data, and returns a
//  value produced by piecewise linear interpolation.
//
        project_value = projection ( sample_node_num, sample_node_xyz, 
          sample_element_order, sample_element_num, sample_element_node, 
          sample_element_neighbor, sample_value_dim, sample_value, 
          project_node_num, project_node_xyz );

        for ( j = 0; j < fem_value_dim; j++ )
        {
          b[ni+j*fem_node_num] = b[ni+j*fem_node_num]
            + volume * ref_weight[quad] * ( project_value[j+0*fem_value_dim] * phi[i] );
        }

        delete [] project_value;
//
//  Consider each basis function in the element.
//
        for ( j = 0; j < 4; j++ )
        {
          nj = fem_element_node[j+element*fem_element_order];

          a[ni+nj*fem_node_num] = a[ni+nj*fem_node_num] 
            + volume * ref_weight[quad] * ( phi[i] * phi[j] );
        }
      }
    }
  }
//
//  SOLVE the linear system A * X = B.
//
  x = r8ge_fss_new ( fem_node_num, a, fem_value_dim, b );
//
//  Copy solution.
//
  fem_value = new double[fem_value_dim*fem_value_num];

  for ( j = 0; j < fem_value_num; j++ )
  {
    for ( i = 0; i < fem_value_dim; i++ )
    {
      fem_value[i+j*fem_value_dim] = x[j+i*fem_value_num];
    }
  }

  delete [] a;
  delete [] b;
  delete [] phi;
  delete [] ref_quad;
  delete [] ref_weight;
  delete [] tet_quad;
  delete [] tet_xyz;
  delete [] x;

  return fem_value;
}
//****************************************************************************80

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words, 
//    separated by spaces.  There may also be some blank lines, and some 
//    comment lines, which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  string text;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) <= 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input >> text;

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
      {
        continue;
      }
      got_one = true;
      break;
    }
  }

  input.close ( );

  if ( !got_one )
  {
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

int file_row_count ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  string line;
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return (-1);
  }

  for ( ; ; )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4i4i4_sort_a ( int i1, int i2, int i3, int *j1, int *j2, int *j3 )

//****************************************************************************80
//
//  Purpose:
//
//    I4I4I4_SORT_A ascending sorts a triple of I4's.
//
//  Discussion:
//
//    The program allows the reasonable call:
//
//      i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
//    and this will return the reasonable result.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, I3, the values to sort.
//
//    Output, int *J1, *J2, *J3, the sorted values.
//
{
  int k1;
  int k2;
  int k3;
//
//  Copy arguments, so that the user can make "reasonable" calls like:
//
//    i4i4i4_sort_a ( i1, i2, i3, &i1, &i2, &i3 );
//
  k1 = i1;
  k2 = i2;
  k3 = i3;

  *j1 = i4_min ( i4_min ( k1, k2 ), i4_min ( k2, k3 ) );
  *j2 = i4_min ( i4_max ( k1, k2 ), 
        i4_min ( i4_max ( k2, k3 ), i4_max ( k3, k1 ) ) );
  *j3 = i4_max ( i4_max ( k1, k2 ), i4_max ( k2, k3 ) );

  return;
}
//****************************************************************************80

int *i4mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int I4MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  int *table;
  int *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new int[m*n];

  x = new int[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_i4vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80
 
void i4mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HEADER_READ reads the header from an I4MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );
 
  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }
 
  *n = file_row_count ( input_filename );
 
  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }
 
  return;
}
//****************************************************************************80

int i4mat_min ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MIN returns the minimum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MIN, the minimum entry of A.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int value;

  value = i4_huge;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

double *projection ( int fem_node_num, double fem_node_xyz[], 
  int fem_element_order, int fem_element_num, int fem_element_node[], 
  int fem_element_neighbor[], int fem_value_dim, double fem_value[], 
  int sample_node_num, double sample_node_xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PROJECTION evaluates an FEM function on a T3 or T6 triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FEM_NODE_NUM, the number of nodes.
//
//    Input, double FEM_NODE_XYZ[3*FEM_NODE_NUM], the coordinates 
//    of the nodes.
//
//    Input, int FEM_ELEMENT_ORDER, the order of the elements.
//
//    Input, int FEM_ELEMENT_NUM, the number of elements.
//
//    Input, int FEM_ELEMENT_NODE(FEM_ELEMENT_ORDER,FEM_ELEMENT_NUM), the
//    nodes that make up each element.
//
//    Input, int FEM_ELEMENT_NEIGHBOR[4*FEM_ELEMENT_NUM], the 
//    index of the neighboring element on each side, or -1 if no neighbor there.
//
//    Input, int FEM_VALUE_DIM, the "dimension" of the values.
//
//    Input, double FEM_VALUE[FEM_VALUE_DIM*FEM_NODE_NUM], the 
//    finite element coefficient values at each node.
//
//    Input, int SAMPLE_NODE_NUM, the number of sample nodes.
//
//    Input, double SAMPLE_NODE_XYZ[3*SAMPLE_NODE_NUM], the sample nodes.
//
//    Output, double PROJECTION[FEM_VALUE_DIM*SAMPLE_NODE_NUM],
//    the sampled values.
//
{
  double *b;
  double dot;
  int face;
  int i;
  int j;
  int k;
  double p_xyz[3];
  double *sample_value;
  int step_num;
  int t;
  int t_node;
  double *t_xyz;

  b = new double[fem_element_order];
  sample_value = new double[fem_value_dim*sample_node_num];
  t_xyz = new double[3*fem_element_order];
//
//  For each sample point: find the element T that contains it,
//  and evaluate the finite element function there.
//
  for ( j = 0; j < sample_node_num; j++ )
  {
    p_xyz[0] = sample_node_xyz[0+3*j];
    p_xyz[1] = sample_node_xyz[1+3*j];
    p_xyz[2] = sample_node_xyz[2+3*j];
//
//  Find the triangle T that contains the point.
//
    t = tet_mesh_search_delaunay ( fem_node_num, fem_node_xyz, 
      fem_element_order, fem_element_num, fem_element_node, 
      fem_element_neighbor, p_xyz, &face, &step_num );

    if ( t == - 1 )
    {
      cerr << "\n";
      cerr << "PROJECTION - Fatal error!\n";
      cerr << "  Search failed.\n";
      exit ( 1 );
    }
//
//  Evaluate the finite element basis functions at the point in T.
//
    for ( i = 0; i < fem_element_order; i++ )
    {
      t_node = fem_element_node[i+t*fem_element_order];
      t_xyz[0+i*3] = fem_node_xyz[0+t_node*3];
      t_xyz[1+i*3] = fem_node_xyz[1+t_node*3];
      t_xyz[2+i*3] = fem_node_xyz[2+t_node*3];
    }

    basis_mn_tet4 ( t_xyz, 1, p_xyz, b );
//
//  Multiply by the finite element values to get the sample values.
//
    for ( i = 0; i < fem_value_dim; i++ )
    {
      dot = 0.0;
      for ( k = 0; k < fem_element_order; k++ )
      {
        t_node = fem_element_node[k+t*fem_element_order];
        dot = dot + fem_value[i+t_node*fem_value_dim] * b[k];
      }
      sample_value[i+j*fem_value_dim] = dot;
    }
  }

  delete [] b;
  delete [] t_xyz;

  return sample_value;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double *r8ge_fss_new ( int n, double a[], int nb, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FSS_NEW factors and solves multiple R8GE systems.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    This routine does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input, double B[N*NB], the right hand sides of the linear systems.
//
//    Output, double R8GE_FSS_NEW[N*NB], the solutions of the linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = new double[n*nb];

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cout << "\n";
      cout << "R8GE_FSS_NEW - Fatal error!\n";
      cout << "  Zero pivot on step " << jcol << "\n";
      return NULL;
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R8MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  double *table;
  double *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80

double r8mat_det_4d ( double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80
 
void r8mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HEADER_READ reads the header from an R8MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

int r8mat_solve ( int n, int rhs_num, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( r8_abs ( apivot ) < r8_abs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//  For greater precision, try
//
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(10) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

double *r8mat_zero_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ZERO_NEW returns a new zeroed R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

bool r8vec_is_nonnegative ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_IS_NONNEGATIVE is true if all entries in an R8VEC are nonnegative.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector to be checked.
//
//    Output, bool R8VEC_IS_NONNEGATIVE is true if all elements of X
//    are nonnegative.
//
{
  int i;

  for ( i = 0; i < n; i++ ) 
  {
    if ( x[i] < 0.0 ) 
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( string s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for ( ; ; ) 
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

bool s_to_i4vec ( string s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

int *tet_mesh_neighbor_tets ( int tetra_order, int tetra_num, 
  int tetra_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
//
//  Discussion:
//
//    A tet mesh of a set of nodes can be completely described by
//    the coordinates of the nodes, and the list of nodes that make up
//    each tetrahedron.  In the most common case, four nodes are used.
//    There is also a 10 node case, where nodes are also placed on
//    the midsides of the tetrahedral edges.
//
//    This routine can handle 4 or 10-node tetrahedral meshes.  The
//    10-node case is handled simply by ignoring the six midside nodes,
//    which are presumed to be listed after the vertices.
//
//    The tetrahedron adjacency information records which tetrahedron
//    is adjacent to a given tetrahedron on a particular face.
//
//    This routine creates a data structure recording this information.
//
//    The primary amount of work occurs in sorting a list of 4 * TETRA_NUM
//    data items.
//
//    The neighbor tetrahedrons are indexed by the face they share with
//    the tetrahedron.
//
//    Each face of the tetrahedron is indexed by the node which is NOT
//    part of the face.  That is:
//
//    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.
//    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
//    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
//    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
//
//    For instance, if the (transposed) TETRA_NODE array was:
//
//    Row       1      2      3      4
//    Col
//
//      1       4      3      5      1
//      2       4      2      5      1
//      3       4      7      3      5
//      4       4      7      8      5
//      5       4      6      2      5
//      6       4      6      8      5
//
//    then the (transposed) TETRA_NEIGHBOR array should be:
//
//    Row       1      2      3      4
//    Col
//
//      1      -1      2     -1      3
//      2      -1      1     -1      5
//      3      -1      1      4     -1
//      4      -1      6      3     -1
//      5      -1      2      6     -1
//      6      -1      4      5     -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TETRA_ORDER, the order of the tetrahedrons.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
//
//    Output, int TET_MESH_NEIGHBORS[4*TETRA_NUM], the four tetrahedrons that
//    are direct neighbors of a given tetrahedron.  If there is no neighbor
//    sharing a given face, the index is set to -1.
//
{
  int a;
  int b;
  int c;
  int face;
  int face1;
  int face2;
  int *faces;
  int i;
  int j;
  int k;
  int l;
  int tetra;
  int *tetra_neighbor;
  int tetra1;
  int tetra2;

  faces = new int[5*(4*tetra_num)];
  tetra_neighbor = new int[4*tetra_num];
//
//  Step 1.
//  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
//  construct the four face relations:
//
//    (J,K,L,1,T)
//    (I,K,L,2,T)
//    (I,J,L,3,T)
//    (I,J,K,4,T)
//
//  In order to make matching easier, we reorder each triple of nodes
//  into ascending order.
//
  for ( tetra = 0; tetra < tetra_num; tetra++ )
  {
    i = tetra_node[0+tetra*tetra_order];
    j = tetra_node[1+tetra*tetra_order];
    k = tetra_node[2+tetra*tetra_order];
    l = tetra_node[3+tetra*tetra_order];

    i4i4i4_sort_a ( j, k, l, &a, &b, &c );

    faces[0+0*5+tetra*5*4] = a;
    faces[1+0*5+tetra*5*4] = b;
    faces[2+0*5+tetra*5*4] = c;
    faces[3+0*5+tetra*5*4] = 0;
    faces[4+0*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, k, l, &a, &b, &c );

    faces[0+1*5+tetra*5*4] = a;
    faces[1+1*5+tetra*5*4] = b;
    faces[2+1*5+tetra*5*4] = c;
    faces[3+1*5+tetra*5*4] = 1;
    faces[4+1*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, j, l, &a, &b, &c );

    faces[0+2*5+tetra*5*4] = a;
    faces[1+2*5+tetra*5*4] = b;
    faces[2+2*5+tetra*5*4] = c;
    faces[3+2*5+tetra*5*4] = 2;
    faces[4+2*5+tetra*5*4] = tetra;

    i4i4i4_sort_a ( i, j, k, &a, &b, &c );

    faces[0+3*5+tetra*5*4] = a;
    faces[1+3*5+tetra*5*4] = b;
    faces[2+3*5+tetra*5*4] = c;
    faces[3+3*5+tetra*5*4] = 3;
    faces[4+3*5+tetra*5*4] = tetra;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:3; the routine we call here
//  sorts on rows 1 through 5 but that won't hurt us.
//
//  What we need is to find cases where two tetrahedrons share a face.
//  By sorting the columns of the FACES array, we will put shared faces
//  next to each other.
//
  i4col_sort_a ( 5, 4*tetra_num, faces );
//
//  Step 3. Neighboring tetrahedrons show up as consecutive columns with
//  identical first three entries.  Whenever you spot this happening,
//  make the appropriate entries in TETRA_NEIGHBOR.
//
  for ( j = 0; j < tetra_num; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      tetra_neighbor[i+j*4] = -1;
    }
  }

  face = 0;

  for ( ; ; )
  {
    if ( 4 * tetra_num - 1 <= face )
    {
      break;
    }

    if ( faces[0+face*5] == faces[0+(face+1)*5] &&
         faces[1+face*5] == faces[1+(face+1)*5] &&
         faces[2+face*5] == faces[2+(face+1)*5] )
    {
      face1 = faces[3+face*5];
      tetra1 = faces[4+face*5];
      face2 = faces[3+(face+1)*5];
      tetra2 = faces[4+(face+1)*5];
      tetra_neighbor[face1+tetra1*4] = tetra2;
      tetra_neighbor[face2+tetra2*4] = tetra1;
      face = face + 2;
    }
    else
    {
      face = face + 1;
    }
  }

  delete [] faces;

  return tetra_neighbor;
}
//****************************************************************************80

int tet_mesh_search_delaunay ( int node_num, double node_xyz[], int tet_order, 
  int tet_num, int tet_node[], int tet_neighbor[], double p[], int *face, 
  int *step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_SEARCH_DELAUNAY searches a Delaunay tet mesh for a point.
//
//  Discussion:
//
//    The algorithm "walks" from one tetrahedron to its neighboring tetrahedron,
//    and so on, until a tetrahedron is found containing point P, or P is found
//    to be outside the convex hull.
//
//    The algorithm computes the barycentric coordinates of the point with
//    respect to the current tetrahedron.  If all 4 quantities are positive,
//    the point is contained in the tetrahedron.  If the I-th coordinate is
//    negative, then P lies on the far side of edge I, which is opposite
//    from vertex I.  This gives a hint as to where to search next.
//
//    For a Delaunay tet mesh, the search is guaranteed to terminate.
//    For other meshes, a continue may occur.
//
//    Note the surprising fact that, even for a Delaunay tet mesh of
//    a set of nodes, the nearest node to P need not be one of the
//    vertices of the tetrahedron containing P.
//
//    The code can be called for tet meshes of any order, but only
//    the first 4 nodes in each tetrahedron are considered.  Thus, if
//    higher order tetrahedrons are used, and the extra nodes are intended
//    to give the tetrahedron a polygonal shape, these will have no effect,
//    and the results obtained here might be misleading.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of 
//    the nodes.
//
//    Input, int TET_ORDER, the order of the tetrahedrons.
//
//    Input, int TET_NUM, the number of tetrahedrons.
//
//    Input, int TET_NODE[TET_ORDER*TET_NUM],
//    the nodes that make up each tetrahedron.
//
//    Input, int TET_NEIGHBOR[4*TET_NUM], the 
//    tetrahedron neighbor list.
//
//    Input, double P[3], the coordinates of a point.
//
//    Output, int *FACE, indicates the position of the point P in
//    face TET_INDEX:
//    0, the interior or boundary of the tetrahedron;
//    -1, outside the convex hull of the tet mesh, past face 1;
//    -2, outside the convex hull of the tet mesh, past face 2;
//    -3, outside the convex hull of the tet mesh, past face 3.
//    -4, outside the convex hull of the tet mesh, past face 4.
//
//    Output, int *STEP_NUM, the number of steps taken.
//
//    Output, int TET_MESH_SEARCH_DELAUNAY, the index of the tetrahedron 
//    where the search ended.  If a cycle occurred, then -1 is returned.
//
{
  double *alpha;
  int i;
  int j;
  int k;
  int tet_index;
  double tet_xyz[3*4];
  static int tet_index_save = -1;
//
//  If possible, start with the previous successful value of TET_INDEX.
//
  if ( tet_index_save < 1 || tet_num < tet_index_save )
  {
    tet_index = ( tet_num + 1 ) / 2;
  }
  else
  {
    tet_index = tet_index_save;
  }

  *step_num = -1;
  *face = 0;

  for ( ; ; )
  {
    *step_num = *step_num + 1;

    if ( tet_num < *step_num )
    {
      cerr << "\n";
      cerr << "TET_MESH_SEARCH_DELAUNAY - Fatal error!\n";
      cerr << "  The algorithm seems to be cycling.\n";
      tet_index = -1;
      *face = -1;
      exit ( 1 );
    }

    for ( j = 0; j < 4; j++ )
    {
      k = tet_node[j+tet_index*4];
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+k*3];
      }
    }

    alpha = tetrahedron_barycentric ( tet_xyz, p );
//
//  If the barycentric coordinates are all positive, then the point
//  is inside the tetrahedron and we're done.
//
    if ( 0.0 <= alpha[0] && 0.0 <= alpha[1] && 0.0 <= alpha[2] && 0.0 <= alpha[3] )
    {
      break;
    }
//
//  At least one barycentric coordinate is negative.
//
//  If there is a negative barycentric coordinate for which there exists an
//  opposing tetrahedron neighbor closer to the point, move to that tetrahedron.
//
    if ( alpha[0] < 0.0 && 0 < tet_neighbor[0+tet_index*4] )
    {
      tet_index = tet_neighbor[0+tet_index*4];
      continue;
    }
    else if ( alpha[1] < 0.0 && 0 < tet_neighbor[1+tet_index*4] )
    {
      tet_index = tet_neighbor[1+tet_index*4];
      continue;
    }
    else if ( alpha[2] < 0.0 && 0 < tet_neighbor[2+tet_index*4] )
    {
      tet_index = tet_neighbor[2+tet_index*4];
      continue;
    }
    else if ( alpha[3] < 0.0 && 0 < tet_neighbor[3+tet_index*4] )
    {
      tet_index = tet_neighbor[3+tet_index*4];
      continue;
    }
//
//  All negative barycentric coordinates correspond to vertices opposite
//  faces on the convex hull.
//
//  Note the face and exit.
//
    if ( alpha[0] < 0.0 )
    {
      *face = -1;
      break;
    }
    else if ( alpha[1] < 0.0 )
    {
      *face = -2;
      break;
    }
    else if ( alpha[2] < 0.0 )
    {
      *face = -3;
      break;
    }
    else if ( alpha[3] < 0.0 )
    {
      *face = -4;
      break;
    }
  }

  tet_index_save = tet_index;

  return tet_index;
}
//****************************************************************************80

int tet_mesh_search_naive ( int node_num, double node_xyz[],
  int tet_order, int tet_num, int tet_node[], double p[], int *step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_SEARCH_NAIVE naively searches a tet mesh.
//
//  Discussion:
//
//    The algorithm simply checks each tetrahedron to see if point P is
//    contained in it.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the coordinates 
//    of the nodes.
//
//    Input, int TET_ORDER, the order of the tetrahedrons.
//
//    Input, int TET_NUM, the number of tetrahedrons in
//    the mesh.
//
//    Input, int TET_NODE[TET_ORDER*TET_NUM], 
//    the nodes that make up each tetrahedron.
//
//    Input, double P[3], the coordinates of a point.
//
//    Output, int TET_MESH_ORDER4_SEARCH_NAIE, the index of the tetrahedron
//    where the search ended, or -1 if no tetrahedron was found containing
//    the point.
//
//    Output, int *STEP_NUM, the number of tetrahedrons examined.
{
  double *alpha;
  int i;
  int j;
  int tet;
  int tet_index;
  double tet_xyz[3*4];

  tet_index = -1;
  *step_num = 0;

  for ( tet = 0; tet < tet_num; tet++ )
  {
    for ( j = 0; j < 4; j++ )
    {
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+tet_node[j+tet*4]*3];
      }
    }
    alpha = tetrahedron_barycentric ( tet_xyz, p );

    if ( r8vec_is_nonnegative ( 4, alpha ) )
    {
      tet_index = tet;
      *step_num = tet;
      return tet_index;
    }

    delete [] alpha;
  }

  return tet_index;
}
//****************************************************************************80

double *tetrahedron_barycentric ( double tetra[3*4], double p[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_BARYCENTRIC returns the barycentric coordinates of a point.
//
//  Discussion:
//
//    The barycentric coordinates of a point P with respect to
//    a tetrahedron are a set of four values C(1:4), each associated
//    with a vertex of the tetrahedron.  The values must sum to 1.
//    If all the values are between 0 and 1, the point is contained
//    within the tetrahedron.
//
//    The barycentric coordinate of point X related to vertex A can be
//    interpreted as the ratio of the volume of the tetrahedron with 
//    vertex A replaced by vertex X to the volume of the original 
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, double C[4], the barycentric coordinates of the point with
//    respect to the tetrahedron.
//
{
# define N 3
# define RHS_NUM 1

  double a[N*(N+RHS_NUM)];
  double *c;
  int info;
//
//  Set up the linear system
//
//    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
//    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
//    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
//
//  which is satisfied by the barycentric coordinates.
//

  a[0+0*N] = tetra[0+1*3] - tetra[0+0*3];
  a[1+0*N] = tetra[1+1*3] - tetra[1+0*3];
  a[2+0*N] = tetra[2+1*3] - tetra[2+0*3];

  a[0+1*N] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*N] = tetra[1+2*3] - tetra[1+0*3];
  a[2+1*N] = tetra[2+2*3] - tetra[2+0*3];

  a[0+2*N] = tetra[0+3*3] - tetra[0+0*3];
  a[1+2*N] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*N] = tetra[2+3*3] - tetra[2+0*3];

  a[0+3*N] = p[0]         - tetra[0+0*3];
  a[1+3*N] = p[1]         - tetra[1+0*3];
  a[2+3*N] = p[2]         - tetra[2+0*3];
//
//  Solve the linear system.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TETRAHEDRON_BARYCENTRIC - Fatal error!\n";
    cout << "  The linear system is singular.\n";
    cout << "  The input data does not form a proper tetrahedron.\n";
    exit ( 1 );
  }

  c = new double[4];

  c[1] = a[0+3*N];
  c[2] = a[1+3*N];
  c[3] = a[2+3*N];

  c[0] = 1.0 - c[1] - c[2] - c[3];

  return c;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

double tetrahedron_volume ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
//
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    { 
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = fabs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
