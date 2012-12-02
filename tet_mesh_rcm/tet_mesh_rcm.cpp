# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int adj_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[] );
int adj_perm_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[], 
  int perm[], int perm_inv[] );
void adj_print ( int node_num, int adj_num, int adj_row[], int adj[], 
  string title );
void adj_print_some ( int node_num, int node_lo, int node_hi, int adj_num, 
  int adj_row[], int adj[], string title );
char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void degree ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int deg[], int *iccsze, int ls[], int node_num );
int file_column_count ( string input_filename );
bool file_exist ( string filename );
int file_row_count ( string input_filename );
int *genrcm ( int node_num, int adj_num, int adj_row[], int adj[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_swap ( int *i, int *j );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_sort2_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_reverse ( int n, int a[] );
void level_set ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num );
bool perm_check ( int n, int p[], int base );
int *perm_inverse3 ( int n, int perm[] );
void r8col_permute ( int m, int n, int p[], int base, double a[] );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void rcm ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int perm[], int *iccsze, int node_num );
void root_find ( int *root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
int tet_mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] );
void tet_mesh_order4_adj_count ( int node_num, int element_num, 
  int element_node[], int *adj_num, int adj_row[] );
int *tet_mesh_order4_adj_set ( int node_num, int element_num, 
  int element_node[], int adj_num, int adj_row[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TET_MESH_RCM.
//
//  Discussion:
//
//    TET_MESH_RCM applies the RCM reordering to a tet mesh.
//
//    The user supplies a node file and a tetrahedron file, containing
//    the coordinates of the nodes, and the indices of the nodes that
//    make up each tetrahedron.  Either 4-node or 10-node tetrahedrons may
//    be used.
//
//    The program reads the data, computes the adjacency information,
//    carries out the RCM algorithm to get the permutation, applies
//    the permutation to the nodes and tetrahedrons, and writes out
//    new node and tetrahedron files that correspond to the RCM permutation.
//
//    Note that node data is normally three dimensional, that is,
//    each node has an X, Y and Z coordinate.  In some applications, it
//    may be desirable to specify more information.  This program
//    will accept node data that includes DIM_NUM entries on each line,
//    as long as DIM_NUM is the same for each entry.  
//
//    Thanks to Xingxing Zhang for pointing out some problems with a
//    previous version of this program, 10 May 2011.
//
//  Usage:
//
//    tet_mesh_rcm prefix
//
//    where prefix is the common file prefix:
//
//    * prefix_nodes.txt,           the node coordinates (input);
//    * prefix_elements.txt,        the element definitions (input).
//    * prefix_rcm_nodes.txt,       the new node coordinates (output);
//    * prefix_rcm_elements.txt,    the new element definitions (output).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *adj;
  int adj_num;
  int *adj_row;
  int bandwidth;
  int base_internal;
  int base_user;
  bool debug = false;
  int dim_num;
  int element;
  int *element_node;
  int element_num;
  int element_order;
  int i;
  string node_filename;
  string element_filename;
  int j;
  int node;
  int node_max;
  int node_min;
  int node_num;
  double *node_xyz;
  int order;
  string node_rcm_filename;
  string element_rcm_filename;
  int *perm;
  int *perm_inv;
  string prefix;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TET_MESH_RCM\n";
  cout << "  C++ version\n";
  cout << "  Read a node dataset of NODE_NUM points in 3 dimensions.\n";
  cout << "  Read an associated tet mesh dataset of TETRA_NUM\n";
  cout << "  tetrahedrons using 4 or 10 nodes.\n";
  cout << "\n";
  cout << "  Apply the RCM reordering (Reverse Cuthill-McKee).\n";
  cout << "\n";
  cout << "  Reorder the data and write it out to files.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  Get the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TET_MESH_RCM:\n";
    cout << "  Please enter the filename prefix.\n";

    cin >> prefix;
  }
  else 
  {
    prefix = argv[1];
  }
//
//  Create the filenames.
//
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";
  node_rcm_filename = prefix + "_rcm_nodes.txt";
  element_rcm_filename = prefix + "_rcm_elements.txt";
//
//  Read the node data.
//
  r8mat_header_read ( node_filename, &dim_num, &node_num );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << dim_num << "\n";
  cout << "  Number of points NODE_NUM  = " << node_num << "\n";

  node_xyz = r8mat_data_read ( node_filename, dim_num, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print_some ( dim_num, node_num, node_xyz, 1, 1, dim_num, 5, 
    "  Coordinates of first 5 nodes:" );
//
//  Read the tet mesh data.
//
  i4mat_header_read ( element_filename, &element_order, &element_num );

  if ( element_order != 4 && element_order != 10 )
  {
    cerr << "\n";
    cerr << "TET_MESH_RCM - Fatal error!\n";
    cerr << "  The tet mesh must have order 4 or order 10.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Read the header of \"" << element_filename << "\".\n";
  cout << "\n";
  cout << "  Tetrahedron order = " << element_order << "\n";
  cout << "  Number of tetras  = " << element_num << "\n";

  element_node = i4mat_data_read ( element_filename, element_order, 
    element_num );

  cout << "\n";
  cout << "  Read the data in \"" << element_filename << "\".\n";

  i4mat_transpose_print_some ( element_order, element_num, 
    element_node, 1, 1, element_order, 5, "  First 5 tetrahedrons:" );
//
//  If the element information is 1-based, make it 0-based.
//
  base_user = tet_mesh_base_zero ( node_num, element_order, element_num, element_node );

  if ( base_user != 0 && base_user != 1 )
  {
    cerr << "\n";
    cerr << "TET_MESH_RCM - Fatal error!\n";
    cerr << "  The input data does not seem to be 0-based or 1-based.\n";
    exit ( 1 );
  }
//
//  Following code depends on the element order.
//
  if ( element_order == 4 )
  {
//
//  Count the number of adjacencies.
//  Set up the ADJ_ROW adjacency pointer array.
//
    adj_row = new int[node_num+1];

    tet_mesh_order4_adj_count ( node_num, element_num, element_node, 
      &adj_num, adj_row );

    if ( debug || node_num < 10 )
    {
      cout << "\n";
      cout << "  ADJ_NUM = " << adj_num << "\n";

      i4vec_print ( node_num + 1, adj_row, "  ADJ_ROW:" );
    }
//
//  Set up the ADJ adjacency array.
//
    adj = tet_mesh_order4_adj_set ( node_num, element_num, element_node,
      adj_num, adj_row );

    if ( node_num < 10 )
    {
      adj_print ( node_num, adj_num, adj_row, adj, "  ADJ" );
    }
  }
  else if ( element_order == 10 )
  {
    cerr << "\n";
    cerr << "TET_MESH_RCM - Fatal error!\n";
    cerr << "  This order 10 code is not ready.\n";
    exit ( 1 );
  }
//
//  Compute the bandwidth.
//
  bandwidth = adj_bandwidth ( node_num, adj_num, adj_row, adj );

  cout << "\n";
  cout << "  ADJ bandwidth = " << bandwidth << "\n";
//
//  GENRCM computes the RCM permutation.
//
  perm = genrcm ( node_num, adj_num, adj_row, adj );
//
//  Compute the inverse permutation.
//
  perm_inv = perm_inverse3 ( node_num, perm );

  if ( node_num < 10 )
  {
    cout << "\n";
    cout << "         I   PERM[I] INVERSE[I]\n";
    cout << "\n";
    for ( i = 0; i < node_num; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(8) << perm[i]
           << "  " << setw(8) << perm_inv[i] << "\n";
    }
  }
//
//  Compute the bandwidth of the permuted array.
//
  bandwidth = adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, 
    perm, perm_inv );

  cout << "\n";
  cout << "  ADJ bandwidth after RCM permutation = " << bandwidth << "\n";
//
//  Permute the nodes in NODE_XYZ.
//
  base_internal = 0;
  r8col_permute ( dim_num, node_num, perm, base_internal, node_xyz );
//
//  Permute the node indices in ELEMENT_NODE.
//
  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < element_order; i++ )
    {
      node = element_node[i+j*element_order];
      element_node[i+j*element_order] = perm_inv[node];
    }
  }
//
//  If the user base was 1, restore it!
//
  if ( base_user == 1 )
  {
    for ( i = 0; i < element_num; i++ )
    {
      for ( j = 0; j < element_order; j++ )
      {
        element_node[i+j*element_order] = element_node[i+j*element_order] + 1;
      }
    }
    cout << "\n";
    cout << "  Output files will use the same 1-based ordering used by the input.\n";
  }
//
//  Write the node and element data.
//
  r8mat_write ( node_rcm_filename, dim_num, node_num, node_xyz );

  cout << "\n";
  cout << "  Created the file \"" << node_rcm_filename << "\".\n";

  i4mat_write ( element_rcm_filename, element_order, element_num, element_node );
    
  cout << "  Created the file \"" << element_rcm_filename << "\".\n";
//
//  Deallocate memory.
//
  delete [] adj;
  delete [] adj_row;
  delete [] node_xyz;
  delete [] perm;
  delete [] perm_inv;
  delete [] element_node;
//
//  Terminate.
//
  cout << "\n";
  cout << "TET_MESH_RCM:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int adj_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int ADJ_BANDWIDTH, the bandwidth of the adjacency
//    matrix.
//
{
  int band_hi;
  int band_lo;
  int col;
  int i;
  int j;
  int value;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j];
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  value = band_lo + 1 + band_hi;

  return value;
}
//****************************************************************************80

int adj_perm_bandwidth ( int node_num, int adj_num, int adj_row[], int adj[], 
  int perm[], int perm_inv[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.  
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
//    and inverse permutation.
//
//    Output, int ADJ_PERM_BANDWIDTH, the bandwidth of the permuted 
//    adjacency matrix.
//
{
  int band_hi;
  int band_lo;
  int bandwidth;
  int col;
  int i;
  int j;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[perm[i]]; j <= adj_row[perm[i]+1]-1; j++ )
    {
      col = perm_inv[adj[j]];
      band_lo = i4_max ( band_lo, i - col );
      band_hi = i4_max ( band_hi, col - i );
    }
  }

  bandwidth = band_lo + 1 + band_hi;

  return bandwidth;
}
//****************************************************************************80

void adj_print ( int node_num, int adj_num, int adj_row[], int adj[], 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PRINT prints adjacency information.
//
//  Discussion:
//
//    The list has the form:
//
//    Row   Nonzeros
//
//    1       2   5   9
//    2       7   8   9   15   78   79   81  86  91  99
//          100 103
//    3      48  49  53
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 December 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
//    into rows.  The entries for row I are in entries ADJ_ROW(I)
//    through ADJ_ROW(I+1)-1.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
//    for each row, the column indices of the nonzero entries.
//
//    Input, string TITLE, a title.
//
{
  adj_print_some ( node_num, 0, node_num - 1, adj_num, adj_row, adj, title );

  return;
}
//****************************************************************************80

void adj_print_some ( int node_num, int node_lo, int node_hi, int adj_num, 
  int adj_row[], int adj[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PRINT_SOME prints some adjacency information.
//
//  Discussion:
//
//    The list has the form:
//
//    Row   Nonzeros
//
//    1       2   5   9
//    2       7   8   9   15   78   79   81  86  91  99
//          100 103
//    3      48  49  53
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_LO, NODE_HI, the first and last nodes for
//    which the adjacency information is to be printed.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1], organizes the adjacency entries
//    into rows.  The entries for row I are in entries ADJ_ROW(I)
//    through ADJ_ROW(I+1)-1.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure, which contains,
//    for each row, the column indices of the nonzero entries.
//
//    Input, string TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;
  int jmax;
  int jmin;

  cout << "\n";
  cout << title << "\n";
  cout << "  Sparse adjacency structure:\n";
  cout << "\n";
  cout << "  Number of nodes       = " << node_num << "\n";;
  cout << "  Number of adjacencies = " << adj_num << "\n";
  cout << "\n";
  cout << "  Node   Min   Max          Nonzeros \n";
  cout << "\n";

  for ( i = node_lo; i <= node_hi; i++ )
  {
    jmin = adj_row[i];
    jmax = adj_row[i+1] - 1;

    if ( jmax < jmin )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(4) << jmin
           << "  " << setw(4) << jmax << "\n";
    }
    else
    {
      for ( jlo = jmin; jlo <= jmax; jlo = jlo + 5 )
      {
        jhi = i4_min ( jlo + 4, jmax );

        if ( jlo == jmin )
        {
          cout << "  " << setw(4) << i
               << "  " << setw(4) << jmin
               << "  " << setw(4) << jmax
               << "   ";
          for ( j = jlo; j <= jhi; j++ )
          {
            cout << setw(8) << adj[j];
          }
          cout << "\n";
        }
        else
        {
          cout << "                     ";
          for ( j = jlo; j <= jhi; j++ )
          {
            cout << setw(8) << adj[j];
          }
          cout << "\n";
        }
      }
    }
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
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the character was
//    'illegal', then DIGIT is -1.
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

void degree ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int deg[], int *iccsze, int ls[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    DEGREE computes the degrees of the nodes in the connected component.
//
//  Discussion:
//
//    The connected component is specified by MASK and ROOT.
//    Nodes for which MASK is zero are ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], is nonzero for those nodes which are
//    to be considered.
//
//    Output, int DEG[NODE_NUM], contains, for each  node in the connected
//    component, its degree.
//
//    Output, int *ICCSIZE, the number of nodes in the connected component.
//
//    Output, int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
//    in the connected component, starting with ROOT, and proceeding 
//    by levels.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int i;
  int ideg;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;
//
//  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
//
  ls[0] = root;
  adj_row[root-1] = - adj_row[root-1];
  lvlend = 0;
  *iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
//
//  Find the degrees of nodes in the current level,
//  and at the same time, generate the next level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = ls[i-1];
      jstrt = - adj_row[node-1];
      jstop = abs ( adj_row[node] ) - 1;
      ideg = 0;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          ideg = ideg + 1;

          if ( 0 <= adj_row[nbr-1] )
          {
            adj_row[nbr-1] = -adj_row[nbr-1];
            *iccsze = *iccsze + 1;
            ls[*iccsze-1] = nbr;
          }
        }
      }
      deg[node-1] = ideg;
    }
//
//  Compute the current level width.
//
    lvsize = *iccsze - lvlend;
//
//  If the current level width is nonzero, generate another level.
//
    if ( lvsize == 0 )
    {
      break;
    }
  }
//
//  Reset ADJ_ROW to its correct sign and return.
//
  for ( i = 0; i < *iccsze; i++ )
  {
    node = ls[i] - 1;
    adj_row[node] = -adj_row[node];
  }

  return;
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

bool file_exist ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_EXIST reports whether a file exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;
  bool value;

  file.open ( filename.c_str ( ), ios::in );

  if ( !file )
  {
    value = false;
  }
  else
  {
    value = true;
  }
  return value;
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

int *genrcm ( int node_num, int adj_num, int adj_row[], int adj[] )

//****************************************************************************80
//
//  Purpose:
//
//    GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
//
//  Discussion:
//
//    For each connected component in the graph, the routine obtains
//    an ordering by calling RCM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int  ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int GENRCM[NODE_NUM], the RCM ordering.
//
//  Local Parameters:
//
//    Local, int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
//    structure.  The level structure is stored in the currently unused
//    spaces in the permutation vector PERM.
//
//    Local, int MASK[NODE_NUM], marks variables that have been numbered.
//
{
  int i;
  int iccsze;
  int level_num;
  int *level_row;
  int *mask;
  int num;
  int *perm;
  int root;
//
//  Assuming the input dat is 0 based, add 1 to ADJ_ROW and ADJ,
//  because GENRCM uses 1-based indexing!
//
  for ( i = 0; i < node_num + 1; i++ )
  {
    adj_row[i] = adj_row[i] + 1;
  }
  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = adj[i] + 1;
  }

  perm = new int[node_num];
  level_row = new int[node_num+1];
  mask = new int[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    mask[i] = 1;
  }

  num = 1;

  for ( i = 0; i < node_num; i++ )
  {
//
//  For each masked connected component...
//
    if ( mask[i] != 0 )
    {
      root = i + 1;
//
//  Find a pseudo-peripheral node ROOT.  The level structure found by
//  ROOT_FIND is stored starting at PERM(NUM).
//
      root_find ( &root, adj_num, adj_row, adj, mask, &level_num,
        level_row, perm+num-1, node_num );
//
//  RCM orders the component using ROOT as the starting node.
//
      rcm ( root, adj_num, adj_row, adj, mask, perm+num-1, &iccsze,
        node_num );

      num = num + iccsze;
    }
//
//  We can stop once every node is in one of the connected components.
//
    if ( node_num < num )
    {
      break;
    }
  }

  delete [] level_row;
  delete [] mask;
//
//  PERM is computed as a 1-based vector.
//  Rewrite it as a 0-based vector.
//
  for ( i = 0; i < node_num; i++ )
  {
    perm[i] = perm[i] - 1;
  }
//
//  Subtract 1 from ADJ_ROW and ADJ because GENRCM used 1-based indexing!
//
  for ( i = 0; i < node_num + 1; i++ )
  {
    adj_row[i] = adj_row[i] - 1;
  }
  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = adj[i] - 1;
  }
  return perm;
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

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;
 
  return;
}
//****************************************************************************80

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
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
    cerr << "\n";
    cerr << "I4COL_COMPARE - Fatal error!\n";
    cerr << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cerr << "\n";
    cerr << "I4COL_COMPARE - Fatal error!\n";
    cerr << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cerr << "\n";
    cerr << "I4COL_COMPARE - Fatal error!\n";
    cerr << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cerr << "\n";
    cerr << "I4COL_COMPARE - Fatal error!\n";
    cerr << "  N = " << n << " is less than column index J = " << j << ".\n";
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
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
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

void i4col_sort2_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A, and the length
//    of a vector of data.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors.
//    On output, the elements of each column of A have been sorted in ascending
//    order.
//
{
  int col;
  int i;
  int indx;
  int isgn;
  int j;
  int row;
  int temp;

  if ( m <= 1 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }
//
//  Initialize.
//
  for ( col = 0; col < n; col++ )
  {
    i = 0;
    indx = 0;
    isgn = 0;
    j = 0;
//
//  Call the external heap sorter.
//
    for ( ; ; )
    {
      sort_heap_external ( m, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
      if ( 0 < indx )
      {
        temp       = a[i-1+col*m];
        a[i-1+col*m] = a[j-1+col*m];
        a[j-1+col*m] = temp;
      }
//
//  Compare the I and J objects.
//
      else if ( indx < 0 )
      {
        if ( a[j-1+col*m] < a[i-1+col*m] )
        {
          isgn = +1;
        }
        else
        {
          isgn = -1;
        }
      }
      else if ( indx == 0 )
      {
        break;
      }
    }
  }

  return;
}
//****************************************************************************80

int i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
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
//    An I4COL is an M by N array of integer values, regarded
//    as an array of N columns of length M.
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
    cerr << "\n";
    cerr << "I4COL_SWAP - Fatal error!\n";
    cerr << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cerr << "\n";
    cerr << "I4COL_SWAP - Fatal error!\n";
    cerr << "  ICOL2 is out of range.\n";
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

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title for the matrix.
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  return;
# undef INCX
}
//****************************************************************************80

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 June 2009
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
//    Input, int TABLE[M*N], the table data.
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
    cerr << "I4MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
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

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void i4vec_reverse ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A(N), the array to be reversed.
//
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}
//****************************************************************************80

void level_set ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_SET generates the connected level structure rooted at a given node.
//
//  Discussion:
//
//    Only nodes for which MASK is nonzero will be considered.
//
//    The root node chosen by the user is assigned level 1, and masked.
//    All (unmasked) nodes reachable from a node in level 1 are
//    assigned level 2 and masked.  The process continues until there
//    are no unmasked nodes adjacent to any node in the current level.
//    The number of levels may vary between 2 and NODE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node at which the level structure
//    is to be rooted.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK[NODE_NUM].  On input, only nodes with nonzero
//    MASK are to be processed.  On output, those nodes which were included
//    in the level set have MASK set to 1.
//
//    Output, int *LEVEL_NUM, the number of levels in the level
//    structure.  ROOT is in level 1.  The neighbors of ROOT
//    are in level 2, and so on.
//
//    Output, int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted 
//    level structure.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int i;
  int iccsze;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  mask[root-1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
//
//  LBEGIN is the pointer to the beginning of the current level, and
//  LVLEND points to the end of this level.
//
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num-1] = lbegin;
//
//  Generate the next level by finding all the masked neighbors of nodes
//  in the current level.
//
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = level[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          iccsze = iccsze + 1;
          level[iccsze-1] = nbr;
          mask[nbr-1] = 0;
        }
      }
    }
//
//  Compute the current level width (the number of nodes encountered.)
//  If it is positive, generate the next level.
//
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
    {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
//
//  Reset MASK to 1 for the nodes in the level structure.
//
  for ( i = 0; i < iccsze; i++ )
  {
    mask[level[i]-1] = 1;
  }

  return;
}
//****************************************************************************80

bool perm_check ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      cerr << "\n";
      cerr << "PERM_CHECK - Fatal error!\n";
      cerr << "  Did not find " << found << "\n";
      exit ( 1 );
    }

  }

  return true;
}
//****************************************************************************80

int *perm_inverse3 ( int n, int perm[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE3 produces the inverse of a given permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items permuted.
//
//    Input, int PERM[N], a permutation.
//
//    Output, int PERM_INVERSE3[N], the inverse permutation.
//
{
  int i;
  int *perm_inv;

  perm_inv = new int[n];

  for ( i = 0; i < n; i++ )
  {
    perm_inv[perm[i]] = i;
  }

  return perm_inv;
}
//****************************************************************************80

void r8col_permute ( int m, int n, int p[], int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_PERMUTE permutes an R8COL in place.
//
//  Discussion:
//
//    An R8COL is an M by N array of R8's, regarded as an array of N columns,
//    each of length M.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      M = 2
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//      BASE = 1
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the length of objects.
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a
//    1-based permutation.
//
//    Input/output, double A[M*N], the array to be permuted.
//
{
  double *a_temp;
  int i;
  int iget;
  int iput;
  int istart;
  int j;

  if ( !perm_check ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "R8COL_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }

  a_temp = new double[m];
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        a_temp[i] = a[i+(istart-1)*m];
      }
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R8COL_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          for ( i = 0; i < m; i++ )
          {
            a[i+(iput-1)*m] = a_temp[i];
          }
          break;
        }
        for ( i = 0; i < m; i++ )
        {
          a[i+(iput-1)*m] = a[i+(iget-1)*m];
        }
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( j = 0; j < n; j++ )
  {
    p[j] = - p[j];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }

  delete [] a_temp;

  return;
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

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
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

void rcm ( int root, int adj_num, int adj_row[], int adj[], int mask[], 
  int perm[], int *iccsze, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
//
//  Discussion:
//
//    The connected component is specified by a node ROOT and a mask.
//    The numbering starts at the root node.
//
//    An outline of the algorithm is as follows:
//
//    X(1) = ROOT.
//
//    for ( I = 1 to N-1)
//      Find all unlabeled neighbors of X(I),
//      assign them the next available labels, in order of increasing degree.
//
//    When done, reverse the ordering.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//    It is used as the starting point for the RCM ordering.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW(NODE_NUM+1).  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ(ADJ_NUM), the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK(NODE_NUM), a mask for the nodes.  Only 
//    those nodes with nonzero input mask values are considered by the 
//    routine.  The nodes numbered by RCM will have their mask values 
//    set to zero.
//
//    Output, int PERM(NODE_NUM), the RCM ordering.
//
//    Output, int ICCSZE, the size of the connected component
//    that has been numbered.
//
//    Input, int NODE_NUM, the number of nodes.
//
//  Local Parameters:
//
//    Workspace, int DEG[NODE_NUM], a temporary vector used to hold 
//    the degree of the nodes in the section graph specified by mask and root.
//
{
  int *deg;
  int fnbr;
  int i;
  int j;
  int jstop;
  int jstrt;
  int k;
  int l;
  int lbegin;
  int lnbr;
  int lperm;
  int lvlend;
  int nbr;
  int node;
//
//  Find the degrees of the nodes in the component specified by MASK and ROOT.
//
  deg = new int[node_num];

  degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num );

  mask[root-1] = 0;

  if ( *iccsze <= 1 )
  {
    delete [] deg;
    return;
  }

  lvlend = 0;
  lnbr = 1;
//
//  LBEGIN and LVLEND point to the beginning and
//  the end of the current level respectively.
//
  while ( lvlend < lnbr )
  {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for ( i = lbegin; i <= lvlend; i++ )
    {
//
//  For each node in the current level...
//
      node = perm[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;
//
//  Find the unnumbered neighbors of NODE.
//
//  FNBR and LNBR point to the first and last neighbors
//  of the current node in PERM.
//
      fnbr = lnbr + 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          lnbr = lnbr + 1;
          mask[nbr-1] = 0;
          perm[lnbr-1] = nbr;
        }
      }
//
//  If no neighbors, skip to next node in this level.
//
      if ( lnbr <= fnbr )
      {
        continue;
      }
//
//  Sort the neighbors of NODE in increasing order by degree.
//  Linear insertion is used.
//
      k = fnbr;

      while ( k < lnbr )
      {
        l = k;
        k = k + 1;
        nbr = perm[k-1];

        while ( fnbr < l )
        {
          lperm = perm[l-1];

          if ( deg[lperm-1] <= deg[nbr-1] )
          {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }
//
//  We now have the Cuthill-McKee ordering.  Reverse it.
//
  i4vec_reverse ( *iccsze, perm );

  delete [] deg;

  return;
}
//****************************************************************************80

void root_find ( int *root, int adj_num, int adj_row[], int adj[], int mask[], 
  int *level_num, int level_row[], int level[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    ROOT_FIND finds a pseudo-peripheral node.
//
//  Discussion:
//
//    The diameter of a graph is the maximum distance (number of edges)
//    between any two nodes of the graph.
//
//    The eccentricity of a node is the maximum distance between that
//    node and any other node of the graph.
//
//    A peripheral node is a node whose eccentricity equals the
//    diameter of the graph.
//
//    A pseudo-peripheral node is an approximation to a peripheral node;
//    it may be a peripheral node, but all we know is that we tried our
//    best.
//
//    The routine is given a graph, and seeks pseudo-peripheral nodes,
//    using a modified version of the scheme of Gibbs, Poole and
//    Stockmeyer.  It determines such a node for the section subgraph
//    specified by MASK and ROOT.
//
//    The routine also determines the level structure associated with
//    the given pseudo-peripheral node; that is, how far each node
//    is from the pseudo-peripheral node.  The level structure is
//    returned as a list of nodes LS, and pointers to the beginning
//    of the list of nodes that are at a distance of 0, 1, 2, ...,
//    NODE_NUM-1 from the pseudo-peripheral node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//    Norman Gibbs, William Poole, Paul Stockmeyer,
//    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
//    SIAM Journal on Numerical Analysis,
//    Volume 13, pages 236-250, 1976.
//
//    Norman Gibbs,
//    Algorithm 509: A Hybrid Profile Reduction Algorithm,
//    ACM Transactions on Mathematical Software,
//    Volume 2, pages 378-387, 1976.
//
//  Parameters:
//
//    Input/output, int *ROOT.  On input, ROOT is a node in the
//    the component of the graph for which a pseudo-peripheral node is
//    sought.  On output, ROOT is the pseudo-peripheral node obtained.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], specifies a section subgraph.  Nodes 
//    for which MASK is zero are ignored by FNROOT.
//
//    Output, int *LEVEL_NUM, is the number of levels in the level structure
//    rooted at the node ROOT.
//
//    Output, int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
//    level structure array pair containing the level structure found.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
  int iccsze;
  int j;
  int jstrt;
  int k;
  int kstop;
  int kstrt;
  int level_num2;
  int mindeg;
  int nabor;
  int ndeg;
  int node;
//
//  Determine the level structure rooted at ROOT.
//
  level_set ( *root, adj_num, adj_row, adj, mask, level_num, 
    level_row, level, node_num );
//
//  Count the number of nodes in this level structure.
//
  iccsze = level_row[*level_num] - 1;
//
//  Extreme case:
//    A complete graph has a level set of only a single level.
//    Every node is equally good (or bad).
//
  if ( *level_num == 1 )
  {
    return;
  }
//
//  Extreme case:
//    A "line graph" 0--0--0--0--0 has every node in its only level.
//    By chance, we've stumbled on the ideal root.
//
  if ( *level_num == iccsze )
  {
    return;
  }
//
//  Pick any node from the last level that has minimum degree
//  as the starting point to generate a new level set.
//
  for ( ; ; )
  {
    mindeg = iccsze;

    jstrt = level_row[*level_num-1];
    *root = level[jstrt-1];

    if ( jstrt < iccsze )
    {
      for ( j = jstrt; j <= iccsze; j++ )
      {
        node = level[j-1];
        ndeg = 0;
        kstrt = adj_row[node-1];
        kstop = adj_row[node] - 1;

        for ( k = kstrt; k <= kstop; k++ )
        {
          nabor = adj[k-1];
          if ( 0 < mask[nabor-1] )
          {
            ndeg = ndeg + 1;
          }
        }

        if ( ndeg < mindeg )
        {
          *root = node;
          mindeg = ndeg;
        }
      }
    }
//
//  Generate the rooted level structure associated with this node.
//
    level_set ( *root, adj_num, adj_row, adj, mask, &level_num2,
      level_row, level, node_num );
//
//  If the number of levels did not increase, accept the new ROOT.
//
    if ( level_num2 <= *level_num )
    {
      break;
    }

    *level_num = level_num2;
//
//  In the unlikely case that ROOT is one endpoint of a line graph,
//  we can exit now.
//
    if ( iccsze <= *level_num )
    {
      break;
    }
  }

  return;
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
//    C++ version by John Burkardt.
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

int tet_mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
//    Output, int TET_MESH_BASE_ZERO, is:
//    0, if the input mesh was based at 0.
//    1, if the input mesh was based at 1.
//   -1, otherwise.
//
{
  int base;
  int element;
  int node;
  int node_max;
  int node_min;
  int order;
//
//  Compute the minimum and maximum node indices encountered.
//
  node_min = node_num + 1;
  node_max = -1;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      node_min = i4_min ( node_min, node );
      node_max = i4_max ( node_max, node );
    }
  }
//
//  Try to determine the base.
//
  if ( node_min == 1 && node_max == node_num )
  {
    base = 1;
    cout << "\n";
    cout << "TET_MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else if ( node_min == 0 && node_max == node_num - 1 )
  {
    base = 0;
    cout << "\n";
    cout << "TET_MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    base = -1;
    cout << "\n";
    cout << "TET_MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
  }
  return base;
}
//****************************************************************************80

void tet_mesh_order4_adj_count ( int node_num, int element_num, 
  int element_node[], int *adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_ADJ_COUNT counts the number of nodal adjacencies.
//
//  Discussion:
//
//    Assuming that the tet mesh is to be used in a finite element
//    computation, we declare that two distinct nodes are "adjacent" if and
//    only if they are both included in some tetrahedron.
//
//    It is the purpose of this routine to determine the number of
//    such adjacency relationships.
//
//    The initial count gets only the (I,J) relationships, for which
//    node I is strictly less than node J.  This value is doubled
//    to account for symmetry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Output, int *ADJ_NUM, the total number of adjacency relationships,
//
//    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
{
  int i;
  int j;
  int k;
  int node;
  int *pair;
  int pair_num;
  int pair_unique_num;
  int tetra;
//
//  Each order 4 tetrahedron defines 6 adjacency pairs.
//
  pair = new int[2*6*element_num];

  for ( tetra = 0; tetra < element_num; tetra++ )
  {
    pair[0+             tetra *2] = element_node[0+tetra*4];
    pair[1+             tetra *2] = element_node[1+tetra*4];

    pair[0+(  element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(  element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(2*element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(2*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(3*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(3*element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(4*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(4*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(5*element_num+tetra)*2] = element_node[2+tetra*4];
    pair[1+(5*element_num+tetra)*2] = element_node[3+tetra*4];
  }
  pair_num = 6 * element_num;
//
//  Force the nodes of each pair to be listed in ascending order.
//
  i4col_sort2_a ( 2, pair_num, pair );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Get the number of unique columns.
//
  pair_unique_num = i4col_sorted_unique_count ( 2, pair_num, pair );
//
//  The number of adjacencies is TWICE this value, plus the number of nodes.
//
  *adj_num = 2 * pair_unique_num;
//
//  Now set up the ADJ_ROW counts.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_row[node] = 0;
  }

  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj_row[i] = adj_row[i] + 1;
    adj_row[j] = adj_row[j] + 1;
  }
//
//  We used ADJ_ROW to count the number of entries in each row.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num-1; 0 <= node; node-- )
  {
    adj_row[node+1] = adj_row[node];
  }

  adj_row[0] = 0;
  for ( node = 1; node <= node_num; node++ )
  {
    adj_row[node] = adj_row[node-1] + adj_row[node];
  }

  delete [] pair;

  return;
}
//****************************************************************************80

int *tet_mesh_order4_adj_set ( int node_num, int element_num, 
  int element_node[], int adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    TET_MESH_ORDER4_ADJ_SET sets the nodal adjacency matrix.
//
//  Discussion:
//
//    A compressed format is used for the nodal adjacency matrix.
//
//    It is assumed that we know ADJ_NUM, the number of adjacency entries
//    and the ADJ_ROW array, which keeps track of the list of slots
//    in ADJ where we can store adjacency information for each row.
//
//    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
//    now we have a place to store the adjacency information.
//
//    A copy of the ADJ_ROW array is useful, as we can use it to keep track
//    of the next available entry in ADJ for adjacencies associated with
//    a given row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TETRA_NUM, the number of tetrahedrons.
//
//    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
//
//    Input, int ADJ_NUM, the total number of adjacency relationships,
//
//    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
//
//    Output, int TET_MESH_ORDER4_ADJ_SET[ADJ_NUM], 
//    the adjacency information.
//
{
  int *adj;
  int *adj_row_copy;
  int i;
  int j;
  int k;
  int node;
  int *pair;
  int pair_num;
  int tetra;
//
//  Each order 4 tetrahedron defines 6 adjacency pairs.
//
  pair = new int[2*6*element_num];

  for ( tetra = 0; tetra < element_num; tetra++ )
  {
    pair[0+             tetra *2] = element_node[0+tetra*4];
    pair[1+             tetra *2] = element_node[1+tetra*4];

    pair[0+(  element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(  element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(2*element_num+tetra)*2] = element_node[0+tetra*4];
    pair[1+(2*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(3*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(3*element_num+tetra)*2] = element_node[2+tetra*4];

    pair[0+(4*element_num+tetra)*2] = element_node[1+tetra*4];
    pair[1+(4*element_num+tetra)*2] = element_node[3+tetra*4];

    pair[0+(5*element_num+tetra)*2] = element_node[2+tetra*4];
    pair[1+(5*element_num+tetra)*2] = element_node[3+tetra*4];
  }
  pair_num = 6 * element_num;
//
//  Force the nodes of each pair to be listed in ascending order.
//
  i4col_sort2_a ( 2, pair_num, pair );
//
//  Rearrange the columns in ascending order.
//
  i4col_sort_a ( 2, pair_num, pair );
//
//  Mark all entries of ADJ so we will know later if we missed one.
//
  adj = new int[adj_num];

  for ( i = 0; i < adj_num; i++ )
  {
    adj[i] = -1;
  }
//
//  Copy the ADJ_ROW array and use it to keep track of the next
//  free entry for each row.
//
  adj_row_copy = new int[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    adj_row_copy[node] = adj_row[node];
  }
//
//  Now set up the ADJ_ROW counts.
//
  for ( k = 0; k < pair_num; k++ )
  {
    if ( 0 < k )
    {
      if ( pair[0+(k-1)*2] == pair[0+k*2] &&
           pair[1+(k-1)*2] == pair[1+k*2] )
      {
        continue;
      }
    }
    i = pair[0+k*2];
    j = pair[1+k*2];

    adj[adj_row_copy[i]] = j;
    adj_row_copy[i] = adj_row_copy[i] + 1;
    adj[adj_row_copy[j]] = i;
    adj_row_copy[j] = adj_row_copy[j] + 1;
  }
  delete [] adj_row_copy;
  delete [] pair;

  return adj;
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
//    03 October 2003
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
