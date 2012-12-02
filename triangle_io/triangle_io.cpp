# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_io.hpp"

//****************************************************************************80

void element_data_example ( int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_DATA_EXAMPLE returns the element information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_ATT_NUM, the number of element attributes.
//
//    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
//    nodes that make up each element.
//
//    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
//    of each element.
//
{
  int element_node_save[3*24] = {
    1,  2,  6,
    7,  6,  2,
    2,  3,  7,
    8,  7,  3,
    3,  4,  8,
    9,  8,  4,
    4,  5,  9,
   10,  9,  5,
    6,  7, 11,
   12, 11,  7,
    7,  8, 12,
   13, 12,  8,
    8,  9, 13,
   14, 13,  9,
    9, 10, 14,
   15, 14, 10,
   11, 12, 16,
   17, 16, 12,
   12, 13, 17,
   18, 17, 13,
   16, 17, 19,
   20, 19, 17,
   17, 18, 20,
   21, 20, 18 };

  i4mat_copy ( element_order, element_num, element_node_save, element_node );

  return;
}
//****************************************************************************80

void element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] )

//*****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_DATA_READ reads the header information from an element file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_ATT_NUM, number of element attributes listed on each 
//    node record.
//
//    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
//    nodes that make up each element.
//
//    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
//    of each element.
//
{
  int element;
  int i;
  int i1;
  int i2;
  int i3;
  ifstream input;
  int ival;
  double value;

  element = -1;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "ELEMENT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
/*
  Read, but ignore, dimension line.
*/
    if ( element == -1 )
    {
      input >> i1 >> i2 >> i3;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < element_order; i++ )
      {
        input >> ival;
        element_node[i+element*element_order] = ival;
      }
      for ( i = 0; i < element_att_num; i++ )
      {
        input >> value;
        element_att[i+element*element_att_num] = value;
      }
    }

    element = element + 1;

    if ( element_num <= element )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void element_size_example ( int *element_num, int *element_order, 
  int *element_att_num )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_SIZE_EXAMPLE returns the element size information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *ELEMENT_NUM, the number of elements.
//
//    Output, int *ELEMENT_ORDER, the order of the elements.
//
//    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
//
{
  *element_num = 24;
  *element_order = 3;
  *element_att_num = 0;

  return;
}
//****************************************************************************80

void element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_SIZE_READ reads the header information from an element file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Output, int *ELEMENT_NUM, the number of elements.
//
//    Output, int *ELEMENT_ORDER, the order of the elements.
//
//    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
//
{
  ifstream input;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "ELEMENT_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  input >> *element_num >> *element_order >> *element_att_num;

  input.close ( );

  return;
}
//****************************************************************************80

void element_write ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_WRITE writes a TRIANGLE ".ele" file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be written.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_ATT_NUM, the number of element attributes.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
//    nodes that make up each element.
//
//    Input, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
//    of each element.
//
{
  int att;
  int dim;
  int element;
  int order;
  ofstream output;
  
  output.open ( element_file.c_str ( ) );

  output << "  " << element_num 
         << "  " << element_order 
         << "  " << element_att_num << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    output << "  " << element + 1;
    for ( order = 0; order < element_order; order++ )
    {
      output << "  " << element_node[order+element*element_order];
    }
    for ( att = 0; att < element_att_num; att++ )
    {
      output << element_att[att+element*element_att_num];
    }
    output << "\n";
  }

  output.close ( );

  return;
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

void i4mat_copy ( int m, int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_COPY copies one I4MAT to another.
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
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A1[M*N], the matrix to be copied.
//
//    Output, int A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
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
//    20 August 2010
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
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
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
      cout << setw(6) << i - 1 << "  ";
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
      cout << setw(5) << j - 1 << ":";
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

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

void node_data_example ( int node_num, int node_dim, int node_att_num, 
  int node_marker_num, double node_coord[], double node_att[], 
  int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_DATA_EXAMPLE returns the node information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
//    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
//
//    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
//
//    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
//
{
  double node_coord_save[2*21] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0 };
  int node_marker_save[21] = {
    1, 1, 1, 1, 1, 1, 0, 0, 0, 1,
    1, 0, 0, 1, 1, 1, 0, 1, 1, 1,
    1 };

  r8mat_copy ( 2, 21, node_coord_save, node_coord );
  i4vec_copy ( 21, node_marker_save, node_marker );

  return;
}
//****************************************************************************80

void node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_HEADER_READ reads the header information from a node file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
//    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
//
//    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
//
//    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
//
{
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  ifstream input;
  int ival;
  int node;
  double value;

  node = -1;

  input.open ( node_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "NODE_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read, but ignore, dimension line.
//
    if ( node == -1 )
    {
      input >> i1 >> i2 >> i3 >> i4;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < node_dim; i++ )
      {
        input >> value;
        node_coord[i+node*node_dim] = value;
      }
      for ( i = 0; i < node_att_num; i++ )
      {
        input >> value;
        node_att[i+node*node_att_num] = value;
      }
      for ( i = 0; i < node_marker_num; i++ )
      {
        input >> ival;
        node_marker[i+node*node_marker_num] = ival;
      }
    }

    node = node + 1;

    if ( node_num <= node )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void node_size_example ( int *node_num, int *node_dim, int *node_att_num, 
  int *node_marker_num )

//****************************************************************************80
//
//  Purpose:
//
//    EXAMPLE_NODE_SIZE returns the sizes of node information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *NODE_DIM, the spatial dimension.
//
//    Output, int *NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Output, int *NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
{
  *node_num = 21;
  *node_dim = 2;
  *node_att_num = 0;
  *node_marker_num = 1;

  return;
}
//****************************************************************************80

void node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_SIZE_READ reads the header information from a node file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *NODE_DIM, the spatial dimension.
//
//    Output, int *NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Output, int *NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
{
  ifstream input;

  input.open ( node_file.c_str ( ) );

  input >> *node_num
        >> *node_dim
        >> *node_att_num
        >> *node_marker_num;

  input.close ( );

  return;
}
//****************************************************************************80

void node_write ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], 
  double node_att[], int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_WRITE writes a TRIANGLE ".node" file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be written.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
//    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
//
//    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
//
//    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
//
{
  int att;
  int dim;
  int node;
  ofstream output;
  
  output.open ( node_file.c_str ( ) );

  output << "  " << node_num
         << "  " << node_dim
         << "  " << node_att_num
         << "  " << node_marker_num << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    output << "  " << node + 1;
    for ( dim = 0; dim < node_dim; dim++ )
    {
      output << "  " << node_coord[dim+node*node_dim];
    }
    for ( att = 0; att < node_att_num; att++ )
    {
      output << "  " << node_att[att+node*node_att_num];
    }
    if ( node_marker_num == 1 )
    {
      output << "  " << node_marker[node];
    }
    output << "\n";
  }

  output.close ( );

  return;
}

//****************************************************************************80

void r8mat_copy ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY copies one R8MAT to another.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
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
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
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
//    Input, string TITLE, a title.
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

  cout << "\n";
  cout << title << "\n";

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
    cout << "\n";

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

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
