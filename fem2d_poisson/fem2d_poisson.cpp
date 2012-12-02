# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void assemble_poisson ( int node_num, double node_xy[],
  int element_num, int element_node[], int quad_num, int ib, double a[],
  double f[] );
int bandwidth ( int element_num, int element_node[] );
void basis_one_t3 ( double t[2*3], int i, double p[2], double *qi,
  double *dqidx, double *dqidy );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
double *dgb_mxv ( int m, int n, int ml, int mu, double a[], double x[] );
void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, string title );
double *dgb_sl ( int n, int ml, int mu, double a_lu[], int pivot[],
  double b[], int job );
void dirichlet_apply ( int node_num, double node_xy[], int node_condition[],
  int ib, double a[], double f[] );
void dirichlet_condition ( int node_num, double node_xy[], double node_rhs[] );
int file_column_count ( string input_filename );
int file_row_count ( string input_filename );
void h_coef ( int node_num, double node_xy[], double node_h[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void k_coef ( int node_num, double node_xy[], double node_k[] );
void lvec_print ( int n, bool a[], string title );
void points_plot ( string filename, int node_num, double node_xy[],
  bool node_label );
void quad_rule ( int quad_num, double quad_w[], double quad_xy[] );
double r8_abs ( double x );
double r8_huge ( void );
int r8_nint ( double x );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double r8vec_amax ( int n, double a[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
void rhs ( int node_num, double node_xy[], double node_rhs[] );
void reference_to_physical_t3 ( double t[2*3], int n, double ref[],
  double phy[] );
double *residual_poisson ( int node_num, double node_xy[], int node_condition[],
  int element_num, int element_node[], int quad_num, int ib, double a[],
  double f[], double node_u[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void solution_evaluate ( double xy[2], double t[2*3], double node_u[3], double *u,
  double *dudx, double *dudy );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( void );
double triangle_area_2d ( double t[2*3] );
bool *triangulation_order3_boundary_node ( int node_num, int element_num,
  int element_node[] );
void triangulation_order3_plot ( string filename, int node_num,
  double node_xy[], int element_num, int element_node[], int node_show,
  int element_show );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_POISSON.
//
//  Discussion:
//
//    FEM2D_POISSON solves the Poisson equation
//
//      -DEL H(X,Y) DEL U(X,Y) + K(X,Y) * U(X,Y) = F(X,Y)
//
//    in a triangulated region in the plane.
//
//    Along the boundary of the region, Dirichlet conditions
//    are imposed:
//
//      U(X,Y) = G(X,Y)
//
//    The code uses continuous piecewise linear basis functions on
//    triangles.
//
//  Problem specification:
//
//    The user defines the geometry by supplying two data files
//    which list the node coordinates, and list the nodes that make up
//    each element.
//
//    The user specifies the right hand side of the Dirichlet boundary
//    conditions by supplying a function
//
//      void dirichlet_condition ( int node_num, double node_xy[2*node_num],
//        double node_bc[node_num] )
//
//    The user specifies the coefficient function K(X,Y) of the Poisson
//    equation by supplying a routine of the form
//
//      void k_coef ( int node_num, double node_xy[2*node_num],
//        double node_k[node_num] )
//
//    The user specifies the right hand side of the Poisson equation
//    by supplying a routine of the form
//
//      void rhs ( int node_num, double node_xy[2*node_num],
//        double node_f[node_num] )
//
//  Usage:
//
//    fem2d_poisson 'prefix'
//
//    where 'prefix' is the common filename prefix so that:
//
//    * prefix_nodes.txt contains the coordinates of the nodes;
//    * prefix_elements.txt contains the indices of nodes forming each element.
//
//    Files created include:
//
//    * prefix_nodes.eps, an image of the nodes;
//    * prefix_elements.eps, an image of the elements;
//    * prefix_solution.txt, the value of the solution at every node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*NODE_NUM], the coefficient matrix.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, integer ELEMENT_ORDER, the element order.
//
//    Local, double F[NODE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, bool NODE_BOUNDARY[NODE_NUM], is TRUE if the node is
//    found to lie on the boundary of the region.
//
//    Local, int NODE_CONDITION[NODE_NUM],
//    indicates the condition used to determine the variable at a node.
//    0, there is no condition (and no variable) at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Local, int NODE_NUM, the number of nodes.
//
//    Local, double NODE_R[NODE_NUM], the residual error.
//
//    Local, double NODE_U[NODE_NUM], the finite element coefficients.
//
//    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Local, integer QUAD_NUM, the number of quadrature points used for
//    assembly.  This is currently set to 3, the lowest reasonable value.
//    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
//    of QUAD_NUM greater than 3 may be appropriate.
//
{
  double *a;
  bool debug = false;
  int dim_num;
  string element_eps_filename;
  string element_filename;
  int *element_node;
  int element_num;
  int element_order;
  int element_show;
  double *f;
  int ib;
  int ierr;
  int job;
  int node;
  bool *node_boundary;
  int *node_condition;
  string node_eps_filename;
  string node_filename;
  bool node_label;
  int node_num;
  int node_show;
  double *node_r;
  double *node_u;
  double *node_xy;
  int *pivot;
  string prefix;
  int quad_num = 7;
  string solution_filename;
  double temp;

  timestamp ( );
  cout << "\n";
  cout << "FEM2D_POISSON:\n";
  cout << "  C++ version:\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Solution of the Poisson equation in an arbitrary region\n";
  cout << "  in 2 dimensions.\n";
  cout << "\n";
  cout << "  - DEL H(x,y) DEL U(x,y) + K(x,y) * U(x,y) = F(x,y) in the region\n";
  cout << "\n";
  cout << "                                     U(x,y) = G(x,y) on the boundary.\n";
  cout << "\n";
  cout << "  The finite element method is used,\n";
  cout << "  with triangular elements,\n";
  cout << "  which must be a 3 node linear triangle.\n";
//
//  Get the filename prefix.
//
  if ( 1 <= argc )
  {
    prefix = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "  Please enter the filename prefix:\n";

    cin >> prefix;
  }
//
//  Create the file names.
//
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";
  node_eps_filename = prefix + "_nodes.eps";
  element_eps_filename = prefix + "_elements.eps";
  solution_filename = prefix + "_solution.txt";

  cout << "\n";
  cout << "  Node file is \"" << node_filename << "\".\n";
  cout << "  Element file is \"" << element_filename << "\".\n";
//
//  Read the node coordinate file.
//
  r8mat_header_read ( node_filename, &dim_num, &node_num );

  cout << "  Number of nodes =          " << node_num << "\n";

  node_condition = new int[node_num];

  node_xy = r8mat_data_read ( node_filename, dim_num, node_num );

  r8mat_transpose_print_some ( dim_num, node_num, node_xy, 1, 1, 2, 10,
    "  First 10 nodes" );
//
//  Read the element description file.
//
  i4mat_header_read ( element_filename, &element_order, &element_num );

  cout << "\n";
  cout << "  Element order =            " << element_order << "\n";
  cout << "  Number of elements =       " << element_num << "\n";

  if ( element_order != 3 )
  {
    cout << "\n";
    cout << "FEM2D_POISSON - Fatal error!\n";
    cout << "  The input triangulation has order " << element_order << "\n";
    cout << "  However, a triangulation of order 3 is required.\n";
    exit ( 1 );
  }
  element_node = i4mat_data_read ( element_filename, element_order,
    element_num );

  i4mat_transpose_print_some ( 3, element_num,
    element_node, 1, 1, 3, 10, "  First 10 elements" );

  cout << "\n";
  cout << "  Quadrature order =          " << quad_num << "\n";
//
//  Determine which nodes are boundary nodes and which have a
//  finite element unknown.  Then set the boundary values.
//
  node_boundary = triangulation_order3_boundary_node ( node_num, element_num,
    element_node );

  if ( debug )
  {
    lvec_print ( node_num, node_boundary, "    Node  Boundary?" );
  }
//
//  Determine the node conditions.
//  For now, we'll just assume all boundary nodes are Dirichlet.
//
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_boundary[node] )
    {
      node_condition[node] = 2;
    }
    else
    {
      node_condition[node] = 1;
    }
  }
//
//  Determine the bandwidth of the coefficient matrix.
//
  ib = bandwidth ( element_num, element_node );

  cout << "\n";
  cout << "  The matrix half bandwidth is " << ib << "\n";
  cout << "  The matrix bandwidth is      " << 2 * ib + 1 << "\n";
  cout << "  The storage bandwidth is     " << 3 * ib + 1 << "\n";
//
//  Make a picture of the nodes.
//
  if ( node_num <= 1000 )
  {
    if ( node_num <= 100 )
    {
      node_label = true;
    }
    else
    {
      node_label = false;
    }

    points_plot ( node_eps_filename, node_num, node_xy, node_label );
  }
//
//  Make a picture of the elements.
//
  if ( node_num <= 1000 )
  {
    if ( node_num <= 100 )
    {
      node_show = 2;
      element_show = 2;
    }
    else
    {
      node_show = 0;
      element_show = 1;
    }

    triangulation_order3_plot ( element_eps_filename, node_num,
      node_xy, element_num, element_node, node_show, element_show );
  }
//
//  Allocate space for the coefficient matrix A and right hand side F.
//
  a = new double[(3*ib+1)*node_num];
  f = new double[node_num];
  node_u = new double[node_num];
  pivot = new int[node_num];
//
//  Assemble the finite element coefficient matrix A and the right-hand side F.
//
  assemble_poisson ( node_num, node_xy, element_num,
    element_node, quad_num, ib, a, f );
//
//  Print a tiny portion of the matrix.
//
  if ( debug )
  {
    dgb_print_some ( node_num, node_num, ib, ib, a, 1, 1, 10, 10,
      "  Initial block of Finite Element matrix A:" );

    r8vec_print_some ( node_num, f, 1, 10,
      "  Part of right hand side vector:" );
  }
//
//  Adjust the linear system to account for Dirichlet boundary conditions.
//
  dirichlet_apply ( node_num, node_xy, node_condition, ib, a, f );

  if ( debug )
  {
    dgb_print_some ( node_num, node_num, ib, ib, a, 1, 1, 10, 10,
      "  Finite Element matrix A after boundary adjustments:" );

    r8vec_print_some ( node_num, f, 1, 10,
      "  Part of right hand side vector:" );
  }
//
//  Solve the linear system using a banded solver.
//
  ierr = dgb_fa ( node_num, ib, ib, a, pivot );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "FEM2D_POISSON - Fatal error!\n";
    cout << "  DGB_FA returned the error condition IERR = " << ierr << ".\n";
    cout << "\n";
    cout << "  The linear system was not factored, and the\n";
    cout << "  algorithm cannot proceed.\n";
    exit ( 1 );
  }

  job = 0;

  delete [] node_u;

  node_u = dgb_sl ( node_num, ib, ib, a, pivot, f, job );

  r8vec_print_some ( node_num, node_u, 1, 10,
    "  Part of the solution vector U:" );
//
//  Write an ASCII file that can be read into MATLAB.
//
  r8mat_write ( solution_filename, 1, node_num, node_u );

  cout << "\n";
  cout << "FEM2D_POISSON:\n";
  cout << "  Wrote an ASCII file\n";
  cout << "    \"" << solution_filename << "\".\n";
  cout << "  of the form\n";
  cout << "    U ( X(I), Y(I) )\n";
  cout << "  which can be used for plotting.\n";

  if ( debug )
  {
    r8vec_print_some ( node_num, node_u, 1, 10,
      "  Part of the solution vector:" );
  }

  node_r = residual_poisson ( node_num, node_xy, node_condition,
    element_num, element_node, quad_num, ib, a, f, node_u );

  temp = r8vec_amax ( node_num, node_r );

  cout << "\n";
  cout << "  Maximum absolute residual = " << temp << "\n";
//
//  Deallocate memory.
//
  delete [] a;
  delete [] f;
  delete [] node_boundary;
  delete [] node_condition;
  delete [] node_r;
  delete [] node_u;
  delete [] node_xy;
  delete [] pivot;
  delete [] element_node;
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_POISSON:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void assemble_poisson ( int node_num, double node_xy[],
  int element_num, int element_node[], int quad_num, int ib, double a[],
  double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_POISSON assembles the system for the Poisson equation.
//
//  Discussion:
//
//    The matrix is known to be banded.  A special matrix storage format
//    is used to reduce the space required.  Details of this format are
//    discussed in the routine DGB_FA.
//
//    Note that a 3 point quadrature rule, which is sometimes used to
//    assemble the matrix and right hand side, is just barely accurate
//    enough for simple problems.  If you want better results, you
//    should use a quadrature rule that is more accurate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Output, double A(3*IB+1,NODE_NUM), the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Output, double F(NODE_NUM), the right hand side.
//
//  Local parameters:
//
//    Local, double BI, DBIDX, DBIDY, the value of some basis function
//    and its first derivatives at a quadrature point.
//
//    Local, double BJ, DBJDX, DBJDY, the value of another basis
//    function and its first derivatives at a quadrature point.
//
{
  double area;
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int j;
  int node;
  double k_value;
  double p[2];
  double *phys_h;
  double *phys_k;
  double *phys_rhs;
  double *phys_xy;
  int quad;
  double *quad_w;
  double *quad_xy;
  double t3[2*3];
  int test;
  double *w;

  phys_h = new double[quad_num];
  phys_k = new double[quad_num];
  phys_rhs = new double[quad_num];
  phys_xy = new double[2*quad_num];
  quad_w = new double[quad_num];
  quad_xy = new double[2*quad_num];
  w = new double[quad_num];
//
//  Initialize the arrays to zero.
//
  for ( node = 0; node < node_num; node++ )
  {
    f[node] = 0.0;
  }
  for ( node = 0; node < node_num; node++ )
  {
    for ( i = 0; i < 3*ib+1; i++ )
    {
      a[i+node*(3*ib+1)] = 0.0;
    }
  }
//
//  Get the quadrature weights and nodes.
//
  quad_rule ( quad_num, quad_w, quad_xy );
//
//  Add up all quantities associated with the ELEMENT-th element.
//
  for ( element = 0; element < element_num; element++ )
  {
//
//  Make a copy of the element.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t3[i+j*2] = node_xy[i+(element_node[j+element*3]-1)*2];
      }
    }
//
//  Map the quadrature points QUAD_XY to points PHYS_XY in the physical element.
//
    reference_to_physical_t3 ( t3, quad_num, quad_xy, phys_xy );

    area = r8_abs ( triangle_area_2d ( t3 ) );

    for ( quad = 0; quad < quad_num; quad++ )
    {
      w[quad] = quad_w[quad] * area;
    }

    rhs ( quad_num, phys_xy, phys_rhs );
    h_coef ( quad_num, phys_xy, phys_h );
    k_coef ( quad_num, phys_xy, phys_k );
//
//  Consider the QUAD-th quadrature point.
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
      p[0] = phys_xy[0+quad*2];
      p[1] = phys_xy[1+quad*2];
//
//  Consider the TEST-th test function.
//
//  We generate an integral for every node associated with an unknown.
//  But if a node is associated with a boundary condition, we do nothing.
//
      for ( test = 1; test <= 3; test++ )
      {
        i = element_node[test-1+element*3];

        basis_one_t3 ( t3, test, p, &bi, &dbidx, &dbidy );

        f[i-1] = f[i-1] + w[quad] * phys_rhs[quad] * bi;
//
//  Consider the BASIS-th basis function, which is used to form the
//  value of the solution function.
//
        for ( basis = 1; basis <= 3; basis++ )
        {
          j = element_node[basis-1+element*3];

          basis_one_t3 ( t3, basis, p, &bj, &dbjdx, &dbjdy );

          a[i-j+2*ib+(j-1)*(3*ib+1)] = a[i-j+2*ib+(j-1)*(3*ib+1)]
            + w[quad] * (
             phys_h[quad] * ( dbidx * dbjdx + dbidy * dbjdy )
            + phys_k[quad] * bj * bi );
        }
      }
    }
  }
  delete [] phys_h;
  delete [] phys_k;
  delete [] phys_rhs;
  delete [] phys_xy;
  delete [] quad_w;
  delete [] quad_xy;
  delete [] w;

  return;
}
//****************************************************************************80

int bandwidth ( int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth of the coefficient matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input,  ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int BANDWIDTH, the half bandwidth of the matrix.
//
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;
  int nhba;

  nhba = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < 3; local_i++ )
    {
      global_i = element_node[local_i+element*3];
      for ( local_j = 1; local_j <= 3; local_j++ )
      {
        global_j = element_node[local_j+element*3];
        nhba = i4_max ( nhba, abs ( global_j - global_i ) );
      }
    }
  }

  return nhba;
}
//****************************************************************************80

void basis_one_t3 ( double t[2*3], int i, double p[2], double *qi,
  double *dqidx, double *dqidy )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_ONE_T3 evaluates basis functions for a linear triangular element.
//
//  Discussion:
//
//    The routine is given the coordinates of the nodes of a triangle.
//
//           3
//          / .
//         /   .
//        /     .
//       1-------2
//
//    It evaluates the linear basis function Q(I)(X,Y) associated with
//    node I, which has the property that it is a linear function
//    which is 1 at node I and zero at the other two nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the nodes.
//
//    Input, int I, the index of the desired basis function.
//    I should be between 1 and 3.
//
//    Input, double P[2], the coordinates of the point where
//    the basis function is to be evaluated.
//
//    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
//    and its X and Y derivatives.
//
{
  double area;
  int ip1;
  int ip2;
  double temp;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] )
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] )
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    cout << "\n";
    cout << "BASIS_ONE_T3 - Fatal error!\n";
    cout << "  Element has zero area.\n";
    cout << "  Area = " << area << "\n";
    cout << "\n";
    cout << "  Node 1: ( " << t[0+0*2] << ", " << t[1+0*2] << " )\n";
    cout << "  Node 2: ( " << t[0+1*2] << ", " << t[1+1*2] << " )\n";
    cout << "  Node 3: ( " << t[0+2*2] << ", " << t[1+2*2] << " )\n";
    exit ( 1 );
  }

  if ( i < 1 || 3 < i )
  {
    cout << "\n";
    cout << "BASIS_ONE_T3 - Fatal error!\n";
    cout << "  Basis index I is not between 1 and 3.\n";
    cout << "  I = " << i << "\n";
    exit ( 1 );
  }

  ip1 = i4_wrap ( i + 1, 1, 3 );
  ip2 = i4_wrap ( i + 2, 1, 3 );

  *qi = ( ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] )
        * ( p[1]           - t[1+(ip1-1)*2] )
        - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] )
        * ( p[0]           - t[0+(ip1-1)*2] ) ) / area;

  *dqidx = - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) / area;
  *dqidy =   ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) / area;

  return;
}
//****************************************************************************80

char ch_cap ( char c )

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
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 )
  {
    c = c - 32;
  }

  return c;
}
//****************************************************************************80*

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80*
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
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 )
  {
    c1 = c1 - 32;
  }
  if ( 97 <= c2 && c2 <= 122 )
  {
    c2 = c2 - 32;
  }

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
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
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
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

int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_FA performs a LINPACK-style PLU factorization of a DGB matrix.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2003
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.
//    On output, A has been overwritten by the LU factors.
//
//    Output, int PIVOT[N], the pivot vector.
//
//    Output, int DGB_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int i0;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  double t;

  m = ml + mu + 1;
//
//  Zero out the initial fill-in columns.
//
  j0 = mu + 2;
  j1 = i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      a[i-1+(jz-1)*col] = 0.0;
    }
  }

  jz = j1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Zero out the next fill-in column.
//
    jz = jz + 1;
    if ( jz <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(jz-1)*col] = 0.0;
      }
    }
//
//  Find L = pivot index.
//
    lm = i4_min ( ml, n-k );
    l = m;

    for ( j = m+1; j <= m + lm; j++ )
    {
      if ( fabs ( a[l-1+(k-1)*col] ) < fabs ( a[j-1+(k-1)*col] ) )
      {
        l = j;
      }
    }

    pivot[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*col] == 0.0 )
    {
      cout << "\n";
      cout << "DGB_FA - Fatal error!\n";
      cout << "  Zero pivot on step " << k << "\n";
      return k;
    }
//
//  Interchange if necessary.
//
    t                = a[l-1+(k-1)*col];
    a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
    a[m-1+(k-1)*col] = t;
//
//  Compute multipliers.
//
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
    }
//
//  Row elimination with column indexing.
//
    ju = i4_max ( ju, mu + pivot[k-1] );
    ju = i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;

      if ( l != mm )
      {
        t                 = a[l-1+(j-1)*col];
        a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
        a[mm-1+(j-1)*col] = t;
      }
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col]
          + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
      }
    }
  }

  pivot[n-1] = n;

  if ( a[m-1+(n-1)*col] == 0.0 )
  {
    cout << "\n";
    cout << "DGB_FA - Fatal error!\n";
    cout << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
}
//****************************************************************************80

double *dgb_mxv ( int m, int n, int ml, int mu, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_MXV multiplies a DGB matrix times a vector.
//
//  Discussion:
//
//    The DGB storage format is for an M by N banded matrix, with lower
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
//    extra superdiagonals, which may be required to store nonzero entries
//    generated during Gaussian elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    LINPACK and LAPACK storage of general band matrices requires
//    an extra ML upper diagonals for possible fill in entries during
//    Gauss elimination.  This routine does not access any entries
//    in the fill in diagonals, because it assumes that the matrix
//    has NOT had Gauss elimination applied to it.  If the matrix
//    has been Gauss eliminated, then the routine DGB_MU must be
//    used instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongarra, Bunch, Cleve Moler, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1.
//
//    Input, double A[(2*ML+MU+1)*N], the DGB matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DGB_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;
  int jhi;
  int jlo;

  b = new double[m];

  for ( i = 1; i <= m; i++ )
  {
    b[i-1] = 0.0;
    jlo = i4_max ( 1, i - ml );
    jhi = i4_min ( n, i + mu );
    for ( j = jlo; j <= jhi; j++ )
    {
      b[i-1] = b[i-1] + a[i-j+ml+mu+(j-1)*(2*ml+mu+1)] * x[j-1];
    }
  }

  return b;
}
//****************************************************************************80

void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_PRINT_SOME prints some of a DGB matrix.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower
//    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
//    superdiagonals, which may be required to store nonzero entries generated
//    during Gaussian elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2006
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
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than min(M,N)-1..
//
//    Input, double A[(2*ML+MU+1)*N], the DGB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int col = 2 * ml + mu + 1;
  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
    cout << "  Col: ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2lo = i4_max ( i2lo, j2lo - mu - ml );

    i2hi = i4_min ( ihi, m );
    i2hi = i4_min ( i2hi, j2hi + ml );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(6) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        if ( i < j - mu - ml || j + ml < i )
        {
          cout << "            ";
        }
        else
        {
          cout << setw(10) << a[i-j+ml+mu+(j-1)*col] << "  ";
        }
      }
      cout << "\n";
    }
  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

double *dgb_sl ( int n, int ml, int mu, double a_lu[], int pivot[],
  double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_SL solves a system factored by DGB_FA.
//
//  Discussion:
//
//    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
//    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
//    which may be required to store nonzero entries generated during Gaussian
//    elimination.
//
//    The original M by N matrix is "collapsed" downward, so that diagonals
//    become rows of the storage array, while columns are preserved.  The
//    collapsed array is logically 2*ML+MU+1 by N.
//
//    The two dimensional array can be further reduced to a one dimensional
//    array, stored by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Bunch, Moler, Stewart.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from DGB_FA.
//
//    Input, int PIVOT[N], the pivot vector from DGB_FA.
//
//    Input, double B[N], the right hand side vector.
//
//    Input, int JOB.
//    0, solve A * x = b.
//    nonzero, solve A' * x = b.
//
//    Output, double DGB_SL[N], the solution.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
//
  m = mu + ml + 1;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    if ( 1 <= ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = i4_min ( ml, n-k );
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
        for ( i = 1; i <= lm; i++ )
        {
          x[k+i-1] = x[k+i-1] + x[k-1] * a_lu[m+i-1+(k-1)*col];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a_lu[la+i-1+(k-1)*col];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - x[lb+i-1] * a_lu[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a_lu[m-1+(k-1)*col];
    }
//
//  Solve L' * X = Y.
//
    if ( 1 <= ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = i4_min ( ml, n-k );
        for ( i = 1; i <= lm; i++ )
        {
          x[k-1] = x[k-1] + x[k+i-1] * a_lu[m+i-1+(k-1)*col];
        }
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

void dirichlet_apply ( int node_num, double node_xy[], int node_condition[],
  int ib, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_APPLY accounts for Dirichlet boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int NODE_CONDITION[NODE_NUM], reports the condition
//    used to set the unknown associated with the node.
//    0, unknown.
//    1, finite element equation.
//    2, Dirichlet condition;
//    3, Neumann condition.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by
//    NODE_NUM coefficient matrix, stored in a compressed format; on output,
//    the matrix has been adjusted for Dirichlet boundary conditions.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//    On output, the right hand side has been adjusted for Dirichlet
//    boundary conditions.
//
{
  int column;
  int column_high;
  int column_low;
  int DIRICHLET = 2;
  int node;
  double *node_bc;

  node_bc = new double[node_num];

  dirichlet_condition ( node_num, node_xy, node_bc );

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_condition[node] == DIRICHLET )
    {
      column_low = i4_max ( node + 1 - ib, 1 );
      column_high = i4_min ( node + 1 + ib, node_num );

      for ( column = column_low; column <= column_high; column++ )
      {
        a[node+1-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
      }
      a[2*ib+node*(3*ib+1)] = 1.0;

      f[node] = node_bc[node];
    }
  }

  delete [] node_bc;

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
    exit ( 1 );
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
    exit ( 1 );
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
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
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
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of integer division.
//
//  Formula:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Discussion:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80*

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80*
//
//  Purpose:
//
//    I4_WRAP forces an integer to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I4_WRAP
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
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

int *i4mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
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
//    Output, int I4MAT_DATA_READ[M*N], the data.
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
    exit ( 1 );
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
//  Discussion:
//
//    An I4MAT is an array of I4's.
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
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    exit ( 1 );
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
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAt, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
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

void lvec_print ( int n, bool a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LVEC_PRINT prints a logical vector.
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
//    Input, int N, the number of components of the vector.
//
//    Input, bool A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i + 1 << "  "
         << setw(1) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void points_plot ( string filename, int node_num, double node_xy[],
  bool node_label )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_PLOT plots a pointset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
//
{
  int circle_size;
  int delta;
  ofstream file_unit;
  int node;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( filename.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "POINTS_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file \"" << filename << "\".\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: points_plot.C\n";
  file_unit << "%%Title: " << filename << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  file_unit << "%\n";
  file_unit << "%  Draw filled dots at each node:\n";
  file_unit << "%\n";
  file_unit << "%  Set the color to blue:\n";
  file_unit << "%\n";
  file_unit << "0.000  0.150  0.750  setrgbcolor\n";
  file_unit << "%\n";

  for ( node = 0; node < node_num; node++ )
  {
    x_ps = ( int ) (
      ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
      / ( x_max                     - x_min ) );

    y_ps = ( int ) (
      ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
      / ( y_max                     - y_min ) );

    file_unit << "newpath  "
      << x_ps << "  "
      << y_ps << "  "
      << circle_size << " 0 360 arc closepath fill\n";
  }
//
//  Label the nodes.
//
  if ( node_label )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps     << "  "
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }

  file_unit << "%\n";
  file_unit << "restore showpage\n";
  file_unit << "%\n";
  file_unit << "% End of page\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80

void quad_rule ( int quad_num, double quad_w[], double quad_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_RULE sets the quadrature rule for assembly.
//
//  Discussion:
//
//    The quadrature rule is given for a reference element.
//
//      0 <= X,
//      0 <= Y, and
//      X + Y <= 1.
//
//      ^
//    1 | *
//      | |.
//    Y | | .
//      | |  .
//    0 | *---*
//      +------->
//        0 X 1
//
//    The rules have the following precision:
//
//    QUAD_NUM  Precision
//
//     1        1
//     3        2
//     4        3
//     6        4
//     7        5
//     9        6
//    13        7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int QUAD_NUM, the number of quadrature nodes.
//
//    Output, double QUAD_W[QUAD_NUM], the quadrature weights.
//
//    Output, double QUAD_XY[2*QUAD_NUM],
//    the coordinates of the quadrature nodes.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  double t;
  double u;
  double v;
  double w;

  if ( quad_num == 1 )
  {
    quad_xy[0+0*2] = 1.0 / 3.0;
    quad_xy[1+0*2] = 1.0 / 3.0;

    quad_w[0] = 1.0;
  }
  else if ( quad_num == 3 )
  {
    quad_xy[0+0*2] = 0.5;
    quad_xy[1+0*2] = 0.0;
    quad_xy[0+1*2] = 0.5;
    quad_xy[1+1*2] = 0.5;
    quad_xy[0+2*2] = 0.0;
    quad_xy[1+2*2] = 0.5;

    quad_w[0] = 1.0 / 3.0;
    quad_w[1] = 1.0 / 3.0;
    quad_w[2] = 1.0 / 3.0;
  }
  else if ( quad_num == 4 )
  {
    a =   6.0 / 30.0;
    b =  10.0 / 30.0;
    c =  18.0 / 30.0;

    d =  25.0 / 48.0;
    e = -27.0 / 48.0;

    quad_xy[0+0*2] = b;
    quad_xy[1+0*2] = b;
    quad_xy[0+1*2] = c;
    quad_xy[1+1*2] = a;
    quad_xy[0+2*2] = a;
    quad_xy[1+2*2] = c;
    quad_xy[0+3*2] = a;
    quad_xy[1+3*2] = a;

    quad_w[0] = e;
    quad_w[1] = d;
    quad_w[2] = d;
    quad_w[3] = d;
  }
  else if ( quad_num == 6 )
  {
    a = 0.816847572980459;
    b = 0.091576213509771;
    c = 0.108103018168070;
    d = 0.445948490915965;
    v = 0.109951743655322;
    w = 0.223381589678011;

    quad_xy[0+0*2] = a;
    quad_xy[1+0*2] = b;
    quad_xy[0+1*2] = b;
    quad_xy[1+1*2] = a;
    quad_xy[0+2*2] = b;
    quad_xy[1+2*2] = b;
    quad_xy[0+3*2] = c;
    quad_xy[1+3*2] = d;
    quad_xy[0+4*2] = d;
    quad_xy[1+4*2] = c;
    quad_xy[0+5*2] = d;
    quad_xy[1+5*2] = d;

    quad_w[0] = v;
    quad_w[1] = v;
    quad_w[2] = v;
    quad_w[3] = w;
    quad_w[4] = w;
    quad_w[5] = w;
  }
  else if ( quad_num == 7 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    u = 0.225;
    v = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    w = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;

    quad_xy[0+0*2] = a;
    quad_xy[1+0*2] = a;
    quad_xy[0+1*2] = b;
    quad_xy[1+1*2] = c;
    quad_xy[0+2*2] = c;
    quad_xy[1+2*2] = b;
    quad_xy[0+3*2] = c;
    quad_xy[1+3*2] = c;
    quad_xy[0+4*2] = d;
    quad_xy[1+4*2] = e;
    quad_xy[0+5*2] = e;
    quad_xy[1+5*2] = d;
    quad_xy[0+6*2] = e;
    quad_xy[1+6*2] = e;

    quad_w[0] = u;
    quad_w[1] = v;
    quad_w[2] = v;
    quad_w[3] = v;
    quad_w[4] = w;
    quad_w[5] = w;
    quad_w[6] = w;
  }
  else if ( quad_num == 9 )
  {
    a = 0.124949503233232;
    b = 0.437525248383384;
    c = 0.797112651860071;
    d = 0.165409927389841;
    e = 0.037477420750088;

    u = 0.205950504760887;
    v = 0.063691414286223;

    quad_xy[0+0*2] = a;
    quad_xy[1+0*2] = b;
    quad_xy[0+1*2] = b;
    quad_xy[1+1*2] = a;
    quad_xy[0+2*2] = b;
    quad_xy[1+2*2] = b;
    quad_xy[0+3*2] = c;
    quad_xy[1+3*2] = d;
    quad_xy[0+4*2] = c;
    quad_xy[1+4*2] = e;
    quad_xy[0+5*2] = d;
    quad_xy[1+5*2] = c;
    quad_xy[0+6*2] = d;
    quad_xy[1+6*2] = e;
    quad_xy[0+7*2] = e;
    quad_xy[1+7*2] = c;
    quad_xy[0+8*2] = e;
    quad_xy[1+8*2] = d;

    quad_w[0] = u;
    quad_w[1] = u;
    quad_w[2] = u;
    quad_w[3] = v;
    quad_w[4] = v;
    quad_w[5] = v;
    quad_w[6] = v;
    quad_w[7] = v;
    quad_w[8] = v;
  }
  else if ( quad_num == 13 )
  {
    h = 1.0 / 3.0;
    a = 0.479308067841923;
    b = 0.260345966079038;
    c = 0.869739794195568;
    d = 0.065130102902216;
    e = 0.638444188569809;
    f = 0.312865496004875;
    g = 0.048690315425316;

    w = -0.149570044467670;
    t =  0.175615257433204;
    u =  0.053347235608839;
    v =  0.077113760890257;

    quad_xy[0+ 0*2] = h;
    quad_xy[1+ 0*2] = h;
    quad_xy[0+ 1*2] = a;
    quad_xy[1+ 1*2] = b;
    quad_xy[0+ 2*2] = b;
    quad_xy[1+ 2*2] = a;
    quad_xy[0+ 3*2] = b;
    quad_xy[1+ 3*2] = b;

    quad_xy[0+ 4*2] = c;
    quad_xy[1+ 4*2] = d;
    quad_xy[0+ 5*2] = d;
    quad_xy[1+ 5*2] = c;
    quad_xy[0+ 6*2] = d;
    quad_xy[1+ 6*2] = d;

    quad_xy[0+ 7*2] = e;
    quad_xy[1+ 7*2] = f;
    quad_xy[0+ 8*2] = e;
    quad_xy[1+ 8*2] = g;
    quad_xy[0+ 9*2] = f;
    quad_xy[1+ 9*2] = e;

    quad_xy[0+10*2] = f;
    quad_xy[1+10*2] = g;
    quad_xy[0+11*2] = g;
    quad_xy[1+11*2] = e;
    quad_xy[0+12*2] = g;
    quad_xy[1+12*2] = f;

    quad_w[ 0] = w;
    quad_w[ 1] = t;
    quad_w[ 2] = t;
    quad_w[ 3] = t;
    quad_w[ 4] = u;
    quad_w[ 5] = u;
    quad_w[ 6] = u;
    quad_w[ 7] = v;
    quad_w[ 8] = v;
    quad_w[ 9] = v;
    quad_w[10] = v;
    quad_w[11] = v;
    quad_w[12] = v;
  }
  else
  {
    cout << "\n";
    cout << "QUAD_RULE - Fatal error!\n";
    cout << "  No rule is available of order QUAD_NUM = " << quad_num << "\n";
    exit ( 1 );
  }
  return;
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
//    02 April 2005
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
  if ( 0.0 <= x )
  {
    return x;
  }
  else
  {
    return ( -x );
  }
}
//****************************************************************************80

double r8_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
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
//    Output, double R8_HUGE, a "huge" R8.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         R8_NINT
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
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
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
//    An R8MAT is an array of R8's.
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
//    Output, double R8MAT_DATA_READ[M*N], the data.
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
    exit ( 1 );
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
//  Discussion:
//
//    An R8MAT is an array of R8's.
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
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    exit ( 1 );
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

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
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
//    Input, double TABLE[M*N], the data.
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
    cerr << "  Could not open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
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

double r8vec_amax ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMAX returns the maximum absolute value in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, double AMAX, the value of the entry
//    of largest magnitude.
//
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }

  return amax;
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       << "  "
         << "  " << setw(14) << a[i-1]  << "\n";
  }

  return;
}
//****************************************************************************80

void reference_to_physical_t3 ( double t[2*3], int n, double ref[],
  double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are literally halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |.
//    |  | .
//    S  |  .
//    |  |   .
//    |  |    .
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of objects to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                 + t[i+1*2] *       + ref[0+j*2]
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

double *residual_poisson ( int node_num, double node_xy[], int node_condition[],
  int element_num, int element_node[], int quad_num, int ib, double a[],
  double f[], double node_u[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESIDUAL_POISSON evaluates the residual for the Poisson equation.
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
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the
//    coordinates of nodes.
//
//    Input, int NODE_CONDITION[NODE_NUM], reports the condition
//    used to set the unknown associated with the node.
//    0, unknown.
//    1, finite element equation.
//    2, Dirichlet condition;
//    3, Neumann condition.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Workspace, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Workspace, double F[NODE_NUM], the right hand side.
//
//    Input, double NODE_U[NODE_NUM], the value of the solution
//    at each node.
//
//    Output, double NODE_R[NODE_NUM], the finite element
//    residual at each node.
//
//  Local parameters:
//
//    Local, double BI, DBIDX, DBIDY, the value of some basis function
//    and its first derivatives at a quadrature point.
//
//    Local, double BJ, DBJDX, DBJDY, the value of another basis
//    function and its first derivatives at a quadrature point.
//
{
  double area;
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int j;
  double k_value;
  int node;
  double *node_r;
  double p[2];
  double *phys_h;
  double *phys_k;
  double *phys_rhs;
  double *phys_xy;
  int quad;
  double *quad_w;
  double *quad_xy;
  double t3[2*3];
  int test;
  double *w;

  phys_h = new double[quad_num];
  phys_k = new double[quad_num];
  phys_rhs = new double[quad_num];
  phys_xy = new double[2*quad_num];

  quad_w = new double[quad_num];
  quad_xy = new double[2*quad_num];

  w = new double[quad_num];
//
//  Initialize the arrays to zero.
//
  for ( node = 0; node < node_num; node++ )
  {
    f[node] = 0.0;
  }
  for ( node = 0; node < node_num; node++ )
  {
    for ( i = 0; i < 3*ib+1; i++ )
    {
      a[i+node*(3*ib+1)] = 0.0;
    }
  }
//
//  Get the quadrature weights and nodes.
//
  quad_rule ( quad_num, quad_w, quad_xy );
//
//  The actual values of A and F are determined by summing up
//  contributions from all the elements.
//
  for ( element = 0; element < element_num; element++ )
  {
//
//  Make a copy of the element.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t3[i+j*2] = node_xy[i+(element_node[j+element*3]-1)*2];
      }
    }
//
//  Map the quadrature points QUAD_XY to points XY in the physical element.
//
    reference_to_physical_t3 ( t3, quad_num, quad_xy, phys_xy );

    area = r8_abs ( triangle_area_2d ( t3 ) );

    for ( quad = 0; quad < quad_num; quad++ )
    {
      w[quad] = quad_w[quad] * area;
    }

    rhs ( quad_num, phys_xy, phys_rhs );
    h_coef ( quad_num, phys_xy, phys_h );
    k_coef ( quad_num, phys_xy, phys_k );
//
//  Consider a quadrature point QUAD, with coordinates (X,Y).
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
      p[0] = phys_xy[0+quad*2];
      p[1] = phys_xy[1+quad*2];
//
//  Consider one of the basis functions, which will play the
//  role of test function in the integral.
//
//  We generate an integral for every node associated with an unknown.
//  But if a node is associated with a boundary condition, we do nothing.
//
      for ( test = 1; test <= 3; test++ )
      {
        i = element_node[test-1+element*3];

        basis_one_t3 ( t3, test, p, &bi, &dbidx, &dbidy );

        f[i-1] = f[i-1] + w[quad] * phys_rhs[quad] * bi;
//
//  Consider another basis function, which is used to form the
//  value of the solution function.
//
        for ( basis = 1; basis <= 3; basis++ )
        {
          j = element_node[basis-1+element*3];

          basis_one_t3 ( t3, basis, p, &bj, &dbjdx, &dbjdy );

          a[i-j+2*ib+(j-1)*(3*ib+1)] = a[i-j+2*ib+(j-1)*(3*ib+1)] + w[quad]
            * ( phys_h[quad] * ( dbidx * dbjdx + dbidy * dbjdy )
            + phys_k[quad] * bj * bi );
        }
      }
    }
  }
//
//  Apply boundary conditions.
//
  dirichlet_apply ( node_num, node_xy, node_condition, ib, a, f );
//
//  Compute A*U.
//
  node_r = dgb_mxv ( node_num, node_num, ib, ib, a, node_u );
//
//  Set RES = A * U - F.
//
  for ( node = 0; node < node_num; node++ )
  {
    node_r[node] = node_r[node] - f[node];
  }

  delete [] phys_h;
  delete [] phys_k;
  delete [] phys_rhs;
  delete [] phys_xy;
  delete [] quad_w;
  delete [] quad_xy;
  delete [] w;

  return node_r;
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

void solution_evaluate ( double xy[2], double t[2*3], double node_u[3],
  double *u, double *dudx, double *dudy )

//****************************************************************************80
//
//  Purpose:
//
//    SOLUTION_EVALUATE evaluates the solution at a point in a element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XY[2], the point where the solution is
//    to be evaluated.
//
//    Input, double T[2*3], the coordinates of the vertices
//    of the element which contains XY.
//
//    Input, double NODE_U[3], the value of the solution
//    at the nodes of the element.
//
//    Output, double U, DUDX, DUDY, the solution and its X and
//    Y derivatives at XY.
//
{
  double b;
  double dbdx;
  double dbdy;
  int i;

  *u = 0.0;
  *dudx = 0.0;
  *dudy = 0.0;

  for ( i = 1; i <= 3; i++ )
  {
    basis_one_t3 ( t, i, xy, &b, &dbdx, &dbdy );

    *u    = *u    + node_u[i-1] * b;
    *dudx = *dudx + node_u[i-1] * dbdx;
    *dudy = *dudy + node_u[i-1] * dbdy;
  }
  return;
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
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis and Herbert Wilf,
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

void timestamp ( void )

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
//****************************************************************************80

double triangle_area_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counterclockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.
//
{
  double area;

  area = 0.5 * (
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) +
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) +
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );

  return area;
}
//****************************************************************************80

bool *triangulation_order3_boundary_node ( int node_num, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_BOUNDARY_NODE indicates nodes on the boundary.
//
//  Discussion:
//
//    This routine is given a triangulation, an abstract list of triples
//    of nodes.  It is assumed that the nodes in each triangle are listed
//    in a counterclockwise order, although the routine should work
//    if the nodes are consistently listed in a clockwise order as well.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
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
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, bool TRIANGULATION_ORDER3_BOUNDARY_NODE[NODE_NUM],
//    is TRUE if the node is on a boundary edge.
//
{
  int e1;
  int e2;
  int *edge;
  bool equal;
  int i;
  int j;
  int m;
  int n;
  bool *node_boundary;

  m = 2;
  n = 3 * element_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < element_num; j++ )
  {
    edge[0+(j               )*m] = element_node[0+j*3];
    edge[1+(j               )*m] = element_node[1+j*3];
    edge[0+(j+  element_num)*m] = element_node[1+j*3];
    edge[1+(j+  element_num)*m] = element_node[2+j*3];
    edge[0+(j+2*element_num)*m] = element_node[2+j*3];
    edge[1+(j+2*element_num)*m] = element_node[0+j*3];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[1+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[1+j*m] );
    edge[0+j*m] = e1;
    edge[1+j*m] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Records which appear twice are internal edges and can be ignored.
//
  node_boundary = new bool[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    node_boundary[i] = false;
  }

  j = 0;

  while ( j < 3 * element_num )
  {
    j = j + 1;

    if ( j == 3 * element_num )
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
      break;
    }

    equal = true;

    for ( i = 0; i < m; i++ )
    {
      if ( edge[i+(j-1)*m] != edge[i+j*m] )
      {
        equal = false;
      }
    }

    if ( equal )
    {
      j = j + 1;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
    }

  }

  return node_boundary;
}
//****************************************************************************80

void triangulation_order3_plot ( string filename, int node_num,
  double node_xy[], int element_num, int element_node[], int node_show,
  int element_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_PLOT plots a triangulation of a set of nodes.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the output file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM], lists, for each triangle,
//    the indices of the nodes that form the vertices of the triangle.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int TRIANGLE_SHOW:
//    0, do not show triangles;
//    1, show triangles;
//    2, show triangles and label them.
//
{
  double ave_x;
  double ave_y;
  int circle_size;
  int delta;
  int e;
  ofstream file_unit;
  int i;
  int node;
  int triangle;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;

//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( filename.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file \"" << filename << "\".\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order3_plot.C\n";
  file_unit << "%%Title: " << filename << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at each node:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750  setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps << "  "
        << y_ps << "  "
        << circle_size << " 0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps     << "  "
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }
//
//  Draw the triangles.
//
  if ( 1 <= element_show )
  {
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to red.\n";
    file_unit << "%\n";
    file_unit << "0.900  0.200  0.100 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Draw the triangles.\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < element_num; triangle++ )
    {
      file_unit << "newpath\n";

      for ( i = 1; i <= 4; i++ )
      {
        e = i4_wrap ( i, 1, 3 );

        node = element_node[e-1+triangle*3] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        if ( i == 1 )
        {
          file_unit << x_ps << "  " << y_ps << "  moveto\n";
        }
        else
        {
          file_unit << x_ps << "  " << y_ps << "  lineto\n";
        }
      }
      file_unit << "stroke\n";
    }
  }
//
//  Label the triangles.
//
  if ( 2 <= element_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the triangles.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker red.\n";
    file_unit << "%\n";
    file_unit << "0.950  0.250  0.150 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < element_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 1; i <= 3; i++ )
      {
        node = element_node[i-1+triangle*3] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / 3.0;
      ave_y = ave_y / 3.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      file_unit << x_ps << "  "
                << y_ps << "  moveto ("
                << triangle+1 << ") show\n";
    }
  }

  file_unit << "%\n";
  file_unit << "restore  showpage\n";
  file_unit << "%\n";
  file_unit << "%  End of page.\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
