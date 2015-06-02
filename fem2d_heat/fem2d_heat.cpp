# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void assemble_backward_euler ( int node_num, double node_xy[],
  int element_order, int element_num, int element_node[], int quad_num,
  int ib, double time, double time_step_size, double u_old[], double a[],
  double f[] );
void assemble_boundary ( int node_num, double node_xy[], int node_condition[],
  int ib, double time, double a[], double f[] );
void assemble_heat ( int node_num, double node_xy[], int node_condition[],
  int element_order, int element_num, int element_node[], int quad_num,
  int ib, double time, double a[], double f[] );
int bandwidth ( int element_order, int element_num, int element_node[] );
void basis_11_t6 ( double t[2*3], int i, double p[2], double *qi, double *dqidx,
  double *dqidy );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int dgb_fa ( int n, int ml, int mu,  double abd[], int ipvt[] );
void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, string title );
double *dgb_sl ( int n, int ml, int mu, double a_lu[], int pivot[],
  double b[], int job );
double *dirichlet_condition ( int node_num, double node_xy[], double time );
int file_column_count ( string input_filename );
void file_name_inc ( string *file_name );
void file_name_specification ( int argc, char *argv[], char *node_file_name,
  char *element_file_name );
int file_row_count ( string input_filename );
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
double *initial_condition ( int node_num, double node_xy[], double time );
void i4vec_print_some ( int n, int a[], int max_print, string title );
double k_coef ( int n, double node_xy[], double time );
void lvec_print ( int n, bool a[], string title );
void quad_rule ( int quad_num, double quad_w[], double quad_xy[] );
double r8_abs ( double x );
double r8_huge ( );
int r8_nint ( double x );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
double rhs ( int n, double node_xy[], double time );
void reference_to_physical_t3 ( double t[2*3], int n, double ref[],
  double phy[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( );
double triangle_area_2d ( double t[2*3] );
bool *triangulation_order6_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] );


//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_HEAT.
//
//  Discussion:
//
//    FEM2D_HEAT solves the heat equation
//
//      dUdT - Laplacian U(X,Y,T) + K(X,Y,T) * U(X,Y,T) = F(X,Y,T)
//
//    in a triangulated region in the plane.
//
//    Along the boundary of the region, Dirichlet conditions
//    are imposed:
//
//      U(X,Y,T) = G(X,Y,T)
//
//    At the initial time T_INIT, the value of U is given
//    at all points in the region:
//
//      U(X,Y,T_INIT) = H(X,Y)
//
//    The code uses continuous piecewise linear basis functions on
//    triangles.
//
//    The backward Euler approximation is used for the time derivatives.
//
//  Problem specification:
//
//    The user defines the geometry by supplying two data files
//    which list the node coordinates, and list the nodes that make up
//    each triangular element..
//
//    The user specifies the coefficient function K(X,Y,T)
//    by supplying a routine of the form
//
//      double k_coef ( int node_num, double node_xy[], double time )
//
//    The user specifies the right hand side
//    by supplying a routine of the form
//
//     double rhs ( int node_num, double node_xy[], double time )
//
//    The user specifies the right hand side of the Dirichlet boundary
//    conditions by supplying a function
//
//      double *dirichlet_condition ( int node_num, double node_xy[],
//        double time )
//
//    The user specifies the initial condition by supplying a function
//
//      double initial_condition ( int node_num, double node_xy[], double time )
//
//  Usage:
//
//    fem2d_heat prefix
//
//    invokes the program:
//
//    * "prefix"_nodes.txt contains the coordinates of the nodes;
//    * "prefix"_elements.txt contains the indices of nodes that make up each
//      triangular element.
//
//    Files created include:
//
//    * "prefix"_u0000.txt, the initial value of the solution;
//    * "prefix"_u0001.txt and so on, the computed solution at later times;
//    * "prefix"_times.txt, the value of time at each step, from the initial to
//      final times.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, double A(3*IB+1,NODE_NUM), the coefficient matrix.
//
//    Local, int DIM_NUM, the spatial dimension, which is 2.
//
//    Local, string ELEMENT_FILE_NAME, the name of the
//    input file containing the element information.
//
//    Local, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, int ELEMENT_ORDER, the order of each element.
//
//    Local, double F[NODE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, int NODE_NUM, the number of nodes.
//
//    Local, logical NODE_BOUNDARY[NODE_NUM], is TRUE if a given
//    node is on the boundary.
//
//    Local, int NODE_CONDITION[NODE_NUM], indicates the type of
//    boundary condition being applied to nodes on the boundary.
//
//    Local, string NODE_FILE_NAME, the name of the
//    input file containing the node coordinate information.
//
//    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Local, integer QUAD_NUM, the number of quadrature points used for
//    assembly.  This is currently set to 3, the lowest reasonable value.
//    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
//    of QUAD_NUM greater than 3 may be appropriate.
//
//    Local, double U[NODE_NUM], the finite element coefficients
//    defining the solution at the current time.
//
//    Local, double U_OLD[NODE_NUM], the finite element coefficients
//    defining the solution at the previous time.
//
{
  double *a;
  bool debug = false;
  int dim_num;
  string element_filename;
  int *element_node;
  int element_num;
  int element_order;
  double *f;
  int i;
  int ib;
  int ierr;
  int job;
  int node;
  bool *node_boundary;
  int *node_condition;
  string node_filename;
  int node_num;
  double *node_xy;
  int *pivot;
  string prefix;
  int quad_num = 7;
  string solution_filename;
  double temp;
  double time;
  string time_filename;
  double time_final;
  double time_init;
  double time_old;
  int time_step;
  int time_step_num;
  double time_step_size;
  ofstream time_unit;
  double *u;
  double *u_old;

  timestamp ( );
  cout << "\n";
  cout << "FEM2D_HEAT\n";
  cout << "  C++ version:\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Solution of the time dependent heat equation\n";
  cout << "  on an arbitrary triangulated region D in 2 dimensions.\n";
  cout << "\n";
  cout << "  Ut - Uxx - Uyy + K(x,y,t) * U = F(x,y,t) in D\n";
  cout << "                              U = G(x,y,t) on the boundary.\n";
  cout << "                              U = H(x,y,t) at initial time.\n";
  cout << "\n";
  cout << "  The finite element method is used with\n";
  cout << "  6 node quadratic triangular elements (\"T6\").\n";
  cout << "\n";
  cout << "  The time derivative is approximated using the\n";
  cout << "  backward Euler method.\n";
//
//  Get the filename prefix.
//
  if ( 1 < argc )
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
  solution_filename = prefix + "_u0000.txt";
  time_filename = prefix + "_times.txt";

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

  if ( element_order != 6 )
  {
    cout << "\n";
    cout << "FEM2D_HEAT - Fatal error!\n";
    cout << "  The input triangulation has order " << element_order << "\n";
    cout << "  However, a triangulation of order 6 is required.\n";
    exit ( 1 );
  }
  element_node = i4mat_data_read ( element_filename, element_order,
    element_num );

  i4mat_transpose_print_some ( element_order, element_num,
    element_node, 1, 1, element_order, 10, "  First 10 elements" );

  cout << "\n";
  cout << "  Quadrature order =          " << quad_num << "\n";
//
//  Determine which nodes are boundary nodes and which have a
//  finite element unknown.  Then set the boundary values.
//
  node_boundary = triangulation_order6_boundary_node ( node_num, element_num,
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
  ib = bandwidth ( element_order, element_num, element_node );

  cout << "\n";
  cout << "  The matrix half bandwidth is " << ib << "\n";
  cout << "  The matrix bandwidth is      " << 2 * ib + 1 << "\n";
  cout << "  The storage bandwidth is     " << 3 * ib + 1 << "\n";
//
//  Set time stepping quantities.
//
  time_init = 0.0;
  time_final = 0.5;
  time_step_num = 10;
  time_step_size = ( time_final - time_init ) / ( double ) ( time_step_num );

  cout << "\n";
  cout << "  Initial time = " << time_init << "\n";
  cout << "  Final time =   " << time_final << "\n";
  cout << "  Step size =    " << time_step_size << "\n";
  cout << "  Number of steps = " << time_step_num << "\n";
//
//  Allocate space for the coefficient matrix A and right hand side F.
//
  a = new double[(3*ib+1)*node_num];
  f = new double[node_num];
  pivot = new int[node_num];
  u_old = new double[node_num];
//
//  Set the value of U at the initial time.
//
  time = time_init;

  u = initial_condition ( node_num, node_xy, time );

  time_unit.open ( time_filename.c_str ( ) );

  if ( !time_unit )
  {
    cout << "\n";
    cout << "FEM2D_HEAT - Fatal error!\n";
    cout << "  Could not open the output time file.\n";
    exit ( 1 );
  }

  time_unit << setw(14) << time << "\n";

  r8mat_write ( solution_filename, 1, node_num, u );
//
//  Time looping.
//
  for ( time_step = 1; time_step <= time_step_num; time_step++ )
  {
    time_old = time;
    for ( node = 0; node < node_num; node++ )
    {
      u_old[node] = u[node];
    }
    time = ( ( double ) ( time_step_num - time_step ) * time_init
           + ( double ) (                 time_step ) * time_final )
           / ( double ) ( time_step_num             );
//
//  Assemble the finite element coefficient matrix A and the right-hand side F.
//
    assemble_heat ( node_num, node_xy, node_condition, element_order,
      element_num, element_node, quad_num, ib, time, a, f );

    if ( debug )
    {
      dgb_print_some ( node_num, node_num, ib, ib, a, 1, 1, 10, 10,
        "  Initial block of Finite Element matrix A:" );

      r8vec_print_some ( node_num, f, 1, 10,
        "  Part of right hand side vector:" );
    }
//
//  Adjust the linear system for the dU/dT term, which we are treating
//  using the backward Euler formula.
//
    assemble_backward_euler ( node_num, node_xy, element_order,
      element_num, element_node, quad_num, ib, time, time_step_size,
      u_old, a, f );

    if ( debug )
    {
      dgb_print_some ( node_num, node_num, ib, ib, a, 1, 1, 10, 10,
        "  A after DT adjustment:" );

      r8vec_print_some ( node_num, f, 1, 10,
        "  F after DT adjustment:" );
    }
//
//  Adjust the linear system to account for Dirichlet boundary conditions.
//
    assemble_boundary ( node_num, node_xy, node_condition, ib, time, a, f );

    if ( debug )
    {
      dgb_print_some ( node_num, node_num, ib, ib, a, 1, 1, 10, 10,
        "  Finite Element matrix A after boundary condition adjustment:" );

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
      cout << "FEM2D_HEAT - Fatal error!\n";
      cout << "  DGB_FA returned the error condition IERR = " << ierr << ".\n";
      cout << "\n";
      cout << "  The linear system was not factored, and the\n";
      cout << "  algorithm cannot proceed.\n";
      exit ( 1 );
    }

    job = 0;

    delete [] u;

    u = dgb_sl ( node_num, ib, ib, a, pivot, f, job );

    if ( debug )
    {
      r8vec_print_some ( node_num, u, 1, 10, "  Part of the solution vector:" );
    }
//
//  Increment the file name, and write the new solution.
//
    time_unit << setw(14) << time << "\n";

    file_name_inc ( &solution_filename );

    r8mat_write ( solution_filename, 1, node_num, u );
  }

  time_unit.close ( );
//
//  Free memory.
//
  delete [] a;
  delete [] element_node;
  delete [] f;
  delete [] node_boundary;
  delete [] node_condition;
  delete [] node_xy;
  delete [] pivot;
  delete [] u;
  delete [] u_old;
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_HEAT:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void assemble_backward_euler ( int node_num, double node_xy[],
  int element_order, int element_num, int element_node[], int quad_num,
  int ib, double time, double time_step_size, double u_old[], double a[],
  double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_BACKWARD_EULER adjusts the system for the backward Euler term.
//
//  Discussion:
//
//    The input linear system
//
//      A * U = F
//
//    is appropriate for the equation
//
//      -Uxx - Uyy - K * U = RHS
//
//    We need to modify the matrix A and the right hand side F to
//    account for the approximation of the time derivative in
//
//      Ut - Uxx - Uyy - K * U = RHS
//
//    by the backward Euler approximation:
//
//      Ut approximately equal to ( U - Uold ) / dT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY(2,NODE_NUM), the coordinates of nodes.
//
//    Input, int ELEMENT_ORDER, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Input, double TIME_STEP_SIZE, the size of the time step.
//
//    Input, double U_OLD(NODE_NUM), the finite element
//    coefficients for the solution at the previous time.
//
//    Input/output, double A(3*IB+1,NODE_NUM), the NODE_NUM
//    by NODE_NUM coefficient matrix, stored in a compressed format.
//
//    Input/output, double F(NODE_NUM), the right hand side.
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
  double p[2];
  double *phys_xy;
  int quad;
  double *quad_w;
  double *quad_xy;
  double t3[2*3];
  double t6[2*6];
  int test;
  double *w;

  phys_xy = new double[2*quad_num];
  quad_w = new double[quad_num];
  quad_xy = new double[2*quad_num];
  w = new double[quad_num];
//
//  Get the quadrature rule weights and nodes.
//
  quad_rule ( quad_num, quad_w, quad_xy );

  for ( element = 0; element < element_num; element++ )
  {
//
//  Make two copies of the triangle.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t3[i+j*2] = node_xy[i+(element_node[j+element*element_order]-1)*2];
      }
    }
    for ( j = 0; j < 6; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t6[i+j*2] = node_xy[i+(element_node[j+element*element_order]-1)*2];
      }
    }
//
//  Map the quadrature points QUAD_XY to points PHYS_XY in the physical triangle.
//
    reference_to_physical_t3 ( t3, quad_num, quad_xy, phys_xy );

    area = r8_abs ( triangle_area_2d ( t3 ) );

    for ( quad = 0; quad < quad_num; quad++ )
    {
      w[quad] = quad_w[quad] * area;
    }

    for ( quad = 0; quad < quad_num; quad++ )
    {
      p[0] = phys_xy[0+quad*2];
      p[1] = phys_xy[1+quad*2];

      for ( test = 1; test <= element_order; test++ )
      {
        node = element_node[test-1+element*element_order];

        basis_11_t6 ( t6, test, p, &bi, &dbidx, &dbidy );
//
//  Carry the U_OLD term to the right hand side.
//
        f[node-1] = f[node-1] + w[quad] * bi * u_old[node-1] / time_step_size;
//
//  Modify the diagonal entries of A.
//
        for ( basis = 1; basis <= element_order; basis++ )
        {
          j = element_node[basis-1+element*element_order];

          basis_11_t6 ( t6, basis, p, &bj, &dbjdx, &dbjdy );

          a[node-j+2*ib+(j-1)*(3*ib+1)] = a[node-j+2*ib+(j-1)*(3*ib+1)]
            + w[quad] * bi * bj / time_step_size;

        }
      }
    }
  }
  delete [] phys_xy;
  delete [] quad_w;
  delete [] quad_xy;
  delete [] w;

  return;
}
//****************************************************************************80

void assemble_boundary ( int node_num, double node_xy[], int node_condition[],
  int ib, double time, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_BOUNDARY modifies the linear system for the boundary conditions.
//
//  Discussion:
//
//    For now, we are only working with Dirichlet boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 January 2007
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
//    Input, double TIME, the current time.
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
  double *bc_value;
  int column;
  int column_high;
  int column_low;
  int DIRICHLET = 2;
  int node;

  bc_value = dirichlet_condition ( node_num, node_xy, time );

  for ( node = 1; node <= node_num; node++ )
  {
    if ( node_condition[node-1] == DIRICHLET )
    {
      column_low = i4_max ( node - ib, 1 );
      column_high = i4_min ( node + ib, node_num );

      for ( column = column_low; column <= column_high; column++ )
      {
        a[node-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
      }
      a[2*ib+(node-1)*(3*ib+1)] = 1.0;

      f[node-1] = bc_value[node-1];
    }
  }

  delete [] bc_value;

  return;
}
//****************************************************************************80

void assemble_heat ( int node_num, double node_xy[], int node_condition[],
  int element_order, int element_num, int element_node[], int quad_num,
  int ib, double time, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_HEAT assembles the finite element system for the heat equation.
//
//  Discussion:
//
//    The matrix is known to be banded.  A special matrix storage format
//    is used to reduce the space required.  Details of this format are
//    discussed in the routine DGB_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 July 2007
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
//    Input, int NODE_CONDITION(NODE_NUM), reports the condition
//    used to set the unknown associated with the node.
//    0, unknown.
//    1, finite element equation.
//    2, Dirichlet condition;
//    3, Neumann condition.
//
//    Input, int ELEMENT_NUM, the number of element.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
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
  double *phys_xy;
  int quad;
  double *quad_w;
  double *quad_xy;
  double rhs_value;
  double t3[2*3];
  double t6[2*6];
  int test;
  double *w;

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
//  Make two copies of the triangle.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t3[i+j*2] = node_xy[i+(element_node[j+element*element_order]-1)*2];
      }
    }
    for ( j = 0; j < 6; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t6[i+j*2] = node_xy[i+(element_node[j+element*element_order]-1)*2];
      }
    }
//
//  Map the quadrature points QUAD_XY to points PHYS_XY in the physical triangle.
//
    reference_to_physical_t3 ( t3, quad_num, quad_xy, phys_xy );

    area = r8_abs ( triangle_area_2d ( t3 ) );

    for ( quad = 0; quad < quad_num; quad++ )
    {
      w[quad] = quad_w[quad] * area;
    }
//
//  Consider the QUAD-th quadrature point.
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
      k_value = k_coef ( 1, phys_xy+quad*2, time );
      rhs_value =  rhs ( 1, phys_xy+quad*2, time );
//
//  Consider the TEST-th test function.
//
//  We generate an integral for every node associated with an unknown.
//  But if a node is associated with a boundary condition, we do nothing.
//
      for ( test = 1; test <= element_order; test++ )
      {
        i = element_node[test-1+element*element_order];

        basis_11_t6 ( t6, test, phys_xy+quad*2, &bi, &dbidx, &dbidy );

        f[i-1] = f[i-1] + w[quad] * rhs_value * bi;
//
//  Consider the BASIS-th basis function, which is used to form the
//  value of the solution function.
//
        for ( basis = 1; basis <= element_order; basis++ )
        {
          j = element_node[basis-1+element*element_order];

          basis_11_t6 ( t6, basis, phys_xy+quad*2, &bj, &dbjdx, &dbjdy );

          a[i-j+2*ib+(j-1)*(3*ib+1)] = a[i-j+2*ib+(j-1)*(3*ib+1)]
            + w[quad] * ( dbidx * dbjdx + dbidy * dbjdy + k_value * bj * bi );
        }
      }
    }
  }
  delete [] phys_xy;
  delete [] quad_w;
  delete [] quad_xy;
  delete [] w;

  return;
}
//****************************************************************************80

int bandwidth ( int element_order, int element_num, int element_node[] )

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
//    02 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
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

  for ( element = 1; element <= element_num; element++ )
  {
    for ( local_i = 1; local_i <= element_order; local_i++ )
    {
      global_i = element_node[local_i-1+(element-1)*element_order];
      for ( local_j = 1; local_j <= element_order; local_j++ )
      {
        global_j = element_node[local_j-1+(element-1)*element_order];
        nhba = i4_max ( nhba, abs ( global_j - global_i ) );
      }
    }
  }

  return nhba;
}
//****************************************************************************80

void basis_11_t6 ( double t[2*6], int i, double p[], double *bi,
  double *dbidx, double *dbidy )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T6: one basis at one point for the T6 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the nodes of a triangle.
//
//           3
//          / \
//         6   5
//        /     \
//       1---4---2
//
//    It evaluates the quadratic basis function B(I)(X,Y) associated with
//    node I, which has the property that it is a quadratic function
//    which is 1 at node I and zero at the other five nodes.
//
//    This routine assumes that the sides of the triangle are straight,
//    so that the midside nodes fall on the line between two vertices.
//
//    This routine relies on the fact that each basis function can be
//    written as the product of two linear factors, which are easily
//    computed and normalized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the coordinates of the nodes.
//
//    Input, int I, the index of the desired basis function.
//    I should be between 1 and 6.
//
//    Input, double P[2], the coordinates of a point at which the basis
//    function is to be evaluated.
//
//    Output, double *BI, *DBIDX, *DBIDY, the values of the basis function
//    and its X and Y derivatives.
//
{
  double gf;
  double gn;
  double hf;
  double hn;
  int j1;
  int j2;
  int k1;
  int k2;

  if ( i < 1 || 6 < i )
  {
    cout << "\n";
    cout << "BASIS_11_T6 - Fatal error!\n";
    cout << "  Basis index I is not between 1 and 6.\n";
    cout << "  I = " << i << "\n";
    exit ( 1 );
  }
//
//  Determine the pairs of nodes.
//
  if ( i <= 3 )
  {
    j1 = i4_wrap ( i + 1, 1, 3 );
    j2 = i4_wrap ( i + 2, 1, 3 );
    k1 = i + 3;
    k2 = i4_wrap ( i + 5, 4, 6 );
  }
  else
  {
    j1 = i - 3;
    j2 = i4_wrap ( i - 3 + 2, 1, 3 );
    k1 = i4_wrap ( i - 3 + 1, 1, 3 );
    k2 = i4_wrap ( i - 3 + 2, 1, 3 );
  }
//
//  For C++ indexing, it is helpful to knock the indices down by one.
//
  i  = i  - 1;
  j1 = j1 - 1;
  j2 = j2 - 1;
  k1 = k1 - 1;
  k2 = k2 - 1;
//
//  Evaluate the two linear factors GF and HF,
//  and their normalizers GN and HN.
//
  gf = ( p[0]      - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] )
     - ( t[0+j2*2] - t[0+j1*2] ) * ( p[1]      - t[1+j1*2] );

  gn = ( t[0+i*2]  - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] )
     - ( t[0+j2*2] - t[0+j1*2] ) * ( t[1+i*2]  - t[1+j1*2] );

  hf = ( p[0]      - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] )
     - ( t[0+k2*2] - t[0+k1*2] ) * ( p[1]      - t[1+k1*2] );

  hn = ( t[0+i*2]  - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] )
     - ( t[0+k2*2] - t[0+k1*2] ) * ( t[1+i*2]  - t[1+k1*2] );
//
//  Construct the basis function and its derivatives.
//
  *bi =        ( gf                      / gn )
             * ( hf                      / hn );

  *dbidx =   ( ( t[1+j2*2] - t[1+j1*2] ) / gn )
           * (   hf                      / hn )
           + (   gf                      / gn )
           * ( ( t[1+k2*2] - t[1+k1*2] ) / hn );

  *dbidy = - ( ( t[0+j2*2] - t[0+j1*2] ) / gn )
           * (   hf                      / hn )
           - (   gf                      / gn )
           * ( ( t[0+k2*2] - t[0+k1*2] ) / hn );

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
//    Input, string TITLE, a title to print.
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

void file_name_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the input string contains no digits or is empty,
//    an error condition results.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILE_NAME_INC - Fatal error!\n";
    cerr << "  Input file name is empty string.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    cerr << "\n";
    cerr << "FILE_NAME_INC - Fatal error!\n";
    cerr << "  Filename contained no digits.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void file_name_specification ( int argc, char *argv[], char *node_file_name,
  char *element_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_SPECIFICATION determines the names of the input files.
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
//    Output, char *NODE_FILE_NAME, the name of the node file.
//
//    Output, char *ELEMENT_FILE_NAME, the name of the element file.
//
{
//
//  Get the number of command line arguments.
//
  if ( 1 <= argc )
  {
    strcpy ( node_file_name, argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "FILE_NAME_SPECIFICATION:\n";
    cout << "  Please enter the name of the node file.\n";

    cin.getline ( node_file_name, sizeof ( node_file_name ) );
  }

//
//  If at least two command line arguments, the second is the element file.
//
  if ( 2 <= argc )
  {
    strcpy ( element_file_name, argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "FILE_NAME_SPECIFICATION:\n";
    cout << "  Please enter the name of the element file.\n";

    cin.getline ( element_file_name, sizeof ( element_file_name ) );
  }

  return;
}
//****************************************************************************80

int file_row_count ( string filename )

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
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int record_num;
  int row_num;
  string text;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the file: \"" << filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( text[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( text ) == 0 )
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
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
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
//        I         J     MOD  I_MODP   I4_MODP Factorization
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
//    I4COL_SWAP swaps two columns of an integer array.
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
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
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
//    Input, string TITLE, a title for the matrix.
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

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_print_some ( int n, int a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT_SOME prints "some" of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4 values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, int I_LO, I_HI, the first and last indices to print.
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
    cout << "  " << setw(8)  << i
         << ": " << setw(12) << a[i-1] << "\n";
  }

  return;
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
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
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

double r8_huge ( )

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
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
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

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
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
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
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
    cerr << "  Could not open the output file.\n";
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

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
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
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i-1]  << "\n";
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
//    02 May 2011
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
  static double ten = 10.0;

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
      rexp = pow ( ten, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( ten, rexp );
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
//    FORTRAN77 original by Nijenhuis and Wilf.
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

bool *triangulation_order6_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_BOUNDARY_NODE indicates nodes on the boundary.
//
//  Discussion:
//
//    This routine is given an order 6 triangulation, an abstract list of
//    sets of six nodes.  The vertices are listed clockwise, then the
//    midside nodes.
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
//    25 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Output, bool TRIANGULATION_ORDER6_BOUNDARY_NODE[NODE_NUM],
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

  m = 3;
  n = 3 * triangle_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < triangle_num; j++ )
  {
    edge[0+(j               )*m] = triangle_node[0+j*6];
    edge[1+(j               )*m] = triangle_node[3+j*6];
    edge[2+(j               )*m] = triangle_node[1+j*6];

    edge[0+(j+  triangle_num)*m] = triangle_node[1+j*6];
    edge[1+(j+  triangle_num)*m] = triangle_node[4+j*6];
    edge[2+(j+  triangle_num)*m] = triangle_node[2+j*6];

    edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*6];
    edge[1+(j+2*triangle_num)*m] = triangle_node[5+j*6];
    edge[2+(j+2*triangle_num)*m] = triangle_node[0+j*6];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[2+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[2+j*m] );
    edge[0+j*m] = e1;
    edge[2+j*m] = e2;
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

  while ( j < 3 * triangle_num )
  {
    j = j + 1;

    if ( j == 3 * triangle_num )
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

  delete [] edge;

  return node_boundary;
}
