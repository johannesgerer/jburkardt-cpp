# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void assemble_stokes ( int node_num, int element_num, int quad_num,
  int variable_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int element_node[],
  double nu, int ib, double a[], double f[] );
int bandwidth ( int element_order, int element_num, int element_node[],
  int node_num, int node_p_variable[], int node_u_variable[],
  int node_v_variable[] );
void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[],
  double dphidx[], double dphidy[] );
void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[],
  double dphidx[], double dphidy[] );
void boundary_type ( int node_num, double node_xy[], bool node_boundary[],
  int node_type[], int node_u_condition[], int node_v_condition[],
  int node_p_condition[] );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, char *title );
double *dgb_sl ( int n, int ml, int mu, double a_lu[], int pivot[],
  double b[], int job );
void dirichlet_apply ( int node_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int node_p_condition[],
  int node_u_condition[], int node_v_condition[], int variable_num, int ib,
  double a[], double f[] );
void dirichlet_condition ( int n, double xy[], double u_bc[], double v_bc[],
  double p_bc[] );
double *dtable_data_read ( char *input_filename, int m, int n );
void dtable_header_read ( char *input_filename, int *m, int *n );
int file_column_count ( char *input_filename );
void file_name_specification ( int argc, char *argv[], char *node_file_name,
  char *element_file_name );
int file_row_count ( char *input_filename );
int i4_huge ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
int *itable_data_read ( char *input_filename, int m, int n );
void itable_header_read ( char *input_filename, int *m, int *n );
void lvec_print ( int n, bool a[], char *title );
void neumann_apply ( int node_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int node_p_condition[],
  int node_u_condition[], int node_v_condition[], int variable_num,
  double f[] );
void neumann_condition ( int n, double xy[], double u_bc[], double v_bc[],
  double p_bc[] );
void nodes3_write ( char *file_name, int node_num, double node_xy[],
  int node_type[] );
void points_plot ( char *file_name, int node_num, double node_xy[],
  bool node_label );
void pressure3_write ( char *file_name, int node_num, int node_p_variable[],
  int variable_num, double node_c[] );
void quad_rule ( int quad_num, double quad_w[], double quad_xy[] );
double r8_huge ( void );
int r8_nint ( double x );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title );
void reference_to_physical_t6 ( double t[2*6], int n, double ref[],
  double phy[] );
void rhs ( int n, double xy[ ], double u_rhs[], double v_rhs[],
  double p_rhs[] );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, bool *error );
bool s_to_i4vec ( char *s, int n, int ivec[] );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
int s_word_count ( char *s );
void solution_evaluate ( double xy[2], double t[2*3], double node_u[3], double *u,
  double *dudx, double *dudy );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void timestamp ( void );
double triangle_area_2d ( double t[2*3] );
void triangles3_write ( char *file_name, int element_num, int element_node[],
  int node_num, int node3_label[] );
bool *triangulation_order6_boundary_node ( int node_num, int element_num,
  int element_node[] );
void triangulation_order6_plot ( char *file_name, int node_num, double node_xy[],
  int element_num, int element_node[], int node_show, int triangle_show );
void velocity6_write ( char *file_name, int node_num, int node_u_variable[],
  int node_v_variable[], int variable_num, double node_c[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_STOKES.
//
//  Discussion:
//
//    FEM2D_STOKES solves the steady incompressible Stokes equations
//    for velocity vector W and scalar pressure P:
//
//      - nu * Laplacian W(X,Y) + Grad P(X,Y) = F(X,Y)
//
//                                 Div W(X,Y) = G(X,Y)
//
//    in an arbitrary triangulated region in the plane.
//
//    Let U and V denote the scalar components of the velocity vector W.
//
//    Along the boundary of the region, the user controls the type of
//    boundary condition to be imposed, if any.  Currently, these
//    conditions may be of Dirichlet form:
//
//      U(X,Y) = U_BC(X,Y)
//      V(X,Y) = V_BC(X,Y)
//      P(X,Y) = P_BC(X,Y)
//
//    or Neumann form with ZERO right hand side:
//
//      dU/dn(X,Y) = 0
//      dV/dn(X,Y) = 0
//      dP/dn(X,Y) = 0
//
//    The code uses the finite element method.  The Taylor-Hood element
//    is used, in which a single reference element is used to define
//    both a piecewise quadratic representation of velocity, and a piecewise
//    linear representation of pressure.
//
//  Geometry specification:
//
//    The user defines the geometry by supplying two data files
//    which list the node coordinates, and list the nodes that make up
//    each element.
//
//  Equation specification:
//
//    The user specifies
//
//    * the kinematic viscosity NU;
//
//    * the type of boundary conditions imposed:
//
//      void boundary_type ( int node_num, double node_xy[],
//        bool node_boundary[], int node_type[], int node_u_condition[],
//        int node_v_condition[], int node_p_condition[] )
//
//    * the right hand side of any Dirichlet boundary conditions:
//
//      void dirichlet_condition ( int node_num, double node_xy[],
//        double u_bc[], double v_bc[], double p_bc[] )
//
//    * the right hand side of any Neumann boundary conditions:
//      (currently, nonzero values will be ignored.)
//
//      void neumann_condition ( int node_num, double node_xy[],
//        double u_bc[], double v_bc[], double p_bc[] )
//
//    * the right hand side of the Stokes equations
//      by supplying a routine of the form
//
//      void rhs ( int node_num, double node_xy[], double u_rhs[],
//        double v_rhs[], double p_rhs[] )
//
//  Usage:
//
//    fem2d_stokes nodes_file element_file
//
//    invokes the program:
//
//    * "nodes_file", the coordinates of the nodes;
//    * "element_file", the indices of nodes that make up each element.
//
//    Graphics files created include:
//
//    * "nodes6.eps", an image of the nodes;
//    * "triangles6.eps", an image of the elements;
//
//    Data files created include:
//
//    * "nodes3.txt", the nodes associated with pressure;
//    * "triangles3.txt", the linear triangles associated with pressure;
//    * "pressure3.txt", the pressure at the pressure nodes;
//    * "velocity6.txt", the velocity at the nodes in "nodes6_file".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Command Line Parameters:
//
//    Command line argument, char *NODE_FILE_NAME,
//    the name of the node file.  If this argument is not supplied,
//    it will be requested.
//
//    Command line argument, char *ELEMENT_FILE_NAME,
//    the name of the element file.  If this argument is not supplied,
//    it will be requested.
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*VARIABLE_NUM], the VARIABLE_NUM by
//    VARIABLE_NUM coefficient matrix, stored in a banded format.
//
//    Local, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, int ELEMENT_ORDER, the element order.
//
//    Local, double F[VARIABLE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, bool NODE_BOUNDARY[NODE_NUM], is TRUE if the node is
//    found to lie on the boundary of the region.
//
//    Local, double NODE_C[VARIABLE_NUM], the finite element coefficients.
//
//    Local, int NODE_NUM, the number of nodes.
//
//    Local, int NODE_P_CONDITION[NODE_NUM],
//    indicates the condition used to determine pressure at a node.
//    0, there is no condition (and no variable) at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Local, int NODE_P_VARIABLE[NODE_NUM],
//    is the index of the pressure variable associated with the node,
//    or -1 if there is no associated pressure variable.
//
//    Local, int NODE_TYPE[NODE_NUM], determines if the node is a
//    vertex or midside node.
//    1, the node is a vertex (P, U, V variables are associated with it).
//    2, the node is a midside node (only U and V variables are associated.)
//
//    Local, int NODE_U_CONDITION[NODE_NUM],
//    indicates the condition used to determine horizontal velocity at a node.
//    0, there is no condition (and no variable) at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Local, int NODE_U_VARIABLE[NODE_NUM],
//    is the index of the horizontal velocity variable associated with the node,
//    or -1 if there is no associated variable.
//
//    Local, int NODE_V_CONDITION[NODE_NUM],
//    indicates the condition used to determine vertical velocity at a node.
//    0, there is no condition (and no variable) at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Local, int NODE_V_VARIABLE[NODE_NUM],
//    is the index of the vertical velocity variable associated with the node,
//    or -1 if there is no associated variable.
//
//    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Local, int NODE3_NUM, the number of order 3 nodes.
//
//    Local, int NODE3_LABEL[NODE_NUM], contains the renumbered
//    label of order3 nodes, and -1 for nodes that are not order3 nodes.
//
//    Local, double NU, the kinematic viscosity.
//
//    Local, integer QUAD_NUM, the number of quadrature points used for
//    assembly.  This is currently set to 3, the lowest reasonable value.
//    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
//    of QUAD_NUM greater than 3 may be appropriate.
//
//    Local, int VARIABLE_NUM, the number of variables.
//
{
  double *a;
  bool debugging = false;
  int dim_num;
  int element;
  char element_file_name[255];
  int *element_node;
  int element_num;
  int element_order;
  double *f;
  char file_name[255];
  int i;
  int ib;
  int ierr;
  int ip;
  int iu;
  int iv;
  int j;
  int job;
  int neumann_num;
  int node;
  bool *node_boundary;
  double *node_c;
  int *node_condition;
  char node_file_name[255];
  bool node_label;
  int node_num;
  int *node_p_condition;
  int *node_p_variable;
  double *node_r;
  int node_show;
  int *node_type;
  int *node_u_condition;
  int *node_u_variable;
  int *node_v_condition;
  int *node_v_variable;
  double *node_xy;
  int node3_num;
  int *node3_label;
  double nu = 1.0;
  int p_node;
  int *pivot;
  int quad_num = 7;
  int triangle_show;
  int variable;
  int variable_num;

  timestamp ( );
  cout << "\n";
  cout << "FEM2D_STOKES\n";
  cout << "  C++ version:\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Finite element solution of the \n";
  cout << "  steady incompressible Stokes equations\n";
  cout << "  on a triangulated region in 2 dimensions.'\n";
  cout << "\n";
  cout << "  - nu * ( Uxx + Uyy ) + dPdx = F1(x,y)\n";
  cout << "  - nu * ( Vxx + Vyy ) + dPdy = F2(x,y)\n";
  cout << "                      Ux + Vy = F3(x,y).\n";
  cout << "\n";
  cout << "  Boundary conditions may be of Dirichlet type:\n";
  cout << "\n";
  cout << "    U(x,y) = U_BC(x,y)\n";
  cout << "    V(x,y) = V_BC(x,y)\n";
  cout << "    P(x,y) = P_BC(x,y)\n";
  cout << "\n";
  cout << "  or of Neumann type with zero right hand side:\n";
  cout << "\n";
  cout << "    dU/dn(x,y) = 0\n";
  cout << "    dV/dn(x,y) = 0\n";
  cout << "    dP/dn(x,y) = 0\n";
  cout << "\n";
  cout << "  The finite element method uses Taylor-Hood\n";
  cout << "  triangular elements which are linear for pressure\n";
  cout << "  and quadratic for velocity.\n";
  cout << "\n";
  cout << "  Quadrature order =           " << quad_num << "\n";
  cout << "  The kinematic viscosity NU = " << nu << "\n";
  cout << "\n";
  cout << "  Current status:\n";
  cout << "\n";
  cout << "  * testing zero Neumann condition option.\n";
//
//  Get the file names.
//
  file_name_specification ( argc, argv, node_file_name, element_file_name );

  cout << "\n";
  cout << "  Node file is \"" << node_file_name << "\".\n";
  cout << "  Element file is \"" << element_file_name << "\".\n";
//
//  Read the node coordinate file.
//
  dtable_header_read ( node_file_name, &dim_num, &node_num );

  cout << "  Number of nodes =          " << node_num << "\n";

  node_p_condition = new int[node_num];
  node_p_variable = new int[node_num];
  node_type = new int[node_num];
  node_u_condition = new int[node_num];
  node_u_variable = new int[node_num];
  node_v_condition = new int[node_num];
  node_v_variable = new int[node_num];

  node_xy = dtable_data_read ( node_file_name, dim_num, node_num );

  r8mat_transpose_print_some ( dim_num, node_num, node_xy, 1, 1, 2, 10,
    "  First 10 nodes" );
//
//  Read the element description file.
//
  itable_header_read ( element_file_name, &element_order, &element_num );

  cout << "\n";
  cout << "  Element order =            " << element_order << "\n";
  cout << "  Number of elements =       " << element_num << "\n";

  if ( element_order != 6 )
  {
    cout << "\n";
    cout << "FEM2D_STOKES - Fatal error!\n";
    cout << "  The input triangulation has order " << element_order << "\n";
    cout << "  However, a triangulation of order 6 is required.\n";
    exit ( 1 );
  }
  element_node = itable_data_read ( element_file_name, element_order,
    element_num );

  i4mat_transpose_print_some ( element_order, element_num,
    element_node, 1, 1, element_order, 10, "  First 10 elements" );
//
//  Determine the "type" of each node.
//  A vertex node, of type 1, has U, V, and P variables.
//  A midside node, of type 2, has U and V only.
//
  for ( node = 0; node < node_num; node++ )
  {
    node_type[node] = 1;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( j = 4; j <= 6; j++ )
    {
      node = element_node[(j-1)+element*6];
      node_type[node-1] = 2;
    }
  }
//
//  Determine which nodes are boundary nodes.
//
  node_boundary = triangulation_order6_boundary_node ( node_num,
    element_num, element_node );

  if ( false )
  {
    lvec_print ( node_num, node_boundary, "    Node  Boundary?" );
  }
//
//  Determine the node conditions:
//  For now, we'll just assume all boundary nodes are Dirichlet.
//
//  All conditions begin as finite element conditions.
//
  for ( node = 0; node < node_num; node++ )
  {
    node_p_condition[node] = 1;
    node_u_condition[node] = 1;
    node_v_condition[node] = 1;
  }
//
//  Conditions on velocities associated with a boundary node are Dirichlet
//  conditions.
//
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_boundary[node] )
    {
      node_u_condition[node] = 2;
      node_v_condition[node] = 2;
    }
  }
//
//  Midside nodes have no associated pressure variable.
//
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_type[node] == 2 )
    {
      node_p_condition[node] = 0;
    }
  }
//
//  Replace a single finite element pressure condition by a Dirichlet
//  condition.
//
  p_node = -1;

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_p_condition[node] == 1 )
    {
      node_p_condition[node] = 2;
      p_node = node + 1;
      break;
    }
  }

  if ( p_node == -1 )
  {
    cout << "\n";
    cout << "FEM2D_STOKES - Fatal error!\n";
    cout << "  Unable to find a finite element pressure condition\n";
    cout << "  suitable for replacement by a Dirichlet condition.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Dirichlet boundary condition on pressure\n";
  cout << "  will be applied at node " << p_node << "\n";
//
//  Allow the user to examine and modify the tentative boundary conditions.
//
  boundary_type ( node_num, node_xy, node_boundary, node_type,
    node_u_condition, node_v_condition, node_p_condition );

  neumann_num = 0;

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_u_condition[node] == 3 )
    {
      neumann_num = neumann_num + 1;
    }

    if ( node_v_condition[node] == 3 )
    {
      neumann_num = neumann_num + 1;
    }

    if ( node_p_condition[node] == 3 )
    {
      neumann_num = neumann_num + 1;
    }
  }

  cout << "\n";
  cout << "  Number of Neumann conditions added = " << neumann_num << "\n";

  cout << "\n";
  cout << "  Boundary conditions per node:\n";
  cout << "\n";
  cout << "      Node    U_cond    V_cond    P_cond\n";
  cout << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(8) << node
         << "  " << setw(8) << node_u_condition[node]
         << "  " << setw(8) << node_v_condition[node]
         << "  " << setw(8) << node_p_condition[node] << "\n";
  }
//
//  Number the variables.
//
  variable_num = 0;

  for ( node = 0; node < node_num; node++ )
  {
    variable_num = variable_num + 1;
    node_u_variable[node] = variable_num;

    variable_num = variable_num + 1;
    node_v_variable[node] = variable_num;

    if ( node_type[node] == 1 )
    {
      variable_num = variable_num + 1;
      node_p_variable[node] = variable_num;
    }
    else
    {
      node_p_variable[node] = -1;
    }
  }

  cout << "\n";
  cout << "  Total number of variables is " << variable_num << "\n";

  cout << "\n";
  cout << "  Variable indices per node:\n";
  cout << "\n";
  cout << "      Node   U_index   V_index   P_index\n";
  cout << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(8) << node+1
         << "  " << setw(8) << node_u_variable[node]
         << "  " << setw(8) << node_v_variable[node]
         << "  " << setw(8) << node_p_variable[node] << "\n";
  }
//
//  Determine the bandwidth of the coefficient matrix.
//
  ib = bandwidth ( element_order, element_num, element_node,
    node_num, node_p_variable, node_u_variable, node_v_variable );

  cout << "\n";
  cout << "  The matrix half bandwidth is " << ib << "\n";
  cout << "  The matrix bandwidth is      " << 2 * ib + 1 << "\n";
  cout << "  The storage bandwidth is     " << 3 * ib + 1 << "\n";
//
//  Plot the nodes.
//
  if ( node_num <= 100 )
  {
    strcpy ( file_name, "nodes6.eps" );
    node_label = true;

    points_plot ( file_name, node_num, node_xy, node_label );

    cout << "\n";
    cout << "  Order 6 nodes plotted in \"" << file_name << "\".\n";
  }
//
//  Plot the triangles.
//
  if ( node_num <= 100 )
  {
    strcpy ( file_name, "triangles6.eps" );
    node_show = 2;
    triangle_show = 2;

    triangulation_order6_plot ( file_name, node_num,
      node_xy, element_num, element_node, node_show, triangle_show );

    cout << "\n";
    cout << "  Order 6 triangles plotted in \"" << file_name << "\".\n";
  }
//
//  Allocate space for the coefficient matrix A and right hand side F.
//
  a = new double[(3*ib+1)*variable_num];
  f = new double[variable_num];
  node_r = new double[variable_num];
  node_c = new double[variable_num];
  pivot = new int[variable_num];
//
//  Assemble the finite element coefficient matrix A and the right-hand side F.
//
  assemble_stokes ( node_num, element_num, quad_num,
    variable_num, node_xy, node_p_variable, node_u_variable,
    node_v_variable, element_node, nu, ib, a, f );
//
//  Print a tiny portion of the matrix.
//
  if ( debugging )
  {
    dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 20, 20,
      "  Part of Stokes matrix::" );

    r8vec_print_some ( variable_num, f, 1, 20,
      "  Part of Stokes right hand side:" );
  }
//
//  Adjust the linear system to account for Dirichlet boundary conditions.
//
  dirichlet_apply ( node_num, node_xy, node_p_variable,
    node_u_variable, node_v_variable, node_p_condition,
    node_u_condition, node_v_condition, variable_num, ib, a, f );

  if ( debugging )
  {
    dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 20, 20,
      "  Part of Stokes matrix, adjusted for Dirichlet BC:" );

    r8vec_print_some ( variable_num, f, 1, 20,
      "  Part of Stokes right hand side, adjusted for Dirichlet BC:" );
  }
//
//  Adjust the linear system to account for Neumann boundary conditions.
//
  neumann_apply ( node_num, node_xy, node_p_variable,
    node_u_variable, node_v_variable, node_p_condition,
    node_u_condition, node_v_condition, variable_num, f );

  if ( false )
  {
    dgb_print_some ( variable_num, variable_num, ib, ib, a, 1, 1, 20, 20,
      "  Part of Stokes matrix, adjusted for Neumann BC:" );

    r8vec_print_some ( variable_num, f, 1, 20,
      "  Part of Stokes right hand side, adjusted for Neumann BC:" );
  }
//
//  Solve the linear system using a banded solver.
//
  ierr = dgb_fa ( variable_num, ib, ib, a, pivot );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "FEM2D_STOKES - Fatal error!\n";
    cout << "  DGB_FA returned the error condition IERR = " << ierr << ".\n";
    cout << "\n";
    cout << "  The linear system was not factored, and the\n";
    cout << "  algorithm cannot proceed.\n";
    exit ( 1 );
  }

  job = 0;

  delete [] node_c;

  node_c = dgb_sl ( variable_num, ib, ib, a, pivot, f, job );

  if ( debugging )
  {
    r8vec_print_some ( variable_num, node_c, 1, 10,
      "  Part of the solution vector:" );
  }
//
//  Print the solution vector based at nodes.
//
  cout << "\n";
  cout << "  Solution values associated with nodes:\n";
  cout << "\n";
  cout << "      Node               U               V               P\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    iu = node_u_variable[node] - 1;
    iv = node_v_variable[node] - 1;
    ip = node_p_variable[node] - 1;

    cout << "  " << setw(8)  << node
         << "  " << setw(14) << node_c[iu]
         << "  " << setw(14) << node_c[iv];

    if ( 0 <= ip )
    {
      cout << "  " << setw(14) << node_c[ip] << "\n";
    }
    else
    {
      cout << "\n";
    }
  }
//
//  Compute a renumbering of the pressure nodes.
//
  node3_num = 0;

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_type[node] == 1 )
    {
      node3_num = node3_num + 1;
    }
  }
  node3_label = new int[node_num];

  node3_num = 0;

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_type[node] == 1 )
    {
      node3_num = node3_num + 1;
      node3_label[node] = node3_num;
    }
    else
    {
      node3_label[node] = -1;
    }
  }
//
//  Write the pressure nodes to a file.
//
  strcpy ( file_name, "nodes3.txt" );

  nodes3_write ( file_name, node_num, node_xy, node_type );

  cout << "\n";
  cout << "  Pressure nodes written to \"" << file_name << "\".\n";
//
//  Write the pressure triangles to a file.
//
  strcpy ( file_name, "triangles3.txt" );

  triangles3_write ( file_name, element_num, element_node,
    node_num, node3_label );

  delete [] node3_label;

  cout << "\n";
  cout << "  Pressure triangles written to \"" << file_name << "\".\n";
//
//  Write the pressures to a file.
//
  strcpy ( file_name, "pressure3.txt" );

  pressure3_write ( file_name, node_num, node_p_variable,
    variable_num, node_c );

  cout << "\n";
  cout << "  Pressures written to \"" << file_name << "\".\n";
//
//  Write the velocities to a file.
//
  strcpy ( file_name, "velocity6.txt" );

  velocity6_write ( file_name, node_num, node_u_variable,
    node_v_variable, variable_num, node_c );

  cout << "\n";
  cout << "  Velocities written to \"" << file_name << "\".\n";
//
//  Deallocate memory.
//
  delete [] a;
  delete [] f;
  delete [] node_boundary;
  delete [] node_c;
  delete [] node_p_condition;
  delete [] node_p_variable;
  delete [] node_r;
  delete [] node_type;
  delete [] node_u_condition;
  delete [] node_u_variable;
  delete [] node_v_condition;
  delete [] node_v_variable;
  delete [] node_xy;
  delete [] pivot;
  delete [] element_node;
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_STOKES:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void assemble_stokes ( int node_num, int element_num, int quad_num,
  int variable_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int element_node[],
  double nu, int ib, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_STOKES assembles the system for the Stokes equations.
//
//  Discussion:
//
//    The matrix is known to be banded.  A special matrix storage format
//    is used to reduce the space required.  Details of this format are
//    discussed in the routine DGB_FA.
//
//    The Stokes equations in weak form are:
//
//      Integral ( nu * ( dBdx(I) * dUdx + dBdy(I) * dUdy )
//        + B(I) * ( dPdx - U_RHS ) ) = 0
//
//      Integral ( nu * ( dBdx(I) * dVdx + dBdy(I) * dVdy )
//        + B(I) * (  dPdy - V_RHS ) ) = 0
//
//      Integral ( Q(I) * ( dUdx + dVdy - P_RHS ) ) = 0
//
//    Once the basic finite element system is set up by this routine, another
//    routine adjusts the system to account for boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, int VARIABLE_NUM, the number of unknowns.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
//
//    Input, int NODE_P_VARIABLE[NODE_NUM], the index of the pressure
//    variable associated with a node, or -1 if there is none.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM], the index of the horizontal
//    velocity variable associated with a node, or -1 if there is none.
//
//    Input, int NODE_V_VARIABLE[NODE_NUM], the index of the vertical
//    velocity variable associated with a node, or -1 if there is none.
//
//    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that form each
//    element.  Nodes 1, 2, and 3 are the vertices.  Node 4 is between 1
//    and 2, and so on.
//
//    Input, double NU, the kinematic viscosity.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Output, double A[(3*IB+1)*VARIABLE_NUM], the VARIABLE_NUM by VARIABLE_NUM
//    coefficient matrix, stored in a banded format.
//
//    Output, double F[VARIABLE_NUM], the right hand side.
//
//  Local Parameters:
//
//    Local, double B[6*QUAD_NUM], DBDX[6*QUAD_NUM],
//    DBDY[6*QUAD_NUM], the values of the quadratic basis functions
//    and first derivatives at the quadrature points.
//
//    Local, double Q[3*QUAD_NUM], DQDX[3*QUAD_NUM],
//    DQDY[3*QUAD_NUM], the values of the linear basis functions
//    and first derivatives at the quadrature points.
//
//    Local, double QUAD_W[QUAD_NUM], quadrature weights.
//
//    Local, double QUAD_XY[2*QUAD_NUM], the quadrature points.
//
{
  double area;
  double *b;
  int bw;
  double *dbdx;
  double *dbdy;
  double *dqdx;
  double *dqdy;
  int element;
  int i;
  int *ip;
  int *iu;
  int *iv;
  int j;
  int node;
  double *p_rhs;
  double *q;
  int quad;
  double *quad_w;
  double *quad_xy;
  double t3[2*3];
  double t6[2*6];
  double *u_rhs;
  double *v_rhs;
  int variable;
  double *w;
  double *xy;
//
//  Allocate temporary arrays.
//
  b = new double[6*quad_num];
  dbdx = new double[6*quad_num];
  dbdy = new double[6*quad_num];
  dqdx = new double[3*quad_num];
  dqdy = new double[3*quad_num];
  ip = new int[3];
  iu = new int[6];
  iv = new int[6];
  p_rhs = new double[quad_num];
  q = new double[3*quad_num];
  quad_w = new double[quad_num];
  quad_xy = new double[2*quad_num];
  u_rhs = new double[quad_num];
  v_rhs = new double[quad_num];
  w = new double[quad_num];
  xy = new double[2*quad_num];

  bw = 3 * ib + 1;
//
//  Initialize the arrays to zero.
//
  for ( variable = 0; variable < variable_num; variable++ )
  {
    f[variable] = 0.0;
  }
  for ( variable = 0; variable < variable_num; variable++ )
  {
    for ( i = 0; i < bw; i++ )
    {
      a[i+variable*bw] = 0.0;
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
//  Make copies of the T3 and T6 triangles.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t3[i+j*2] = node_xy[i+(element_node[j+element*6]-1)*2];
      }
    }
    for ( j = 0; j < 6; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        t6[i+j*2] = node_xy[i+(element_node[j+element*6]-1)*2];
      }
    }
//
//  Map the quadrature points QUAD_XY to points XY in the physical triangle.
//
    reference_to_physical_t6 ( t6, quad_num, quad_xy, xy );
    area = triangle_area_2d ( t3 );
    for ( quad = 0; quad < quad_num; quad++ )
    {
      w[quad] = quad_w[quad] * area;
    }
    rhs ( quad_num, xy, u_rhs, v_rhs, p_rhs );
//
//  Evaluate the basis functions at the quadrature points.
//
    basis_mn_t6 ( t6, quad_num, xy, b, dbdx, dbdy );
    basis_mn_t3 ( t3, quad_num, xy, q, dqdx, dqdy );
//
//  Extract the indices of the finite element coefficients for this element.
//
    for ( j = 0; j < 6; j++ )
    {
      iu[j] = node_u_variable[element_node[j+element*6]-1];
      iv[j] = node_v_variable[element_node[j+element*6]-1];
    }
    for ( j = 0; j < 3; j++ )
    {
      ip[j] = node_p_variable[element_node[j+element*6]-1];
    }
//
//  C++ doesn't offer us a chance to vectorize over QUAD,
//  so we have to handle each component separately.
//
//  This means that much of the effort made above to extract
//  all the local information into arrays can't actually be
//  used effectively to speed up the calculation.
//
//  However, we keep to the plan, in order to synchronize with
//  the FORTRAN90 and MATLAB codes, which can vectorize.
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
//
//  The horizontal momentum equation.
//
      for ( i = 0; i < 6; i++ )
      {
        f[iu[i]-1] = f[iu[i]-1] + w[quad] * u_rhs[quad] * b[i+quad*6];

        for ( j = 0; j < 6; j++ )
        {
          a[iu[i]-iu[j]+2*ib+(iu[j]-1)*bw]
        = a[iu[i]-iu[j]+2*ib+(iu[j]-1)*bw] + w[quad] *
          nu * ( dbdx[j+quad*6] * dbdx[i+quad*6]
               + dbdy[j+quad*6] * dbdy[i+quad*6] );
        }

        for ( j = 0; j < 3; j++ )
        {
          a[iu[i]-ip[j]+2*ib+(ip[j]-1)*bw]
        = a[iu[i]-ip[j]+2*ib+(ip[j]-1)*bw] + w[quad] *
          dqdx[j+quad*3] * b[i+quad*6];
        }
      }
//
//  The vertical momementum equation.
//
      for ( i = 0; i < 6; i++ )
      {
        f[iv[i]-1] = f[iv[i]-1] + w[quad] * v_rhs[quad] * b[i+quad*6];

        for ( j = 0; j < 6; j++ )
        {
          a[iv[i]-iv[j]+2*ib+(iv[j]-1)*bw]
        = a[iv[i]-iv[j]+2*ib+(iv[j]-1)*bw] + w[quad] *
          nu * ( dbdx[j+quad*6] * dbdx[i+quad*6]
               + dbdy[j+quad*6] * dbdy[i+quad*6] );
        }

        for ( j = 0; j < 3; j++ )
        {
          a[iv[i]-ip[j]+2*ib+(ip[j]-1)*bw]
        = a[iv[i]-ip[j]+2*ib+(ip[j]-1)*bw] + w[quad] *
          dqdy[j+quad*3] * b[i+quad*6];
        }
      }
//
//  The continuity equation.
//
      for ( i = 0; i < 3; i++ )
      {
        f[ip[i]-1] = f[ip[i]-1] + w[quad] * p_rhs[quad] * q[i+quad*3];

        for ( j = 0; j < 6; j++ )
        {
          a[ip[i]-iu[j]+2*ib+(iu[j]-1)*bw]
        = a[ip[i]-iu[j]+2*ib+(iu[j]-1)*bw] + w[quad] *
          dbdx[j+quad*6] * q[i+quad*3];

          a[ip[i]-iv[j]+2*ib+(iv[j]-1)*bw]
        = a[ip[i]-iv[j]+2*ib+(iv[j]-1)*bw] + w[quad] *
          dbdy[j+quad*6] * q[i+quad*3];
        }
      }
    }
  }
//
//  Delete temporary arrays.
//
  delete [] b;
  delete [] dbdx;
  delete [] dbdy;
  delete [] dqdx;
  delete [] dqdy;
  delete [] ip;
  delete [] iu;
  delete [] iv;
  delete [] p_rhs;
  delete [] q;
  delete [] quad_w;
  delete [] quad_xy;
  delete [] u_rhs;
  delete [] v_rhs;
  delete [] w;
  delete [] xy;

  return;
}
//****************************************************************************80

int bandwidth ( int element_order, int element_num, int element_node[],
  int node_num, int node_p_variable[], int node_u_variable[],
  int node_v_variable[] )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth of the coefficient matrix.
//
//  Discussion:
//
//    We take the bandwidth to be the maximum difference between the
//    indices of two variables associated with nodes that share a element.
//
//    Therefore, we can compute the bandwidth by examining each element,
//    and finding the maximum difference in indices of any two variables
//    associated with nodes in that element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the number of nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_P_VARIABLE[NODE_NUM],
//    is the index of the pressure variable associated with the node,
//    or -1 if there is no associated pressure variable.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM],
//    is the index of the horizontal velocity variable associated with the node.
//
//    Input, int NODE_V_VARIABLE[NODE_NUM],
//    is the index of the vertical velocity variable associated with the node.
//
//    Output, integer BANDWIDTH, the half bandwidth of the matrix.
//
{
  int element;
  int ib;
  int local;
  int node;
  int v;
  int v_max;
  int v_min;

  ib = 0;

  for ( element = 0; element < element_num; element++ )
  {
    v_max = -i4_huge ( );
    v_min =  i4_huge ( );

    for ( local = 0; local < element_order; local++ )
    {
      node = element_node[local+element*6];

      v = node_u_variable[node-1];
      v_max = i4_max ( v_max, v );
      v_min = i4_min ( v_min, v );

      v = node_v_variable[node-1];
      v_max = i4_max ( v_max, v );
      v_min = i4_min ( v_min, v );

      if ( 0 < node_p_variable[node-1] )
      {
        v = node_p_variable[node-1];
        v_max = i4_max ( v_max, v );
        v_min = i4_min ( v_min, v );
      }
    }
    ib = i4_max ( ib, v_max - v_min );
  }
  return ib;
}
//****************************************************************************80

void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[],
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T3: all bases at N points for a T3 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a triangle.
//    It works directly with these coordinates, and does not refer to a
//    reference element.
//
//    The sides of the triangle DO NOT have to lie along a coordinate
//    axis.
//
//    The routine evaluates the basis functions associated with each vertex,
//    and their derivatives with respect to X and Y.
//
//  Physical Element T3:
//
//            3
//           / \
//          /   \
//         /     \
//        /       \
//       1---------2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices
//    of the triangle.  It is common to list these points in counter clockwise
//    order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the points where the basis functions
//    are to be evaluated.
//
//    Output, double PHI[3*N], the value of the basis functions
//    at the evaluation points.
//
//    Output, double DPHIDX[3*N], DPHIDY[3*N], the value of the
//    derivatives at the evaluation points.
//
//  Local parameters:
//
//    Local, double AREA, is (twice) the area of the triangle.
//
{
  double area;
  int i;
  int j;
  double temp;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] )
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] )
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    cout << "\n";
    cout << "BASIS_MN_T3 - Fatal error!\n";
    cout << "  Element has zero area.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*3] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*3] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*3] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*3] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*3] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*3] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*3] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*3] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*3] =      ( t[0+1*2] - t[0+0*2] );
  }
//
//  Normalize.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*3]    = phi[i+j*3]    / area;
      dphidx[i+j*3] = dphidx[i+j*3] / area;
      dphidy[i+j*3] = dphidy[i+j*3] / area;
    }
  }

  return;
}
//****************************************************************************80

void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[],
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T6: all bases at N points for a T6 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices and midside
//    nodes of a triangle.  It works directly with these coordinates, and does
//    not refer to a reference element.
//
//    This routine requires that the midside nodes be "in line"
//    with the vertices, that is, that the sides of the triangle be
//    straight.  However, the midside nodes do not actually have to
//    be halfway along the side of the triangle.
//
//  The physical element T6:
//
//    This picture indicates the assumed ordering of the six nodes
//    of the triangle.
//
//    |
//    |
//    |        3
//    |       / \
//    |      /   \
//    Y     6     5
//    |    /       \
//    |   /         \
//    |  1-----4-----2
//    |
//    +--------X-------->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the nodal oordinates of the element.
//    It is common to list these points in counter clockwise order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the coordinates of the point where
//    the basis functions are to be evaluated.
//
//    Output, double PHI[6*N], the value of the basis functions at P.
//
//    Output, double DPHIDX[6*N], DPHIDY[6*N], the value of the X
//    and Y derivatives of the basis functions at P.
//
{
  double gn;
  double gx;
  double hn;
  double hx;
  int j;

  for ( j = 0; j < n; j++ )
  {
//
//  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] )
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] )
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] )
       - ( t[0+5*2] - t[0+3*2] ) * ( p[1+j*2] - t[1+3*2] );

    hn = ( t[0+0*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] )
       - ( t[0+5*2] - t[0+3*2] ) * ( t[1+0*2] - t[1+3*2] );

    phi[0+j*6] =     ( gx * hx ) / ( gn * hn );
    dphidx[0+j*6] =  (      ( t[1+2*2] - t[1+1*2] ) * hx
                     + gx * ( t[1+5*2] - t[1+3*2] ) ) / ( gn * hn );
    dphidy[0+j*6] = -(      ( t[0+2*2] - t[0+1*2] ) * hx
                     + gx * ( t[0+5*2] - t[0+3*2] ) ) / ( gn * hn );
//
//  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
//
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] )
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] )
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] )
       - ( t[0+3*2] - t[0+4*2] ) * ( p[1+j*2] - t[1+4*2] );

    hn = ( t[0+1*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] )
       - ( t[0+3*2] - t[0+4*2] ) * ( t[1+1*2] - t[1+4*2] );

    phi[1+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[1+j*6] =  (      ( t[1+2*2] - t[1+0*2] ) * hx
                     + gx * ( t[1+3*2] - t[1+4*2] ) ) / ( gn * hn );
    dphidy[1+j*6] = -(      ( t[0+2*2] - t[0+0*2] ) * hx
                     + gx * ( t[0+3*2] - t[0+4*2] ) ) / ( gn * hn );
//
//  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] )
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] )
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] )
       - ( t[0+4*2] - t[0+5*2] ) * ( p[1+j*2] - t[1+5*2] );

    hn = ( t[0+2*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] )
       - ( t[0+4*2] - t[0+5*2] ) * ( t[1+2*2] - t[1+5*2] );

    phi[2+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[2+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx
                     + gx * ( t[1+4*2] - t[1+5*2] ) ) / ( gn * hn );
    dphidy[2+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx
                     + gx * ( t[0+4*2] - t[0+5*2] ) ) / ( gn * hn );
//
//  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
//
    gx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] )
       - ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    gn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] )
       - ( t[0+0*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    hx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] )
       - ( t[0+1*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    hn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] )
       - ( t[0+1*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    phi[3+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[3+j*6] =  (      ( t[1+0*2] - t[1+2*2] ) * hx
                     + gx * ( t[1+1*2] - t[1+2*2] ) ) / ( gn * hn );
    dphidy[3+j*6] = -(      ( t[0+0*2] - t[0+2*2] ) * hx
                     + gx * ( t[0+1*2] - t[0+2*2] ) ) / ( gn * hn );
//
//  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
//
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] )
       - ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] )
       - ( t[0+1*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] )
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    hn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] )
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    phi[4+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[4+j*6] =  (      ( t[1+1*2] - t[1+0*2] ) * hx
                     + gx * ( t[1+2*2] - t[1+0*2] ) ) / ( gn * hn );
    dphidy[4+j*6] = -(      ( t[0+1*2] - t[0+0*2] ) * hx
                     + gx * ( t[0+2*2] - t[0+0*2] ) ) / ( gn * hn );
//
//  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] )
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] )
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] )
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    hn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] )
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    phi[5+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[5+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx
                     + gx * ( t[1+2*2] - t[1+1*2] ) ) / ( gn * hn );
    dphidy[5+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx
                     + gx * ( t[0+2*2] - t[0+1*2] ) ) / ( gn * hn );
  }

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
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
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
  int jlo, int ihi, int jhi, char *title )

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
//    Input, char *TITLE, a title to print.
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

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
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
//    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
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

void dirichlet_apply ( int node_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int node_p_condition[],
  int node_u_condition[], int node_v_condition[], int variable_num, int ib,
  double a[], double f[] )

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
//    06 February 2006
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
//    Input, int NODE_P_VARIABLE[NODE_NUM],
//    is the index of the pressure variable associated with the node,
//    or -1 if there is no associated pressure variable.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM],
//    is the index of the horizontal velocity variable associated with the node.
//
//    Input, int NODE_V_VARIABLE[NODE_NUM],
//    is the index of the vertical velocity variable associated with the node.
//
//    Input, int NODE_P_CONDITION[NODE_NUM],
//    indicates the condition used to determine pressure at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int NODE_U_CONDITION[NODE_NUM],
//    indicates the condition used to determine horizontal velocity at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int NODE_V_CONDITION[NODE_NUM],
//    indicates the condition used to determine vertical velocity at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input/output, double A[(3*IB+1)*VARIABLE_NUM], the VARIABLE_NUM by
//    VARIABLE_NUM coefficient matrix, stored in a compressed format; on output,
//    the matrix has been adjusted for Dirichlet boundary conditions.
//
//    Input/output, double F[VARIABLE_NUM], the right hand side.
//    On output, the right hand side has been adjusted for Dirichlet
//    boundary conditions.
//
{
  int column;
  int column_high;
  int column_low;
  int DIRICHLET = 2;
  int ip;
  int iu;
  int iv;
  int node;
  double *p_bc;
  double *u_bc;
  double *v_bc;
  double value;

  u_bc = new double[node_num];
  v_bc = new double[node_num];
  p_bc = new double[node_num];

  dirichlet_condition ( node_num, node_xy, u_bc, v_bc, p_bc );

  for ( node = 0; node < node_num; node++ )
  {
    iu = node_u_variable[node];
    iv = node_v_variable[node];
    ip = node_p_variable[node];

    if ( node_u_condition[node] == DIRICHLET )
    {
      column_low = i4_max ( iu - ib, 1 );
      column_high = i4_min ( iu + ib, variable_num );

      for ( column = column_low; column <= column_high; column++ )
      {
        a[iu-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
      }
      a[2*ib+(iu-1)*(3*ib+1)] = 1.0;

      f[iu-1] = u_bc[node];
    }
    if ( node_v_condition[node] == DIRICHLET )
    {
      column_low = i4_max ( iv - ib, 1 );
      column_high = i4_min ( iv + ib, variable_num );

      for ( column = column_low; column <= column_high; column++ )
      {
        a[iv-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
      }
      a[2*ib+(iv-1)*(3*ib+1)] = 1.0;

      f[iv-1] = v_bc[node];
    }
    if ( 0 < ip )
    {
      if ( node_p_condition[node] == DIRICHLET )
      {
        column_low = i4_max ( ip - ib, 1 );
        column_high = i4_min ( ip + ib, variable_num );

        for ( column = column_low; column <= column_high; column++ )
        {
          a[ip-column+2*ib+(column-1)*(3*ib+1)] = 0.0;
        }
        a[2*ib+(ip-1)*(3*ib+1)] = 1.0;

        f[ip-1] = p_bc[node];
      }
    }
  }
  delete [] u_bc;
  delete [] v_bc;
  delete [] p_bc;

  return;
}
//****************************************************************************80

double *dtable_data_read ( char *input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_READ reads the data from a real TABLE file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
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
//    27 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double DTABLE_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "DTABLE_DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    input.getline ( line, sizeof ( line ) );

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

void dtable_header_read ( char *input_filename, int *m, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_HEADER_READ reads the header from a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

int file_column_count ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
//    by spaces.  There may also be some blank lines, and some comment lines,
//    which have a "#" in column 1.
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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char line[256];
//
//  Open the file.
//
  input.open ( input_filename );

  if ( !input )
  {
    column_num = -1;
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Fatal error!\n";
    cout << "  Could not open the file:\n";
    cout << "  \"" << input_filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = true;
    break;

  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( input_filename );

    for ( ; ; )
    {
      input.getline ( line, sizeof ( line ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
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
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Warning!\n";
    cout << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( line );

  return column_num;
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
//  If at least two command line arguments, the second is the triangulation file.
//
  if ( 2 <= argc )
  {
    strcpy ( element_file_name, argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "FILE_NAME_SPECIFICATION:\n";
    cout << "  Please enter the name of the triangulation file.\n";

    cin.getline ( element_file_name, sizeof ( element_file_name ) );
  }

  return;
}
//****************************************************************************80

int file_row_count ( char *input_filename )

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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  char line[100];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "FILE_ROW_COUNT - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return (-1);
  }

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

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

int i4_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
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

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

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
//    Input, char *TITLE, a title for the matrix.
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
      cout << setw(8) << i << "  ";
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
        cout << setw(8) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

int *itable_data_read ( char *input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ITABLE_DATA_READ reads data from an integer TABLE file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
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
//    11 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int ITABLE_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  char line[255];
  int *table;
  int *x;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "ITABLE_DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new int[m*n];

  x = new int[m];

  j = 0;

  while ( j < n )
  {
    input.getline ( line, sizeof ( line ) );

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

void itable_header_read ( char *input_filename, int *m, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    ITABLE_HEADER_READ reads the header from an integer TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cout << "\n";
    cout << "ITABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cout << "\n";
    cout << "ITABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

void lvec_print ( int n, bool a[], char *title )

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
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(8) << i + 1 << "  "
         << setw(1) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void neumann_apply ( int node_num, double node_xy[], int node_p_variable[],
  int node_u_variable[], int node_v_variable[], int node_p_condition[],
  int node_u_condition[], int node_v_condition[], int variable_num,
  double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    NEUMANN_APPLY accounts for Neumann boundary conditions.
//
//  Discussion:
//
//    At the moment, this program only allows Neumann boundary conditions
//    of the form
//
//      dU/dn = 0
//      dV/dn = 0
//      dP/dn = 0
//
//    For such conditions, there is NO change necessary to the linear system.
//    So this routine actually does nothing.  It is here as preparation
//    for later treatment of nonzero Neumann conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2006
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
//    Input, int NODE_P_VARIABLE[NODE_NUM],
//    is the index of the pressure variable associated with the node,
//    or -1 if there is no associated pressure variable.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM],
//    is the index of the horizontal velocity variable associated with the node.
//
//    Input, int NODE_V_VARIABLE[NODE_NUM],
//    is the index of the vertical velocity variable associated with the node.
//
//    Input, int NODE_P_CONDITION[NODE_NUM],
//    indicates the condition used to determine pressure at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int NODE_U_CONDITION[NODE_NUM],
//    indicates the condition used to determine horizontal velocity at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int NODE_V_CONDITION[NODE_NUM],
//    indicates the condition used to determine vertical velocity at a node.
//    0, there is no condition at this node.
//    1, a finite element equation is used;
//    2, a Dirichlet condition is used.
//    3, a Neumann condition is used.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input/output, double F[VARIABLE_NUM], the right hand side.
//    On output, the right hand side has been adjusted for Neumann
//    boundary conditions.
//
{
  int ip;
  int iu;
  int iv;
  int NEUMANN = 3;
  int node;
  double *p_bc;
  double *u_bc;
  double *v_bc;
  double value;
//
//  The user routine supplies a right hand side value for a possible
//  Neumann condition at EVERY node.
//
  u_bc = new double[node_num];
  v_bc = new double[node_num];
  p_bc = new double[node_num];

  neumann_condition ( node_num, node_xy, u_bc, v_bc, p_bc );

  for ( node = 0; node < node_num; node++ )
  {
    iu = node_u_variable[node];
    iv = node_v_variable[node];
    ip = node_p_variable[node];

    if ( node_u_condition[node] == NEUMANN )
    {
//    f[iu-1] = f[iu-1] + line integral;
    }
    if ( node_v_condition[node] == NEUMANN )
    {
//    f[iv-1] = f[iv-1] + line integral;
    }
    if ( 0 < ip )
    {
      if ( node_p_condition[node] == NEUMANN )
      {
//      f[ip-1] = f[ip-1] + line integral;
      }
    }
  }
  delete [] u_bc;
  delete [] v_bc;
  delete [] p_bc;

  return;
}
//****************************************************************************80

void nodes3_write ( char *file_name, int node_num, double node_xy[],
  int node_type[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODES3_WRITE writes the pressure nodes to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the file name.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int NODE_TYPE[NODE_NUM], determines if the node is a
//    vertex or midside node.
//    1, the node is a vertex (P, U, V variables are associated with it).
//    2, the node is a midside node (only U and V variables are associated.)
//
{
  ofstream file_unit;
  int node;

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "NODES3_WRITE - Warning!\n";
    cout << "  Could not write the file \"" << file_name << "\".\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_type[node] == 1 )
    {
      file_unit << "  " << setw(14) << node_xy[0+node*2]
                << "  " << setw(14) << node_xy[1+node*2] << "\n";
    }
  }

  file_unit.close ( );

  return;
}
//****************************************************************************80

void points_plot ( char *file_name, int node_num, double node_xy[],
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
//    Input, char *FILE_NAME, the name of the file to create.
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

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "POINTS_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: points_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";
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

void pressure3_write ( char *file_name, int node_num, int node_p_variable[],
  int variable_num, double node_c[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRESSURE3_WRITE writes the pressures to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the file name.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_P_VARIABLE[NODE_NUM],
//    is the index of the pressure variable associated with the node,
//    or -1 if there is no associated pressure variable.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, double NODE_C[VARIABLE_NUM], the finite element coefficients.
//
{
  ofstream file_unit;
  int node;
  int variable;

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "PRESSURE3_WRITE - Warning!\n";
    cout << "  Could not write the file \"" << file_name << "\".\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    variable = node_p_variable[node];

    if ( 0 < variable )
    {
      file_unit << "  " << setw(14) << node_c[variable-1] << "\n";
    }
  }

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
//  Examples:
//
//        X       Value
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
  int ihi, int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    A DMAT is a doubly dimensioned array of double precision values, which
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
//    Input, char *TITLE, an optional title.
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
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title )

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
//    Input, char *TITLE, an optional title.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       << "  "
         << "  " << setw(14) << a[i-1]  << "\n";
  }

  return;
}
//****************************************************************************80

void reference_to_physical_t6 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 6 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    The mapping from (XSI,ETA) to (X,Y) has the form:
//
//      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
//                 + D1 * XSI    + E1 * ETA     + F1
//
//      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
//                 + D2 * XSI    + E2 * ETA     + F2
//
//  Reference Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0),
//    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
//
//    Input, integer N, the number of objects to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double e[2];
  double f[2];
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    a[i] =   2.0 * t[i+0*2] + 2.0 * t[i+1*2]
           - 4.0 * t[i+3*2];

    b[i] =   4.0 * t[i+0*2]
           - 4.0 * t[i+3*2] + 4.0 * t[i+4*2] - 4.0 * t[i+5*2];

    c[i] =   2.0 * t[i+0*2]                  + 2.0 * t[i+2*2]
                                             - 4.0 * t[i+5*2];

    d[i] = - 3.0 * t[i+0*2] -       t[i+1*2]
           + 4.0 * t[i+3*2];

    e[i] = - 3.0 * t[i+0*2]                  -       t[i+2*2]
                                             + 4.0 * t[i+5*2];
    f[i] =         t[i+0*2];

  }

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = a[i] * ref[0+j*2] * ref[0+j*2]
                 + b[i] * ref[0+j*2] * ref[1+j*2]
                 + c[i] * ref[1+j*2] * ref[1+j*2]
                 + d[i] * ref[0+j*2]
                 + e[i] * ref[1+j*2]
                 + f[i];
    }
  }

  return;
}
//****************************************************************************80

int s_len_trim ( char *s )

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
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( char *s, int *last, bool *error )

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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a string to be examined.
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

  while ( *s )
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

bool s_to_i4vec ( char *s, int n, int ivec[] )

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
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  error = false;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s, &lchar, &error );

    if ( error )
    {
      cout << "\n";
      cout << "S_TO_I4VEC - Fatal error!\n";
      cout << "  S_TO_I4 returned error while reading item " << i << "\n";
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( char *s, int *lchar, bool *error )

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
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
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

bool s_to_r8vec ( char *s, int n, double rvec[] )

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
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

int s_word_count ( char *s )

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
//    08 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int nword;

  nword = 0;
  blank = true;

  while ( *s )
  {
    if ( *s == ' ' )
    {
      blank = true;
    }
    else if ( blank )
    {
      nword = nword + 1;
      blank = false;
    }
    *s++;
  }

  return nword;
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
//    Original FORTRAN77 version by Nijenhuis and Wilf.
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

void triangles3_write ( char *file_name, int element_num, int element_node[],
  int node_num, int node3_label[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLES3_WRITE writes the pressure triangles to a file.
//
//  Discussion:
//
//    The first three rows of the array ELEMENT_NODE(6,NODE) contain
//    exactly the nodes that make up the pressure triangles.
//
//    However, we must relabel the nodes!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the file name.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that
//    make up each triangle.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE3_LABEL[NODE_NUM], contains the renumbered
//    label of order3 nodes, and -1 for nodes that are not order3 nodes.
//
{
  ofstream file_unit;
  int i;
  int triangle;

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGLES3_WRITE - Warning!\n";
    cout << "  Could not write the file \"" << file_name << "\".\n";
    return;
  }
  for ( triangle = 0; triangle < element_num; triangle++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      file_unit << "  " << setw(8)
                << node3_label [ element_node[i+triangle*6] - 1 ];
    }
    file_unit << "\n";
  }

  file_unit.close ( );

  return;
}
//****************************************************************************80

bool *triangulation_order6_boundary_node ( int node_num, int element_num,
  int element_node[] )

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
//    14 June 2005
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
//    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that make up the
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
  n = 3 * element_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < element_num; j++ )
  {
    edge[0+(j               )*m] = element_node[0+j*6];
    edge[1+(j               )*m] = element_node[3+j*6];
    edge[2+(j               )*m] = element_node[1+j*6];

    edge[0+(j+  element_num)*m] = element_node[1+j*6];
    edge[1+(j+  element_num)*m] = element_node[4+j*6];
    edge[2+(j+  element_num)*m] = element_node[2+j*6];

    edge[0+(j+2*element_num)*m] = element_node[2+j*6];
    edge[1+(j+2*element_num)*m] = element_node[5+j*6];
    edge[2+(j+2*element_num)*m] = element_node[0+j*6];
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

void triangulation_order6_plot ( char *file_name, int node_num, double node_xy[],
  int element_num, int element_node[], int node_show, int triangle_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//    This routine has been specialized to deal correctly ONLY with
//    a mesh of 6 node elements, with the property that starting
//    from local node 1 and traversing the edges of the element will
//    result in encountering local nodes 1, 4, 2, 5, 3, 6 in that order.
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
//    Input, char *FILE_NAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[6*ELEMENT_NUM], lists, for each triangle,
//    the indices of the nodes that form the vertices and midsides
//    of the triangle.
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
  int ip1;
  int node;
  int order[6] = { 1, 4, 2, 5, 3, 6 };
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

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER6_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order6_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";
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
  if ( 1 <= triangle_show )
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
      node = element_node[order[0]-1+triangle*6] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  " << x_ps << "  " << y_ps << "  moveto\n";

      for ( i = 1; i <= 6; i++ )
      {
        ip1 = ( i % 6 ) + 1;
        node = element_node[order[ip1-1]-1+triangle*6] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        file_unit << x_ps << "  " << y_ps << "  lineto\n";
      }
      file_unit << "stroke\n";
    }
  }
//
//  Label the triangles.
//
  if ( 2 <= triangle_show )
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

      for ( i = 0; i < 6; i++ )
      {
        node = element_node[i+triangle*6] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }

      ave_x = ave_x / 6.0;
      ave_y = ave_y / 6.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      file_unit << setw(4) << x_ps << "  "
                << setw(4) << y_ps << "  "
                << "moveto (" << triangle+1 << ") show\n";
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

void velocity6_write ( char *file_name, int node_num, int node_u_variable[],
  int node_v_variable[], int variable_num, double node_c[] )

//****************************************************************************80
//
//  Purpose:
//
//    VELOCITY6_WRITE writes the velocities to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_NAME, the file name.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
//    the indices of the horizontal and vertical velocity variables
//    associated with the node, or -1 if there is none.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, double NODE_C[VARIABLE_NUM], the finite element coefficients.
//
{
  ofstream file_unit;
  int node;
  double u;
  int u_index;
  double v;
  int v_index;
  int variable;

  file_unit.open ( file_name );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "VELOCITY6_WRITE - Warning!\n";
    cout << "  Could not write the file \"" << file_name << "\"\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    u_index = node_u_variable[node];

    if ( 0 < u_index )
    {
      u = node_c[u_index-1];
    }
    else
    {
      u = 0.0;
    }

    v_index = node_v_variable[node];

    if ( 0 < v_index )
    {
      v = node_c[v_index-1];
    }
    else
    {
      v = 0.0;
    }
    file_unit << "  " << setw(14) << u
              << "  " << setw(14) << v << "\n";
  }

  file_unit.close ( );

  return;
}
