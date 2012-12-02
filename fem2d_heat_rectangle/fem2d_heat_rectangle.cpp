# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( void );
void adjust_backward_euler ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num, double wq[],
  double xq[], double yq[], double element_area[], int ib, double time,
  double time_step_size, double u_old[], double a[], double f[] );
void adjust_boundary ( int node_num, double node_xy[], int node_boundary[],
  int ib, double time, double a[], double f[] );
void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] );
void assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num,
  double wq[], double xq[], double yq[], double element_area[],
  int ib, double time, double a[], double f[] );
int bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num );
void compare ( int node_num, double node_xy[], double time, double u[] );
int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
void dgb_print_some ( int m, int n, int ml, int mu, double a[], int ilo,
  int jlo, int ihi, int jhi, char *title );
double *dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job );
void element_write ( int nnodes, int element_num, int element_node[],
  char *triangulation_txt_file_name );
void errors ( double element_area[], int element_node[], double node_xy[],
  double u[], int element_num, int nnodes,
  int node_num, double time, double *el2, double *eh1 );
void exact_u ( int node_num, double node_xy[], double time, double u_exact[],
  double dudx_exact[], double dudy_exact[] );
void file_name_inc ( char *file_name );
void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_print_some ( int n, int a[], int max_print, char *title );
int *node_boundary_set ( int nx, int ny, int node_num );
void nodes_plot ( char *file_name, int node_num, double node_xy[],
  bool node_label );
void nodes_write ( int node_num, double node_xy[], char *output_filename );
void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *bb, double *bx, double *by );
void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] );
void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] );
double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title );
double rhs ( double x, double y, double time );
int s_len_trim ( char *s );
void solution_write ( int node_num, double u[], char *u_file_name ) ;
void timestamp ( void );
void triangulation_order6_plot ( char *file_name, int node_num, double node_xy[],
  int tri_num, int triangle_node[], int node_show, int triangle_show );
void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_HEAT_RECTANGLE.
//
//  Discussion:
//
//    FEM2D_HEAT_RECTANGLE solves
//
//      dUdT - Laplacian U(X,Y,T) = F(X,Y,T)
//
//    in a rectangular region in the plane.
//
//    Along the boundary of the region, Dirichlet conditions
//    are imposed:
//
//      U(X,Y,T) = G(X,Y,T)
//
//    At the initial time T_INIT, the value of U is given at all points
//    in the region:
//
//      U(X,Y,T) = H(X,Y,T)
//
//    The code uses continuous piecewise quadratic basis functions on
//    triangles determined by a uniform grid of NX by NY points.
//
//    The backward Euler approximation is used for the time derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*NODE_NUM], the coefficient matrix.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Local, int ELEMENT_NODE[ELEMENT_NUM*NNODES]; ELEMENT_NODE(I,J) is the
//    global node index of the local node J in element I.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, double F[NODE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, integer NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
//    Local, int NNODES, the number of nodes used to form one element.
//
//    Local, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Local, int NX, the number of points in the X direction.
//
//    Local, int NY, the number of points in the Y direction.
//
//    Local, int QUAD_NUM, the number of quadrature points used for assembly.
//
//    Local, double TIME, the current time.
//
//    Local, double TIME_FINAL, the final time.
//
//    Local, double TIME_INIT, the initial time.
//
//    Local, double TIME_OLD, the time at the previous time step.
//
//    Local, int TIME_STEP_NUM, the number of time steps to take.
//
//    Local, double TIME_STEP_SIZE, the size of the time steps.
//
//    Local, double U[NODE_NUM], the finite element coefficients
//    defining the solution at the current time.
//
//    Local, double WQ[QUAD_NUM], quadrature weights.
//
//    Local, double XL, XR, YB, YT, the X coordinates of
//    the left and right sides of the rectangle, and the Y coordinates
//    of the bottom and top of the rectangle.
//
//    Local, double XQ[QUAD_NUM*ELEMENT_NUM], YQ[QUAD_NUM*ELEMENT_NUM], the X and Y
//    coordinates of the quadrature points in each element.
//
{
# define NNODES 6
# define QUAD_NUM 3
# define NX 5
# define NY 5

# define ELEMENT_NUM ( NX - 1 ) * ( NY - 1 ) * 2
# define NODE_NUM ( 2 * NX - 1 ) * ( 2 * NY - 1 )

  double *a;
  double dt;
  double *dudx_exact;
  double *dudy_exact;
  double eh1;
  double el2;
  int element;
  double element_area[ELEMENT_NUM];
  bool *element_mask;
  int element_node[NNODES*ELEMENT_NUM];
  double *f;
  int i;
  int ib;
  int ierr;
  int job;
  int local;
  int node;
  int *node_boundary;
  char *node_eps_file_name;
  char *node_txt_file_name = "rectangle_nodes.txt";
  bool node_label;
  int node_show;
  double node_xy[2*NODE_NUM];
  int *pivot;
  double time;
  char *time_file_name = "rectangle_time.txt";
  double time_final;
  double time_init;
  double time_old;
  int time_step;
  int time_step_num;
  double time_step_size;
  ofstream time_unit;
  int triangle_show;
  char *triangulation_eps_file_name = "rectangle_elements.eps";
  char *triangulation_txt_file_name = "rectangle_elements.txt";
  double *u;
  double *u_exact;
  char u_file_name[100];
  double *u_old;
  double wq[QUAD_NUM];
  double xl = 0.0;
  double xq[QUAD_NUM*ELEMENT_NUM];
  double xr = 1.0;
  double yb = 0.0;
  double yq[QUAD_NUM*ELEMENT_NUM];
  double yt = 1.0;

  strcpy ( node_eps_file_name, "rectangle_nodes.eps" );
  strcpy ( u_file_name, "rectangle_u0000.txt" );

  timestamp ( );

  cout << "\n";
  cout << "FEM2D_HEAT_RECTANGLE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Solution of the time-dependent heat equation\n";
  cout << "  on a unit box in 2 dimensions.\n";
  cout << "\n";
  cout << "  Ut - Uxx - Uyy = F(x,y,t) in the box\n";
  cout << "        U(x,y,t) = G(x,y,t) for (x,y) on the boundary.\n";
  cout << "        U(x,y,t) = H(x,y,t) for t = T_INIT.\n";
  cout << "\n";
  cout << "  The finite element method is used, with piecewise\n";
  cout << "  quadratic basis functions on 6 node triangular\n";
  cout << "  elements.\n";
  cout << "\n";
  cout << "  The corner nodes of the triangles are generated by an\n";
  cout << "  underlying grid whose dimensions are\n";
  cout << "\n";
  cout << "  NX =                 " << NX << "\n";
  cout << "  NY =                 " << NY << "\n";
  cout << "\n";
  cout << "  Number of nodes    = " << NODE_NUM << "\n";
  cout << "  Number of elements = " << ELEMENT_NUM << "\n";
//
//  Set the coordinates of the nodes.
//
  xy_set ( NX, NY, NODE_NUM, xl, xr, yb, yt, node_xy );
//
//  Organize the nodes into a grid of 6-node triangles.
//
  grid_t6 ( NX, NY, NNODES, ELEMENT_NUM, element_node );
//
//  Set the quadrature rule for assembly.
//
  quad_a ( node_xy, element_node, ELEMENT_NUM, NODE_NUM,
    NNODES, wq, xq, yq );
//
//  Determine the areas of the elements.
//
  area_set ( NODE_NUM, node_xy, NNODES, ELEMENT_NUM,
    element_node, element_area );
//
//  Determine which nodes are boundary nodes and which have a
//  finite element unknown.  Then set the boundary values.
//
  node_boundary = node_boundary_set ( NX, NY, NODE_NUM );

  if ( false )
  {
    i4vec_print_some ( NODE_NUM, node_boundary, 10, "  NODE_BOUNDARY:" );
  }
//
//  Determine the bandwidth of the coefficient matrix.
//
  ib = bandwidth ( NNODES, ELEMENT_NUM, element_node, NODE_NUM );

  cout << "\n";
  cout << "  The matrix half bandwidth is " << ib << "\n";
  cout << "  The matrix row size is       " << 3 * ib + 1 << "\n";
//
//  Make an EPS picture of the nodes.
//
  if ( NX <= 10 && NY <= 10 )
  {
    node_label = true;
    nodes_plot ( node_eps_file_name, NODE_NUM, node_xy, node_label );

    cout << "\n";
    cout << "FEM2D_HEAT_RECTANGLE:\n";
    cout << "  Wrote an EPS file\n";
    cout << "    \"" << node_eps_file_name << "\".\n";
    cout << "  containing a picture of the nodes.\n";
  }
//
//  Write the nodes to an ASCII file that can be read into MATLAB.
//
  nodes_write ( NODE_NUM, node_xy, node_txt_file_name );

  cout << "\n";
  cout << "FEM2D_HEAT_RECTANGLE:\n";
  cout << "  Wrote an ASCII node file\n";
  cout << "    " << node_txt_file_name << "\n";
  cout << "  of the form\n";
  cout << "    X(I), Y(I)\n";
  cout << "  which can be used for plotting.\n";
//
//  Make a picture of the elements.
//
  if ( NX <= 10 && NY <= 10 )
  {
    node_show = 2;
    triangle_show = 2;

    triangulation_order6_plot ( triangulation_eps_file_name, NODE_NUM,
      node_xy, ELEMENT_NUM, element_node, node_show, triangle_show );

    cout << "\n";
    cout << "FEM2D_HEAT_RECTANGLE:\n";
    cout << "  Wrote an EPS file\n";
    cout << "    \"" << triangulation_eps_file_name << "\".\n";
    cout << "  containing a picture of the elements.\n";
  }
//
//  Write the elements to a file that can be read into MATLAB.
//
  element_write ( NNODES, ELEMENT_NUM, element_node,
    triangulation_txt_file_name );

  cout << "\n";
  cout << "FEM2D_HEAT_RECTANGLE:\n";
  cout << "  Wrote an ASCII element file\n";
  cout << "    \"" << triangulation_txt_file_name << "\".\n";
  cout << "  of the form\n";
  cout << "    Node(1) Node(2) Node(3) Node(4) Node(5) Node(6)\n";
  cout << "  which can be used for plotting.\n";
//
//  Set time stepping quantities.
//
  time_init = 0.0;
  time_final = 0.5;
  time_step_num = 10;
  time_step_size = ( time_final - time_init ) / ( double ) ( time_step_num );
//
//  Allocate space.
//
  a = new double[(3*ib+1)*NODE_NUM];
  dudx_exact = new double[NODE_NUM];
  dudy_exact = new double[NODE_NUM];
  f = new double[NODE_NUM];
  pivot = new int[NODE_NUM];
  u = new double[NODE_NUM];
  u_exact = new double[NODE_NUM];
  u_old = new double[NODE_NUM];
//
//  Set the value of U at the initial time.
//
  time = time_init;
  exact_u ( NODE_NUM, node_xy, time, u_exact, dudx_exact, dudy_exact );

  for ( node = 0; node < NODE_NUM; node++ )
  {
    u[node] = u_exact[node];
  }

  time_unit.open ( time_file_name );

  if ( !time_unit )
  {
    cout << "\n";
    cout << "FEM2D_HEAT_RECTANGLE- Warning!\n";
    cout << "  Could not write the time file \""
         << time_file_name << "\".\n";
    exit ( 1 );
  }

  time_unit << "  " << setw(14) << time << "\n";

  solution_write ( NODE_NUM, u, u_file_name );
//
//  Time looping.
//
  cout << "\n";
  cout << "     Time        L2 Error       H1 Error\n";
  cout << "\n";

  for ( time_step = 1; time_step <= time_step_num; time_step++ )
  {
    time_old = time;
    for ( node = 0; node < NODE_NUM; node++ )
    {
      u_old[node] = u[node];
    }
    delete [] u;

    time = ( ( double ) ( time_step_num - time_step ) * time_init
           + ( double ) (                 time_step ) * time_final )
           / ( double ) ( time_step_num             );
//
//  Assemble the coefficient matrix A and the right-hand side F of the
//  finite element equations.
//
    assemble ( NODE_NUM, node_xy, NNODES,
      ELEMENT_NUM, element_node, QUAD_NUM,
      wq, xq, yq, element_area, ib, time, a, f );

    if ( false )
    {
      dgb_print_some ( NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25,
        "  Initial block of coefficient matrix A:" );

      r8vec_print_some ( NODE_NUM, f, 1, 10,
        "  Part of the right hand side F:" );
    }
//
//  Modify the coefficient matrix and right hand side to account for the dU/dt
//  term, which we are treating using the backward Euler formula.
//
    adjust_backward_euler ( NODE_NUM, node_xy, NNODES, ELEMENT_NUM,
      element_node, QUAD_NUM, wq, xq, yq, element_area, ib, time,
      time_step_size, u_old, a, f );

    if ( false )
    {
      dgb_print_some ( NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25,
        "  A after DT adjustment:" );

      r8vec_print_some ( NODE_NUM, f, 1, 10,
        "  F after DT adjustment:" );
    }
//
//  Modify the coefficient matrix and right hand side to account for
//  boundary conditions.
//
    adjust_boundary ( NODE_NUM, node_xy, node_boundary, ib, time, a, f );

    if ( false )
    {
      dgb_print_some ( NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25,
        "  A after BC adjustment:" );

      r8vec_print_some ( NODE_NUM, f, 1, 10,
        "  F after BC adjustment:" );
    }
//
//  Solve the linear system using a banded solver.
//
    ierr = dgb_fa ( NODE_NUM, ib, ib, a, pivot );

    if ( ierr != 0 )
    {
      cout << "\n";
      cout << "FEM2D_HEAT_RECTANGLE - Error!\n";
      cout << "  DGB_FA returned an error condition.\n";
      cout << "\n";
      cout << "  The linear system was not factored, and the\n";
      cout << "  algorithm cannot proceed.\n";
      exit ( 1 );
    }

    job = 0;
    u = dgb_sl ( NODE_NUM, ib, ib, a, pivot, f, job );

    if ( false )
    {
      r8vec_print_some ( NODE_NUM, u, 1, 10,
        "  Part of the solution vector:" );
    }
//
//  Calculate error using 13 point quadrature rule.
//
    errors ( element_area, element_node, node_xy, u,
      ELEMENT_NUM, NNODES, NODE_NUM, time, &el2, &eh1 );
//
//  Compare the exact and computed solutions just at the nodes.
//
    if ( false )
    {
      compare ( NODE_NUM, node_xy, time, u );
    }
//
//  Increment the file name, and write the new solution.
//
    time_unit << setw(14) << time << "\n";

    file_name_inc ( u_file_name );

    solution_write ( NODE_NUM, u, u_file_name );
  }
//
//  Deallocate memory.
//
  delete [] a;
  delete [] dudx_exact;
  delete [] dudy_exact;
  delete [] f;
  delete [] node_boundary;
  delete [] pivot;
  delete [] u;
  delete [] u_exact;
  delete [] u_old;

  time_unit.close ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_HEAT_RECTANGLE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
# undef ELEMENT_NUM
# undef NNODES
# undef NODE_NUM
# undef NX
# undef NY
# undef QUAD_NUM
}
//****************************************************************************80

void adjust_backward_euler ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int quad_num, double wq[],
  double xq[], double yq[], double element_area[], int ib, double time,
  double time_step_size, double u_old[], double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJUST_BACKWARD_EULER adjusts the system for the backward Euler term.
//
//  Discussion:
//
//    The input linear system
//
//      A * U = F
//
//    is appropriate for the equation
//
//      -Uxx - Uyy = RHS
//
//    We need to modify the matrix A and the right hand side F to
//    account for the approximation of the time derivative in
//
//      Ut - Uxx - Uyy = RHS
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
//    01 May 2006
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
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, double WQ[QUAD_NUM], quadrature weights.
//
//    Input, double XQ[QUAD_NUM*ELEMENT_NUM],
//    YQ[QUAD_NUM*ELEMENT_NUM], the
//    coordinates of the quadrature points in each element.
//
//    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Input, double TIME_STEP_SIZE, the size of the time step.
//
//    Input, double U_OLD[NODE_NUM], the finite element
//    coefficients for the solution at the previous time.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM
//    by NODE_NUM coefficient matrix, stored in a compressed format.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//
{
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int ip;
  int ipp;
  int j;
  int node;
  int quad;
  int test;
  double w;
  double x;
  double y;

  for ( element = 0; element < element_num; element++ )
  {
    for ( quad = 0; quad < quad_num; quad++ )
    {
      x = xq[quad+element*quad_num];
      y = yq[quad+element*quad_num];
      w = element_area[element] * wq[quad];

      for ( test = 0; test < nnodes; test++ )
      {
        node = element_node[test+element*nnodes];

        qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );
//
//  Carry the U_OLD term to the right hand side.
//
        f[node] = f[node] + w * bi * u_old[node] / time_step_size;
//
//  Modify the diagonal entries of A.
//
        for ( basis = 0; basis < nnodes; basis++ )
        {
          j = element_node[basis+element*nnodes];

          qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          a[node-j+2*ib+j*(3*ib+1)] = a[node-j+2*ib+j*(3*ib+1)]
            + w * bi * bj / time_step_size;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

void adjust_boundary ( int node_num, double node_xy[], int node_boundary[],
  int ib, double time, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJUST_BOUNDARY modifies the linear system for boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2006
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
//    Input, int NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//    On output, A has been adjusted for boundary conditions.
//
//    Input/output, double F[NODE_NUM], the right hand side.
//    On output, F has been adjusted for boundary conditions.
//
{
  double *dudx_exact;
  double *dudy_exact;
  int i;
  int j;
  int jhi;
  int jlo;
  int node;
  double *u_exact;
//
//  Get the exact solution at every node.
//
  u_exact = new double[node_num];
  dudx_exact = new double[node_num];
  dudy_exact = new double[node_num];

  exact_u ( node_num, node_xy, time, u_exact, dudx_exact, dudy_exact );

  for ( node = 0; node < node_num; node++ )
  {
    if ( node_boundary[node] != 0 )
    {
      jlo = i4_max ( node - ib, 0 );
      jhi = i4_min ( node + ib, node_num - 1 );

      for ( j = jlo; j <= jhi; j++ )
      {
        a[node-j+2*ib+j*(3*ib+1)] = 0.0;
      }
      a[node-node+2*ib+node*(3*ib+1)] = 1.0;

      f[node] = u_exact[node];
    }
  }

  delete [] u_exact;
  delete [] dudx_exact;
  delete [] dudy_exact;

  return;
}
//****************************************************************************80

void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_SET sets the area of each element.
//
//  Discussion:
//
//    The areas of the elements are needed in order to adjust
//    the integral estimates produced by the quadrature formulas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the
//    coordinates of the nodes.
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
//
{
  int element;
  int i1;
  int i2;
  int i3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  for ( element = 0; element < element_num; element++ )
  {
    i1 = element_node[0+element*nnodes];
    x1 = node_xy[0+i1*2];
    y1 = node_xy[1+i1*2];

    i2 = element_node[1+element*nnodes];
    x2 = node_xy[0+i2*2];
    y2 = node_xy[1+i2*2];

    i3 = element_node[2+element*nnodes];
    x3 = node_xy[0+i3*2];
    y3 = node_xy[1+i3*2];

    element_area[element] = 0.5 * fabs
      ( y1 * ( x2 - x3 )
      + y2 * ( x3 - x1 )
      + y3 * ( x1 - x2 ) );
  }

  return;
}
//****************************************************************************80*

void assemble ( int node_num, double node_xy[], int nnodes, int element_num,
  int element_node[], int quad_num, double wq[], double xq[], double yq[],
  double element_area[], int ib, double time, double a[], double f[] )

//****************************************************************************80*
//
//  Purpose:
//
//    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.
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
//    14 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int QUAD_NUM, the number of quadrature points used in assembly.
//
//    Input, double WQ[QUAD_NUM], quadrature weights.
//
//    Input, double XQ[QUAD_NUM*ELEMENT_NUM], YQ[QUAD_NUM*ELEMENT_NUM], the X and Y
//    coordinates of the quadrature points in each element.
//
//    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Input, int IB, the half-bandwidth of the matrix.
//
//    Input, double TIME, the current time.
//
//    Output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
//    coefficient matrix, stored in a compressed format.
//
//    Output, double F[NODE_NUM], the right hand side.
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
  double aij;
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
  int quad;
  int test;
  double w;
  double x;
  double y;
//
//  Initialize the arrays to zero.
//
  for ( i = 0; i < node_num; i++ )
  {
    f[i] = 0.0;
  }

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < 3*ib + 1; i++ )
    {
      a[i+j*(3*ib+1)] = 0.0;
    }
  }
//
//  The actual values of A and F are determined by summing up
//  contributions from all the elements.
//
  for ( element = 0; element < element_num; element++ )
  {
    for ( quad = 0; quad < quad_num; quad++ )
    {
      x = xq[quad+element*quad_num];
      y = yq[quad+element*quad_num];
      w = element_area[element] * wq[quad];

      for ( test = 0; test < nnodes; test++ )
      {
        node = element_node[test+element*nnodes];

        qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        f[node] = f[node] + w * rhs ( x, y, time ) * bi;
//
//  We are about to compute a contribution associated with the
//  I-th test function and the J-th basis function, and add this
//  to the entry A(I,J).
//
//  Because of the compressed storage of the matrix, the element
//  will actually be stored in A(I-J+2*IB+1,J).
//
//  An extra complication: we are storing the array as a vector.
//
//  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + J * (3*IB+1)];
//
        for ( basis = 0; basis < nnodes; basis++ )
        {
          j = element_node[basis+element*nnodes];

          qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          aij = dbidx * dbjdx + dbidy * dbjdy;

          a[node-j+2*ib+j*(3*ib+1)] = a[node-j+2*ib+j*(3*ib+1)] + w * aij;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

int bandwidth ( int nnodes, int element_num, int element_node[], int node_num )

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
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int BANDWIDTH, the half bandwidth of the matrix.
//
{
  int element;
  int i;
  int iln;
  int in;
  int j;
  int jln;
  int jn;
  int nhba;

  nhba = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( iln = 0; iln < nnodes; iln++ )
    {
      i = element_node[iln+element*nnodes];
      for ( jln = 0; jln < nnodes; jln++ )
      {
        j = element_node[jln+element*nnodes];
        nhba = i4_max ( nhba, j - i );
      }
    }
  }

  return nhba;
}
//****************************************************************************80

void compare ( int node_num, double node_xy[], double time, double u[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMPARE compares the exact and computed solution at the nodes.
//
//  Discussion:
//
//    This is a rough comparison, done only at the nodes.  Such a pointwise
//    comparison is easy, because the value of the finite element
//    solution is exactly the value of the finite element coefficient
//    associated with that node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, double TIME, the current time.
//
//    Input, double F[NUNK], the solution vector of the finite
//    element system.
//
{
  double *dudx_exact;
  double *dudy_exact;
  int node;
  double *u_exact;

  u_exact = new double[node_num];
  dudx_exact = new double[node_num];
  dudy_exact = new double[node_num];

  exact_u ( node_num, node_xy, time, u_exact, dudx_exact, dudy_exact );

  cout << "\n";
  cout << "COMPARE:\n";
  cout << "  Compare computed and exact solutions at the nodes.\n";
  cout << "\n";
  cout << "         X           Y          U           U\n";
  cout << "                              exact       computed\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << setw(12) << node_xy[0+node*2] << "  "
         << setw(12) << node_xy[1+node*2] << "  "
         << setw(12) << u_exact[node]     << "  "
         << setw(12) << u[node]           << "\n";
  }

  delete [] u_exact;
  delete [] dudx_exact;
  delete [] dudy_exact;

  return;
}
//****************************************************************************80

int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
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
//    20 February 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
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
//    Output, int SGB_FA, singularity flag.
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
//    Input, double A[(2*ML+MU+1)*N], the SGB matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title to print.
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
    i2lo = i4_max ( i2lo, j2lo - mu );

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
        if ( ml < i-j || mu < j-i )
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

double *dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
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
//    20 February 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, int ML, MU, the lower and upper bandwidths.
//    ML and MU must be nonnegative, and no greater than N-1.
//
//    Input, double A[(2*ML+MU+1)*N], the LU factors from DGB_FA.
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
          x[k+i-1] = x[k+i-1] + x[k-1] * a[m+i-1+(k-1)*col];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
      lm = i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a[la+i-1+(k-1)*col];
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
        x[k-1] = x[k-1] - x[lb+i-1] * a[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
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
          x[k-1] = x[k-1] + x[k+i-1] * a[m+i-1+(k-1)*col];
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
//****************************************************************************80*

void element_write ( int nnodes, int element_num, int element_node[],
 char *output_filename )

//****************************************************************************80*
//
//  Purpose:
//
//    ELEMENT_WRITE writes the elements to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, char *OUTPUT_FILENAME, the name of the file
//    in which the data should be stored.
//
{
  int element;
  int i;
  ofstream output;

  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "ELEMENT_WRITE - Warning!\n";
    cout << "  Could not write the node file.\n";
    return;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( i = 0; i < nnodes; i++ )
    {
      output << setw(8)  << element_node[i+element*nnodes] << "  ";
    }
    output << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80*

void errors ( double element_area[], int element_node[], double node_xy[],
  double u[], int element_num, int nnodes, int node_num, double time,
  double *el2, double *eh1 )

//****************************************************************************80*
//
//  Purpose:
//
//    ERRORS calculates the error in the L2 and H1-seminorm.
//
//  Discussion:
//
//    This routine uses a 13 point quadrature rule in each element,
//    in order to estimate the values of
//
//      EL2 = Sqrt ( Integral ( U(x,y) - Uh(x,y) )**2 dx dy )
//
//      EH1 = Sqrt ( Integral ( Ux(x,y) - Uhx(x,y) )**2 +
//                            ( Uy(x,y) - Uhy(x,y) )**2 dx dy )
//
//    Here U is the exact solution, and Ux and Uy its spatial derivatives,
//    as evaluated by a user-supplied routine.
//
//    Uh, Uhx and Uhy are the computed solution and its spatial derivatives,
//    as specified by the computed finite element solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, double U[NUNK], the coefficients of the solution.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, real ( kind = 8 ) TIME, the current time.
//
//    Output, double precision *EL2, the L2 error.
//
//    Output, double precision *EH1, the H1 seminorm error.
//
//  Local Parameters:
//
//    Local, double AR, the weight for a given quadrature point
//    in a given element.
//
//    Local, double BI, DBIDX, DBIDY, a basis function and its first
//    derivatives evaluated at a particular quadrature point.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, int NQE, the number of points in the quadrature rule.
//    This is actually fixed at 13.
//
//    Local, double UEX, UEXX, UEXY, the exact solution and its first
//    derivatives evaluated at a particular quadrature point.
//
//    Local, double UH, UHX, UHY, the computed solution and its first
//    derivatives evaluated at a particular quadrature point.
//
//    Local, double WQE[NQE], stores the quadrature weights.
//
//    Local, double X, Y, the coordinates of a particular
//    quadrature point.
//
//    Local, double XQE[NQE], YQE[NQE], stores the location
//    of quadrature points in a given element.
//
{
# define NQE 13

  double ar;
  double bi;
  double dbidx;
  double dbidy;
  double dudx_exact[1];
  double dudxh;
  double dudy_exact[1];
  double dudyh;
  int element;
  int i;
  int in1;
  int ip;
  int quad;
  double u_exact[1];
  double uh;
  double wqe[NQE];
  double x;
  double x1;
  double xqe[NQE];
  double xy[2];
  double y;
  double y1;
  double yqe[NQE];

  *el2 = 0.0;
  *eh1 = 0.0;
//
//  For each element, retrieve the nodes, area, quadrature weights,
//  and quadrature points.
//
  for ( element = 0; element < element_num; element++ )
  {
    quad_e ( node_xy, element_node, element, element_num,
      nnodes, node_num, NQE, wqe, xqe, yqe );
//
//  For each quadrature point, evaluate the computed solution and its X and
//  Y derivatives.
//
    for ( quad = 0; quad < NQE; quad++ )
    {
      ar = element_area[element] * wqe[quad];
      x = xqe[quad];
      y = yqe[quad];

      uh = 0.0;
      dudxh = 0.0;
      dudyh = 0.0;

      for ( in1 = 0; in1 < nnodes; in1++ )
      {
        i = element_node[in1+element*nnodes];

        qbf ( x, y, element, in1, node_xy,
          element_node, element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        uh    = uh    + bi    * u[i];
        dudxh = dudxh + dbidx * u[i];
        dudyh = dudyh + dbidy * u[i];
      }
//
//  Evaluate the exact solution and its X and Y derivatives.
//
      xy[0] = x;
      xy[1] = y;

      exact_u ( 1, xy, time, u_exact, dudx_exact, dudy_exact );
//
//  Add the weighted value at this quadrature point to the quadrature sum.
//
      *el2 = *el2 + ar *   pow ( ( uh    - u_exact[0]  ), 2 );

      *eh1 = *eh1 + ar * ( pow ( ( dudxh - dudx_exact[0] ), 2 )
                         + pow ( ( dudyh - dudy_exact[0] ), 2 ) );
    }
  }

  *el2 = sqrt ( *el2 );
  *eh1 = sqrt ( *eh1 );

  cout << setw(14) << time
       << setw(14) << *el2
       << setw(14) << *eh1 << "\n";

  return;
# undef NQE
}
//****************************************************************************80

void exact_u ( int node_num, double node_xy[], double time, double u[],
  double dudx[], double dudy[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT_U calculates the exact solution and its first derivatives.
//
//  Discussion:
//
//    It is assumed that the user knows the exact solution and its
//    derivatives.  This, of course, is NOT true for a real computation.
//    But for this code, we are interested in studying the convergence
//    behavior of the approximations, and so we really need to assume
//    we know the correct solution.
//
//    As a convenience, this single routine is used for several purposes:\
//
//    * it supplies the initial value function H(X,Y,T);
//    * it supplies the boundary value function G(X,Y,T);
//    * it is used by the COMPARE routine to make a node-wise comparison
//      of the exact and approximate solutions.
//    * it is used by the ERRORS routine to estimate the integrals of
//      the L2 and H1 errors of approximation.
//
//    DUDX and DUDY are only needed for the ERRORS calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes at which
//    a value is desired.
//
//    Input, double NODE_XY(2,NODE_NUM), the coordinates of
//    the points where a value is desired.
//
//    Input, double TIME, the current time.
//
//    Output, double U[NODE_NUM], the exact solution.
//
//    Output, double DUDX[NODE_NUM], DUDY[NODE_NUM],
//    the X and Y derivatives of the exact solution.
//
{
# define PI 3.141592653589793

  int node;
  double x;
  double y;

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    u[node]    =      sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
    dudx[node] = PI * cos ( PI * x ) * sin ( PI * y ) * exp ( - time );
    dudy[node] = PI * sin ( PI * x ) * cos ( PI * y ) * exp ( - time );
  }

  return;
# undef PI
}
//****************************************************************************80

void file_name_inc ( char *file_name )

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
//    If the input string contains no digits, a blank string is returned.
//
//    If a blank string is input, then an error condition results.
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
//    14 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, character *FILE_NAME, (a pointer to) the character string
//    to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = s_len_trim ( file_name );

  if ( lens <= 0 )
  {
    cout << "\n";
    cout << "FILE_NAME_INC - Fatal error!\n";
    cout << "  Input file name is blank.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = *(file_name+i);

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      if ( c == '9' )
      {
        c = '0';
        *(file_name+i) = c;
      }
      else
      {
        c = c + 1;
        *(file_name+i) = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    strcpy ( file_name, " " );
  }

  return;
}
//****************************************************************************80

void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T6 produces a grid of pairs of 6 node triangles.
//
//  Example:
//
//    Input:
//
//      NX = 4, NY = 3
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  3, 15,  2,  9,  8;
//        17, 15,  3, 16,  9, 10;
//         3,  5, 17,  4, 11, 10;
//        19, 17,  5, 18, 11, 12;
//         5,  7, 19,  6, 13, 12;
//        21, 19,  7, 20, 13, 14;
//        15, 17, 29, 16, 23, 22;
//        31, 29, 17, 30, 23, 24;
//        17, 19, 31, 18, 25, 24;
//        33, 31, 19, 32, 25, 26;
//        19, 21, 33, 20, 27, 26;
//        35, 33, 21, 34, 27, 28.
//
//  Diagram:
//
//   29-30-31-32-33-34-35
//    |\ 8  |\10  |\12  |
//    | \   | \   | \   |
//   22 23 24 25 26 27 28
//    |   \ |   \ |   \ |
//    |  7 \|  9 \| 11 \|
//   15-16-17-18-19-20-21
//    |\ 2  |\ 4  |\ 6  |
//    | \   | \   | \   |
//    8  9 10 11 12 13 14
//    |   \ |   \ |   \ |
//    |  1 \|  3 \|  5 \|
//    1--2--3--4--5--6--7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, controls the number of elements along the
//    X and Y directions.  The number of elements will be
//    2 * ( NX - 1 ) * ( NY - 1 ).
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
//
{
  int c;
  int e;
  int element;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element = 0;

  for ( j = 1; j <= ny - 1; j++ )
  {
    for ( i = 1; i <= nx - 1; i++ )
    {
      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 2;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s + 1;
      n  = s + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element_node[0+element*nnodes] = sw;
      element_node[1+element*nnodes] = se;
      element_node[2+element*nnodes] = nw;
      element_node[3+element*nnodes] = s;
      element_node[4+element*nnodes] = c;
      element_node[5+element*nnodes] = w;
      element = element + 1;

      element_node[0+element*nnodes] = ne;
      element_node[1+element*nnodes] = nw;
      element_node[2+element*nnodes] = se;
      element_node[3+element*nnodes] = n;
      element_node[4+element*nnodes] = c;
      element_node[5+element*nnodes] = e;
      element = element + 1;
    }
  }

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
//    Input, int I1, I2, are two ints to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
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
//    Input, int I1, I2, two ints to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
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

void i4vec_print_some ( int n, int a[], int max_print, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT_SOME prints "some" of an I4VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2003
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
//    Input, int MAX_PRINT, the maximum number of lines to print.
//
//    Input, char *TITLE, an optional title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
    cout << "\n";
  }

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i]  << "\n";
    }
    cout << "......  ..............\n";
    i = n - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(10) << a[i]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print-1; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i]  << "\n";
    }
    i = max_print - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(10) << a[i]  << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

int *node_boundary_set ( int nx, int ny, int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_BOUNDARY_SET assigns an unknown value index at each node.
//
//  Discussion:
//
//    Every node is assigned a value which indicates whether it is
//    an interior node, or a boundary node.
//
//  Example:
//
//    On a simple 5 by 5 grid, where the nodes are numbered starting
//    at the lower left, and increasing in X first, we would have the
//    following values of NODE_BOUNDARY:
//
//       1  1  1  1  1
//       1  0  0  0  1
//       1  0  0  0  1
//       1  0  0  0  1
//       1  1  1  1  1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of elements in the X and Y directions.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
{
  int i;
  int j;
  int node;
  int *node_boundary;

  node_boundary = new int[node_num];

  node = 0;

  for ( j = 1; j <= 2 * ny - 1; j++ )
  {
    for ( i = 1; i <= 2 * nx - 1; i++ )
    {
      if ( j == 1 ||
           j == 2 * ny - 1 ||
           i == 1 ||
           i == 2 * nx - 1 )
      {
        node_boundary[node] = 1;
      }
      else
      {
        node_boundary[node] = 0;
      }

      node = node + 1;
    }
  }
  return node_boundary;
}
//****************************************************************************80

void nodes_plot ( char *file_name, int node_num, double node_xy[],
  bool node_label )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_PLOT plots a pointset.
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
  int i;
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
  file_unit << "%%Creator: nodes_plot.C\n";
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
      << node     << ") show\n";
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
//****************************************************************************80*

void nodes_write ( int node_num, double node_xy[], char *output_filename )

//****************************************************************************80*
//
//  Purpose:
//
//    NODES_WRITE writes the nodes to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, char *OUTPUT_FILENAME, the name of the file
//    in which the data should be stored.
//
{
  int node;
  ofstream output;
  double x;
  double y;

  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "NODES_WRITE - Warning!\n";
    cout << "  Could not write the node file.\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    output << setw(8)  << x << "  "
           << setw(8)  << y << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *b, double *dbdx, double *dbdy )

//****************************************************************************80
//
//  Purpose:
//
//    QBF evaluates the quadratic basis functions.
//
//  Discussion:
//
//    This routine assumes that the "midpoint" nodes are, in fact,
//    exactly the average of the two extreme nodes.  This is NOT true
//    for a general quadratic triangular element.
//
//    Assuming this property of the midpoint nodes makes it easy to
//    determine the values of (R,S) in the reference element that
//    correspond to (X,Y) in the physical element.
//
//    Once we know the (R,S) coordinates, it's easy to evaluate the
//    basis functions and derivatives.
//
//  The physical element T6:
//
//    In this picture, we don't mean to suggest that the bottom of
//    the physical triangle is horizontal.  However, we do assume that
//    each of the sides is a straight line, and that the intermediate
//    points are exactly halfway on each side.
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
//  Reference element T6:
//
//    In this picture of the reference element, we really do assume
//    that one side is vertical, one horizontal, of length 1.
//
//    |
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-------->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the (global) coordinates of the point
//    at which the basis function is to be evaluated.
//
//    Input, int ELEMENT, the index of the element which contains the point.
//
//    Input, int INODE, the local index, between 0 and 5, that
//    specifies which basis function is to be evaluated.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, double *B, *DBDX, *DBDY, the value of the basis function
//    and its X and Y derivatives at (X,Y).
//
{
  double dbdr;
  double dbds;
  double det;
  double drdx;
  double drdy;
  double dsdx;
  double dsdy;
  int i;
  double r;
  double s;
  double xn[6];
  double yn[6];

  for ( i = 0; i < 6; i++ )
  {
    xn[i] = node_xy[0+element_node[i+element*nnodes]*2];
    yn[i] = node_xy[1+element_node[i+element*nnodes]*2];
  }
//
//  Determine the (R,S) coordinates corresponding to (X,Y).
//
//  What is happening here is that we are solving the linear system:
//
//    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
//    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
//
//  by computing the inverse of the coefficient matrix and multiplying
//  it by the right hand side to get R and S.
//
//  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
//  for R and S.
//
  det =   ( xn[1] - xn[0] ) * ( yn[2] - yn[0] )
        - ( xn[2] - xn[0] ) * ( yn[1] - yn[0] );

  r = ( ( yn[2] - yn[0] ) * ( x     - xn[0] )
      + ( xn[0] - xn[2] ) * ( y     - yn[0] ) ) / det;

  drdx = ( yn[2] - yn[0] ) / det;
  drdy = ( xn[0] - xn[2] ) / det;

  s = ( ( yn[0] - yn[1] ) * ( x     - xn[0] )
      + ( xn[1] - xn[0] ) * ( y     - yn[0] ) ) / det;

  dsdx = ( yn[0] - yn[1] ) / det;
  dsdy = ( xn[1] - xn[0] ) / det;
//
//  The basis functions can now be evaluated in terms of the
//  reference coordinates R and S.  It's also easy to determine
//  the values of the derivatives with respect to R and S.
//
  if ( inode == 0 )
  {
    *b   =   2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
    dbdr = - 3.0 + 4.0 * r + 4.0 * s;
    dbds = - 3.0 + 4.0 * r + 4.0 * s;
  }
  else if ( inode == 1 )
  {
    *b   =   2.0 * r * ( r - 0.5 );
    dbdr = - 1.0 + 4.0 * r;
    dbds =   0.0;
  }
  else if ( inode == 2 )
  {
    *b   =   2.0 * s * ( s - 0.5 );
    dbdr =   0.0;
    dbds = - 1.0               + 4.0 * s;
  }
  else if ( inode == 3 )
  {
    *b   =   4.0 * r * ( 1.0 - r - s );
    dbdr =   4.0 - 8.0 * r - 4.0 * s;
    dbds =           - 4.0 * r;
  }
  else if ( inode == 4 )
  {
    *b   =   4.0 * r * s;
    dbdr =                           4.0 * s;
    dbds =             4.0 * r;
  }
  else if ( inode == 5 )
  {
    *b   =   4.0 * s * ( 1.0 - r - s );
    dbdr =                         - 4.0 * s;
    dbds =   4.0 - 4.0 * r - 8.0 * s;
  }
  else
  {
    cout << "\n";
    cout << "QBF - Fatal error!\n";
    cout << "  Request for local basis function INODE = " << inode << "\n";
    exit ( 1 );
  }
//
//  We need to convert the derivative information from (R(X,Y),S(X,Y))
//  to (X,Y) using the chain rule.
//
  *dbdx = dbdr * drdx + dbds * dsdx;
  *dbdy = dbdr * drdy + dbds * dsdy;

  return;
}
//****************************************************************************80

void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_A sets the quadrature rule for assembly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Output, double WQ[3], quadrature weights.
//
//    Output, double XQ[3*ELEMENT_NUM], YQ[3*ELEMENT_NUM], the
//    coordinates of the quadrature points in each element.
//
{
  int element;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  wq[0] = 1.0 / 3.0;
  wq[1] = wq[0];
  wq[2] = wq[0];

  for ( element = 0; element < element_num; element++ )
  {
    ip1 = element_node[0+element*nnodes];
    ip2 = element_node[1+element*nnodes];
    ip3 = element_node[2+element*nnodes];

    x1 = node_xy[0+ip1*2];
    x2 = node_xy[0+ip2*2];
    x3 = node_xy[0+ip3*2];

    y1 = node_xy[1+ip1*2];
    y2 = node_xy[1+ip2*2];
    y3 = node_xy[1+ip3*2];

    xq[0+element*3] = 0.5 * ( x1 + x2 );
    xq[1+element*3] = 0.5 * ( x2 + x3 );
    xq[2+element*3] = 0.5 * ( x1 + x3 );

    yq[0+element*3] = 0.5 * ( y1 + y2 );
    yq[1+element*3] = 0.5 * ( y2 + y3 );
    yq[2+element*3] = 0.5 * ( y1 + y3 );
  }
  return;
}
//****************************************************************************80*

void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] )

//****************************************************************************80*
//
//  Purpose:
//
//    QUAD_E sets a quadrature rule for the error calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
//    index of local node I in element J.
//
//    Input, int ELEMENT, the index of the element for which the quadrature
//    points are to be computed.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int NNODES, the number of nodes used to form one element.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NQE, the number of points in the quadrature rule.
//    This is actually fixed at 13.
//
//    Output, double WQE[NQE], the quadrature weights.
//
//    Output, double XQE[NQE], YQE[NQE], the X and Y coordinates
//    of the quadrature points.
//
{
  int i;
  int ii;
  int iii;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;
  double z4;
  double z5;
  double z6;
  double z7;

  for ( i = 0; i < 3; i++ )
  {
    wqe[i] = 0.175615257433204;
    ii = i + 3;
    wqe[ii] = 0.053347235608839;
    ii = i + 6;
    iii = ii + 3;
    wqe[ii] = 0.077113760890257;
    wqe[iii] = wqe[ii];
  }

  wqe[12] = -0.14957004446767;

  z1 = 0.479308067841923;
  z2 = 0.260345966079038;
  z3 = 0.869739794195568;
  z4 = 0.065130102902216;
  z5 = 0.638444188569809;
  z6 = 0.312865496004875;
  z7 = 0.048690315425316;

  ip1 = element_node[0+element*nnodes];
  ip2 = element_node[1+element*nnodes];
  ip3 = element_node[2+element*nnodes];

  x1 = node_xy[0+ip1*2];
  x2 = node_xy[0+ip2*2];
  x3 = node_xy[0+ip3*2];

  y1 = node_xy[1+ip1*2];
  y2 = node_xy[1+ip2*2];
  y3 = node_xy[1+ip3*2];

  xqe[ 0] = z1 * x1 + z2 * x2 + z2 * x3;
  yqe[ 0] = z1 * y1 + z2 * y2 + z2 * y3;
  xqe[ 1] = z2 * x1 + z1 * x2 + z2 * x3;
  yqe[ 1] = z2 * y1 + z1 * y2 + z2 * y3;
  xqe[ 2] = z2 * x1 + z2 * x2 + z1 * x3;
  yqe[ 2] = z2 * y1 + z2 * y2 + z1 * y3;
  xqe[ 3] = z3 * x1 + z4 * x2 + z4 * x3;
  yqe[ 3] = z3 * y1 + z4 * y2 + z4 * y3;
  xqe[ 4] = z4 * x1 + z3 * x2 + z4 * x3;
  yqe[ 4] = z4 * y1 + z3 * y2 + z4 * y3;
  xqe[ 5] = z4 * x1 + z4 * x2 + z3 * x3;
  yqe[ 5] = z4 * y1 + z4 * y2 + z3 * y3;
  xqe[ 6] = z5 * x1 + z6 * x2 + z7 * x3;
  yqe[ 6] = z5 * y1 + z6 * y2 + z7 * y3;
  xqe[ 7] = z5 * x1 + z7 * x2 + z6 * x3;
  yqe[ 7] = z5 * y1 + z7 * y2 + z6 * y3;
  xqe[ 8] = z6 * x1 + z5 * x2 + z7 * x3;
  yqe[ 8] = z6 * y1 + z5 * y2 + z7 * y3;
  xqe[ 9] = z6 * x1 + z7 * x2 + z5 * x3;
  yqe[ 9] = z6 * y1 + z7 * y2 + z5 * y3;
  xqe[10] = z7 * x1 + z5 * x2 + z6 * x3;
  yqe[10] = z7 * y1 + z5 * y2 + z6 * y3;
  xqe[11] = z7 * x1 + z6 * x2 + z5 * x3;
  yqe[11] = z7 * y1 + z6 * y2 + z5 * y3;
  xqe[12] = ( x1 + x2 + x3 ) / 3.0;
  yqe[12] = ( y1 + y2 + y3 ) / 3.0;

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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" value.
//
{
  return ( double ) HUGE_VAL;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
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
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
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
//        X        R8_NINT
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

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
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
//****************************************************************************80*

double rhs ( double x, double y, double time )

//****************************************************************************80*
//
//  Purpose:
//
//    RHS gives the right-hand side of the differential equation.
//
//  Discussion:
//
//    The function specified here depends on the problem being
//    solved.  This is one of the routines that a user will
//    normally want to change.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point
//    in the region, at which the right hand side of the
//    differential equation is to be evaluated.
//
//    Input, double TIME, the current time.
//
//    Output, double RHS, the value of the right
//    hand side of the differential equation at (X,Y).
//
{
# define PI 3.141592653589793

  double ut;
  double uxx;
  double uyy;
  double value;

  ut =            - sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
  uxx = - PI * PI * sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
  uyy = - PI * PI * sin ( PI * x ) * sin ( PI * y ) * exp ( - time );

  value = ut - uxx - uyy;

  return value;
# undef PI
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
  char* t;

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
//****************************************************************************80*

void solution_write ( int node_num, double u[], char *u_file_name )

//****************************************************************************80*
//
//  Purpose:
//
//    SOLUTION_WRITE writes the solution to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double U[NODE_NUM], the coefficients of the solution.
//
//    Input, char *U_FILE_NAME, the name of the file
//    in which the data should be stored.
//
{
  int node;
  ofstream u_file;

  u_file.open ( u_file_name );

  if ( !u_file )
  {
    cout << "\n";
    cout << "SOLUTION_WRITE - Warning!\n";
    cout << "  Could not write the solution file \""
         << u_file_name << "\".\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    u_file << setw(14) << u[node] << "\n";
  }

  u_file.close ( );

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
//    24 September 2003
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

void triangulation_order6_plot ( char *file_name, int node_num, double node_xy[],
  int tri_num, int triangle_node[], int node_show, int triangle_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a pointset.
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
//    Input, double precision NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRI_NUM], lists, for each triangle,
//    the indices of the points that form the vertices and midsides
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
  int order[6] = { 0, 3, 1, 4, 2, 5 };
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
        << node     << ") show\n";
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

    for ( triangle = 0; triangle < tri_num; triangle++ )
    {
      node = triangle_node[order[0]+triangle*6];

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
        ip1 = ( i % 6 );
        node = triangle_node[order[ip1]+triangle*6];

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

    for ( triangle = 0; triangle < tri_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 0; i < 6; i++ )
      {
        node = triangle_node[i+triangle*6];
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
                << "moveto (" << triangle << ") show\n";
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

void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    XY_SET sets the XY coordinates of the nodes.
//
//  Discussion:
//
//    The nodes are laid out in an evenly spaced grid, in the unit square.
//
//    The first node is at the origin.  More nodes are created to the
//    right until the value of X = 1 is reached, at which point
//    the next layer is generated starting back at X = 0, and an
//    increased value of Y.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2004
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of elements in the X and
//    Y direction.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double XL, XR, YB, YT, the X coordinates of
//    the left and right sides of the rectangle, and the Y coordinates
//    of the bottom and top of the rectangle.
//
//    Output, double NODE_XY[2*NODE_NUM], the nodes.
//
{
  int i;
  int j;

  for ( j = 0; j < 2*ny - 1; j++ )
  {
    for ( i = 0; i < 2*nx - 1; i++ )
    {
      node_xy[0+(i+j*(2*nx-1))*2] =
        ( double ( 2 * nx - i - 2 ) * xl
        + double (          i     ) * xr )
        / double ( 2 * nx     - 2 );

      node_xy[1+(i+j*(2*nx-1))*2] =
        ( double ( 2 * ny - j - 2 ) * yb
        + double (          j     ) * yt )
        / double ( 2 * ny     - 2 );
    }
  }
  return;
}
