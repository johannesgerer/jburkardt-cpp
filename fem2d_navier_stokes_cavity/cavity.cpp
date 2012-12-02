# include <cstdlib>
# include <cmath>

using namespace std;

void dirichlet_condition ( int n, double xy[], double u_bc[], double v_bc[], 
  double p_bc[] );
void rhs ( int n, double xy[], double u_rhs[], double v_rhs[], double p_rhs[] );

//******************************************************************************

void dirichlet_condition ( int n, double xy[], double u_bc[], double v_bc[], 
  double p_bc[] )

//******************************************************************************
//
//  Purpose:
//
//    DIRICHLET_CONDITION sets the value of a Dirichlet boundary condition.
//
//  Discussion:
//
//
//                           U = 1  V = 0
//
//                       1 +---------------+
//                         |               |
//                         |               |
//                         |               |
//            U = V = 0    |               | U = V = 0
//                         |               |
//                         |               |
//                         |               |
//                       0 +---------------+
//
//                         0               1
//
//                            U = V = 0
//
//  Modified:
//
//    08 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double XY[2*N], the coordinates of the points.
//
//    Output, double U_BC[N], V_BC[N], P_BC[N], the values of the
//    Dirichlet boundary conditions on horizontal velocity, vertical velocity,
//    and pressure.
//
{
  int node;
  double tol = 0.0001;
  double x;
  double y;

  for ( node = 0; node < n; node++ )
  {
    x = xy[0+node*2];
    y = xy[1+node*2];

    if ( fabs ( y - 1.0 ) <= tol )
    {
      u_bc[node] = 1.0;
    }
    else
    {
      u_bc[node] = 0.0;
    }

    v_bc[node] = 0.0;
    p_bc[node] = 0.0;
  }

  return;
}
//******************************************************************************

void rhs ( int n, double xy[], double u_rhs[], double v_rhs[], double p_rhs[] )

//******************************************************************************
//
//  Purpose:
//
//    RHS gives the right-hand side of the differential equation.
//
//  Modified:
//
//    08 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double XY[2*N], the coordinates of the points.
//
//    Output, double U_RHS[N], V_RHS[N], P_RHS[N], the right
//    hand sides of the differential equations at the points.
//
{
  int node;

  for ( node = 0; node < n; node++ )
  {
    p_rhs[node] = 0.0;
    u_rhs[node] = 0.0;
    v_rhs[node] = 0.0;
  }

  return;
}
