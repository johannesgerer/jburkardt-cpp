# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[], 
  double mass, double f[], double &pot, double &kin );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
void timestamp ( );
void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double dt );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MD.
//
//  Discussion:
//
//    MD implements a simple molecular dynamics simulation.
//
//    The velocity Verlet time integration scheme is used. 
//
//    The particles interact with a central pair potential.
//
//    This program is based on a FORTRAN90 program by Bill Magro.
//
//  Usage:
//
//    md nd np step_num dt
//
//    where
//
//    * nd is the spatial dimension (2 or 3);
//    * np is the number of particles (500, for instance);
//    * step_num is the number of time steps (500, for instance).
//    * dt is the time step (0.1 for instance)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 December 2014
//
//  Author:
//
//    John Burkardt.
//
{
  double *acc;
  double ctime;
  double dt;
  double e0;
  double *force;
  int i;
  int id;
  double kinetic;
  double mass = 1.0;
  int nd;
  int np;
  double *pos;
  double potential;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;

  timestamp ( );
  cout << "\n";
  cout << "MD\n";
  cout << "  C++ version\n";
  cout << "  A molecular dynamics program.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    nd = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter ND, the spatial dimension (2 or 3).\n";
    cin >> nd;
  }
//
//  Get the number of particles.
//
  if ( 2 < argc )
  {
    np = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter NP, the number of particles (500, for instance).\n";
    cin >> np;
  }
//
//  Get the number of time steps.
//
  if ( 3 < argc )
  {
    step_num = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter STEP_NUM, the number of time steps (500 or 1000, for instance).\n";
    cin >> step_num;
  }
//
//  Get the time step.
//
  if ( 4 < argc )
  {
    dt = atof ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter DT, the time step size (0.1, for instance).\n";
    cin >> dt;
  }
//
//  Report.
//
  cout << "\n";
  cout << "  ND, the spatial dimension, is " << nd << "\n";
  cout << "  NP, the number of particles in the simulation is " << np << "\n";
  cout << "  STEP_NUM, the number of time steps, is " << step_num << "\n";
  cout << "  DT, the size of each time step, is " << dt << "\n";
//
//  Allocate memory.
//
  acc = new double[nd*np];
  force = new double[nd*np];
  pos = new double[nd*np];
  vel = new double[nd*np];
//
//  This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//
  cout << "\n";
  cout << "  At each step, we report the potential and kinetic energies.\n";
  cout << "  The sum of these energies should be a constant.\n";
  cout << "  As an accuracy check, we also print the relative error\n";
  cout << "  in the total energy.\n";
  cout << "\n";
  cout << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
  cout << "                Energy P        Energy K       Relative Energy Error\n";
  cout << "\n";

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;
  
  ctime = cpu_time ( );

  for ( step = 0; step <= step_num; step++ )
  {
    if ( step == 0 )
    {
      initialize ( np, nd, pos, vel, acc );
    }
    else
    {
      update ( np, nd, pos, vel, force, acc, mass, dt );
    }

    compute ( np, nd, pos, vel, mass, force, potential, kinetic );

    if ( step == 0 )
    {
      e0 = potential + kinetic;
    }

    if ( step == step_print )
    {
      cout << "  " << setw(8) << step
           << "  " << setw(14) << potential
           << "  " << setw(14) << kinetic
           << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 << "\n";
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }

  }
//
//  Report timing.
//
  ctime = cpu_time ( ) - ctime;
  cout << "\n";
  cout << "  Elapsed cpu time " << ctime << " seconds.\n";
//
//  Free memory.
//
  delete [] acc;
  delete [] force;
  delete [] pos;
  delete [] vel;
//
//  Terminate.
//
  cout << "\n";
  cout << "MD\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void compute ( int np, int nd, double pos[], double vel[], double mass, 
  double f[], double &pot, double &kin )

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces and energies.
//
//  Discussion:
//
//    The computation of forces and energies is fully parallel.
//
//    The potential function V(X) is a harmonic well which smoothly
//    saturates to a maximum value at PI/2:
//
//      v(x) = ( sin ( min ( x, PI2 ) ) )**2
//
//    The derivative of the potential is:
//
//      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
//            = sin ( 2.0 * min ( x, PI2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2008
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double POS[ND*NP], the position of each particle.
//
//    Input, double VEL[ND*NP], the velocity of each particle.
//
//    Input, double MASS, the mass of each particle.
//
//    Output, double F[ND*NP], the forces.
//
//    Output, double &POT, the total potential energy.
//
//    Output, double &KIN, the total kinetic energy.
//
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pot = 0.0;
  kin = 0.0;

  for ( k = 0; k < np; k++ )
  {
//
//  Compute the potential energy and forces.
//
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist ( nd, pos+k*nd, pos+j*nd, rij );
//
//  Attribute half of the potential energy to particle J.
//
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pot = pot + 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < nd; i++ )
        {
          f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
//
//  Compute the kinetic energy.
//
    for ( i = 0; i < nd; i++ )
    {
      kin = kin + vel[i+k*nd] * vel[i+k*nd];
    }
  }

  kin = kin * 0.5 * mass;
  
  return;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double dist ( int nd, double r1[], double r2[], double dr[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIST computes the displacement (and its norm) between two particles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double R1[ND], R2[ND], the positions.
//
//    Output, double DR[ND], the displacement vector.
//
//    Output, double D, the Euclidean norm of the displacement.
//
{
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
//****************************************************************************80

void initialize ( int np, int nd, double pos[], double vel[], double acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes the positions, velocities, and accelerations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 December 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Output, double POS[ND*NP], the positions.
//
//    Output, double VEL[ND*NP], the velocities.
//
//    Output, double ACC[ND*NP], the accelerations.
//
{
  int i;
  int j;
  int seed;
//
//  Set the positions.
//
  seed = 123456789;
  r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed, pos );
//
//  Set the velocities.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
//
//  Set the accelerations.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }
  return;
}
//****************************************************************************80

void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
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
//    31 May 2001 09:45:54 AM
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

void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double dt )

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE updates positions, velocities and accelerations.
//
//  Discussion:
//
//    The time integration is fully parallel.
//
//    A velocity Verlet algorithm is used for the updating.
//
//    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
//    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
//    a(t+dt) = f(t) / m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input/output, double POS[ND*NP], the positions.
//
//    Input/output, double VEL[ND*NP], the velocities.
//
//    Input, double F[ND*NP], the forces.
//
//    Input/output, double ACC[ND*NP], the accelerations.
//
//    Input, double MASS, the mass of each particle.
//
//    Input, double DT, the time step.
//
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}
