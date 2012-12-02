# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "netcdf_mpas.hpp"
# include "netcdf.hpp"

int main ( );
void test01 ( string filename );
void test02 ( string filename );
void test03 ( string filename );
void test04 ( string filename );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_PRB tests the NETCDF_MPAS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
{
  string filename = "x1.642.grid.nc";

  timestamp ( );
  cout << "\n";
  cout << "NETCDF_MPAS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NETCDF_MPAS library.\n";
//
//  TEST01 won't run on my home computer because of difficulties
//  with the version of NETCDF installed there.
//
  if ( false )
  {
    test01 ( filename );
  }
  test02 ( filename );
  test03 ( filename );
  test04 ( filename );
//
//  Terminate.
//
  cout << "\n";
  cout << "NETCDF_MPAS_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 reads an MPAS NETCDF grid file and reports some information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be examined.
//
{
  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Open an MPAS NETCDF grid file,\n";
  cout << "  use inquire functions to examine information,\n";
  cout << "  and make a report.\n";

  netcdf_mpas_report ( filename );

  return;
}
//****************************************************************************80

void test02 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 extracts the cell center information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be examined.
//
{
  int i;
  int ncells;
  double *xcell;
  double *ycell;
  double *zcell;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Open an MPAS NETCDF grid file,\n";
  cout << "  find the number of cells,\n";
  cout << "  and get the cell center coordinates.\n";

  ncells = netcdf_mpas_read_ncells ( filename );

  xcell = new double[ncells];
  ycell = new double[ncells];
  zcell = new double[ncells];

  netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );

  cout << "\n";
  cout << "  First 10 cell centers:\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << setw(2)  << i
         << "  " << setw(10) << xcell[i]
         << "  " << setw(10) << ycell[i]
         << "  " << setw(10) << zcell[i] << "\n";
  }
  delete [] xcell;
  delete [] ycell;
  delete [] zcell;

  return;
}
//****************************************************************************80

void test03 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 extracts the cell centers and cells on edges.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be examined.
//
{
  int *cellsonedge;
  int i;
  int j;
  int ncells;
  int nedges;
  double *xcell;
  double *ycell;
  double *zcell;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Open an MPAS NETCDF grid file,\n";
  cout << "  find the number of cells and edges,\n";
  cout << "  get the cell center coordinates,\n";
  cout << "  and the cellsOnEdge array.\n";

  ncells = netcdf_mpas_read_ncells ( filename );

  xcell = new double[ncells];
  ycell = new double[ncells];
  zcell = new double[ncells];

  netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );

  nedges = netcdf_mpas_read_nedges ( filename );
  cout << "\n";
  cout << "  Number of edges NEDGES = " << nedges << "\n";

  cellsonedge = new int[2*nedges];

  netcdf_mpas_read_cellsonedge ( filename, nedges, cellsonedge );

  cout << "\n";
  cout << "  First ten entries of CELLSONEDGE:\n";
  cout << "\n";
  for ( j = 0; j < 10; j++ )
  {
    cout << "  " << setw(2) << j
         << "  " << setw(4) << cellsonedge[0+2*j]
         << "  " << setw(4) << cellsonedge[1+2*j] << "\n";
  }
  delete [] cellsonedge;
  delete [] xcell;
  delete [] ycell;
  delete [] zcell;

  return;
}
//****************************************************************************80

void test04 ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 extracts the cell centers and cells on vertices.
//
//  Discussion:
//
//    The data retrieved in this example should be sufficient to plot
//    the triangulation, including the triangular faces.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User"s Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file to be examined.
//
{
  int *cellsonvertex;
  int i;
  int j;
  int ncells;
  int nvertices;
  double *xcell;
  double *ycell;
  double *zcell;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Open an MPAS NETCDF grid file,\n";
  cout << "  find the number of cells and edges,\n";
  cout << "  get the cell center coordinates,\n";
  cout << "  and the cellsOnVertex array.\n";

  ncells = netcdf_mpas_read_ncells ( filename );

  xcell = new double[ncells];
  ycell = new double[ncells];
  zcell = new double[ncells];

  netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );

  nvertices = netcdf_mpas_read_nvertices ( filename );
  cout << "\n";
  cout << "  Number of verties NVERTICES = " << nvertices << "\n";

  cellsonvertex = new int[3*nvertices];

  netcdf_mpas_read_cellsonvertex ( filename, nvertices, cellsonvertex );

  cout << "\n";
  cout << "  First ten entries of CELLSONVERTEX:\n";
  cout << "\n";
  for ( j = 0; j < 10; j++ )
  {
    cout << "  " << setw(2) << j
         << "  " << setw(4) << cellsonvertex[0+3*j]
         << "  " << setw(4) << cellsonvertex[1+3*j]
         << "  " << setw(4) << cellsonvertex[2+3*j] << "\n";
  }
  delete [] cellsonvertex;
  delete [] xcell;
  delete [] ycell;
  delete [] zcell;

  return;
}


