# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "netcdf_mpas.h"
# include "netcdf.hpp"

//****************************************************************************80

void netcdf_mpas_read_cellsonedge ( string filename, int nedges,
  int cellsonedge[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_CELLSONEDGE gets the cellsOnEdge information.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NEDGES, the number of edges.
//
//    Output, int CELLSONEDGE[2*NEDGES];
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "cellsOnEdge" );
  (*var_id).get ( &cellsonedge[0], nedges, 2 );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices,
  int cellsonvertex[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_CELLSONVERTEX gets the cellsOnVertex information.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NEDGES, the number of edges.
//
//    Output, int CELLSONVERTEX[3*NVERTICES];
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "cellsOnVertex" );
  (*var_id).get ( &cellsonvertex[0], nvertices, 3 );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

int netcdf_mpas_read_maxedges ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_MAXEDGES gets MAXEDGES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_MAXEDGES, the value of MAXEDGES.
//
{
  int maxedges;
  long int maxedges_size;
  NcDim *maxedges_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get MAXEDGES, which is a NETCDF dimension.
//
  maxedges_id = ncid.get_dim ( "maxEdges" );

  maxedges_size = (*maxedges_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  maxedges = ( int ) maxedges_size;

  return maxedges;
}
//****************************************************************************80

int netcdf_mpas_read_ncells ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NCELLS gets the number of cells.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_NCELLS, the value of NCELLS.
//
{
  int ncells;
  long int ncells_size;
  NcDim *ncells_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get NCELLS, which is a NETCDF dimension.
//
  ncells_id = ncid.get_dim ( "nCells" );

  ncells_size = (*ncells_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  ncells = ( int ) ncells_size;

  return ncells;
}
//****************************************************************************80

int netcdf_mpas_read_nedges ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NEDGES gets the number of edges.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_NEDGES, the value of NEDGES.
//
{
  int nedges;
  long int nedges_size;
  NcDim *nedges_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get NCELLS, which is a NETCDF dimension.
//
  nedges_id = ncid.get_dim ( "nEdges" );

  nedges_size = (*nedges_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  nedges = ( int ) nedges_size;

  return nedges;
}
//****************************************************************************80

void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, 
  int nedgesoncell[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NEDGESONCELLS gets nedgesOnCells.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NCELLS, the number of cells.
//
//    Output, int NEDGESONCELLS[NCELLS];
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "nEdgesOnCell" );
  (*var_id).get ( &nedgesoncell[0], ncells );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

int netcdf_mpas_read_nvertices ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NVERTICES gets the number of vertices.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_NVERTICES, the value of NVERTICES.
//
{
  int nvertices;
  long int nvertices_size;
  NcDim *nvertices_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get NCELLS, which is a NETCDF dimension.
//
  nvertices_id = ncid.get_dim ( "nVertices" );

  nvertices_size = (*nvertices_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  nvertices = ( int ) nvertices_size;

  return nvertices;
}
//****************************************************************************80

void netcdf_mpas_read_verticesoncell ( string filename, int maxedges,
  int ncells, int verticesoncell[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_VERTICESONCELLS gets verticesOnCells.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int MAXEDGES, the maximum number of edges for a cell.
//
//    Input, int NCELLS, the number of cells.
//
//    Output, int VERTICESONCELLS[MAXEDGES*NCELLS];
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "verticesOnCell" );
  (*var_id).get ( &verticesoncell[0], ncells, maxedges );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[],
  double ycell[], double zcell[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_XYZCELL reads xCell, yCell, zCell.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NCELLS, the number of nodes.
//
//    Output, double XCELL[NCELLS], YCELL[NCELLS], ZCELL[NCELLS], the
//    coordinates of the nodes.
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "xCell" );
  (*var_id).get ( &xcell[0], ncells );
  var_id = ncid.get_var ( "yCell" );
  (*var_id).get ( &ycell[0], ncells );
  var_id = ncid.get_var ( "zCell" );
  (*var_id).get ( &zcell[0], ncells );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

void netcdf_mpas_read_xyzvertex ( string filename, int nvertices, 
  double xvertex[], double yvertex[], double zvertex[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_CELLS gets the cell center coordinates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NVERTICES, the number of vertices.
//
//    Output, double XVERTEX[NVERTICES], YVERTEXL[NVERTICES], 
//    ZVERTEX[NVERTICES], the coordinates of the nodes.
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "xVertex" );
  (*var_id).get ( &xvertex[0], nvertices );
  var_id = ncid.get_var ( "yVertex" );
  (*var_id).get ( &yvertex[0], nvertices );
  var_id = ncid.get_var ( "zVertex" );
  (*var_id).get ( &zvertex[0], nvertices );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

void netcdf_mpas_report ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_REPORT reads an MPAS NETCDF grid file and reports.
//
//  Discussion:
//
//    In this example, we want to extract all the information from a file
//    of unknown type.
//
//    Here, we are pretending we have no idea what's in the file, so we
//    have to go step by step, making inquiries to NETCDF.
//
//    As we go, we print out what we have discovered.  We don't attempt
//    to return any of the data.
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
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the NETCDF file to examine.
//
{
  NcAtt *att;
  NcDim *dim;
  int i;
  int j;
  long int len;

  NcToken name;
  int num_atts;
  int num_dims;
  int num_vars;
  NcType type;
  NcDim *unlimdimid;
  NcVar *var;

  cout << "\n";
  cout << "NETCDF_MPAS_REPORT:\n";
  cout << "  Report the information stored in a NETCDF\n";
  cout << "  file.  Although we wish to examine a file containing\n";
  cout << "  grid data generated by MPAS, we will assume we do not\n";
  cout << "  have any idea of what is in the file.  So we just read,\n";
  cout << "  inquire, and print out.\n";

  cout << "\n";
  cout << "  The name of the file is \"" << filename << "\"\n";
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Return information about the NETCDF file.
//
  num_dims = ncid.num_dims ( );
  num_vars = ncid.num_vars ( );
  num_atts = ncid.num_atts ( );
  unlimdimid = ncid.rec_dim ( );

  cout << "\n";
  cout << "PRIMARY PARAMETERS:\n";
  cout << "\n";
  cout << "  The number of dimensions         NUM_DIMS   = " << num_dims << "\n";
  cout << "  The number of variables          NUM_VARS   = " << num_vars << "\n";
  cout << "  The number of global attributes  NUM_ATTS   = " << num_atts << "\n";
  cout << "  The unlimited dimension (if any) UNLIMDIMID = \"" << (*unlimdimid).name ( ) << "\"\n";
//
//  Retrieve global attributes.
//  First, we must evaluate the constant "NC_GLOBAL".
//
//  nc_global = netcdf.getConstant ( 'NC_GLOBAL' );

  cout << "\n";
  cout << "GLOBAL ATTRIBUTES:\n";
  cout << " Att  --------Name--------  Type   Len\n";
  cout << "\n";
  for ( i = 0; i < num_atts; i++ )
  {
    att = ncid.get_att ( i );
    name = (*att).name ( );
    type = (*att).type ( );
    len = (*att).num_vals ( );
    cout << "  " << setw(2) << i
         << "  \"" << setw(18) << name << "\""
         << "  " << setw(4) << type
         << "  " << setw(4) << len << "\n";
  }
//
//  Determine names and extents of dimensions.
//  Since each NAME is a char array, we make a cell array to hold them.
//
  cout << "\n";
  cout << "DIMENSIONS:\n";
  cout << " Dim  --------Name--------  Extent\n";
  cout << "\n";

  for ( i = 0; i < num_dims; i++ )
  {
    dim = ncid.get_dim ( i );
    name = (*dim).name ( );
    len = (*dim).size ( );
    cout << "  " << setw(2) << i
         << "  \"" << setw(18) << name << "\""
         << "  " << setw(6) << len << "\n";
  }
//
//  Get variable names, types, dimensions, number of attributes.
//
  cout << "\n";
  cout << "VARIABLES:\n";
  cout << " Var  --------Name--------  Type Natts Ndims  Dims\n";
  cout << "\n";

  for ( i = 0; i < num_vars; i++ )
  {
    var = ncid.get_var ( i );
    name = (*var).name ( );
    type = (*var).type ( );
    num_dims = (*var).num_dims ( );
    num_atts = (*var).num_atts ( );
    cout << "  " << setw(2) << i
         << "  \"" << setw(18) << name << "\""
         << "  " << setw(4) << type
         << "  " << setw(4) << num_atts
         << "  " << setw(4) << num_dims
         << "  ";
    for ( j = 0; j < num_dims; j++ )
    {
      dim = (*var).get_dim ( j );
      if ( j == 0 )
      {
        cout << "[";
      }
      cout << (*dim).name ( );
      if ( j < num_dims - 1 )
      {
        cout <<",";
      }
      else
      {
        cout << "]";
      }
    }
    cout << "\n";
  }
//
//  Close the file.
//
  ncid.close ( );

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
