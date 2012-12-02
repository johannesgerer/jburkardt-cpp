#include <string>
#include <netcdfcpp.H>
#include <cstdlib>
#include <vector>

using namespace std;


int netcdf_mpas_read_dim ( string filename, string dim_name ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_DIM gets the size of the dimension with name dim_name
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
	//    John Burkardt, Doug Jacobsen
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
	//    Input, string DIM_NAME, the name of the dimension in question
	//
	//    Output, int NETCDF_MPAS_READ_DIM, the value of the dimension.
	//
	int ntime;
	long int dim_size;
	NcDim *dim_id;
	bool valid;
	string tmp_name;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	valid = false;

	for(int i = 0; i < ncid.num_dims() && !valid; i++){
		tmp_name = ncid.get_dim(i)->name();

		if(dim_name == tmp_name){
			valid = true;
		}
	}

	if(valid){
		dim_id = ncid.get_dim(dim_name.c_str());
		dim_size = (*dim_id).size ( );
	} else {
		dim_size = 1;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	ntime = ( int ) dim_size;

	return ntime;
}/*}}}*/
//****************************************************************************80
int netcdf_mpas_read_num_vars(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_NUM_VARS gets the number of variables in the netcdf file
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
	//    John Burkardt, Doug Jacobsen
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
	//    Output, int NETCDF_MPAS_READ_NUM_VARS, the number of variables
	//
	long int var_size;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get Ntime, which is a NETCDF dimension.
	//
	var_size = ncid.num_vars();

	ncid.close();

	return var_size;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[], double ycell[], double zcell[] ){/*{{{*/

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
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_xyzvertex ( string filename, int nvertices, double xvertex[], double yvertex[], double zvertex[] ){/*{{{*/
	//****************************************************************************80

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
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_xyzedge ( string filename, int nedges, double xedge[], double yedge[], double zedge[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_EDGES gets the edge midpoint coordinates.
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
	//    John Burkardt, Doug Jacobsen
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
	//    Output, double XEDGE[NEDGES], YEDGE[NEDGES],
	//    ZEDGE[NEDGES], the coordinates of the edge midpoints.
	//
	NcVar *var_id;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var ( "xEdge" );
	(*var_id).get ( &xedge[0], nedges );
	var_id = ncid.get_var ( "yEdge" );
	(*var_id).get ( &yedge[0], nedges );
	var_id = ncid.get_var ( "zEdge" );
	(*var_id).get ( &zedge[0], nedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_verticesoncell ( string filename, int maxedges, int ncells, int verticesOnCell[] ){/*{{{*/

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
	NcVar *var_id;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var ( "verticesOnCell" );
	(*var_id).get ( &verticesOnCell[0], ncells, maxedges );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_verticesonedge ( string filename, int nedges, int verticesonedge[] ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_VERTICESONEDGE gets the verticesOnEdge information.
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
	//    John Burkardt, Doug Jacobsen
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
	//    Output, int VERTICESONEDGE[2*NEDGES];
	//
	NcVar *var_id;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var ( "verticesOnEdge" );
	(*var_id).get ( &verticesonedge[0], nedges, 2 );
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices, int cellsonvertex[] ){/*{{{*/
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
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_cellsonedge ( string filename, int nedges, int cellsonedge[] ){/*{{{*/
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
	//    30 December 2010
	//
	//  Author:
	//
	//    John Burkardt, Doug Jacobsen
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
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, int nedgesoncell[] ){/*{{{*/

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
}/*}}}*/
//****************************************************************************80
int netcdf_mpas_list_ncell_fields(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_LIST_NCELL_FIELDS lists all fields that are printable based on nCells
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
	//    Doug Jacobsen
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
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	NcVar *var_id;
	NcDim *ncell_dim, *nvert_dim, *nedge_dim;
	NcToken ncell_tok, nvert_tok, nedge_tok;
	NcToken dim_name;
	vector<NcToken> var_names;
	vector<int> var_ids;
	string name;
	int included;
	int has_ncells;
	int num_vars;
	int num_dims;
	int i;
	int j;
	int valid = 0;
	int choice = 0;

	vector<string> excluded_vars;

	excluded_vars.push_back("xCell");
	excluded_vars.push_back("yCell");
	excluded_vars.push_back("zCell");
	excluded_vars.push_back("latCell");
	excluded_vars.push_back("lonCell");
	excluded_vars.push_back("indexToCellID");
	excluded_vars.push_back("edgesOnCell");
	excluded_vars.push_back("cellsOnCell");
	excluded_vars.push_back("localVerticalUnitVectors");
	excluded_vars.push_back("verticesOnCell");

	num_vars = ncid.num_vars();
	ncell_dim = ncid.get_dim("nCells");
	ncell_tok = ncell_dim->name();
	nvert_dim = ncid.get_dim("nVertices");
	nvert_tok = nvert_dim->name();
	nedge_dim = ncid.get_dim("nEdges");
	nedge_tok = nedge_dim->name();

	for(i = 0; i < num_vars; i++){
		var_id = ncid.get_var(i);
		num_dims = var_id->num_dims();
		included = 1;
		name = var_id->name();
		for(j = 0; j < excluded_vars.size(); j++){
			if(excluded_vars.at(j) == name){
				included = 0;
			}
		}

		if(included){
			has_ncells = 0;

			for(j = 0; j < num_dims; j++){
				dim_name = (var_id->get_dim(j))->name();
				if(ncell_tok == dim_name){
					has_ncells = 1;
				}
			}

			if(has_ncells){
				var_names.push_back(var_id->name());
				var_ids.push_back(i);
			}
		}
	}

	valid = 0;

	if(var_names.empty()){
		cout << endl;

		cout << "No variables indexed by nCells to color with." << endl;
	}

	while(!valid && !var_names.empty()){
		cout << endl;
		cout << "Available variables indexed by nCells." << endl;
		for(i = 0; i < var_names.size(); i++){
			cout << var_ids.at(i) << "\t" << var_names.at(i) << ":\t";
			num_dims = ncid.get_var(var_ids.at(i))->num_dims();
			for(j = 0; j < num_dims-1; j++){
				cout << ncid.get_var(var_ids.at(i))->get_dim(j)->name() << "*";
			}
			cout << ncid.get_var(var_ids.at(i))->get_dim(num_dims-1)->name() << endl;
		}
		cout << endl << endl;

		cout << "Enter the number of the field you would like to print" << endl;
		cout << "Choose from list above, or enter 0 for default:" << endl;
		cin >> choice;

		for(i = 0; i < var_ids.size(); i++){
			if(var_ids.at(i) == choice){
				valid = 1;
			}
		}

		if(choice == 0){
			valid = 1;
		}
	}

	return choice;
}/*}}}*/
//****************************************************************************80
int netcdf_mpas_list_nvertex_fields(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_LIST_NVERTEX_FIELDS lists all fields that are printable based on nVertices
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
	//    Doug Jacobsen
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
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	NcVar *var_id;
	NcDim *ncell_dim, *nvert_dim, *nedge_dim;
	NcToken ncell_tok, nvert_tok, nedge_tok;
	NcToken dim_name;
	string name;
	vector<NcToken> var_names;
	vector<int> var_ids;
	int included;
	int has_nverts;
	int num_vars;
	int num_dims;
	int i;
	int j;
	int choice = 0;
	int valid = 0;

	vector<string> excluded_vars;

	excluded_vars.push_back("xVertex");
	excluded_vars.push_back("yVertex");
	excluded_vars.push_back("zVertex");
	excluded_vars.push_back("latVertex");
	excluded_vars.push_back("lonVertex");
	excluded_vars.push_back("indexToVertexID");
	excluded_vars.push_back("edgesOnVertex");
	excluded_vars.push_back("cellsOnVertex");
	excluded_vars.push_back("kiteAreasOnVertex");

	excluded_vars.push_back("verticesOnVertex");
	excluded_vars.push_back("kiteAreaVertex");

	num_vars = ncid.num_vars();
	ncell_dim = ncid.get_dim("nCells");
	ncell_tok = ncell_dim->name();
	nvert_dim = ncid.get_dim("nVertices");
	nvert_tok = nvert_dim->name();
	nedge_dim = ncid.get_dim("nEdges");
	nedge_tok = nedge_dim->name();

	for(i = 0; i < num_vars; i++){
		var_id = ncid.get_var(i);
		num_dims = var_id->num_dims();
		included = 1;
		name = var_id->name();
		for(j = 0; j < excluded_vars.size(); j++){
			if(excluded_vars.at(j) == name){
				included = 0;
			}
		}

		if(included){
			has_nverts = 0;

			for(j = 0; j < num_dims; j++){
				dim_name = (var_id->get_dim(j))->name();
				if(nvert_tok == dim_name){
					has_nverts = 1;
				}
			}

			if(has_nverts){
				var_names.push_back(var_id->name());
				var_ids.push_back(i);
			}
		}
	}

	valid = 0;

	if(var_names.empty()){
		cout << endl;
		cout << "No variables indexed by nVertices to color with." << endl;
	}

	while(!valid && !var_names.empty()){
		cout << endl;
		cout << "Available variables indexed by nVertices." << endl;
		for(i = 0; i < var_names.size(); i++){
			cout << var_ids.at(i) << "\t" << var_names.at(i) << ":\t";
			num_dims = ncid.get_var(var_ids.at(i))->num_dims();
			for(j = 0; j < num_dims-1; j++){
				cout << ncid.get_var(var_ids.at(i))->get_dim(j)->name() << "*";
			}
			cout << ncid.get_var(var_ids.at(i))->get_dim(num_dims-1)->name() << endl;
		}
		cout << endl << endl;

		cout << "Enter the number of the field you would like to print" << endl;
		cout << "Choose from list above, or enter 0 for default:" << endl;
		cin >> choice;

		for(i = 0; i < var_ids.size(); i++){
			if(var_ids.at(i) == choice){
				valid = 1;
			}
		}

		if(choice == 0){
			valid = 1;
		}
	}

	return choice;

}/*}}}*/
//****************************************************************************80
int netcdf_mpas_list_nedge_fields(string filename){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_LIST_NEDGE_FIELDS lists all fields that are printable based on nEdges
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
	//    Doug Jacobsen
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
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	NcVar *var_id;
	NcDim *ncell_dim, *nvert_dim, *nedge_dim;
	NcToken ncell_tok, nvert_tok, nedge_tok;
	NcToken dim_name;
	string name;
	vector<NcToken> var_names;
	vector<int> var_ids;
	int included;
	int has_nedges;
	int num_vars;
	int num_dims;
	int i;
	int j;
	int choice = 0;
	int valid = 0;

	vector<string> excluded_vars;

	excluded_vars.push_back("xEdge");
	excluded_vars.push_back("yEdge");
	excluded_vars.push_back("zEdge");
	excluded_vars.push_back("latEdge");
	excluded_vars.push_back("lonEdge");
	excluded_vars.push_back("indexToEdgeID");
	excluded_vars.push_back("cellsOnEdge");
	excluded_vars.push_back("edgesOnEdge");
	excluded_vars.push_back("edgeNormalVectors");
	excluded_vars.push_back("cellTangentPlane");
	excluded_vars.push_back("weightsOnEdge");
	excluded_vars.push_back("verticesOnEdge");

	num_vars = ncid.num_vars();
	ncell_dim = ncid.get_dim("nCells");
	ncell_tok = ncell_dim->name();
	nvert_dim = ncid.get_dim("nVertices");
	nvert_tok = nvert_dim->name();
	nedge_dim = ncid.get_dim("nEdges");
	nedge_tok = nedge_dim->name();

	for(i = 0; i < num_vars; i++){
		var_id = ncid.get_var(i);
		num_dims = var_id->num_dims();
		included = 1;
		name = var_id->name();
		for(j = 0; j < excluded_vars.size(); j++){
			if(excluded_vars.at(j) == name){
				included = 0;
			}
		}

		if(included){
			has_nedges = 0;

			for(j = 0; j < num_dims; j++){
				dim_name = (var_id->get_dim(j))->name();
				if(nedge_tok == dim_name){
					has_nedges = 1;
				}
			}

			if(has_nedges){
				var_names.push_back(var_id->name());
				var_ids.push_back(i);
			}
		}
	}

	if(var_names.empty()){
		cout << endl;
		cout << "No variables indexed by nEdges to color with." << endl;
	}

	while(!valid && !var_names.empty()){
		cout << endl;
		cout << "Available variables indexed by nEdges." << endl;
		for(i = 0; i < var_names.size(); i++){
			cout << var_ids.at(i) << "\t" << var_names.at(i) << ":\t";
			num_dims = ncid.get_var(var_ids.at(i))->num_dims();
			for(j = 0; j < num_dims-1; j++){
				cout << ncid.get_var(var_ids.at(i))->get_dim(j)->name() << "*";
			}
			cout << ncid.get_var(var_ids.at(i))->get_dim(num_dims-1)->name() << endl;
		}
		cout << endl << endl;
		cout << "Enter the number of the field you would like to print" << endl;
		cout << "Choose from list above, or enter 0 for default:" << endl;
		cin >> choice;

		for(i = 0; i < var_ids.size(); i++){
			if(var_ids.at(i) == choice){
				valid = 1;
			}
		}

		if(choice == 0){
			valid = 1;
		}
	}
	return choice;
}/*}}}*/
//****************************************************************************80
int netcdf_mpas_field_num_dims(string filename, int id){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_FIELD_NUM_DIMS gets the number of dimensions for a field
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
	//    John Burkardt, Doug Jacobsen
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
	//    Input, int id, the variable id
	//
	//    Output, double dims, the number of dimensions for the field
	//
	NcVar *var_id;
	int num_dims;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var (id);
	if(!var_id->is_valid()){
		cout << "Field " << id << " doesn't exist." << endl;
		cout << "Exiting." << endl;
		exit(1);
	}
	num_dims = var_id->num_dims();
	ncid.close ( );

	return num_dims;
}/*}}}*/
//****************************************************************************80
int netcdf_mpas_field_num_items(string filename, int id){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_FIELD_NUM_items gets the number of items for a field
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
	//    John Burkardt, Doug Jacobsen
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
	//    Input, int id, the variable id
	//
	//    Output, int num_items, the number of items for the field
	//
	NcVar *var_id;
	int num_dims;
	int num_items;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var (id);
	if(!var_id->is_valid()){
		cout << "Field " << id << " doesn't exist." << endl;
		cout << "Exiting." << endl;
		exit(1);
	}
	num_dims = var_id->num_dims();

	num_items = var_id->get_dim(0)->size();

	for(int i = 1; i < num_dims; i++){
		num_items *= var_id->get_dim(i)->size();
	}
	ncid.close ( );

	return num_items;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_print_field_info(string filename, int id){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_PRINT_FIELD_INFO prints a fields information, like dimensions and name
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
	//    John Burkardt, Doug Jacobsen
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
	//    Input, int id, the variable id
	//
	NcVar *var_id;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var (id);

	if(!var_id->is_valid()){
		cout << "netcdf_mpas_print_field_info: field " << id << " doesn't exist." << endl;
		cout << "Exiting." << endl;
		exit(1);
	}

	cout << endl << endl;
	cout << "Field " << id << " is " << var_id->name() << endl;
	cout << "Dimensions: ";
	for(int i = 0; i < var_id->num_dims()-1; i++){
		cout << var_id->get_dim(i)->name() << "*";
	}
	cout << var_id->get_dim(var_id->num_dims()-1)->name() << endl;
	cout << "Variable type: " << var_id->type() << endl;

	//
	//  Close the file.
	//
	ncid.close ( );

	return;

}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_field(string filename, int id, double values[], int cur_time, int cur_level){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_1D_FIELD gets a field based on one dimension.
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
	//    Doug Jacobsen
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
	//    Input, int id, the variable id
	//
	//    Output, double values[] are the values of the field from the netcdf file
	//
	//	  Output, long dims[] the dimensions of the field from the netcdf file
	NcVar *var_id;
	int num_dims;
	int type;
	long *dims;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var (id);
	num_dims = var_id->num_dims();
	dims = new long[num_dims];
	for(int i = 0; i < num_dims; i++){
		dims[i] = var_id->get_dim(i)->size();
	}

	type = var_id->type();
	if(type == 6){ // NCDOUBLE
		switch(num_dims){
			case 2:
				(*var_id).set_cur(0, cur_level);
				(*var_id).get ( &values[0], dims[0], 1);
				break;
			case 3:
				(*var_id).set_cur(cur_time, 0, cur_level);
				(*var_id).get ( &values[0], 1, dims[1], 1);
				break;
			default: //1 dim
				(*var_id).get ( &values[0], dims[0]);
				break;
		}
	} else if(type == 4){ //NCINT
		int *tmp_values;
		int num_items;

		switch(num_dims){
			case 2:
				num_items = dims[0];
				tmp_values = new int[num_items];
				(*var_id).set_cur(0, cur_level);
				(*var_id).get ( &tmp_values[0], dims[0], 1);
				break;
			case 3:
				num_items = dims[1];
				tmp_values = new int[num_items];
				(*var_id).set_cur(cur_time, 0, cur_level);
				(*var_id).get ( &tmp_values[0], 1, dims[1], 1);
				break;
			default: //1 dim
				num_items = dims[0];
				tmp_values = new int[num_items];
				(*var_id).get ( &tmp_values[0], dims[0]);
				break;
		}

		for(int i = 0; i < num_items; i++){
			values[i] = tmp_values[i]*1.0;
		}

		delete [] tmp_values;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
//****************************************************************************80
void netcdf_mpas_read_full_field(string filename, int id, double values[]){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    NETCDF_MPAS_READ_FULL_FIELD gets a full field.
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
	//    Doug Jacobsen
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
	//    Input, int id, the variable id
	//
	//    Output, double values[] are the values of the field from the netcdf file
	//
	//	  Output, long dims[] the dimensions of the field from the netcdf file
	NcVar *var_id;
	int num_dims;
	int type;
	long *dims;
	//
	//  Open the file.
	//
	NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
	//
	//  Get the variable values.
	//
	var_id = ncid.get_var (id);
	num_dims = var_id->num_dims();
	dims = new long[num_dims];
	for(int i = 0; i < num_dims; i++){
		dims[i] = var_id->get_dim(i)->size();
	}

	type = var_id->type();
	if(type == 6){ // NCDOUBLE
		switch(num_dims){
			case 2:
				(*var_id).get ( &values[0], dims[0], dims[1]);
				break;
			case 3:
				(*var_id).get ( &values[0], dims[0], dims[1], dims[2]);
				break;
			default: //1 dim
				(*var_id).get ( &values[0], dims[0]);
				break;
		}
	} else if(type == 4){ //NCINT
		int *tmp_values;
		int num_items;

		switch(num_dims){
			case 2:
				num_items = dims[0]*dims[1];
				tmp_values = new int[num_items];
				(*var_id).get ( &tmp_values[0], dims[0], dims[1]);
				break;
			case 3:
				num_items = dims[0]*dims[1]*dims[2];
				tmp_values = new int[num_items];
				(*var_id).get ( &tmp_values[0], dims[0], dims[1], dims[2]);
				break;
			default: //1 dim
				num_items = dims[0];
				tmp_values = new int[num_items];
				(*var_id).get ( &tmp_values[0], dims[0]);
				break;
		}

		for(int i = 0; i < num_items; i++){
			values[i] = tmp_values[i]*1.0;
		}

		delete [] tmp_values;
	}
	//
	//  Close the file.
	//
	ncid.close ( );

	return;
}/*}}}*/
