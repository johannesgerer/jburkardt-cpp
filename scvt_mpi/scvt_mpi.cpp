//**************************************************
// scvt-mpi.cpp
//
//  Purpose:
//   
//   mpi-scvt.cpp is used to compute spherical centroidal Voronoi tessellations using a modified Lloyd's algorithm
//   in parallel.
//
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license. 
// 
//  Modified:
//
//    02 February 2012
//
//  Author:
//
//   Doug Jacobsen
//   Geoff Womeldorff
//
//**************************************************

// Enable _DEBUG for output of routines. Useful for debugging any issues in mpi.
// All messages from turning this flag on are written to cerr
// #define _DEBUG

#include <stdlib.h>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <tr1/unordered_set>
#include <vector>
#include <math.h>
#include <assert.h>

#ifdef USE_NETCDF
#include <netcdfcpp.h>
#endif

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

# include "triangle.h"
# include "triangulation.h"

#define SEED	3729

// Uses boost mpi to make use of serialization routines for packing and unpacking of classes
namespace mpi = boost::mpi;

typedef boost::optional<mpi::status> optional;

// Uses namespaces std and tr1. tr1 is used for unordered_set which gives unique triangulation at the end.
using namespace std;
using namespace tr1;

class bpt {/*{{{*/
	public:
struct bisect_hasher {/*{{{*/
	size_t operator()(const pair<int,int> &p) const {
		uint32_t hash; 
		size_t i, key[2] = { (size_t)p.first, (size_t)p.second};
		for(hash = i = 0; i < sizeof(key); ++i) {
			hash += ((uint8_t *)key)[i];
			hash += (hash << 10);
			hash ^= (hash >> 6);
		}
		hash += (hash << 3);
		hash ^= (hash >> 11);
		hash += (hash << 15);
		return hash;
	}
};/*}}}*/
};/*}}}*/
class region{/*{{{*/
	private:
		friend class boost::serialization::access;
		template<class Archive>
			void serialize(Archive & ar, const unsigned int version)
			{
				ar & center;
				ar & radius;
				ar & input_radius;
				ar & triangles;
				ar & neighbors;
				ar & neighbors1;
				ar & neighbors2;
				ar & boundary_points;
				ar & loop_start;
				ar & loop_stop;
			}

	public:
		pnt center;
		double radius;
		double input_radius;
		vector<pnt> points;
		vector<tri> triangles;
		vector<int> neighbors; // First Level of Neighbors
		vector<int> neighbors1; // First Level of Neighbors + Self
		vector<int> neighbors2; // Second Level of Neighbors + First Level of Neighbors + Self
		vector<pnt> boundary_points;
		vector<int> loop_start; // beginning point in loop
		vector<int> loop_stop; // ending point in loop
};/*}}}*/

struct int_hasher {/*{{{*/
	  size_t operator()(const int v) const { return v; }
};/*}}}*/
class mpi_timer{/*{{{*/
	public:
		mpi::timer my_timer;
		double total_time;
		int num_calls;
		string name;

		mpi_timer() : total_time(0), num_calls(0), name("Default") { };
		mpi_timer(string in_name) : total_time(0), num_calls(0), name(in_name) { };

		mpi_timer operator+(const mpi_timer &t) const {/*{{{*/
			mpi_timer nt;
			nt.init("Sum");
			nt.num_calls = t.num_calls;
			nt.total_time = t.total_time + total_time;

			return nt;
		}/*}}}*/
		void init(string in_name){/*{{{*/
			total_time = 0.0;
			num_calls = 0;
			name = in_name;
		}/*}}}*/
		void init(){/*{{{*/
			total_time = 0;
			num_calls = 0;
		}/*}}}*/
		void start(){/*{{{*/
			my_timer.restart();
		}/*}}}*/
		void stop(){/*{{{*/
			total_time += my_timer.elapsed();
			num_calls++;
		}/*}}}*/
};/*}}}*/

inline std::ostream & operator<<(std::ostream &os, const mpi_timer &t){/*{{{*/
	if(t.num_calls > 0){
	os << t.name << ": " << t.total_time*1e3 << " (ms), " << (t.total_time/t.num_calls)*1e3 << " (ms). Called " << t.num_calls << " times." << endl;
	} else {
	os << t.name << ": Never called. Time = 0.0 (ms)" << endl;
	}
	return os;
}/*}}}*/

// Triangles input flags
string flags_str = "QBPIOYYiz";
char * flags;

//Processor information and global communicator
const int master = 0;
int id, num_procs;
mpi::communicator world;

// Message tags for sending and receiving non-blocking messages
enum {msg_points, msg_tri_print , msg_restart, msg_ave, msg_max, msg_l1};

// Sort types
enum {sort_dot, sort_vor};

// Global constants
int points_begin = 0;
int num_pts = 162;
int num_bdry = 0;
int max_it = 100;
int max_it_no_proj = 100;
int max_it_scale_alpha = 0;
int div_levs = 1;
int num_bisections = 0;
int conv = 0;
int restart = 0;
int sort_method = sort_dot;
double min_bdry_angle = 1.0;
double eps = 1.0E-10;
double proj_alpha;
double max_resolution = 4.0;

//gw: restart mode type and variable (move to a header?)
enum restart_mode_type { RESTART_OVERWRITE, RESTART_RETAIN };
restart_mode_type restart_mode = RESTART_OVERWRITE;
enum fileio_mode_type { FILEIO_TXT, FILEIO_NETCDF, FILEIO_BOTH };
fileio_mode_type fileio_mode = FILEIO_TXT;

//Define variable for quadrature rules
int quad_rule = 2;
string quad_names[6] = {"Centroid Rule", "Vertex Rule", "Midpoint Rule", "7 Point Rule", "13 Point Rule", "19 Point Rule"};

double *wq, *q;
double wq7[7], q7[4];
double wq13[13], q13[8];
double wq19[19], q19[12];

//Define timers for performance studies
const int num_timers = 8;
const int num_global_timers = 4;

mpi_timer my_timers[num_timers];
mpi_timer global_timers[num_global_timers];
string names[num_timers] = {"Total", "Iteration", "Triangulation", "Integration", "Metrics", "Communication", "Convergence Check", "Sort"};
string global_names[num_global_timers] = {"Global Time", "Final Gather", "Final Triangulation", "Final Bisection"};

//Each processor has a copy of all points, and has a n_points (new points) vector for which points it should update.
vector<pnt> points;
vector<pnt> n_points;
vector<pnt>::iterator point_itr;

vector<pnt> boundary_points;
vector<pnt>::iterator boundary_itr;

vector<int> loop_start;
vector<int> loop_stop;

//Each processor has a list of all regions, as well as it's own regions (only one per processor currently)
vector<region> my_regions;
vector<region> regions;
vector<region>::iterator region_itr;

vector<tri> all_triangles;

//Region neighbors hold the connectivity of neighbors, for communication purposes.
vector<unordered_set<int, int_hasher> > region_neighbors;
unordered_set<int, int_hasher>::iterator region_neigh_itr;

//Iterators for neighbors and triangles.
vector<int>::iterator neighbor_itr;
vector<tri>::iterator tri_itr;

/* ***** Setup Routines *****{{{*/
void readParamsFile();
void buildRegions();
void printRegions();
void readBoundaries();
/*}}}*/
/* ***** Bisect Edges Routines *****{{{*/
void bisectEdges(int end);
void bisectTriangulation(int output);
/*}}}*/
/* ***** Point Init Routines ***** {{{*/
void readPoints();
void makeMCPoints(int n);
void makeGeneralizedSpiralPoints(int n);
/*}}}*/
/* ***** Integration Routines ***** {{{*/
void divideIntegrate(const int levs, const pnt &A, const pnt &B, const pnt &C, pnt &Top, double &bot);
void init_quadrature();
void quadrature(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadratureCR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadratureVR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadratureMP(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadrature7P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadrature13P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
void quadrature19P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot);
/*}}}*/
/* ***** Generic Region Routines *****{{{ */
void sortPoints(int sort_type, vector<region> &region_vec);
void sortBoundaryPoints(vector<region> &region_vec);
void triangulateRegions(vector<region> &region_vec);
void integrateRegions(vector<region> &region_vec);
void computeMetrics(double &ave, double &max, double &l1);
void clearRegions(vector<region> &region_vec);
void makeFinalTriangulations(vector<region> &region_vec);
void projectToBoundary(vector<region> &region_vec);
/*}}}*/
/* ***** Specific Region Routines ***** {{{ */
void printAllFinalTriangulation();
void printMyFinalTriangulation();
void storeMyFinalTriangulation();
/*}}}*/
/* ***** Communication Routines ***** {{{*/
void transferUpdatedPoints();
void gatherAllUpdatedPoints();
/*}}}*/
/* ***** Routines for Points *****{{{*/
void writePointsAsRestart(const int it, const restart_mode_type restart_mode, const fileio_mode_type fileio_mode);
#ifdef USE_NETCDF
int writeRestartFileOverwriteNC( const int it, const vector<pnt> &points );
int writeRestartFileRetainNC( const int it, const vector<pnt> &points );
#endif
int writeRestartFileOverwriteTXT( const int it );
int writeRestartFileRetainTXT( const int it );
double density(const pnt &p);
double ellipse_density(const pnt &p, double lat_c, double lon_c, double lat_width, double lon_width);
/*}}}*/

int main(int argc, char **argv){
	int bisection;
	int it, i;
	int stop, do_proj;
	int ave_points, my_points;
	mpi::request *ave_comms, *max_comms, *l1_comms;
	double *my_ave, *my_max, *my_l1;
	double glob_ave, glob_max, glob_l1;
	optional ave_opti, max_opti, l1_opti;
	pnt p;

	flags = new char[flags_str.size()+1];
	strcpy(flags,flags_str.c_str());

	for(i = 0; i < num_global_timers; i++){
		global_timers[i] = mpi_timer(global_names[i]);
	}
	global_timers[0].start(); // Global Time Timer

	mpi::environment env(argc, argv);	
	id = world.rank();
	num_procs = world.size();

	my_ave = new double[num_procs];
	my_max = new double[num_procs];
	my_l1 = new double[num_procs];

	ave_comms = new mpi::request[num_procs];
	max_comms = new mpi::request[num_procs];
	l1_comms = new mpi::request[num_procs];

	// Read in parameters and regions. Setup initial point set
	if(id == master){
		readParamsFile();

		cout << "Using " << quad_names[quad_rule] << " for integration." << endl;
		cout << "Writing restart files every " << restart << " steps." << endl;

		if (points_begin == 0){
			cout <<"Points being read in from SaveVertices." << endl;
			readPoints();
		} else if(points_begin == 1){
			cout << "Points being created with Monte Carlo." << endl;
			makeMCPoints(num_pts);
		} else if(points_begin == 2){
			cout << "Points being created with Generalized Spiral." << endl;
			makeGeneralizedSpiralPoints(num_pts);
		}
		readBoundaries();
		buildRegions();

		ofstream pts_out("point_initial.dat");
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			pts_out << (*point_itr) << endl;
		}
		pts_out.close();
	}

	// Each processor needs to setup the quadrature rules.
	init_quadrature();

	// Broadcast parameters, regions, and initial point set to each processor.
	mpi::broadcast(world,num_pts,master);
	mpi::broadcast(world,max_it,master);
	mpi::broadcast(world,restart,master);
	mpi::broadcast(world,max_it_no_proj,master);
	mpi::broadcast(world,max_it_scale_alpha,master);
	mpi::broadcast(world,num_bisections,master);
	mpi::broadcast(world,div_levs,master);
	mpi::broadcast(world,conv,master);
	mpi::broadcast(world,eps,master);
	mpi::broadcast(world,quad_rule,master);
	mpi::broadcast(world,regions,master);
	mpi::broadcast(world,points,master);
	mpi::broadcast(world,boundary_points,master);

	//printRegions();

	// Setup timers
	for(i = 0; i < num_timers; i++){
		my_timers[i] = mpi_timer(names[i]);
	}

	// Each processor clears my_regions (to make sure it's empty) and add it's own region into it's list.
	// If there is only 1 processor, that processor takes all regions, which is a serial computation.
	my_regions.clear();
	if(num_procs > 1){
		my_regions.push_back(regions.at(id));
		if(num_procs != regions.size()){
			cout << "Region error ---- Region size must equal number of processors." << endl;
			cout << " Or you must only use 1 processor" << endl;
			assert(num_procs == regions.size());
		}
	} else {
		for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
			my_regions.push_back((*region_itr));
		}
	}

	sortBoundaryPoints(my_regions);

	// Loop over bisections
	for(bisection = 0; bisection <= num_bisections; bisection++){
		// Clear loop timers
		for(i = 0; i < num_timers; i++){
			my_timers[i].init();
		}

		// Start loop timer
		my_timers[0].start();

		stop = 0;
		do_proj = 1;
		for(it = 0; it < max_it && !stop; it++){
			glob_ave = 0.0;
			glob_max = 0.0;
			glob_l1 = 0.0;
			for(i = 0; i < num_procs; i++){
				my_ave[i] = 0.0;
				my_max[i] = 0.0;
				my_l1[i] = 0.0;
			}
			my_timers[1].start(); // Iteration Timer

			clearRegions(my_regions);

			my_timers[2].start(); // Triangulation Timer

			my_timers[7].start(); // Sort Timer

			sortPoints(sort_method, my_regions);

			my_timers[7].stop(); // Sort Timer

			triangulateRegions(my_regions);

			my_timers[2].stop();

			my_timers[3].start(); // Integration Timer

			integrateRegions(my_regions);

			my_timers[3].stop();

			if(it > max_it_no_proj){
				proj_alpha = max((double)(it-max_it_no_proj), 0.0)/max((double)max_it_scale_alpha, 1.0);
				projectToBoundary(my_regions);
			}

			my_timers[4].start(); // Metrics Timer

			computeMetrics(my_ave[id],my_max[id], my_l1[id]);

			// Start non-blocking sends and receives of metrics
			if(id == master){
				for(i = 1; i < num_procs; i++){
					ave_comms[i] = world.irecv(i,msg_ave,my_ave[i]);
					max_comms[i] = world.irecv(i,msg_max,my_max[i]);
					l1_comms[i] = world.irecv(i,msg_l1,my_l1[i]);
				}
			} else {
				ave_comms[id] = world.isend(master,msg_ave,my_ave[id]);
				max_comms[id] = world.isend(master,msg_max,my_max[id]);
				l1_comms[id] = world.isend(master,msg_l1,my_l1[id]);
			}

			my_timers[4].stop();
			
			my_timers[5].start(); // Communication Timer
			transferUpdatedPoints();
			my_timers[5].stop();

			my_timers[6].start();
			// Finish metrics sends and receives, check for convergence and broadcast a stop request to all processors.
			if(id == master){
				glob_ave = my_ave[id];
				glob_max = my_max[id];
				glob_l1 = my_l1[id];

				for(i = 1; i < num_procs; i++){
					ave_opti = ave_comms[i].test();
					max_opti = max_comms[i].test();
					l1_opti = l1_comms[i].test();

					if(!ave_opti) ave_comms[i].wait();
					if(!max_opti) max_comms[i].wait();
					if(!l1_opti) l1_comms[i].wait();

					glob_ave += my_ave[i];
					glob_max = std::max(glob_max, my_max[i]);
					glob_l1 += my_l1[i];
				}
				glob_ave = sqrt(glob_ave)/points.size();
				glob_l1 = glob_l1/points.size();

				cout << it << " " << glob_ave << " " << glob_l1 << " " << glob_max << endl;

				if(conv == 1 && glob_ave < eps){
					cout << "Converged on average movement." << endl;
					stop = 1;
				} else if(conv == 2 && glob_max < eps){
					cout << "Converged on maximum movement." << endl;
					stop = 1;
				}
			}

			mpi::broadcast(world,stop,master);
			my_timers[6].stop();
			my_timers[1].stop();
			
			if(restart > 0 && it > 0){
				if(it%restart == 0){
					writePointsAsRestart(it, restart_mode, fileio_mode);
				}
			}
		}
		// Finish loop timer
		my_timers[0].stop();

		// Bisect if needed
		if(bisection < num_bisections){
			bisectTriangulation(0);
		} else {
			if(id == master){
				cout << "No more bisections for convergence" << endl;
			}
		}
		//Print out loop timers
		if(id == master){
			cout << endl << endl;
			for(i = 0; i < num_timers; i++){
				cout << my_timers[i];
			}
			cout << endl;
		}
	}

	//Compute average points per region for diagnostics
	ave_points = 0;
	my_points = 0;
	clearRegions(my_regions);
	sortPoints(sort_method, my_regions);
	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		my_points += (*region_itr).points.size();
	}

	mpi::reduce(world, my_points, ave_points, std::plus<int>(), master);

	ave_points = ave_points / regions.size();	
	if(id == master){
		cout << endl;
		cout << "Average points per region: " << ave_points << endl;
		cout << endl;
	}

	global_timers[0].stop();
	//Gather all updated points onto master processor, for printing to end_points.dat
	global_timers[1].start(); // Global Gather Timer
	gatherAllUpdatedPoints();
	global_timers[1].stop();

	// Compute final triangulation by merging all triangulations from each processor into an
	// unordered_set, and then ordering them ccw before printing them out.
	// write triangles to triangles.dat
	global_timers[2].start(); // Final Triangulation Timer
	clearRegions(my_regions);
	sortPoints(sort_dot, my_regions);
	triangulateRegions(my_regions);
	makeFinalTriangulations(my_regions);
	printMyFinalTriangulation();
	global_timers[2].stop();

	if(id == master){
		ofstream end_pts("end_points.dat");
		ofstream pt_dens("point_density.dat");
		ofstream bdry_pts("boundary_points.dat");
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			end_pts << (*point_itr) << endl;
			pt_dens << density((*point_itr)) << endl;
		}
		for(boundary_itr = boundary_points.begin(); boundary_itr != boundary_points.end(); boundary_itr++){
			bdry_pts << (*boundary_itr) << endl;
		}
		end_pts.close();
		pt_dens.close();
		bdry_pts.close();
	}

	//Bisect all edges of all triangles to give an extra point set at the end, bisected_points.dat
	global_timers[3].start(); // Final Bisection Timer
	bisectTriangulation(1);
	global_timers[3].stop();

	//Print out final timers, for global times.
	if(id == master){
		cout << endl << " ---- Final Timers ---- " << endl;
		for(i = 0; i < num_global_timers; i++){
			cout << global_timers[i];
		}

	}

	return 0;
}

/* ***** Setup Routines ***** {{{*/
void readParamsFile(){/*{{{*/
	//Read in parameters from Params.
	//If Params doesn't exist, write out Params with a default set of parameters
	string junk;
	ifstream params("Params");
	int temp_restart_mode;
	int temp_fileio_mode;

	if(!params){
		cout << "Error opening Params file." << endl;
		cout << "Writing a default Params file." << endl;
		cout << "Exiting, please set up Params, and rerun." << endl;
		ofstream pout("Params");
		pout << "How do you want the points created? (0 - Read from SaveVertices.dat, 1 - Monte Carlo, 2 - Generalized Spiral)" << endl;
		pout << "0" << endl;
		pout << "If you want them generated, how many points do you want?" << endl;
		pout << "162" << endl;
		pout << "How many iterations do you want to run for, if convergence isn't reached?" << endl;
		pout << "1000" << endl;
		pout << "How often, in iterations, do you want the point set written to a file? (Longer is better)" << endl << 500 << endl;
		pout << "How many iterations do you want to run without projection onto the boundary?" << endl;
		pout << "10000" << endl;
		pout << "How many iterations do you want with a variable projection distance? (Minimum of 1)" << endl;
		pout << "0" << endl;
		pout << "How many sub-triangle divisions would you like? (Minimum of 1, Causes every triangle to be divided into 4^n triangles)" << endl;
		pout << "1" << endl;
		pout << "How many bisections do you want until your final point set? (Each bisection maps n -> 4*n-6)" << endl;
		pout << "0" << endl;
		pout << "What do you want to check convergence on? (0 - Max Iterations, "
			 << "1 - Average Generator Movement, 2 - Maximum Generator Movement)" << endl;
		pout << "0" << endl;
		pout << "What do you want your convergence criteria to be? (eps = 1E-10)" << endl;
		pout << "1E-10" << endl;
		pout << "What Quadrature Rule do you want to use? (0 - Centroid, 1 - Vertex, 2 - Midpoint, 3 - 7 Point, 4 - 13 Point, 5 - 19 Point)" << endl;
		pout << "2" << endl;
		pout << "What sorting method do you want to use? (0 - dot product, 1 - voronoi)" << endl;
		pout << "0" << endl;
		pout << "What is the maximum allowable distance between boundary points? (Given in km)" << endl;
		pout << "4.0" << endl;
		pout << "Which format for restart files would you like? (0 - text, 1 - netcdf, 2 - both, 0 will be selected if netcdf is not linked)" << endl;
		pout << "0" << endl;
		pout << "Would you like one restart file, or a series? (0 - overwrite, 1 - retain, ignored if restart files are disabled above)" << endl;
		pout << "0" << endl;

		pout.close();

		exit(1);
	}
	
	getline(params,junk);
	params >> points_begin;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> num_pts;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> max_it;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> restart;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> max_it_no_proj;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> max_it_scale_alpha;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> div_levs;
	div_levs = std::max(1,div_levs);
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> num_bisections;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> conv;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> eps;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> quad_rule;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> sort_method;
	params.ignore(10000,'\n');
	getline(params,junk);
// 	params >> min_bdry_angle;
	params >> max_resolution;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> temp_fileio_mode;
	params.ignore(10000,'\n');
	getline(params,junk);
	params >> temp_restart_mode;
	params.ignore(10000,'\n');
	
	switch (temp_fileio_mode) {
		case 0:
			fileio_mode = FILEIO_TXT;
			break;
		case 1:
			fileio_mode = FILEIO_NETCDF;
			break;
		case 2:
			fileio_mode = FILEIO_BOTH;
			break;
		default:
			cout << "Restart file format has incorrect value in params file. Should be one of {0,1,2}.";
			exit(1);
	}

	switch (temp_restart_mode) {
		case 0:
			restart_mode = RESTART_OVERWRITE;
			break;
		case 1:
			restart_mode = RESTART_RETAIN;
			break;
		default:
			cout << "Restart mode (overwrite/retain) has incorrect value in params file. Should be one of {0,1}";
			exit(1);
	}

// 	min_bdry_angle = min_bdry_angle * M_PI/180.0;

	params.close();
}/*}}}*/
void readBoundaries(){/*{{{*/
	int i, j, n_pts;
	int count_count, bdry_count, fill_count, bdry_total;
	int add_count;
	int count_start, count_stop;
	int cur_loop_start, cur_loop_length;
	int p0_idx, p1_idx;
	//pnt p;
	//pnt p_b, p_e;
	double bdry_lon, bdry_lat;
	double dtr = M_PI/180.0;
	double point_delta;
	double add_spacing;
	double denom;
	double c0, c1;
	double t, omega;
	double r_earth = 6371.0; //avg radius in km
	
	// gw: read boundary points file
	bdry_count = 0;
	ifstream bdry_in("SaveBoundaries");
	if(!bdry_in)
		return;

	while(!bdry_in.eof()){
		bdry_in >> bdry_lon >> bdry_lat;
		bdry_in.ignore(10000,'\n');
	
		if(bdry_in.good()){
			
			bdry_lon *= dtr;
			bdry_lat *= dtr;
			pnt p = pntFromLatLon(bdry_lat, bdry_lon);
			
			p.normalize();
			p.idx = j;
			p.isBdry = 0;
			bdry_count++;
			boundary_points.push_back(p);
		} // end if input good
			
	} // end while not eof
	bdry_in.close();

	// gw: read loop counts file
	count_count = 0;
	ifstream count_in("SaveLoopCounts");
	if(!count_in)
		return;

	while(!count_in.eof()){
		count_in >> count_start >> count_stop;
		count_in.ignore(10000,'\n');
		
		if(count_in.good()){
			loop_start.push_back(count_start);
			loop_stop.push_back(count_stop);
			count_count++;
			
		} // end if input good
	} // end while not eof
	count_in.close();
	// gw: loop over loops
	fill_count = 0;
	for(int cur_loop = 0; cur_loop < count_count; cur_loop++){
		cur_loop_start = loop_start.at(cur_loop) - 1;
        cur_loop_length = loop_stop.at(cur_loop) - cur_loop_start;
	
		// gw: loop over point pairs in current loop
		for(int cur_pair = 0; cur_pair < cur_loop_length; cur_pair++){
		
            p0_idx = cur_loop_start + cur_pair;
            p1_idx = cur_loop_start + (cur_pair+1)%cur_loop_length;
            pnt p0 = boundary_points.at(p0_idx);
            pnt p1 = boundary_points.at(p1_idx);
            point_delta = p1.dotForAngle(p0);
            
			// gw: if distance between pair is greater than allowed amount
			//     then add some additional points
			if ( (point_delta * r_earth) > max_resolution ) {
			
				// gw: figure out how many points to added
				add_count = (int)ceil( (point_delta * r_earth) ) / max_resolution;
                add_spacing = 1.0 / ((double)add_count + 1);
				denom = sin( point_delta );
				bdry_lat = p1.getLat() - p0.getLat();
				bdry_lon = p1.getLon() - p0.getLon();
			
				// gw: loop for adding point(s)
				for(int cur_add = 0; cur_add <= add_count; cur_add++){
					
					/* Slerp - Doesn't work for constant lat, lon curves
					// http://en.wikipedia.org/wiki/Slerp
					t = add_spacing * ( (double)cur_add + 1 );
					c0 = sin( (1.0 - t) * point_delta ) / denom;
					c1 = sin( t * point_delta ) / denom;
					pnt temp_point = c0 * p0 + c1 * p1; // */

					// Linear interpolation, using Lat Lon.
					pnt temp_point = pntFromLatLon( p0.getLat() + cur_add * add_spacing * bdry_lat, 
							                        p0.getLon() + cur_add * add_spacing * bdry_lon);

					temp_point.normalize();
					temp_point.idx = bdry_count + fill_count;
					boundary_points.push_back(temp_point);
					
					fill_count++;
				}
				
			} // end if need to add fill points
			
		} // end loop over point pairs in current loop
		
	} // end loop over loops
	cout << "Read in " << bdry_count << " boundary points." << endl;
	cout << "Made " << fill_count << " fill points." << endl;
	cout << "There are " << boundary_points.size() << " boundary points total." << endl;

	num_bdry = boundary_points.size(); 

}/*}}}*/
void buildRegions(){/*{{{*/
	//Read in region centers, and connectivity (triangulation) from files
	//RegionList, and RegionTriangulation
	//From these, regions are setup to have a center and a radius.
	ifstream region_list("RegionList");
	ifstream region_connectivity("RegionTriangulation");
	unordered_set<int, int_hasher> neighbors1, neighbors2;
	unordered_set<int, int_hasher>::iterator neigh_itr;
	vector<int>::iterator cur_neigh_itr;
	region r;
	pnt p;
	double loc_radius, max_radius;
	double alpha, beta;
	int min_connectivity;
	int t1, t2, t3;
	int i;

	alpha = 2.0/3.0;
	beta = 4;

#ifdef _DEBUG
	cerr << "Building regions " << id << endl;
#endif


	if(!region_list){
		cout << "Failed to open file RegionList." << endl;
		exit(1);	
	}

	// Read in region centers
	i = 0;
	while(!region_list.eof()){
		region_list >> p;
		p.idx = i;
		i++;
		p.normalize();
		r.center = p;
		r.radius = 0.0;

		if(region_list.good()){
			regions.push_back(r);
		}
	}
	region_list.close();

	region_neighbors.resize(regions.size());

	if(!region_connectivity){
		cout << "Failed to open file RegionTriangulation" << endl;
		exit(1);
	}

	//Read in region triangulation, and insert into hash table to get unique connectivity.
	min_connectivity = regions.size();
	while(!region_connectivity.eof()){
		region_connectivity >> t1 >> t2 >> t3;

		if(t1 < min_connectivity)
			min_connectivity = t1;
		if(t2 < min_connectivity)
			min_connectivity = t2;
		if(t3 < min_connectivity)
			min_connectivity = t3;
	}
	region_connectivity.close();

	region_connectivity.open("RegionTriangulation");
	while(!region_connectivity.eof()){
		region_connectivity >> t1 >> t2 >> t3;

		t1 = t1 - min_connectivity;
		t2 = t2 - min_connectivity;
		t3 = t3 - min_connectivity;

		if(region_connectivity.good()){
			region_neighbors[t1].insert(t2);
			region_neighbors[t1].insert(t3);

			region_neighbors[t2].insert(t1);
			region_neighbors[t2].insert(t3);

			region_neighbors[t3].insert(t1);
			region_neighbors[t3].insert(t2);
		}
	}
	region_connectivity.close();

	//Compute region radii by dotting with each neighbor, and taking the max distance.
	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		max_radius = 0.0;
		loc_radius = 0.0;

		neighbors1.insert((*region_itr).center.idx);

		for(region_neigh_itr = region_neighbors[(*region_itr).center.idx].begin(); 
					region_neigh_itr != region_neighbors[(*region_itr).center.idx].end();
					++region_neigh_itr){

			loc_radius = (*region_itr).center.dotForAngle(regions[(*region_neigh_itr)].center);
			if(loc_radius > max_radius){
				max_radius = loc_radius;
			}
			(*region_itr).neighbors.push_back((*region_neigh_itr));
		}

		(*region_itr).radius = std::min(max_radius,M_PI);
		(*region_itr).input_radius = std::min(max_radius,M_PI);

		(*region_itr).points.clear();
		(*region_itr).triangles.clear();
		(*region_itr).boundary_points.clear();
	}

	// Build first and second levels of neighbors, for use in more complicated sort method.
	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		neighbors1.clear();
		neighbors2.clear();

		neighbors1.insert((*region_itr).center.idx);
		neighbors2.insert((*region_itr).center.idx);
		for(cur_neigh_itr = (*region_itr).neighbors.begin(); cur_neigh_itr != (*region_itr).neighbors.end(); ++cur_neigh_itr){
			neighbors1.insert((*cur_neigh_itr));
			neighbors2.insert((*cur_neigh_itr));

			for(neighbor_itr = regions.at((*cur_neigh_itr)).neighbors.begin();
				neighbor_itr != regions.at((*cur_neigh_itr)).neighbors.end(); ++neighbor_itr){
				
				neighbors2.insert((*neighbor_itr));
			}
		}

		for(neigh_itr = neighbors1.begin(); neigh_itr != neighbors1.end(); ++neigh_itr){
			(*region_itr).neighbors1.push_back((*neigh_itr));
		}

		for(neigh_itr = neighbors2.begin(); neigh_itr != neighbors2.end(); ++neigh_itr){
			(*region_itr).neighbors2.push_back((*neigh_itr));
		}
	}

	region_neighbors.clear();


}/*}}}*/
void printRegions(){/*{{{*/
	// This function is only for debugging purposes.
	// It's used to print out the information associated with each region
	// to verify the region information gets read correctly.
	ofstream r_out("RegionCenters.dat");
	ofstream rr_out("RegionRadii.dat");
	ofstream rp_out("RegionProperties.dat");
	double radius;

	clearRegions(regions);
	sortPoints(sort_method, regions);

	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		radius = (*region_itr).radius;

		radius = acos(-radius + 1.0);
		r_out << (*region_itr).center << endl;
		rr_out << radius << endl;
		rp_out << (*region_itr).points.size() << endl;
	}

	r_out.close();
	rr_out.close();

	rp_out.close();
	clearRegions(regions);
}/*}}}*/
/* }}} */
/* ***** Bisect Edges Routines ***** {{{ */
void bisectEdges(int end){/*{{{*/
	//void bisectEdges(int end)
	// This function actually performs the bisection of the point set.
	// All bisections are performed on the master process
	// INPUT:
	// 		int end - detemines if the routine prints or not.
	//
	// 	No output.
	// 	All new points are added directly into the points vector
	unordered_set<pair<int,int>, bpt::bisect_hasher> bisection_pts;
	unordered_set<pair<int,int>, bpt::bisect_hasher>::iterator bi_pts_itr;
	pair<int, int> in_pair;
	tri t;
	int vi1, vi2, vi3;
	int i;

	world.barrier();
	if(id == master){
		if(!end)
			cout << "Bisecting from " << points.size() << " points to ";

		i = points.size();
		for(tri_itr = all_triangles.begin(); tri_itr != all_triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;
			pnt new_p;

			if(vi1 < vi2){
				in_pair.first = vi1;
				in_pair.second = vi2;
				bisection_pts.insert(in_pair);
				i++;
			}

			if(vi2 < vi3){
				in_pair.first = vi2;
				in_pair.second = vi3;
				bisection_pts.insert(in_pair);
				i++;
			}

			if(vi3 < vi1){
				in_pair.first = vi3;
				in_pair.second = vi1;
				bisection_pts.insert(in_pair);
				i++;
			}
		}

		i = points.size();
		for(bi_pts_itr = bisection_pts.begin(); bi_pts_itr != bisection_pts.end(); ++bi_pts_itr){
			pnt p = 0.5*(points.at((*bi_pts_itr).first) + points.at((*bi_pts_itr).second));
			p.idx = i;
			p.isBdry = 0;
			p.normalize();
			points.push_back(p);
			i++;
		}
		if(!end)
			cout << points.size() << " points." << endl;
	} else {
		points.clear();
	}

	clearRegions(regions);
	world.barrier();
	mpi::broadcast(world, points, master);
	return;
}/*}}}*/
void bisectTriangulation(int output){/*{{{*/
	// void bisectTriangulation(int end)
	//
	// Sets up the required datastructures to perform the bisection with
	// bisectEdges
	//
	// Input
	if(!output)
		makeFinalTriangulations(my_regions);

	storeMyFinalTriangulation();
	bisectEdges(output);

	if(id == master && output){
		ofstream pts_out("bisected_points.dat");
		for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
			pts_out << (*point_itr) << endl;
		}
		pts_out.close();
	}

	return;
}/*}}}*/
/* }}} */
/* ***** Point Init Routines ***** {{{*/
void readPoints(){/*{{{*/
	//Read in initial point set from SaveVertices
	ifstream points_in("SaveVertices");
	pnt p;
	int i;

	i = 0;
	while(!points_in.eof()){
		points_in >> p;
		points_in.ignore(10000,'\n');
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		if(points_in.good()){
			points.push_back(p);
		}
		i++;
	}

	num_pts = points.size();

	cout << "Read in " << num_pts << " points from SaveVertices" << endl;

}/*}}}*/
void makeMCPoints(int n){/*{{{*/
	//Create Monte Carlo random point set
	int i, j; 
	srand48(time(NULL));
	double dlon, dlat;
	double dtr;
	double lat, lon;
	double x, y, z;
	double dens_check, dens_comp;
	pnt p;
	
	dtr = M_PI/180.0;

	// Uniform Spherical Points
	for(i = 0; i < n; i++){
		x = drand48()*2.0-1.0;
		y = drand48()*2.0-1.0;
		z = drand48()*2.0-1.0;

		p = pnt(x,y,z);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		points.push_back(p);
	}

	// */
	/* Uniform points on equator between -20 and +20 latitude
	lat = 20.0 * dtr;
	dlon = 2.0*M_PI/(n/10);
	j = 0;
	for(i = 0; i < (n / 10); i++){
		lon = i*dlon;
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
	}
	lat = -1.0*lat;
	for(i = n/10; i < (2 * n / 10); i++){
		lon = (i-(n/10))*dlon;
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
	}
	for(i = (2*n/10); i < n; i++){
		lon = drand48()*2*M_PI;
//		do{
			lat = (drand48()*2.0 - 1.0)*(15.0*dtr);
//		}while(!(fabs(lat) < 15.0*dtr));
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 0;
		points.push_back(p);
		j++;
	}
	// */
	/* Uniform points in rectangle, from 0 to 45 north lat, and 0 to 90 lon
	lat = 0.0;
	dlon = 90.0*dtr/(n/10);
	lon = 0.0;
	j = 0;
	for(i = 0; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lon += dlon;
	}
	lat = 45.0*dtr;
//	dlon = 90.0*dtr/(2*n/30);
	lon = 0.0;
	for(i = 0; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lon += dlon;
	}
	lon = 0.0;
	dlat = 45.0*dtr/(n/20);
	lat = dlat;
	for(i = 1; i < (n/20); i++){
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lat += dlat;
	}
	lon = M_PI/2.0;
	lat = dlat;
	for(i = 1; i < (n/20); i++){
		p = pntFromLatLon(lat,lon);	
		p.normalize();
		p.idx = j;
		p.isBdry = 1;
		points.push_back(p);
		j++;
		lat += dlat;
	}
	for(i = j; i < n; i++){
		lon = (drand48()*0.85 + 0.05)*M_PI/2.0;
		lat = (drand48()*0.85 + 0.05)*M_PI/4.0;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}
	// */
	/* Uniform points in a triangle from 1,0,0, 0,1,0, and 0,0,1
	lat = 0.0;	
	dlon = (M_PI/2.0)/(n/10);
	j = 0;
	for(i = 0; i < (n/10)+1; i++){
		lon = i * dlon;
		p = pntFromLatLon(lat,lon);
		p.idx = j;
		p.isBdry = 2;
		p.normalize();
		points.push_back(p);
		j++;
	}
	lon = 0.0;
	dlat = (M_PI/2.0)/(n/10);
	lat = dlat;
	for(i = 1; i < (n/10)+1; i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 2;
		points.push_back(p);
		lat += dlat;
		j++;
	}
	lon = M_PI/2.0;
	lat = dlat;
	for(i = 1; i < (n/10); i++){
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = j;
		p.isBdry = 2;
		points.push_back(p);
		lat += dlat;
		j++;
	}
	for(i = j; i < n; i++){
		lat = (drand48()*0.85 + 0.05)*M_PI/2.0;
		lon = (drand48()*0.85 + 0.05)*M_PI/2.0;
		p = pntFromLatLon(lat,lon);
		p.normalize();
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}
	// */
	cout << "Created " << points.size() << " points using monte carlo." << endl;
}/*}}}*/
void makeGeneralizedSpiralPoints(int n){/*{{{*/
	//Create Generalize Spiral point set
	int i;
	int idx;
	pnt p;
	double phi_curr, h, theta, aa, bb, cc;
	double gsC = 3.809;
	double twopi_dp;

	twopi_dp = 2.0*M_PI;

	p = pnt(0.0,0.0,0.0,0,0);

	//	first pt, loop primer
	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	points.push_back(p);

	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}

	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	points.push_back(p);

	cout << "Created " << points.size() << " points using generalized spiral." << endl;
}/*}}}*/
/*}}}*/
/* ***** Integration Routines *****  {{{*/
void divideIntegrate(const int levs, const pnt &A, const pnt &B, const pnt &C, pnt &Top, double &bot){/*{{{*/
	//divideIntegrate integrates a triangle, with an arbitrary number of subdivisions.
	//Each subdivision produces 4 triangles from 1 triangle.
	//Master triangle vertices are A, B, C.
	//Number of divisions is levs
	//Top and bot are returned, being the top and bottom of the centroid calculation.
	pnt z1, z2, z3;
	pnt temp;
	double xlam, ylam, ylam2, ylam3, zlam;
	double dt, bot_temp;
	double area;
	int zn, xn;
	int i, j;

	Top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	zn = (1 << levs);
	dt = 1.0/zn;

	for(i = 0; i < zn; i++){
		xn = zn-i;
		zlam = i*dt;
		for(j = 0; j < xn-1; j++){
			xlam = 1.0 - i*dt - j*dt;
			ylam = 1.0 - (xlam + zlam);
			z1 =     xlam    * A +    ylam     * B +    zlam     * C;
			z2 = (xlam - dt) * A + (ylam + dt) * B +    zlam     * C;
			z3 = (xlam - dt) * A +    ylam     * B + (zlam + dt) * C;
			z1.normalize();
			z2.normalize();
			z3.normalize();

			quadrature(z1,z2,z3,temp,bot_temp);
			area = triArea(z1,z2,z3);
			Top += temp*area;
			bot += bot_temp*area;

			z1 =   (xlam - dt)   * A +     ylam    * B + (zlam + dt) * C;
			z2 =   (xlam - dt)   * A + (ylam + dt) * B +    zlam     * C;
			z3 = (xlam - 2.0*dt) * A + (ylam + dt) * B + (zlam + dt) * C;
			z1.normalize();
			z2.normalize();
			z3.normalize();

			quadrature(z1,z2,z3,temp,bot_temp);
			area = triArea(z1,z2,z3);
			Top += temp * area;
			bot += bot_temp*area;
		}

		xlam = dt;
		ylam = 1.0 - (xlam + zlam);

		z1 =    xlam     * A +    ylam     * B +    zlam     * C;
		z2 = (xlam - dt) * A + (ylam + dt) * B +    zlam     * C;
		z3 = (xlam - dt) * A +    ylam     * B + (zlam + dt) * C;
		z1.normalize();
		z2.normalize();
		z3.normalize();

		quadrature(z1,z2,z3,temp,bot_temp);
		area = triArea(z1,z2,z3);
		Top += temp*area;
		bot += bot_temp*area;
	}
}/*}}}*/
void init_quadrature(){/*{{{*/
	/****************************************************************************************
	 *   -This function initializes the quadrature weights (wq) and combination factors (q)
	 *		for a 7 point quadrature rule on triangles.
	 *
	 * 	 - No input or output
	 *****************************************************************************************/
	q7[0] = (6.0 - sqrt(15.0))/21.0;
	q7[1] = (9.0 + 2.0 * sqrt(15.0))/21.0;
	q7[2] = (6.0 + sqrt(15.0))/21.0;
	q7[3] = (9.0 - 2.0 * sqrt(15.0))/21.0;

	wq7[0] = (155.0 - sqrt(15.0))/1200.0;
	wq7[1] = (155.0 - sqrt(15.0))/1200.0;
	wq7[2] = (155.0 - sqrt(15.0))/1200.0;
	wq7[3] = (155.0 + sqrt(15.0))/1200.0;
	wq7[4] = (155.0 + sqrt(15.0))/1200.0;
	wq7[5] = (155.0 + sqrt(15.0))/1200.0;
	wq7[6] = 9.0/40.0;

	q13[0] = 1.0/3.0;
	q13[1] = 0.479308067841923;
	q13[2] = 0.260345966079038;
	q13[3] = 0.869739794195568;
	q13[4] = 0.065130102902216;
	q13[5] = 0.638444188569809;
	q13[6] = 0.312865496004875;
	q13[7] = 0.048690315425316;

	wq13[0] = -0.149570044467670;
	wq13[1] = 0.175615257433204;
	wq13[2] = 0.175615257433204;
	wq13[3] = 0.175615257433204;
	wq13[4] = 0.053347235608839;
	wq13[5] = 0.053347235608839;
	wq13[6] = 0.053347235608839;
	wq13[7] = 0.077113760890257;
	wq13[8] = 0.077113760890257;
	wq13[9] = 0.077113760890257;
	wq13[10] = 0.077113760890257;
	wq13[11] = 0.077113760890257;
	wq13[12] = 0.077113760890257;

	q19[0] = 1.0 / 3.0;
	q19[1] = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
	q19[2] = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
	q19[3] = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
	q19[4] = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
	q19[5] = ( 40.0 - 10.0 * sqrt ( 15.0 )+ 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
	q19[6] = ( 25.0 +  5.0 * sqrt ( 15.0 )-  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
	q19[7] = ( 40.0 + 10.0 * sqrt ( 15.0 )+ 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
	q19[8] = ( 25.0 -  5.0 * sqrt ( 15.0 )-  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
	q19[9] = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
	q19[10] = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 )- sqrt ( 105.0 ) ) / 90.0;
	q19[11] = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 )+ sqrt ( 105.0 ) ) / 90.0;

	wq19[0] = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
	wq19[1] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[2] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[3] = (-9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		+ 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[4] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[5] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[6] = (-9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 )
		/ 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0
		- 198763.0 * sqrt ( 105.0 ) / 939008.0)/3.0;
	wq19[7] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
			/3.0;
	wq19[8] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
			/3.0;
	wq19[9] = (( 102791225.0 - 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	+ 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0)
			/3.0;
	wq19[10] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
	wq19[11] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
	wq19[12] = (( 102791225.0 + 23876225.0 * sqrt ( 15.0 )
			- 34500875.0 * sqrt ( 7.0 )	- 9914825 * sqrt ( 105.0 ) ) / 59157504.0)/3.0;
	wq19[13] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
	wq19[14] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
	wq19[15] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
	wq19[16] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
	wq19[17] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
	wq19[18] = (( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0)/ 6.0;
}/*}}}*/
void quadrature(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	//Call the appropriate quadrature rule that has been requested
	switch(quad_rule){
		case 0:
			quadratureCR(A,B,C,top,bot);
			break;;
		case 1:
			quadratureVR(A,B,C,top,bot);
			break;;
		case 2:
			quadratureMP(A,B,C,top,bot);
			break;;
		case 3:
			quadrature7P(A,B,C,top,bot);
			break;;
		case 4:
			quadrature13P(A,B,C,top,bot);
			break;;
		case 5:
			quadrature19P(A,B,C,top,bot);
			break;;
		default:
			quadratureMP(A,B,C,top,bot);
			break;;
	}
}/*}}}*/
void quadratureCR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	//Integrate a triangle using the Circumcenter rule
	pnt cent;
	double d;

	circumcenter(A,B,C,cent);
	cent.normalize();

	d = density(cent);
	bot = d;

	top = d*cent;
}/*}}}*/
void quadratureVR(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	//Integrate a triangle using the vertex rule
	double d1, d2, d3;

	d1 = density(A)/3.0;
	d2 = density(B)/3.0;
	d3 = density(C)/3.0;

	bot = d1+d2+d3;

	top = d1*A + d2*B + d3*C;
}/*}}}*/
void quadratureMP(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle 
	 *      v defined by the points p1, p2, and p3.
	 *
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: Cent and dens
	 *		  Cent: Cent contains the three components of the first integral
	 *        dens: Dens contains the value of the bottom integral.
	 *****************************************************************************/
	pnt p1, p2, p3;
	double da, db, dc;
	double d1, d2, d3;
	int i, j; // Loop Index

	top = pnt(0.0,0.0,0.0,0,0);
	bot = 0.0;
	
	p1 = (A+B)/2.0;
	p2 = (B+C)/2.0;
	p3 = (C+A)/2.0;

	p1.normalize();
	p2.normalize();
	p3.normalize();

	da = density(A);
	db = density(B);
	dc = density(C);

	d1 = (da+db)/2.0;
	d2 = (db+dc)/2.0;
	d3 = (dc+da)/2.0;

	bot = (d1+d2+d3)/3.0;;

	top = (p1*d1 + p2*d2 + p3*d3)/3.0;
}/*}}}*/
void quadrature7P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle 
	 *      v defined by the points p1, p2, and p3.
	 *
	 *    - Input: 3 points on a sphere, p1, p2, p3.
	 *
	 *    - Output: Cent and dens
	 *		  Cent: Cent contains the three components of the first integral
	 *        dens: Dens contains the value of the bottom integral.
	 *****************************************************************************/
	pnt c_vals; // Temp aray for values
	double d;
	pnt z[7];
	double norm;
	int i, j; // Loop Index

	wq = wq7;
	q = q7;
	top = pnt(0.0,0.0,0.0,0,0);	
	bot = 0.0;

	z[0] = A*q[0] + B*q[0] + C*q[1];
	z[1] = A*q[0] + B*q[1] + C*q[0];
	z[2] = A*q[1] + B*q[0] + C*q[0];
	z[3] = A*q[2] + B*q[2] + C*q[3];
	z[4] = A*q[2] + B*q[3] + C*q[2];
	z[5] = A*q[3] + B*q[2] + C*q[2];
	z[6] = (A + B + C)/3.0;

	for(i = 0; i < 7; i ++){
		z[i].normalize();

		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}
}/*}}}*/
void quadrature13P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle 
	 *      v defined by the points A, B, and C.
	 *
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: top and bot
	 *		  top: top contains the three components of the first integral
	 *        bot: bot contains the value of the bottom integral.
	 *****************************************************************************/
	double d;
	pnt z[13];
	int i;

	wq = wq13;
	q = q13;
	
	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[1]*C;
	z[2] = q[2]*A+q[1]*B+q[1]*C;
	z[3] = q[1]*A+q[1]*B+q[2]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[7]*C;
	z[8] = q[5]*A+q[7]*B+q[6]*C;
	z[9] = q[6]*A+q[5]*B+q[7]*C;
	z[10] = q[6]*A+q[7]*B+q[5]*C;
	z[11] = q[7]*A+q[5]*B+q[6]*C;
	z[12] = q[7]*A+q[6]*B+q[5]*C;

	for(i = 0; i < 13; i ++){
		z[i].normalize();

		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}
}/*}}}*/
void quadrature19P(const pnt &A, const pnt &B, const pnt &C, pnt &top, double &bot){/*{{{*/
	/****************************************************************************
	 *   - This function calculates the integral int_v(rho(x)*x*dx) and int_v(rho(x)dx) over the triangle 
	 *      v defined by the points A, B, and C.
	 *
	 *    - Input: 3 points on a sphere, A, B, C.
	 *
	 *    - Output: top and bot
	 *		  top: top contains the three components of the first integral
	 *        bot: bot contains the value of the bottom integral.
	 *****************************************************************************/
	double d;
	pnt z[19];
	int i;

	wq = wq19;
	q = q19;

	top = pnt(0.0,0.0,0.0);
	bot = 0.0;

	z[0] = q[0]*A+q[0]*B+q[0]*C;
	z[1] = q[1]*A+q[2]*B+q[2]*C;
	z[2] = q[2]*A+q[1]*B+q[2]*C;
	z[3] = q[2]*A+q[2]*B+q[1]*C;
	z[4] = q[3]*A+q[4]*B+q[4]*C;
	z[5] = q[4]*A+q[3]*B+q[4]*C;
	z[6] = q[4]*A+q[4]*B+q[3]*C;
	z[7] = q[5]*A+q[6]*B+q[6]*C;
	z[8] = q[6]*A+q[5]*B+q[6]*C;
	z[9] = q[6]*A+q[6]*B+q[5]*C;
	z[10] = q[7]*A+q[8]*B+q[8]*C;
	z[11] = q[8]*A+q[7]*B+q[8]*C;
	z[12] = q[8]*A+q[8]*B+q[7]*C;
	z[13] = q[9]*A+q[10]*B+q[11]*C;
	z[14] = q[9]*A+q[11]*B+q[10]*C;
	z[15] = q[10]*A+q[9]*B+q[11]*C;
	z[16] = q[10]*A+q[11]*B+q[9]*C;
	z[17] = q[11]*A+q[9]*B+q[10]*C;
	z[18] = q[11]*A+q[10]*B+q[9]*C;

	for(i = 0; i < 13; i ++){
		z[i].normalize();

		d = density(z[i])*wq[i];

		top += z[i]*d;
		bot += d;
	}
}/*}}}*/
/* }}} */
/* ***** Generic Region Routines ***** {{{ */
void sortPoints(int sort_type, vector<region> &region_vec){/*{{{*/
	//Sort points into my region(s).
	//This is done using a dot product and checking if the dot product is inside of the region radius
	double val;
#ifdef _DEBUG
	cerr << "Sorting Points " << id << endl;
#endif
	if(sort_type == sort_dot){
		//Simple Dot Product Sort using the radius based decomposition.
		for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
			(*region_itr).radius = (*region_itr).input_radius;
			for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
				val = (*point_itr).dotForAngle((*region_itr).center);	

				if(val < (*region_itr).radius){
					(*region_itr).points.push_back((*point_itr));
				} 
			}
		}
	} else if (sort_type == sort_vor){
		//More complicated sort, that sorts by Voronoi cells keeping current regions points, as well as neighboring regions points.
		//Should handle variable resolution meshes better than the more simple dot product sorting.
		double my_val;
		int added;
		double min_val;
		double max_dist;
		int min_region;
		vector<int>::iterator cur_neigh_itr;

		for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
			max_dist = 0.0;
			for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
				min_val = M_PI;
				my_val = (*point_itr).dotForAngle((*region_itr).center);

				if(my_val < (*region_itr).input_radius){
					for(neighbor_itr = (*region_itr).neighbors2.begin(); 
							neighbor_itr != (*region_itr).neighbors2.end(); ++neighbor_itr){

						val = (*point_itr).dotForAngle(regions.at((*neighbor_itr)).center);

						if(val < min_val){
							min_region = (*neighbor_itr);
							min_val = val;
						}
					}

					added = 0;
					for(neighbor_itr = (*region_itr).neighbors1.begin(); 
							neighbor_itr != (*region_itr).neighbors1.end() && added == 0; 
							++neighbor_itr){
						if(min_region == (*neighbor_itr)){
							val = (*region_itr).center.dotForAngle(regions[min_region].center);

							if(min_region == (*region_itr).center.idx){
								(*region_itr).points.push_back((*point_itr));
							} else if(my_val < val) {
								(*region_itr).points.push_back((*point_itr));
							}

							added = 1;
						}
					}

					if(my_val > max_dist){
						max_dist = my_val;
					}
				}

				(*region_itr).radius = max_dist;
			}
		}
	}
#ifdef _DEBUG
	cerr << "Done Sorting Points (Local) " << id << endl;
#endif
	return;
}/*}}}*/
void sortBoundaryPoints(vector<region> &region_vec){/*{{{*/
	//Sort points into my region(s).
	//This is done using a dot product and checking if the dot product is inside of the region radius
	double val, min_val;
	int index;
	//vector<region>::iterator region_itr;
#ifdef _DEBUG
	cerr << "Sorting Boundary Points " << id << endl;
#endif
	//simple dot product sort using the radius based decomposition.
	for(point_itr = boundary_points.begin(); point_itr != boundary_points.end(); point_itr++){
		min_val = M_PI;
		for(region_itr = regions.begin(); region_itr != regions.end(); region_itr++){
			val = (*point_itr).dotForAngle((*region_itr).center);

			if(val < min_val){
				min_val = val;
				index = (*region_itr).center.idx;
			}
		}

		for(region_itr = region_vec.begin(); region_itr != region_vec.end(); region_itr++){
			if((*region_itr).center.idx == index){
				(*region_itr).boundary_points.push_back((*point_itr));
			}
		}
	}
#ifdef _DEBUG
	cerr << "Done Sorting Boundary Points (Local) " << id << endl;
#endif
	return;
}/*}}}*/
void triangulateRegions(vector<region> &region_vec){/*{{{*/
	//Triangulate my region(s) points
	//Points are first stereographically projected into a plane tanget to region center
	//Projected point set is then triangulated using Triangle
	//Indexing is mapped from local (projected) point set to global (spherical) point set
	//Triangles whose circumcircles extend outside of the region radius are considered "bad" and are not kept
	struct triangulateio in, out, vorout;
	pnt ccenter;
	pnt a, b, c;
	pnt ac, bc;
	pnt x_hat, y_hat, axis;
	pnt Q;
	double cradius;
	double scale;
	double x, y, z;
	double criteria;
	double s;
	double min_dir;
	tri t;
	int i;
	int vi1, vi2, vi3;
	int swp_vi;

#ifdef _DEBUG
	cerr << "Triangulating points (local) " << id << endl;
#endif

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		in.numberofpoints = (*region_itr).points.size();
		in.numberofpointattributes = 0;
		in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = (double *) NULL;
		in.pointmarkerlist = (int *) NULL;

		axis = (*region_itr).center;
		min_dir = min(fabs(axis.x),min(fabs(axis.y),fabs(axis.z)));
		if(min_dir == fabs(axis.x)){
			axis.x = 1.0;
			axis.y = 0.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.y)){
			axis.x = 0.0;
			axis.y = 1.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.z)){
			axis.x = 0.0;
			axis.y = 0.0;
			axis.z = 1.0;
		}

		x_hat = (*region_itr).center.cross(axis);
		x_hat.normalize();
		y_hat = (*region_itr).center.cross(x_hat);
		y_hat.normalize();

		i = 0;

		for(point_itr = (*region_itr).points.begin(); point_itr != (*region_itr).points.end(); ++point_itr){
			s = 2.0/(*region_itr).center.dot((*point_itr) + (*region_itr).center);
			Q = s*(*point_itr) + (s - 1.0) * (*region_itr).center;
			in.pointlist[2*i] = x_hat.dot(Q);
			in.pointlist[2*i+1] = y_hat.dot(Q);

			i++;
		}

		out.pointlist = (double *)NULL;
		out.trianglelist = (int *)NULL;

		triangulate(flags,&in,&out,&vorout);

		for(i = 0; i < out.numberoftriangles; i++){
			vi1 = out.trianglelist[3*i];
			vi2 = out.trianglelist[3*i+1];
			vi3 = out.trianglelist[3*i+2];

			a = (*region_itr).points.at(vi1);
			b = (*region_itr).points.at(vi2);
			c = (*region_itr).points.at(vi3);

			if(!isCcw(a,b,c)){
				b = (*region_itr).points.at(vi3);
				c = (*region_itr).points.at(vi2);
			}

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			circumcenter(a,b,c,ccenter);
			ccenter.normalize();
			cradius = circumradius(a, b, c);

			criteria = (*region_itr).center.dotForAngle(ccenter);
			criteria += cradius;

			if(criteria < (*region_itr).radius){
				t = tri(vi1, vi2, vi3);
				(*region_itr).triangles.push_back(t);
			}
		}
		free(in.pointlist);
		free(in.regionlist);
		free(in.pointmarkerlist);
		free(out.pointlist);
		free(out.trianglelist);
	}
#ifdef _DEBUG
	cerr << "Done triangulating points (local) " << id << endl;
#endif
	return;
}/*}}}*/
void integrateRegions(vector<region> &region_vec){/*{{{*/
	// Integrate Voronoi cells inside of my region
	// Every region updates all points that are closer to their region center than any other region center.
	// This ensures that each point is only updated once.
	pnt *tops;
	double *bots;
	int vi1, vi2, vi3;
	pnt a, b, c;
	pnt ab, bc, ca;
	pnt ccenter;
	pnt np;
	double dist_temp;
	double a_dist_to_region, a_min_dist;
	double b_dist_to_region, b_min_dist;
	double c_dist_to_region, c_min_dist;
	int a_min_region, b_min_region, c_min_region;
	pnt top_val;
	double bot_val;
	int i;

	#ifdef _DEBUG
		cerr << "Integrating regions " << id << endl;
	#endif

	tops = new pnt[points.size()];
	bots = new double[points.size()];

	for(i = 0; i < points.size(); i++){
		tops[i] = pnt(0.0,0.0,0.0,0,i);
		bots[i] = 0.0;
	}

	i = 0;

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			vi1 = (*tri_itr).vi1;
			vi2 = (*tri_itr).vi2;
			vi3 = (*tri_itr).vi3;

			a = points[vi1];
			b = points[vi2];
			c = points[vi3];

			ab = (a+b)/2.0;
			bc = (b+c)/2.0;
			ca = (c+a)/2.0;

			ab.normalize();
			bc.normalize();
			ca.normalize();

			circumcenter(a,b,c,ccenter);
			ccenter.normalize();

			a_dist_to_region = a.dotForAngle((*region_itr).center);
			b_dist_to_region = b.dotForAngle((*region_itr).center);
			c_dist_to_region = c.dotForAngle((*region_itr).center);

			a_min_dist = 10.0;
			b_min_dist = 10.0;
			c_min_dist = 10.0;

			for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
				if(regions[(*neighbor_itr)].center.idx != (*region_itr).center.idx){
					dist_temp = a.dotForAngle(regions[(*neighbor_itr)].center);
					if(a_min_dist > dist_temp){
						a_min_dist = dist_temp;
						a_min_region = (*neighbor_itr);
					}
					dist_temp = b.dotForAngle(regions[(*neighbor_itr)].center);
					if(b_min_dist > dist_temp){
						b_min_dist = dist_temp;
						b_min_region = (*neighbor_itr);
					}
					dist_temp = c.dotForAngle(regions[(*neighbor_itr)].center);
					if(c_min_dist > dist_temp){
						c_min_dist = dist_temp;
						c_min_region = (*neighbor_itr);
					}
				}
			}

			if(a_dist_to_region < a_min_dist){
				//Triangle 1 - a ab ccenter
				divideIntegrate(div_levs,a,ab,ccenter,top_val,bot_val);
				tops[a.idx] += top_val;
				bots[a.idx] += bot_val;

				//Triangle 2 - a ccenter ca
				divideIntegrate(div_levs,a,ccenter,ca,top_val,bot_val);
				tops[a.idx] += top_val;
				bots[a.idx] += bot_val;
			} else if (a_dist_to_region == a_min_dist && (*region_itr).center.idx < a_min_region){
				//Triangle 1 - a ab ccenter
				divideIntegrate(div_levs,a,ab,ccenter,top_val,bot_val);
				tops[a.idx] += top_val;
				bots[a.idx] += bot_val;

				//Triangle 2 - a ccenter ca
				divideIntegrate(div_levs,a,ccenter,ca,top_val,bot_val);
				tops[a.idx] += top_val;
				bots[a.idx] += bot_val;
			}

			if(b_dist_to_region < b_min_dist){
				//Triangle 1 - b bc ccenter
				divideIntegrate(div_levs,b,bc,ccenter,top_val,bot_val);
				tops[b.idx] += top_val;
				bots[b.idx] += bot_val;

				//Triangle 2 - b ccenter ab
				divideIntegrate(div_levs,b,ccenter,ab,top_val,bot_val);
				tops[b.idx] += top_val;
				bots[b.idx] += bot_val;
			} else if (b_dist_to_region == b_min_dist && (*region_itr).center.idx < b_min_region){
				//Triangle 1 - b bc ccenter
				divideIntegrate(div_levs,b,bc,ccenter,top_val,bot_val);
				tops[b.idx] += top_val;
				bots[b.idx] += bot_val;

				//Triangle 2 - b ccenter ab
				divideIntegrate(div_levs,b,ccenter,ab,top_val,bot_val);
				tops[b.idx] += top_val;
				bots[b.idx] += bot_val;
			}

			if(c_dist_to_region < c_min_dist){
				//Triangle 1 - c ca ccenter
				divideIntegrate(div_levs,c,ca,ccenter,top_val,bot_val);
				tops[c.idx] += top_val;
				bots[c.idx] += bot_val;

				//Triangle 2 - c ccenter bc
				divideIntegrate(div_levs,c,ccenter,bc,top_val,bot_val);
				tops[c.idx] += top_val;
				bots[c.idx] += bot_val;
			} else if (c_dist_to_region == c_min_dist && (*region_itr).center.idx < c_min_region){
				//Triangle 1 - c ca ccenter
				divideIntegrate(div_levs,c,ca,ccenter,top_val,bot_val);
				tops[c.idx] += top_val;
				bots[c.idx] += bot_val;

				//Triangle 2 - c ccenter bc
				divideIntegrate(div_levs,c,ccenter,bc,top_val,bot_val);
				tops[c.idx] += top_val;
				bots[c.idx] += bot_val;
			}
		}
	}

	n_points.clear();
	for(point_itr = points.begin(); point_itr != points.end(); ++point_itr){
		if(bots[(*point_itr).idx] != 0.0){
			if(!(*point_itr).isBdry){
				np = tops[(*point_itr).idx]/bots[(*point_itr).idx];
				np.idx = (*point_itr).idx;
				np.isBdry = (*point_itr).isBdry;
				np.normalize();
				n_points.push_back(np);
			} else if((*point_itr).isBdry){
				if((*point_itr).isBdry == 2){
					(*point_itr).isBdry = 0;
				}

				n_points.push_back((*point_itr));
			}
		}
	}

	delete(tops);
	delete(bots);

	#ifdef _DEBUG
		cerr << "Done Integrating regions " << id << endl;
	#endif
	return;
}/*}}}*/
void computeMetrics(double &ave, double &max, double &l1){/*{{{*/
	//Metrics are computed for my updated points
	pnt norm_pt;
	double val;

	max = 0.0;
	ave = 0.0;
	l1 = 0.0;
#ifdef _DEBUG
	cerr << "Computing local metrics " << id << endl;
#endif

	for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
		norm_pt = points.at((*point_itr).idx) - (*point_itr);
		val = norm_pt.magnitude();

		ave += val*val;
		max = std::max(val,max);
		l1 += val;
	}
#ifdef _DEBUG
	cerr << "Done Computing local metrics " << id << endl;
#endif
	return;
}/*}}}*/
void clearRegions(vector<region> &region_vec){/*{{{*/
	//Clear all of my regions
#ifdef _DEBUG
	cerr << "Clearing local regions " << id << endl;
#endif
	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		(*region_itr).points.clear();
		(*region_itr).triangles.clear();

		assert((*region_itr).points.empty());
		assert((*region_itr).triangles.empty());
	}
#ifdef _DEBUG
	cerr << "Done Clearing local regions " << id << endl;
#endif
}/*}}}*/
void makeFinalTriangulations(vector<region> &region_vec){/*{{{*/
	//Make final triangulations, triangulates regions as in the function triangulateRegions, however
	//here the triangle vertices are ordered from smallest index to largest index
	int i;
	int vi1, vi2, vi3;
	int swp_vi;
	struct triangulateio in, out, vorout;
	tri t;
	double cradius;
	pnt ccenter;
	pnt a, b, c;
	pnt ac, bc;
	pnt x_hat, y_hat, axis;
	pnt Q;
	double s, min_dir;
	double a_dist, b_dist, c_dist;
	double a_dist_min, b_dist_min, c_dist_min;
	double dist_temp;
	double x, y, z;
	double criteria;

#ifdef _DEBUG
	cerr << "Triangulating points (local) " << id << endl;
#endif

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); ++region_itr){
		in.numberofpoints = (*region_itr).points.size();
		in.numberofpointattributes = 0;
		in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
		in.numberofsegments = 0;
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.regionlist = (double *) NULL;
		in.pointmarkerlist = (int *) NULL;

		axis = (*region_itr).center;
		min_dir = min(fabs(axis.x),min(fabs(axis.y),fabs(axis.z)));
		if(min_dir == fabs(axis.x)){
			axis.x = 1.0;
			axis.y = 0.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.y)){
			axis.x = 0.0;
			axis.y = 1.0;
			axis.z = 0.0;
		} else if (min_dir == fabs(axis.z)){
			axis.x = 0.0;
			axis.y = 0.0;
			axis.z = 1.0;
		}

		x_hat = (*region_itr).center.cross(axis);
		x_hat.normalize();
		y_hat = (*region_itr).center.cross(x_hat);
		y_hat.normalize();

		i = 0;
		for(point_itr = (*region_itr).points.begin(); point_itr != (*region_itr).points.end(); ++point_itr){
			s = 2.0/(*region_itr).center.dot((*point_itr) + (*region_itr).center);
			Q = s*(*point_itr) + (s - 1.0) * (*region_itr).center;
			in.pointlist[2*i] = x_hat.dot(Q);
			in.pointlist[2*i+1] = y_hat.dot(Q);
			i++;
		}

		out.pointlist = (double *)NULL;
		out.trianglelist = (int *)NULL;

		triangulate(flags,&in,&out,&vorout);

		for(i = 0; i < out.numberoftriangles; i++){
			vi1 = out.trianglelist[3*i];
			vi2 = out.trianglelist[3*i+1];
			vi3 = out.trianglelist[3*i+2];

			a = (*region_itr).points.at(vi1);
			b = (*region_itr).points.at(vi2);
			c = (*region_itr).points.at(vi3);

			if(!isCcw(a,b,c)){
				b = (*region_itr).points.at(vi3);
				c = (*region_itr).points.at(vi2);
			}

			vi1 = a.idx;
			vi2 = b.idx;
			vi3 = c.idx;

			circumcenter(a,b,c,ccenter);
			ccenter.normalize();
			cradius = circumradius(a, b, c);

			criteria = (*region_itr).center.dotForAngle(ccenter);
			criteria += cradius;

			if(criteria < (*region_itr).radius){
				a_dist = a.dotForAngle((*region_itr).center);
				b_dist = b.dotForAngle((*region_itr).center);
				c_dist = c.dotForAngle((*region_itr).center);

				a_dist_min = 10;
				b_dist_min = 10;
				c_dist_min = 10;

				for(neighbor_itr = (*region_itr).neighbors.begin(); neighbor_itr != (*region_itr).neighbors.end(); ++neighbor_itr){
					dist_temp = a.dotForAngle(regions.at((*neighbor_itr)).center);
					a_dist_min = min(a_dist_min,dist_temp);
					dist_temp = b.dotForAngle(regions.at((*neighbor_itr)).center);
					b_dist_min = min(b_dist_min,dist_temp);
					dist_temp = c.dotForAngle(regions.at((*neighbor_itr)).center);
					c_dist_min = min(c_dist_min,dist_temp);
				}

				if((a_dist < a_dist_min) || (b_dist < b_dist_min) || (c_dist < c_dist_min)){
					t = tri(vi1, vi2, vi3);
					t = t.sortedTri();
					(*region_itr).triangles.push_back(t);
				}
			}
		}
		free(in.pointlist);
		free(in.regionlist);
		free(in.pointmarkerlist);
		free(out.pointlist);
		free(out.trianglelist);
	}
#ifdef _DEBUG
	cerr << "Done triangulating points (local) " << id << endl;
#endif
	return;
}/*}}}*/
void projectToBoundary(vector<region> &region_vec){/*{{{*/
	double new_min_dist, min_dist, dist, alpha, beta;	
	pnt p, p_n;
	pnt a, b, c;
	pnt v1, v2, v3;
	vector<int> closest_cell;
	vector<int> proj_1_point;
	vector<int> proj_2_point;
	vector<int> search_flags;
	int i, index1, index2;

	for(i = 0; i < boundary_points.size(); i++){
		closest_cell.push_back(-1);
	}
	for(i = 0; i < points.size(); i++){
		proj_1_point.push_back(-1);
		proj_2_point.push_back(-1);
		search_flags.push_back(1);
	}

	for(region_itr = region_vec.begin(); region_itr != region_vec.end(); region_itr++){
		for(boundary_itr = (*region_itr).boundary_points.begin(); boundary_itr != (*region_itr).boundary_points.end(); boundary_itr++){

			min_dist = M_PI;
			for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
				dist = (*boundary_itr).dotForAngle((*point_itr));

				if(dist < min_dist){
					min_dist = dist;
					closest_cell.at((*boundary_itr).idx) = (*point_itr).idx;
				}
			}

		}
	}

	for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
		min_dist = M_PI;

		for(i = 0; i < closest_cell.size(); i++){
			if(closest_cell.at(i) == (*point_itr).idx){
				dist = (*point_itr).dotForAngle(boundary_points.at(i));

				if(dist < min_dist){
					min_dist = dist;
					proj_1_point.at((*point_itr).idx) = i;
				}
			}
		}

		if(proj_1_point.at((*point_itr).idx) > -1){
			closest_cell.at(proj_1_point.at((*point_itr).idx)) = -1;
		}

		min_dist = M_PI;

		for(i = 0; i < closest_cell.size(); i++){
			if(closest_cell.at(i) == (*point_itr).idx){
				dist = (*point_itr).dotForAngle(boundary_points.at(i));

				if(dist < min_dist){
					min_dist = dist;
					proj_2_point.at((*point_itr).idx) = i;
				}
			}
		}

	}

	for(point_itr = n_points.begin(); point_itr != n_points.end(); point_itr++){
		if(proj_1_point.at((*point_itr).idx) > -1){
			index1 = proj_1_point.at((*point_itr).idx);
			index2 = proj_2_point.at((*point_itr).idx);
			if(index2 > -1){
				a = boundary_points.at(index1);
				b = boundary_points.at(index2);

				v1 = a - (*point_itr);
				v2 = b - a;

				alpha = v2.magnitude();
				alpha = alpha*alpha;
				alpha = -v1.dot(v2)/alpha;

				if(alpha < 0.0){
					p = boundary_points.at(index1);
					p.idx = (*point_itr).idx;
					p.isBdry = 0;
					p.normalize();
				} else {
					p = a + alpha * (b - a);

					p.idx = (*point_itr).idx;
					p.isBdry = 0;
					p.normalize();
				}
			} else {
				p = boundary_points.at(index1);
				p.idx = (*point_itr).idx;
				p.isBdry = 0;
				p.normalize();
			}

			alpha = std::min(1.0, proj_alpha);
			p_n = alpha*p + (1.0 - alpha) * (*point_itr);
			p_n.normalize();
			p_n.idx = (*point_itr).idx;
			p_n.isBdry = 0;


			(*point_itr).x = p_n.x;
			(*point_itr).y = p_n.y;
			(*point_itr).z = p_n.z;
			(*point_itr).idx = p_n.idx;
			(*point_itr).isBdry = p_n.isBdry;
		}
	}

}/*}}}*/
/* }}} */
/* ***** Specifc Region Routines ***** {{{*/
void printAllFinalTriangulation(){/*{{{*/
	//Single processor version of printMyFinalTriangulation
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;

	for(region_itr = regions.begin(); region_itr != regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			temp_tris_out.push_back((*tri_itr));
			unique_tris.insert((*tri_itr).sortedTri());
		}
	}

	ofstream utris_out("triangles.dat");

	for(utri_itr = unique_tris.begin(); utri_itr != unique_tris.end(); ++utri_itr){
		t = (*utri_itr);

		if(!isCcw(points.at(t.vi1),points.at(t.vi2),points.at(t.vi3))){
			int swp_v;

			swp_v = t.vi2;
			t.vi2 = t.vi3;
			t.vi3 = swp_v;
		}
		utris_out << t << endl;
	}
	utris_out.close();
}/*}}}*/
void printMyFinalTriangulation_t(){/*{{{*/
	//Merge all finalTriangulations made with makeFinalTriangulations onto master processor.
	//Master processor inserts all triangles into an unordered_set to create a list of unique triangles
	//Triangles are then written into a file triangles.dat in ccw order
	mpi::request mycomm;
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;

	#ifdef _DEBUG
		cerr << " Print my final triangulation " << id << endl;
	#endif
	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			temp_tris_out.push_back((*tri_itr));
			unique_tris.insert((*tri_itr).sortedTri());
		}
	}

	for(int i = 1; i < num_procs; i++){
		if(id == i){
			mycomm = world.isend(master,msg_tri_print,temp_tris_out);	
		} else if (id == master){
			temp_tris_in.clear();
			world.recv(i, msg_tri_print, temp_tris_in);
			
			for(tri_itr = temp_tris_in.begin(); tri_itr != temp_tris_in.end(); ++tri_itr){
				unique_tris.insert((*tri_itr).sortedTri());
			}
		}
	}

	if(id != master){
		mycomm.wait();
		temp_tris_out.clear();
	} else {
		ofstream utris_out("triangles.dat");

		for(utri_itr = unique_tris.begin(); utri_itr != unique_tris.end(); ++utri_itr){
			t = (*utri_itr);

			if(!isCcw(points.at(t.vi1),points.at(t.vi2),points.at(t.vi3))){
				int swp_v;

				swp_v = t.vi2;
				t.vi2 = t.vi3;
				t.vi3 = swp_v;
			}
//			utris_out << (*utri_itr) << endl;
			utris_out << t << endl;
		}

		utris_out.close();
	}

	#ifdef _DEBUG
		cerr << " Print my final triangulation done " << id << endl;
	#endif
}/*}}}*/
void printMyFinalTriangulation(){/*{{{*/
	//Merge all finalTriangulations made with makeFinalTriangulations onto master processor.
	//Master processor inserts all triangles into an unordered_set to create a list of unique triangles
	//Triangles are then written into a file triangles.dat in ccw order
	mpi::request mycomm;
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;


	#ifdef _DEBUG
		cerr << " Print my final triangulation " << id << endl;
	#endif

	storeMyFinalTriangulation();
	if(id == master){
		ofstream tris_out("triangles.dat");

		for(tri_itr = all_triangles.begin(); tri_itr != all_triangles.end(); ++tri_itr){
			t = (*tri_itr);

			if(!isCcw(points.at(t.vi1),points.at(t.vi2),points.at(t.vi3))){
				int swp_v;

				swp_v = t.vi2;
				t.vi2 = t.vi3;
				t.vi3 = swp_v;
			}
			tris_out << t << endl;
		}

		tris_out.close();
	}

	all_triangles.clear();

	#ifdef _DEBUG
		cerr << " Print my final triangulation done " << id << endl;
	#endif
}/*}}}*/
void storeMyFinalTriangulation(){/*{{{*/
	//Merge all finalTriangulations made with makeMyFinalTriangulations onto master processor.
	//Master processor inserts all triangles into an unordered_set to create a list of unique triangles
	//Triangles are then written into a file triangles.dat in ccw order
	mpi::request mycomm;
	vector<tri> temp_tris_out;
	vector<tri> temp_tris_in;
	unordered_set<tri, tri::hasher> unique_tris;
	unordered_set<tri, tri::hasher>::iterator utri_itr;
	tri t;

	#ifdef _DEBUG
		cerr << " Store my final triangulation " << id << endl;
	#endif

	for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
		for(tri_itr = (*region_itr).triangles.begin(); tri_itr != (*region_itr).triangles.end(); ++tri_itr){
			temp_tris_out.push_back((*tri_itr));
			unique_tris.insert((*tri_itr).sortedTri());
		}
	}

	for(int i = 1; i < num_procs; i++){
		if(id == i){
			mycomm = world.isend(master,i,temp_tris_out);	
		} else if (id == master){
			world.recv(i, i, temp_tris_in);
			
			for(tri_itr = temp_tris_in.begin(); tri_itr != temp_tris_in.end(); ++tri_itr){
				unique_tris.insert((*tri_itr).sortedTri());
			}
		}
	}

	if(id != master){
		mycomm.wait();
		temp_tris_out.clear();
	} else {

		all_triangles.clear();
		for(utri_itr = unique_tris.begin(); utri_itr != unique_tris.end(); ++utri_itr){
			t = (*utri_itr);

			if(!isCcw(points.at(t.vi1),points.at(t.vi2),points.at(t.vi3))){
				int swp_v;

				swp_v = t.vi2;
				t.vi2 = t.vi3;
				t.vi3 = swp_v;
			}
			all_triangles.push_back(t);
		}
	}

	#ifdef _DEBUG
		cerr << " Store my final triangulation done " << id << endl;
	#endif
}/*}}}*/
/*}}}*/
/* ***** Communication Routines ***** {{{*/
void transferUpdatedPoints(){/*{{{*/
	//Each processor transfers it's updated point set (stored in n_points) to it's region neighbors
	//as defined in RegionTriangulation
	//
	//This keeps communications minimal and only updates the points that need to be updated for each region.
	vector<pnt> temp_points_in;
	vector<pnt> temp_points_out;
	vector<mpi::request> comms;
	optional options;

#ifdef _DEBUG
	cerr << "Transfering updated points " << id << endl;
#endif

	if(num_procs > 1){
		temp_points_out.clear();
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			points.at((*point_itr).idx) = (*point_itr);
			temp_points_out.push_back((*point_itr));
		}

		for(region_itr = my_regions.begin(); region_itr != my_regions.end(); ++region_itr){
			comms.resize((*region_itr).neighbors.size());

			for(int i = 0; i < (*region_itr).neighbors.size(); i++){
				comms[i] = world.isend((*region_itr).neighbors.at(i), msg_points, temp_points_out);
			}

			for(int i = 0; i < (*region_itr).neighbors.size(); i++){
				temp_points_in.clear();
				world.recv((*region_itr).neighbors.at(i), msg_points, temp_points_in);

				for(point_itr = temp_points_in.begin(); point_itr != temp_points_in.end(); ++point_itr){
					points.at((*point_itr).idx) = (*point_itr);
				}
			}
			
			mpi::wait_all(&comms[0],&(comms[(*region_itr).neighbors.size()]));

			comms.clear();
		}
	} else {
		points.swap(n_points);
	}

	temp_points_in.clear();
	temp_points_out.clear();
#ifdef _DEBUG
	cerr << "Done Transfering updated points " << id << endl;
#endif
}/*}}}*/
void gatherAllUpdatedPoints(){/*{{{*/
	//All updated points (stored in n_points) are gathered onto the master processor and broadcasted to every processor
	vector<pnt> temp_points;
	mpi::request mycomm;
	optional options;
#ifdef _DEBUG
	cerr << "Gatering all updated points " << id << endl;
#endif

	if(num_procs > 1){
		temp_points.clear();
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			temp_points.push_back((*point_itr));
			if(id == master){
				points.at((*point_itr).idx) = (*point_itr);
			}
		}
		for(int i = 1; i < num_procs; i++){
			if(id == i){
				mycomm = world.isend(master, msg_points, temp_points);
				points.clear();
			}else if(id == master){
				world.recv(i, msg_points, temp_points);
				for(point_itr = temp_points.begin(); point_itr != temp_points.end(); ++point_itr){
					points.at((*point_itr).idx) = (*point_itr);
				}
				temp_points.clear();
			}
		}
		if(id != master){
			points.clear();
			options = mycomm.test();
			if(!options){
				mycomm.wait();
			}
			temp_points.clear();
		}
		mpi::broadcast(world,points,master);
	} else {
		for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
			points.at((*point_itr).idx) = (*point_itr);
		}
	}
#ifdef _DEBUG
	cerr << "Done Gatering all updated points " << id << endl;
#endif
}/*}}}*/
/*}}}*/
/* ***** Routines for Points ***** {{{*/
void writePointsAsRestart(const int it, const restart_mode_type restart_mode, const fileio_mode_type fileio_mode){/*{{{*/

	#ifdef _DEBUG
		cerr << "Writing restart file " << id << endl;
	#endif

	if(it > 0){
		gatherAllUpdatedPoints();

		if(id == master){

		#ifdef USE_NETCDF
			cout << "restart mode " << restart_mode << " fileio_mode " << fileio_mode << endl;
			if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_NETCDF ) {
				writeRestartFileOverwriteNC( it, points );
			} else if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_TXT ) {
				writeRestartFileOverwriteTXT( it );
		    } else if ( restart_mode == RESTART_OVERWRITE && fileio_mode == FILEIO_BOTH ) {
		    	writeRestartFileOverwriteNC( it, points );
		    	writeRestartFileOverwriteTXT( it );
		    } else if ( restart_mode == RESTART_RETAIN  && fileio_mode == FILEIO_NETCDF ) {
		    	writeRestartFileRetainNC( it, points );
			} else if ( restart_mode == RESTART_RETAIN && fileio_mode == FILEIO_TXT ) {
				writeRestartFileRetainTXT( it );
			} else if ( restart_mode == RESTART_RETAIN && fileio_mode == FILEIO_BOTH ) {
				writeRestartFileRetainNC( it, points );
				writeRestartFileRetainTXT( it );
			}
		#else
			if ( restart_mode == RESTART_OVERWRITE ) {
				writeRestartFileOverwriteTXT( it );
			} else if ( restart_mode == RESTART_RETAIN ) {
				writeRestartFileRetainTXT( it );
		    }
		#endif
		}
		
		/*
		if(id == master){
			ofstream pts_out(temp);
			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}
			pts_out.close();

			if(num_procs > 1){
				world.isend(id+1, msg_restart, junk);
			}
		} else if(id == num_procs-1){
			world.recv(id-1, msg_restart, junk);
			ofstream pts_out(temp, fstream::app);

			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}

			pts_out.close();
		} else {
			world.recv(id-1, msg_restart, junk);
			ofstream pts_out(temp, fstream::app);

			for(point_itr = n_points.begin(); point_itr != n_points.end(); ++point_itr){
				pts_out << (*point_itr) << endl;
			}

			pts_out.close();
			world.isend(id+1, msg_restart, junk);
		}
		*/
	} else {
		return;
	}

	#ifdef _DEBUG
		cerr << "Done Writing restart file " << id << endl;
	#endif
}/*}}}*/

#ifdef USE_NETCDF
int writeRestartFileOverwriteNC( const int it, const vector<pnt> &points ) {

	// set up the file and create
	static const int NC_ERR = 2;
	NcError err(NcError::verbose_nonfatal);
	NcFile grid("point_restart.nc", NcFile::Replace);
	if(!grid.is_valid()) return NC_ERR;
	
	// define dimensions
	NcDim *nPointsDim;
	NcDim *itDim;
	
	// store dimensions
	if (!(nPointsDim = grid.add_dim( "nPoints",	points.size() ))) return NC_ERR;
	if (!(itDim = grid.add_dim( "it", it ))) return NC_ERR;
		
	// create/populate point arrays
	double *xPoint, *yPoint, *zPoint;
	xPoint = new double[nPointsDim->size()];
	yPoint = new double[nPointsDim->size()];
	zPoint = new double[nPointsDim->size()];
	
	for ( int ii = 0; ii < nPointsDim->size(); ii++ ) {
		xPoint[ii] = points.at(ii).x;
		yPoint[ii] = points.at(ii).y;
		zPoint[ii] = points.at(ii).z;
	}
	
	// define coordinate ncvars
	NcVar *xPointVar, *yPointVar, *zPointVar;
	if (!(xPointVar = grid.add_var("xPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(yPointVar = grid.add_var("yPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(zPointVar = grid.add_var("zPoint", ncDouble, nPointsDim))) return NC_ERR;
	
	// write coordinate data
	if (!xPointVar->put(xPoint,nPointsDim->size())) return NC_ERR;
	if (!yPointVar->put(yPoint,nPointsDim->size())) return NC_ERR;
	if (!zPointVar->put(zPoint,nPointsDim->size())) return NC_ERR;
	
	// clear memory
	delete xPoint;
	delete yPoint;
	delete zPoint;
	
	// scope closure destroys NC objs
	return 0;
	
}

int writeRestartFileRetainNC( const int it, const vector<pnt> &points ) {

	// set up the file and create
	static const int NC_ERR = 2;
	NcError err(NcError::verbose_nonfatal);
	std::ostringstream int_hole;
	int_hole << it;
	std::string restart_name = "point_restart_" + int_hole.str() + ".nc";
	NcFile grid(restart_name.c_str(), NcFile::Replace);
	if(!grid.is_valid()) return NC_ERR;

	// define dimensions
	NcDim *nPointsDim;
	NcDim *itDim;
	
	// store dimensions
	if (!(nPointsDim = grid.add_dim( "nPoints",	points.size() ))) return NC_ERR;
	if (!(itDim = grid.add_dim( "it", it ))) return NC_ERR;
		
	// create/populate point arrays
	double *xPoint, *yPoint, *zPoint;
	xPoint = new double[nPointsDim->size()];
	yPoint = new double[nPointsDim->size()];
	zPoint = new double[nPointsDim->size()];
	
	for ( int ii = 0; ii < nPointsDim->size(); ii++ ) {
		xPoint[ii] = points.at(ii).x;
		yPoint[ii] = points.at(ii).y;
		zPoint[ii] = points.at(ii).z;
	}
	
	// define coordinate ncvars
	NcVar *xPointVar, *yPointVar, *zPointVar;
	if (!(xPointVar = grid.add_var("xPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(yPointVar = grid.add_var("yPoint", ncDouble, nPointsDim))) return NC_ERR;
	if (!(zPointVar = grid.add_var("zPoint", ncDouble, nPointsDim))) return NC_ERR;
	
	// write coordinate data
	if (!xPointVar->put(xPoint,nPointsDim->size())) return NC_ERR;
	if (!yPointVar->put(yPoint,nPointsDim->size())) return NC_ERR;
	if (!zPointVar->put(zPoint,nPointsDim->size())) return NC_ERR;
	
	// clear memory
	delete xPoint;
	delete yPoint;
	delete zPoint;
	
	// scope closure destroys NC objs
	return 0;

}
#endif

int writeRestartFileOverwriteTXT ( const int it ) 
{
  char temp[32];

//sprintf(temp,"point_restart.dat\0");
  sprintf(temp,"point_restart.dat");

  ofstream pts_out(temp);

  for ( point_itr = points.begin(); point_itr != points.end(); ++point_itr )
  {
    pts_out << (*point_itr) << endl;
  }	
  pts_out.close();

  return 0;	
}

int writeRestartFileRetainTXT ( const int it ) 
{	
  char temp[32];

//sprintf ( temp, "point_restart_%d.dat\0", it );
  sprintf ( temp, "point_restart_%d.dat", it );
  ofstream pts_out ( temp );

  for ( point_itr = points.begin(); point_itr != points.end(); ++point_itr )
  {
    pts_out << ( *point_itr ) << endl;
  }	
  pts_out.close ( );

  return 0;
}

double density(const pnt &p){/*{{{*/
	//density returns the value of the density function at point p
	return 1.0; // Uniform density

	/* Density function for Shallow Water Test Case 5 
	pnt cent;
	double r;
	double norm;
	double density;
	double min_val;
	double width, trans_cent;

	cent = pnt(0.0, -0.866025403784439, 0.5);
	cent.normalize();

	width = 0.15;
	min_val = 1.0/8.0;
	min_val = pow(min_val,4);
	trans_cent = 30.0*M_PI/180.0;
	norm = 1.0/(1.0-min_val);

	r = cent.dotForAngle(p);

	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2.0)/norm + min_val;

	return density;
	// */
	
	// /* Ellipse density function.
	
	return ellipse_density(p, 40.0, 0.0, 1.0, 0.5);
	// */

}/*}}}*/
double ellipse_density(const pnt &p, double lat_c, double lon_c, double lat_width, double lon_width){/*{{{*/
	//density returns the value of the density function at point p
	//	return 1.0; // Uniform density
	//	/* Ellipse Density function
	pnt work;
	double r1, r2, r;
	double dtr;
	double width, trans_cent, min_val, norm, density;

	dtr = M_PI/180.0;

	work = pntFromLatLon(p.getLat(), lon_c*dtr);
	r1 = work.dotForAngle(p);

	work = pntFromLatLon(lat_c*dtr, p.getLon());
	r2 = work.dotForAngle(p);

	r1 = r1/lon_width;
	r2 = r2/lat_width;
	r = sqrt (r1*r1 + r2*r2);

	width = 0.15;
	trans_cent = 30.0 * dtr;
	min_val = 1.0/12.0;
	min_val = pow(min_val,4);
	norm = 1.0/(1.0-min_val);
	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2)/norm + min_val;

	return density;
	// */

}/*}}}*/
/*}}}*/


