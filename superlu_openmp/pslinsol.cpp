# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "pssp_defs.h"

int main ( int argc, char *argv[] );
void parse_command_line ( int argc, char *argv[], int *procs, int *n,
  int *b, int *w, int *r, int *maxsup );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PSLINSOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2014
//
//  Author:
//
//    Xiaoye Li
//
{
  SuperMatrix   A;
  NCformat *Astore;
  float   *a;
  int      *asub, *xa;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  SuperMatrix   L;       /* factor L */
  SCPformat *Lstore;
  SuperMatrix   U;       /* factor U */
  NCPformat *Ustore;
  SuperMatrix   B;
  int      nrhs, ldx, info, m, n, nnz, b;
  int      nprocs; /* maximum number of processors to use. */
  int      panel_size, relax, maxsup;
  int      permc_spec;
  trans_t  trans;
  float   *xact, *rhs;
  superlu_memusage_t   superlu_memusage;

  timestamp ( );
  cout << "\n";
  cout << "PSLINSOL:\n";
  cout << "  C++/OpenMP version\n";
  cout << "  Call the OpenMP version of SuperLU to solve a linear system.\n";

  nrhs              = 1;
  trans             = NOTRANS;
  nprocs             = 1;
  n                 = 1000;
  b                 = 1;
  panel_size        = sp_ienv(1);
  relax             = sp_ienv(2);
  maxsup            = sp_ienv(3);
//
//  Check for any commandline input.
//
  parse_command_line ( argc, argv, &nprocs, &n, &b, &panel_size, &relax, 
    &maxsup );

#if ( PRNTlevel>=1 || DEBUGlevel>=1 )
  cpp_defs();
#endif

#define HB

#if defined( DEN )
    m = n;
    nnz = n * n;
    sband(n, n, nnz, &a, &asub, &xa);
#elif defined( BAND )
    m = n;
    nnz = (2*b+1) * n;
    sband(n, b, nnz, &a, &asub, &xa);
#elif defined( BD )
    nb = 5;
    bs = 200;
    m = n = bs * nb;
    nnz = bs * bs * nb;
    sblockdiag(nb, bs, nnz, &a, &asub, &xa);
#elif defined( HB )
    sreadhb(&m, &n, &nnz, &a, &asub, &xa);
#else    
    sreadmt(&m, &n, &nnz, &a, &asub, &xa);
#endif

    sCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_S, SLU_GE);
//
//  I had to add the ( NCformat* )...
//
    Astore = ( NCformat* ) A.Store;

    cout << "Dimension " << A.nrow << "x" << A.ncol << "; # nonzeros " << Astore->nnz << "\n";
    
    if (!(rhs = floatMalloc(m * nrhs))) 
    {
      SUPERLU_ABORT("Malloc fails for rhs[].");
    }

    sCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_S, SLU_GE);
    xact = floatMalloc(n * nrhs);
    ldx = n;
    sGenXtrue(n, nrhs, xact, ldx);
    sFillRHS(trans, nrhs, xact, ldx, &A, &B);

    if (!(perm_r = intMalloc(m))) 
    {
      SUPERLU_ABORT("Malloc fails for perm_r[].");
    }

    if (!(perm_c = intMalloc(n))) 
    {
      SUPERLU_ABORT("Malloc fails for perm_c[].");
    }
//
//  Get column permutation vector perm_c[], according to permc_spec:
//  0: natural ordering 
//  1: minimum degree ordering on structure of A'*A
//  2: minimum degree ordering on structure of A'+A
//  3: approximate minimum degree for unsymmetric matrices
//	
    permc_spec = 1;
    get_perm_c(permc_spec, &A, perm_c);

    psgssv ( nprocs, &A, perm_c, perm_r, &L, &U, &B, &info );
 /* 
  Inf. norm of the error 
*/
    if ( info == 0 ) 
    {
	  sinf_norm_error ( nrhs, &B, xact ); 

	  Lstore = ( SCPformat * ) L.Store;
	  Ustore = ( NCPformat * ) U.Store;

      cout << "#NZ in factor L = " << Lstore->nnz << "\n";
      cout << "#NZ in factor U = " << Ustore->nnz << "\n";
      cout << "#NZ in L+U = " << Lstore->nnz + Ustore->nnz - L.ncol << "\n";
	
	  superlu_sQuerySpace(nprocs, &L, &U, panel_size, &superlu_memusage);

	  cout << "L\\U MB " << superlu_memusage.for_lu/1024/1024
           << "  total MB needed " << superlu_memusage.total_needed/1024/1024
           << "  expansions " << superlu_memusage.expansions << "\n";
    }

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
//  
//  Terminate.
//  
  cout << "\n";
  cout << "PSLINSOL:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void parse_command_line ( int argc, char *argv[], int *procs, int *n,
  int *b, int *w, int *r, int *maxsup )

//****************************************************************************80
//
//  Purpose:
//
//    PARSE_COMMAND_LINE parses the command line.
//
//  Discussion:
//
//    The user can include command line arguments to get relaxed snode 
//    size, panel size, etc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 March 2014
//
//  Author:
//
//    Xiaoye Li
//
{
    register int c;
    extern char *optarg;

    while ( (c = getopt(argc, argv, "ht:p:n:b:w:x:s:")) != EOF ) {
	switch (c) {
	  case 'h':
	    cout << "Options: (default values are in parenthesis)\n";
	    cout << "-p <int> - number of processes     ( " << *procs << " )\n";
	    cout << "-n <int> - dimension               ( " << *n << " )\n"; 
	    cout << "-b <int> - semi-bandwidth          ( " << *b << " )\n";
	    cout << "-w <int> - panel size              ( " << *w << " )\n";
	    cout << "-x <int> - relax                   ( " << *r << " )\n";
	    cout << "-s <int> - maximum supernode size  ( " << *maxsup << " )\n";
	    exit(1);
	    break;
	  case 'p': *procs = atoi(optarg); 
	            break;
	  case 'n': *n = atoi(optarg);
	            break;
	  case 'b': *b = atoi(optarg);
	            break;
	  case 'w': *w = atoi(optarg);
	            break;
	  case 'x': *r = atoi(optarg);
	            break;
	  case 's': *maxsup = atoi(optarg);
	            break;
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
