// Author : Eric Thiemard 
// Version : 20.10.99       
// This program is an implementation of the method proposed in Eric 
// Thiemard's paper "An algorithm to compute bounds for the star    
// discrepancy" submitted to Journal of Complexity. 

# include <cstdlib>
# include <iostream>
# include <fstream>
# include <cmath>
# include <ctime>
# include <iomanip>
# include <cstring>

using namespace std;
//
//  GLOBAL DATA
//
typedef struct som sommet;
struct som
{
  int id;
  sommet *pt;
};
static int s;
static int n;
static int num;
static int den;
static int num2;
static int subtotala = 0;
static int subtotalb = 0;
static int *lexsizealpha;
static int *maxsizealpha;
static int *lexsizebeta;
static int *maxsizebeta;
static double epsilon;
static double borne_sup = 0.0;
static double borne_inf = 0.0;
static double *points;
static sommet *superarbre;
static sommet **lexalpha;
static sommet **lexbeta;
//
//  Declaration of routines.
//
int main ( int argc, char *argv[] );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double *data_read ( char *file_in_name, int m, int n );
void decomposition ( double *alpha, double *beta, int min, double value );
int explore ( sommet *liste, double *pave, int dim );
int fastexplore ( double *pave, int range, int *maxsize, int *lexsize,
  sommet **lex, int *subtotal );
int file_column_count ( char *file_in_name );
int file_row_count ( char *file_in_name );
static void file_usage ( );
void freetree ( sommet *noeud );
void initlex ( );
void initparameters ( int argc, char *argv[] );
void insertlex ( sommet *noeud, int range, int *maxsize, int *lexsize, 
  sommet **lex );
double lowbound ( int npoints, double volume, double *pave );
static void memory ( );
void quicksort ( sommet *liste, int dim, int l, int r );
void readfile ( char *filename );
int s_len_trim ( char *s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double dvec[] );
int s_word_count ( char *s );
sommet *subtree ( sommet *liste, int min, int next, int dim );
void supertree ( );
void timestamp ( );
void traiter ( double *outputalpha, double *outputbeta, int range );
static void usage ( char *exec_name );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
// 
//    MAIN is the main program for the star discrepancy bound computation.
//
//  Modified:
//
//    30 September 2003
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  int j;
  int n_max;
  double *oalpha;
  double *obeta;

  timestamp ( );
  cout << "\n";
  cout << "TABLE_DISCREPANCY:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  A program to estimate the discrepancy\n";
  cout << "  of a set of N points in M dimensions.\n";
//
//  Get the values of EPSILON, N, and FILENAME from the command line.
//
  initparameters ( argc, argv );
//
//  Get the spatial dimension by reading the file.
//
  s = file_column_count ( argv[3] );

  if ( s < 2 ) 
  {
    cout << "\n";
    cout << "DISCRBOUND - Fatal error!\n";
    cout << "  The spatial dimension S must be at least 2.\n";
    exit ( 1 );
  }
//
//  Determine the number of points avaiable in the file.
//
  n_max = file_row_count ( argv[3] );

  if ( n_max < n )
  {
    cout << "\n";
    cout << "DISCRBOUND - Fatal error!\n";
    cout << "  The number of data points requested is N = " << n << ".\n";
    cout << "  The data file only contains N_MAX = " << n_max << ".\n";
    exit ( 1 );
  }
//
//  Read the data.
//
  points = data_read ( argv[3], s, n );
//
//  Make sure the points are in the unit hypercube.
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < s; j++ )
    {
      if ( points[j+i*s] < 0.0 || 1.0 < points[j+i*s] )
      {
        cout << "\n";
        cout << "DISCRBOUND - Fatal error!\n";
        cout << "  Every coordinate of every point must be between 0 and 1.\n";
        cout << "  However, coordinate " << j << " of point " << i << "\n";
        cout << "  is equal to " << points[j+i*s] << "\n";
        exit ( 1 );
      }
    }
  }
//
//  Print the points, if there aren't too many.
//
  if ( n <= 20 )
  {
    cout << "\n";
    cout << "x = \n";
    cout << "{\n";
    for ( i = 0; i < n; i++ )
    {
      cout << "  {";
      for ( j = 0; j < s; j++ )
      {
        cout << setw(10) << points[j+i*s] << "  ";
      }
      cout << "  }\n";
    }
    cout << "}\n";
    cout << "\n";
  }

  supertree ( );

  initlex ( );
  
  oalpha = (double *) calloc ( (unsigned) s, sizeof(double) );
  obeta  = (double *) calloc ( (unsigned) s, sizeof(double) );

  for ( i = 0; i < s; i++ )
  {
    obeta[i] = 1.0;
  }

  decomposition ( oalpha, obeta, 0, 1.0 );
//
//  Report the results.
//
  cout << "\n";
  cout << "  S = " << s << "\n";
  cout << "  Epsilon = " << epsilon << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Estimate of Discrepancy D_n^*(x):\n";
  cout << "\n";
  cout << "   Lower         Upper\n";
  cout << "   bound         bound\n";
  cout << "------------  ------------\n";
  cout << "\n";
  cout << setw(12) << borne_inf << "  " 
       << setw(12) << borne_sup << "\n";

  cout << "\n";
  cout << "DISCRBOUND:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80*

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80*
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
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
//
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

double *data_read ( char *file_in_name, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_READ reads coordinate data from a file.
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
//  Modified:
//
//    30 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double DATA_READ[M*N], the point coordinates.
//
{
  double *coord;
  bool error;
  ifstream file_in;
  int i;
  int j;
  char line[255];
  double *x;

  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << file_in_name << "\"\n";
    exit ( 1 );
  }

  coord = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof ( ) )
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
      coord[i+j*m] = x[i];
    }
    j = j + 1;

  }

  file_in.close ( );

  delete [] x;

  return coord;
}
//****************************************************************************80

void decomposition ( double *alpha, double *beta, int min, double value )

//****************************************************************************80
//
//  Purpose:
// 
//    DECOMPOSITION carries out the decomposition of a subinterval.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  double pbetaminp = 1.0;
  double palpha = 1.0;
  double pbeta;
  double delta;
  double *subalpha;
  double *subbeta;
  double *gamma;

  subalpha = (double *) calloc ( (unsigned) s, sizeof(double) );
  subbeta  = (double *) calloc ( (unsigned) s+1, sizeof(double) );
  gamma  = (double *) calloc ( (unsigned) s+1, sizeof(double) );

  for ( i = min; i < s; i++ )
  {
    pbetaminp = pbetaminp * beta[i];
  }
  pbeta = pbetaminp;

  for ( i = 0; i < min; i++ )
  {
    pbetaminp = pbetaminp * beta[i];
    palpha = palpha * alpha[i];
  }

  pbetaminp = pbetaminp - epsilon;
  delta = pow ( pbetaminp / ( pbeta * palpha ), 1.0 / ( s - min ) );

  for ( i = 0; i < min; i++ )
  {
    gamma[i] = alpha[i];
    subalpha[i] = gamma[i];
    subbeta[i] = beta[i];
  }

  for ( i = min; i < s; i++ )
  {
    gamma[i] = delta * beta[i];
    subalpha[i] = alpha[i];
    subbeta[i] = beta[i];
  }

  subbeta[min] = gamma[min];

  value = value * delta;
  if ( epsilon < value )
  {
    for ( i = min; i < s; i++ )
    {
      decomposition ( subalpha, subbeta, i, value );
      subalpha[i]  = gamma[i];
      subbeta[i]   = beta[i];
      subbeta[i+1] = gamma[i+1];
    }
  }
  else
  {
    for ( i = min; i < s; i++ )
    {
      traiter ( subalpha, subbeta, (i==0)?0:i-1 );
      subalpha[i]  = gamma[i];
      subbeta[i]   = beta[i];
      subbeta[i+1] = gamma[i+1];
    }
  }

  traiter ( gamma, beta, s-1 );

  free ( gamma );
  free ( subalpha );
  free ( subbeta );

  return;
}
//****************************************************************************80

int explore ( sommet *liste, double *pave, int dim )

//****************************************************************************80
//
//  Purpose:
// 
//    EXPLORE ???
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  int min = 1;
  int max;
  int next;
  int total;

  if ( pave[dim] <= points[dim+(liste[1].id)*s] )
  {
    return 0;
  }

  if ( liste[0].id == 1 )
  {
    total = 1;
    next = liste[1].id; 
    for ( i = dim; i < s; i++ )
    {
      if ( pave[i] <= points[i+next*s] )
      {
        total = 0;
        break;
      }
    }
  }
  else
  {
    total = 0;
    max = liste[0].id;

    if ( dim == s-1 )
    {
      if ( points[dim+(liste[max].id)*s] < pave[dim] )
      {
        total = max;
      }
      else
      {
        while ( min <= max )
        {
          next = ( min + max + 1 ) / 2;
          if ( points[dim+(liste[next].id)*s] < pave[dim] )
          {    
            total = total + next - min + 1;
            min = next + 1;
          }
          else
          {
            max = next - 1;
          }
        }
      }
    }
    else
    {
      while ( min <= max )
      {
        next = ( ( 1 + min ) * num2 + max * num ) / den;
        if ( points[dim+(liste[next].id)*s] < pave[dim] )
        {
          if ( liste[next].pt == NULL )
          {
            liste[next].pt = subtree ( liste, min, next, dim+1 );
          }
          total = total + explore ( liste[next].pt, pave, dim+1 );
          min = next + 1;
        }
        else
        {
          max = next - 1;
        }
      }
    }
  }
  return total;
}
//****************************************************************************80

int fastexplore ( double *pave, int range, int *maxsize, int *lexsize,
  sommet **lex, int *subtotal )

//****************************************************************************80
//
//  Purpose:
// 
//    FASTEXPLORE ???
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int j;
  int i;
  int min;
  int max;
  int next;
  int start;
  int size = lexsize[range];
  int right;
  int total = 0;
  double seuil = pave[range];
  sommet refnoeud;
  sommet *noeud;

  if ( range == s - 1 )
  {
    for ( i = size-1; 0 <= i; i-- )
    {
      refnoeud = lex[range][i];
      noeud = refnoeud.pt;
      min = refnoeud.id;
      max = noeud[0].id;
      if ( max < min )
      {
        lexsize[range]--;
        lex[range][i] = lex[range][lexsize[range]];
        *subtotal = *subtotal + min - 1;
      }
      else
      {
        total = total + min - 1;
        right = 1;
        while ( min <= max )
        {
          next = ( min + max + 1 ) / 2;
          if ( points[range+(noeud[next].id)*s] < seuil )
          {    
            total = total + next - min + 1;
            min = next + 1;
            if ( right == 1 )
            {
              lex[range][i].id = min;
            }
          }
          else
          {
            right = 0;
            max = next - 1;
          }
        }
      }
    }
    total = total + *subtotal;
  }
  else
  {
    *subtotal = 0;
    lexsize[range+1] = 0;
    for ( i = 0; i < size; i++ )
    {
      refnoeud = lex[range][i];
      noeud = refnoeud.pt;
      start = refnoeud.id;
      min = 1;
      max = noeud[0].id;
      while ( min != start )
      {
        next = ( ( 1 + min ) * num2 + max * num ) / den;
        insertlex ( noeud[next].pt, range+1, maxsize, lexsize, lex );
        total = total + explore ( noeud[next].pt, pave, range+1 );
        min = next + 1;
      }
      right = 1;
      while ( min <= max )
      {
        next = ( ( 1 + min ) * num2 + max * num ) / den;
        if ( points[range+(noeud[next].id)*s] < seuil )
        {
          if ( noeud[next].pt == NULL )
          {
            noeud[next].pt = subtree ( noeud, min, next, range+1 );
          }
          insertlex ( noeud[next].pt, range+1, maxsize, lexsize, lex );
          total = total + explore ( noeud[next].pt, pave, range+1 );
          min = next + 1;
          if ( right == 1 )
          {
            if ( range == 0 )
            {
              if ( lex == lexalpha )
              {
                for ( j = lex[range][i].id; j < next; j++ )
                {
                  if ( noeud[j].pt != NULL )
                  {
                    freetree ( noeud[j].pt );
                  }
                }
              }
            }
            lex[range][i].id = min;
          }
        }
        else
        {
          right = 0;
          max = next - 1;
        }
      }
    }
  }
  return total;
}
//****************************************************************************80

int file_column_count ( char *file_in_name )

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
//    Most lines of the file is presumed to consist of NCOLUMN words, separated
//    by spaces.  There may also be some blank lines, and some comment lines,
//    which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
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
//    Input, char *FILE_IN_NAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  ifstream file_in;
  bool got_one;
  char line[256];
  int ncolumn;
//
//  Open the file.
//
  file_in.open ( file_in_name );

  if ( !file_in )
  {
    ncolumn = -1;
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Fatal error!\n";
    cout << "  Could not open the file:\n";
    cout << "  \"" << file_in_name << "\"\n";
    return ncolumn;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof ( ) )
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
    file_in.close ( );

    file_in.open ( file_in_name );

    for ( ; ; )
    {
      file_in.getline ( line, sizeof ( line ) );

      if ( file_in.eof ( ) )
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

  file_in.close ( );

  if ( !got_one )
  {
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Warning!\n";
    cout << "  The file does not seem to contain any data.\n";
    return 0;
  }

  ncolumn = s_word_count ( line );

  return ncolumn;
}
//****************************************************************************80

int file_row_count ( char *file_in_name )

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
//    Input, char *FILE_IN_NAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream file_in;
  int i;
  char line[100];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "FILE_ROW_COUNT - Fatal error!\n";
    cout << "  Could not open the input file: \"" << file_in_name << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof ( ) )
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

  file_in.close ( );

  return row_num;
}
//****************************************************************************80

static void file_usage ( )

//****************************************************************************80
//
//  Purpose:
// 
//    FILE_USAGE lists the legal input file format.
//
//  Modified:
//
//    30 September 2003.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
//  Parameters:
//
//    None.
//
{
  cout << "\n";
  cout << "FILE_USAGE:\n";
  cout << "  This program requires an input file that defines\n";
  cout << "  the point set.\n";
  cout << "  We assume that the points are S-dimensional, and \n";
  cout << "  that there are N_MAX of them.\n";
  cout << "\n";
  cout << "  We also require that the points be in the unit hypercube.\n";
  cout << "\n";
  cout << "  Then the input file contains N_MAX lines or records.\n";
  cout << "  Each record contains the S coordinate values of one point.\n";
  cout << "\n";
  cout << "  Comment lines begin with the character \"#\",\n";
  cout << "  and may occur anywhere in the file.\n";
  cout << "\n";
  cout << "  At runtime, the user specifies N, the number of points\n";
  cout << "  to be read from the file, which must be less than or\n";
  cout << "  equal to N_MAX.\n";
  cout << "\n";
  cout << "  Here is an example of a file containing N_MAX = 5\n";
  cout << "  points in S = 3 dimensions:\n";
  cout << "\n";
  cout << "#  input_3x5.txt\n";
  cout << "#\n";
  cout << "#  Suitable for input to DISCRBOUND.\n";
  cout << "#  Spatial dimension = 3\n";
  cout << "#  Number of points = 5\n";
  cout << "#\n";
  cout << "  0.00  0.20  0.90\n";
  cout << "  0.80  0.40  0.30\n";
  cout << "  0.20  1.00  0.50\n";
  cout << "  0.40  0.60  0.10\n";
  cout << "  0.60  0.80  0.70\n";
  cout << "\n";

  return;
}
//****************************************************************************80

void freetree ( sommet *noeud )

//****************************************************************************80
//
//  Purpose:
// 
//    FREETREE frees storage associated with a tree.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  int max = noeud[0].id;

  for ( i = 1; i <= max; i++ )
  {
    if ( noeud[i].pt != NULL )
    {
      freetree ( noeud[i].pt );
    }
  }

  free ( noeud );

  return;
}
//****************************************************************************80

void initlex ( )

//****************************************************************************80
//
//  Purpose:
// 
//    INITLEX initializes the lexicon.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;

  maxsizealpha = (int *) calloc ( (unsigned) s, sizeof(int) );

  for ( i = 0; i < s; i++ )
  {
    maxsizealpha[i] = 1;
  }

  lexsizealpha = (int *) calloc ( (unsigned) s, sizeof(int) );
  lexsizealpha[0] = 1;
  lexalpha = (sommet **) calloc ( (unsigned) s, sizeof(sommet *) );

  for ( i = 0; i < s; i++ )
  {
    lexalpha[i] = (sommet *) calloc ( (unsigned) maxsizealpha[i], sizeof(sommet) );
  }

  lexalpha[0][0].id = 1;
  lexalpha[0][0].pt = superarbre;

  maxsizebeta = (int *) calloc ( (unsigned) s, sizeof(int) );

  for ( i = 0; i < s; i++ )
  {
    maxsizebeta[i] = 1;
  }

  lexsizebeta = (int *) calloc ( (unsigned) s, sizeof(int) );
  lexsizebeta[0] = 1;
  lexbeta = (sommet **) calloc ( (unsigned) s, sizeof(sommet *) );
 
  for ( i = 0; i < s; i++ )
  {
    lexbeta[i] = (sommet *) calloc ( (unsigned) maxsizebeta[i], sizeof(sommet) );
  }

  lexbeta[0][0].id = 1;
  lexbeta[0][0].pt = superarbre;

  return;
}
//****************************************************************************80

void initparameters ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
// 
//    INITPARAMETERS sets program parameters based on user input and defaults.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  if ( ( argc != 4 ) && ( argc != 6 ) )
  {
    usage ( argv[0] );
    exit ( 1 );
  }

  epsilon = atof ( argv[1] );

  if ( epsilon <= 0.0 || 1.0 <= epsilon )
  {
    cout << "\n";
    cout << "INITPARAMETERS: Fatal error!\n";
    cout << "  The error tolerance EPSILON = " << epsilon << ".\n";
    cout << "  But EPSILON must be between 0 and 1.\n";
    exit ( 1 );
  }

  n = atoi ( argv[2] );

  if ( n < 1 )
  {
    usage ( argv[0] );
    exit ( 1 );
  }

  if ( argc == 6 )
  {
    num = atoi ( argv[4] );
    den = atoi ( argv[5] );
    if ( num < 1 || den <= num )
    {
      usage ( argv[0] );
      exit ( 1 );
    }
  }
  else
  {
    num = 1;
    den = 2;
  }

  num2 = den - num;

  return;
}
//****************************************************************************80

void insertlex ( sommet *noeud, int range, int *maxsize, int *lexsize, 
  sommet **lex )

//****************************************************************************80
//
//  Purpose:
// 
//    INSERTLEX inserts an item into the lexicon.
//
//  Modified:
//
//    30 September 2003
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i = lexsize[range];

  if ( i == maxsize[range] )
  {
    maxsize[range] *= 2;

    lex[range] = ( sommet * ) realloc ( lex[range], 
      (unsigned) maxsize[range]*sizeof(sommet) );

    if ( lex[range] == NULL)
    {
      memory();
    }
  }

  lex[range][i].pt = noeud;
  lex[range][i].id = 1;
  lexsize[range] = ++i;

  return;
}
//****************************************************************************80

double lowbound ( int npoints, double volume, double *pave )

//****************************************************************************80
//
//  Purpose:
// 
//    LOWBOUND computes the lower bound.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  int j;
  double tmp;

  if ( borne_inf < fabs ( volume - ( (double) npoints / n ) ) )
  {
    if ( volume < ( (double) npoints / n ) )
    {
      volume = 1.0;
      for ( j = 0; j < s; j++ )
      {
        tmp = 0.0;
        for ( i = 0; i < n; i++ )
        {
          if ( tmp < points[j+i*s] && points[j+i*s] <= pave[j] )
          {
            tmp = points[j+i*s];
          }
        }
        volume = volume * tmp;
      }
    }
    else
    {
      volume = 1.0;
      for ( j = 0; j < s; j++ )
      {
        tmp = 1.0;
        for ( i = 0; i < n; i++ )
        {
          if ( points[j+i*s] < tmp && pave[j] <= points[j+i*s] )
          {
            tmp = points[j+i*s];
          }
        }
        volume = volume * tmp;
      }
    }
    return fabs ( volume - ( (double) npoints / n ) );
  }
  else
  {
    return borne_inf;
  }
}
//****************************************************************************80

static void memory ( )

//****************************************************************************80
//
//  Purpose:
// 
//    MEMORY prints a message and terminates on memory allocation errors.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  cout << "\n";
  cout << "MEMORY - Fatal error!\n";
  cout << "  Memory allocation problem\n";

  exit ( 1 );
}
//****************************************************************************80

void quicksort ( sommet *liste, int dim, int l, int r )

//****************************************************************************80
//
//  Purpose:
// 
//    QUICKSORT uses Quicksort to sort an array.
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i = l;
  int j = r+1;
  int tmp;
  double pivot = points[dim+(liste[l].id)*s];

  while ( i < j )
  {
    do
    i++;
    while ( i < r && points[dim+(liste[i].id)*s] < pivot );

    do
    j--;
    while ( pivot < points[dim+(liste[j].id)*s] );

    if ( i < j )
    {
      tmp = liste[i].id;
      liste[i].id = liste[j].id;
      liste[j].id = tmp;
    }

  }

  tmp = liste[l].id;
  liste[l].id = liste[j].id;
  liste[j].id = tmp;

  if ( l < j-1 )
  {
    quicksort ( liste, dim, l, j-1 );
  }

  if ( j+1 < r )
  {
    quicksort ( liste, dim, j+1, r );
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
//    part of the number.
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
//  Examples:
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
//    characters can be read to form a legal value.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the value that was read from the string.
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

bool s_to_r8vec ( char *s, int n, double dvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
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
//    Output, double DVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  for ( i = 0; i < n; i++ )
  {
    dvec[i] = s_to_r8 ( s, &lchar, &error );

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
//    Input, char *S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int i;
  int nword;

  nword = 0;
  blank = true;

  while ( *s != NULL ) 
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

sommet *subtree ( sommet *liste, int min, int next, int dim )

//****************************************************************************80
//
//  Purpose:
// 
//    SUBTREE ???
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  int aux;
  int n2;
  sommet *newarbre;

  aux = min - 1;
  n2 = next - min + 1;
  newarbre = (sommet *) calloc ( (unsigned) n2+1, sizeof(sommet) );

  if ( newarbre == NULL )
  {
    memory();
  }

  for ( i = 1; i <= n2; i++ )
  {
    newarbre[i].id = liste[i+aux].id;
  }

  newarbre[0].id = n2;

  if ( 1 < n2 )
  {
    quicksort ( newarbre, dim, 1, n2 );
  }

  return newarbre;
}
//****************************************************************************80

void supertree ( )

//****************************************************************************80
//
//  Purpose:
// 
//    SUPERTREE ???
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  superarbre = (sommet *) calloc ( (unsigned) n+1, sizeof(sommet) );

  if ( superarbre == NULL )
  {
    memory();
  }

  for ( i = 1; i <= n; i++ )
  {
    superarbre[i].id = i - 1;
  }

  superarbre[0].id = n;

  quicksort ( superarbre, 0, 1, n );

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

void traiter ( double *outputalpha, double *outputbeta, int range )

//****************************************************************************80
//
//  Purpose:
// 
//    TRAITER ???
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
{
  int i;
  double valpha = 1.0;
  double vbeta = 1.0;
  double newborne;
  int nalpha;
  int nbeta;

  for ( i = 0; i < s; i++ )
  {
    valpha = valpha * outputalpha[i];
    vbeta = vbeta * outputbeta[i];
  }

  nalpha = fastexplore ( outputalpha, range, maxsizealpha, lexsizealpha,
   lexalpha, &subtotala );

  nbeta = fastexplore ( outputbeta, range, maxsizebeta, lexsizebeta, 
    lexbeta, &subtotalb );

  newborne = ( (double) nbeta / n ) - valpha;

  if ( borne_sup < newborne )
  {
    borne_sup = newborne;
  }

  newborne = vbeta - ( (double) nalpha / n );

  if ( borne_sup < newborne )
  {
    borne_sup = newborne;
  }

  borne_inf = lowbound ( nalpha, valpha, outputalpha );
  borne_inf = lowbound ( nbeta,  vbeta,  outputbeta );

  return;
}
//****************************************************************************80

static void usage ( char *exec_name )

//****************************************************************************80
//
//  Purpose:
// 
//    USAGE prints usage information for the program.
//
//  Modified:
//
//    30 September 2003
//
//  Reference:
// 
//    Eric Thiemard,
//    An Algorithm to Compute Bounds for the Star Discrepancy,
//    Journal of Complexity,
//    Volume 17, pages 850-880, 2001.
//
//  Parameters:
//
//    Input, char *EXEC_NAME, the name of the executable version of this program.
//
{
  cout << "\n";
  cout << "Usage:\n";
  cout << "\n";
  cout << "  " << exec_name << " EPSILON N FILENAME\n";
  cout << "\n";
  cout << "  EPSILON,  accuracy parameter, 0.0 < epsilon < 1.0\n";
  cout << "  N,        number of points, integer, at least 1.\n";
  cout << "  FILENAME, containing the pointset.\n";
  cout << "\n";
  cout << "  " << exec_name << " EPSILON N FILENAME NUM DEN\n";
  cout << "\n";
  cout << "  Same as above, plus:\n";
  cout << "\n";
  cout << "  NUM,	optional balance numerator,\n";
  cout << "  DEN,       optional balance denominator,\n";
  cout << "		such that 0.0 < NUM / DEN < 1.0\n";
  cout << "             Default is NUM=1, DEN=2.\n"; 

  return;
}
