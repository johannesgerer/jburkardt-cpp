# include <cstdlib>
# include <iostream>
# include <fstream>
# include <cstring>
# include <ctime>

using namespace std;

typedef struct dictent
 {
  struct dictent *alph;
  struct dictent *len;
  struct dictent *gl;
  char *word;
  int le;
} dent;

# define DE_SIZE    ( sizeof ( dent ) )

dent *aptr[26];
dent *cptr[26][40];
dent *gptr[40];
char  my_argv[10][50];
int   nlets[26];
int   reql[10];
int   rl;
dent *sptr[26][40];
int   stack[50];
int   tl;
int   trl;

int main ( int argc, char *argv[] );
void check_dict ( );
void do_search ( char *es, char *sw );
void find_valid_words ( );
int  get_args ( );
void initialize ( );
int  parse_args ( int nargs );
void read_dict ( );
void search ( int suggp, int sugg );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    ANAGRAM finds anagrams of a given set of letters.
//
//  Modified:
//
//    06 February 2003
//
//  Author:
//
//    James Cherry
//
{
  int nargs;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "ANAGRAM\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Written by James Cherry.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
  }

  initialize ( );

  read_dict ( );

  if ( argc == 1 )
  {
    for ( ;; )
    {
      nargs = get_args ( );

      if ( nargs == -1 )
      {
        exit( 0 );
      }
      search ( parse_args ( nargs + 2 ), nargs + 2 );
    }

  }
  else
  {

    for ( nargs = 1; nargs < argc; nargs++ )
    {
      strcpy ( my_argv[nargs], argv[nargs] );
    }

    search ( parse_args ( argc ), argc );

  }

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "ANAGRAM\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void check_dict ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHECK_DICT ??
//
{
  int i;
  dent *pt;

  for ( i = 0; i < 26; i++ )
  {
    cout << aptr[i]->word << " ";
  }

  cout << "\n";

  for ( i = 1; i < 39; i++ )
  {
    if ( sptr[0][i] != NULL )
    {
      cout << sptr[0][i]->word << " ";
    }
  }

  cout << "\n";
  pt = sptr[3][12];

  while ( pt != NULL )
  {
    cout << pt->word << " ";
    pt = pt->len;
  }

  cout << "\n";
  return;
}
//****************************************************************************80

void do_search ( char *es, char *sw )

//****************************************************************************80
//
//  Purpose:
//
//    DO_SEARCH ??
//
{
  int cl;
  int co;
  int cp;
  int i;
  int j;
  long na;
  int rs;
  int rs2;
  int slev;
  dent *sp[50];
  bool VERBOSE = false;
  int tn[26];

  na = 0;
  co = 0;

  if ( VERBOSE )
  {
    cout << "\n";
    cout << es << "  Going through the dictionary...\n";
  }

  find_valid_words();

  if ( VERBOSE )
  {
    cout << "\n";
    cout << es << "  Starting search...\n";
    cout << "\n";
  }

  for ( i = 0; i < 50; i++ )
  {
    stack[i] = 0;
    sp[i] = NULL;
  }

  cl = tl;
  slev = 0;
  stack[0] = 1;
  sp[0] = gptr[1];

  for ( ; ; )
  {
    for ( i = 0, j = 0, rs = tl, rs2 = trl; i < rl && j <= slev; )
    {
      if ( reql[i] <= stack[j] )
      {
        rs -= stack[j];
        if ( stack[j] == reql[i] )
        {
          rs2 -= stack[j];
          i++;
        }
        if ( rs < rs2 )
        {
          j = -1;
          break;
        }
        j++;
      }
      else
      {
        j = -1;
        break;
      }
    }

    if ( j == -1 )
    {
      sp[slev] = NULL;
    }

    while ( sp[slev] != NULL )
    {
      for ( i = 0; i < 26; i++ )
      {
        tn[i] = nlets[i];
      }

      for ( i = 0; i < stack[slev]; i++ )
      {
        if ( --tn[( (int) sp[slev]->word[i] ) - 'a'] < 0 )
        {
          break;
        }
      }

      if ( i == stack[slev] )
      {
        break;
      }
      else
      {
        sp[slev] = sp[slev]->gl;
      }

    }

    if ( sp[slev] == NULL )
    {
      if ( slev != 0 )
      {
        if ( stack[slev] < stack[slev - 1] && stack[slev] < cl )
        {
          stack[slev]++;
          if ( stack[slev] != stack[slev - 1] )
          {
            sp[slev] = gptr[stack[slev]];
          }
          else
          {
            sp[slev] = sp[slev - 1];
          }
        }
        else
        {
          stack[slev] = 0;
          sp[slev] = NULL;
          slev--;
          for ( i = 0; i < stack[slev]; i++ )
          {
            nlets[( (int) sp[slev]->word[i] ) - 'a']++;
          }
          sp[slev] = sp[slev]->gl;
          cl = cl + stack[slev];
        }
      }
      else
      {
        if ( stack[0] < tl )
        {
          stack[0]++;
          sp[0] = gptr[stack[0]];
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      if ( stack[slev] < cl )
      {
        for ( i = 0; i < 26; i++ )
        {
          nlets[i] = tn[i];
        }
        cl -= stack[slev];
        slev++;
        stack[slev] = 1;
        sp[slev] = gptr[1];
      }
      else
      {
        cp = tl + slev + 3;

        if ( sw != NULL )
        {
          cp = cp + 1 + strlen( sw );
        }

        if ( 79 < co + cp )
        {
          co = cp;
          cout << "\n";
        }
        else
        {
          co = co + cp;
        }

        cout << "(";

        if ( sw != NULL )
        {
          cout << sw << " ";
        }

        for ( i = 0; i <= slev; i++ )
        {
          cout << sp[i]->word;
          if ( i != slev )
          {
            cout << " ";
          }
          else
          {
            cout << ") ";
          }
        }
        na++;
        sp[slev] = sp[slev]->gl;
      }
    }
  }

  cout << "\n";

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "  Number of anagrams found was " << na << "\n";
  }
  return;
}
//****************************************************************************80

void find_valid_words ( )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_VALID_WORDS ??
//
{
  int c;
  long checked;
  dent *dptr;
  long found;
  int i;
  int j;
  int k;
  long perc;
  int tn[26];
  dent *tp[40];
  bool VERBOSE = false;

  checked = 0;
  found = 0;

  for ( i = 0; i < 40; i++ )
  {
    gptr[i] = NULL;
  }

  for ( i = 0; i < 26; i++ )
  {
    for ( j = 1; j < 40; j++ )
    {
      dptr = sptr[i][j];
      while ( dptr != NULL )
      {
        for ( k = 0; k < 26; k++ )
        {
          tn[k] = nlets[k];
        }

        for ( k = 0; k < j; k++ )
        {
          c = dptr->word[k] - 'a';
          if ( --tn[c] < 0 )
          {
            break;
          }
        }

        if ( k == j )
        {
          if ( gptr[j] == NULL )
          {
            gptr[j] = dptr;
            tp[j] = dptr;
          }
          else
          {
            tp[j]->gl = dptr;
            tp[j] = dptr;
          }
          dptr->gl = NULL;
          found++;
        }
        dptr = dptr->len;
        checked++;
      }
    }
  }

  perc = 1000 - ( ( found * 10000 + 5 ) / 10 ) / checked;

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "FIND_VALID_WORDS:\n";
    cout << "  Eliminated " << perc / 10 << " percent of words\n";
  }
  return;
}
//****************************************************************************80

int get_args ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_ARGS prompts for the arguments if they weren't on the command line.
//
{
  int i;
  char inp;
  int j;

  cout << "\n";
  cout << "Type in a list of letters, optionally followed by word lengths,\n";
  cout << "optionally followed by word suggestions.  A blank line exits.\n";
  cout << "-> ";

  i = 0;
  j = 0;

  do
  {
    scanf( "%c", &inp );
    if ( ' ' < inp )
    {
      if ( j == 0 )
      {
        i++;
      }
      if ( 'A' <= inp && inp <= 'Z' )
      {
        inp = inp + 'a' - 'A';
      }
      my_argv[i][j++] = inp;
    }
    else
    {
      if ( 0 < i )
      {
        my_argv[i][j] = '\0';
      }
      j = 0;
    }
  } while ( inp != '\n' );

  return ( i - 1 );
}
//****************************************************************************80

void initialize ( )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes certain data.
//
{
  int i;
  int j;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "INITIALIZE:\n";
    cout << "  Initializing data...\n";
  }

  for ( i = 0; i < 26; i++ )
  {
    aptr[i] = NULL;
    for ( j = 1; j < 40; j++ )
    {
      sptr[i][j] = NULL;
      cptr[i][j] = NULL;
    }
  }
  return;
}
//****************************************************************************80

int parse_args ( int nargs )

//****************************************************************************80
//
//  Purpose:
//
//    PARSE_ARGS parses the command line arguments.
//
//  Parameters:
//
//    Output, int parse_args, the number of command line arguments.
{
  int c;
  int i;
  int j;
  int t;
//
//  Replace the arguments by a simple count of each character.
//
  for ( i = 0; i < 26; i++ )
  {
    nlets[i] = 0;
  }

  tl = strlen ( my_argv[1] );

  for ( i = 0; i < tl; i++ )
  {
    c = my_argv[1][i];
    if ( c < 'a' || 'z' < c )
    {
      break;
    }
    nlets[c - 'a']++;
  }

  if ( i != tl )
  {
    cout << "\n";
    cout << "PARSE_ARGS - Fatal error!\n";
    cout << "  Non-alphabetic characters in the input set.\n";
    exit ( 1 );
  }

  i = 2;
  rl = 0;
//
//  The optional numbers field forces the program to use words of specific lengths.
//
  while ( i < nargs )
  {
    if ( my_argv[i][0] < '1' || '9' < my_argv[i][0] )
    {
      break;
    }
    reql[rl++] = atoi( my_argv[i++] );
  }
//
//  Descending sort the numbers list.
//
  for ( j = rl - 1; 0 < j; j-- )
  {
    for ( c = 0; c < j; c++ )
    {
      if ( reql[c] < reql[c + 1] )
      {
        t = reql[c];
        reql[c] = reql[c + 1];
        reql[c + 1] = t;
      }
    }
  }
//
//  Sum the numbers.
//
  for ( j = 0, trl = 0; j < rl; j++ )
  {
    trl = trl + reql[j];
  }
  return ( i );
}
//****************************************************************************80

void read_dict ( )

//****************************************************************************80
//
//  Purpose:
//
//    READ_DICT reads the word list from the dictionary file.
//
//  Modified:
//
//    14 January 2000
//
{
  int ai;
  ifstream file_in;
  char *file_in_name = "anagram_dictionary.txt";
  int i;
  int le;
  dent *tdent;
  char *tword;
  bool VERBOSE = false;
  char word[40];
  dent *wptr;
  int wi;
//
//  Open the dictionary file.
//
  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "READ_DICT - Fatal error!\n";
    cout << "  Could not find the dictionary file:\n";
    cout << file_in_name << "\n";
    exit ( 1 );
  }

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "READ_DICT:\n";
    cout << "  Reading the dictionary file...\n";
  }

  wptr = NULL;
  ai = 0;

  while ( 1 )
  {

    file_in.getline ( word, sizeof ( word ) );

    if ( file_in.eof ( ) )
    {
      break;
    }

    le = strlen ( word );
//
//  Move any uppercase characters to lowercase.
//
    for ( i = 0; i < le; i++ )
    {
      if ( 'A' <= word[i] && word[i] <= 'Z' )
      {
        word[i] = word[i] + ( 'a' - 'A' );
      }
      if ( word[i] < 'a' || 'z' < word[i] )
      {
        break;
      }
    }

    if ( i != le )
    {
      continue;
    }

    tdent = ( dent * ) malloc( DE_SIZE );
    tword = ( char * ) malloc( le + 1 );

    if ( tdent == NULL || tword == NULL )
    {
      cout << "\n";
      cout << "READ_DICT - Fatal error!\n";
      cout << "  malloc() failed\n";
      exit ( 1 );
    }

    wi = word[0] - 'a';
    strcpy( tword, word );
    tdent->word = tword;
    tdent->le = le;
    tdent->alph = NULL;
    tdent->len = NULL;

    if ( wptr == NULL )
    {
      aptr[ai++] = tdent;
    }
    else
    {
      if ( wptr->word[0] != word[0] )
      {
        aptr[ai++] = tdent;
      }
      wptr->alph = tdent;
    }

    if ( cptr[wi][le] == NULL )
    {
      sptr[wi][le] = tdent;
      cptr[wi][le] = tdent;
    }
    else
    {
      cptr[wi][le]->len = tdent;
      cptr[wi][le] = tdent;
    }

    wptr = tdent;
  }
  file_in.close ( );
  return;
}
//****************************************************************************80

void search ( int suggp, int sugg )

//****************************************************************************80
//
//  Purpose:
//
//    SEARCH ??
//
{
  int c;
  int i;
  int j;
  int tn[26];

  if ( suggp == sugg )
  {
    do_search ( "", NULL );
    return;
  }

  for ( i = suggp; i < sugg; i++ )
  {
    cout << "Target letters: " << my_argv[i] << "\n";

    for ( j = 0; j < 26; j++ )
    {
      tn[j] = nlets[j];
    }

    for ( j = 0; j < strlen( my_argv[i] ); j++ )
    {
      c = my_argv[i][j];
      if ( c < 'a' || 'z' < c )
      {
        break;
      }
      if ( --tn[c - 'a'] < 0 )
      {
        break;
      }
    }

    if ( j == strlen ( my_argv[i] ) )
    {
      for ( j = 0; j < 26; j++ )
      {
        nlets[j] = tn[j];
      }
      tl -= strlen( my_argv[i] );
      do_search( "*** ", my_argv[i] );
      for ( j = 0; j < strlen( my_argv[i] ); j++ )
      {
        nlets[( (int) my_argv[i][j] ) - 'a']++;
      }
      tl = tl + strlen( my_argv[i] );
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
//    May 31 2001 09:45:54 AM
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
