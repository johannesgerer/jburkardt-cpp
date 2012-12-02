# include <cstdlib>
# include <iostream>
# include <cstring>

using namespace std;

int main ( int argc, char **argv );
void badchk ( char *buf );
void blankfix ( register char *s );
void checkchar ( register char *s );
void checkline ( FILE *f );
void colout ( char *head, long nz, int what );
void early_eof ( void );
char *exform ( char *s0, char **Z );
int exindx ( char **s );
void namfetch ( int, char s[8] );
void namstore ( int, char s[8] );
FILE *newfile ( int n );
void newofile ( char *buf );
void process ( FILE *, char *infile1 );
char *rdline ( char s[77] );
void scream ( char *, char * );
void usage ( char **, int rc );
//
//  Define trtab to make this source self-contained...
//
char trtab[] = "!\"#$%&'()*+,-./0123456789;<=>?@\
ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^_`abcdefghijklmnopqrstuvwxyz{|}~";

char *bigstore;
int blanksubst;
unsigned BSsize;
int canend;
char chkbuf[76];
int cn;
long fline[72];
char *fname[72];
FILE *inf;
char *infile;
char invtrtab[256];
int just1;
int keepmyst = 1;
int kmax;
char *lastl;
int ncs;
char *progname;
int sflag;
long nline;
long nrow;
char *ss;
char ***xargv;

# define tr(x) Tr[x]

//****************************************************************************80

int main ( int argc, char **argv )

//****************************************************************************80
//
//  Purpose:
//
//    EMPS expands a compressed linear programming file to MPS format.
//
//  Discussion:
//
//    Expand compressed LP programs (in netlib format) to MPS format.
//    This is similar to the Fortran program emps.f , except that it
//    understands command-line arguments, including the -m option,
//    which causes "mystery lines" to be included in the output.
//    ("Mystery lines" allow extensions to MPS format.  At the time of
//    this writing, however, none of the netlib problems contain
//    mystery lines.)
//
//  Author:
//
//    David Gay
//    AT&T Bell Laboratories.
//
{
  char *s;
  char *se;
  FILE *f;
  static char *options[] =
  {
    "-1  {output only 1 nonzero per line}",
    "-b  {replace blanks within names by _'s}",
    "-m  {skip mystery lines}",
    "-S  {split output by problems: put each problem in the file",
    "\tnamed by the first word after \"NAME\" on the NAME line}",
    "-s  {like -S, but force file names to lower case}",
    0
  };

  for ( s = invtrtab, se = s + sizeof(invtrtab); s < se; s++ )
  {
    *s = 92;
  }

  for (s = se = trtab; *s; s++)
  {
    invtrtab[*s] = s - se;
  }

  *chkbuf = ' ';

  progname = *argv;

  xargv = &argv;

  while ( s = *++argv )
  {
    if (*s == '-') switch(s[1])
    {
      case 0: process(stdin, "<stdin>"); break;

      case '1': just1 = 1; break;

      case 'b': blanksubst = '_'; break;

      case 'm': keepmyst = 0; break;

      case 'S': sflag = 1; break;
      case 's': sflag = 2; break;

      case '?': usage(options,0);

      default:
        cerr << progname << ": invalid option " << s << "\n";
        usage ( options, 1 );
    }
    else
    {
      f = fopen ( s, "r" );
      if ( !f )
      {
        cerr << progname << ": can't open " << s << "\n";
        exit ( 1 );
      }
      process(f, s);
    }
  }

  if ( !infile )
  {
    process ( stdin, "<stdin>" );
  }

  return 0;
}
//****************************************************************************80

void badchk ( char *buf )

//****************************************************************************80
//
//  Purpose:
//
//    BADCHK handles the case when a bad checksum is encountered.
//
{
  int i;
  static char csl[] = "Check sum line =";
  char *mb = csl, msgbuf[64];

  cerr << progname << ": Check sum error: expected\n";
  cerr << chkbuf << "\n";
  cerr << "but got:\n";
  cerr << buf << "\n";

  lastl = buf;

  if ( *buf == ' ' )
  {
    for(i = 1; chkbuf[i] == buf[i]; i++);
    sprintf(msgbuf, "Bad check sum for line %ld of %s\n%%s",fline[i], fname[i]);
    mb = msgbuf;
  }

  scream ( mb, csl );
}
//****************************************************************************80

void blankfix ( register char *s )

//****************************************************************************80
//
//  Purpose:
//
//    BLANKFIX replaces blank characters by a substitute character.
//
{
  for ( ; *s; s++ )
  {
    if ( *s == ' ' )
    {
      *s = blanksubst;
    }
  }
}
//****************************************************************************80

void checkchar ( register char *s )

//****************************************************************************80
//
//  Purpose:
//
//    CHECKCHAR ???
{
  register int c;
  register unsigned x;
  register char *Tr = invtrtab;

  for ( x = 0; c = *s; s++ )
  {

    if ( c == '\n' )
    {
      *s = 0;
      break;
    }

    c = tr ( c );

    if ( x & 1 )
    {
      x = ( x >> 1 ) + c + 16384;
    }
    else
    {
      x = ( x >> 1 ) + c;
    }

  }

  fname[ncs] = infile;
  fline[ncs] = nline;
  chkbuf[ncs++] = trtab[x % 92];
}
//****************************************************************************80

void checkline ( FILE *f )

//****************************************************************************80
//
//  Purpose:
//
//    CHECKLINE ???
//
{
  char chklin[76];

  canend = 0;

 again:
  chkbuf[ncs++] = '\n';
  chkbuf[ncs] = 0;
  nline++;

  while ( !fgets ( chklin, 76, f ) )
  {
    if ( !( f = newfile(1) ) )
    {
      early_eof();
    }
  }

  if ( strcmp ( chklin,chkbuf ) )
  {

    if ( *chklin == ':' && ncs <= 72 )
    {
      ncs--;
      checkchar(chklin);
      if (keepmyst)
        cout << chklin+1 << "\n";
      goto again;
    }

    badchk(chklin);

  }

  ncs = 1;
}
//****************************************************************************80

void colout ( char *head, long nz, int what )

//****************************************************************************80
//
//  Purpose:
//
//    COLOUT ???
//
{
  static char *bt[] = {"UP", "LO", "FX", "FR", "MI", "PL"},
    fmt2[] = "    %-8.8s  %-8.8s  %-15.15s%-8.8s  %.15s\n";
  char buf[80];
  char curcol[8];
  char msgbuf[32];
  char *rc1;
  char *rc2;
  char rcbuf1[16], rcbuf2[16], rownm[2][8], *z;
  int first;
  int k;
  int n;

  if (!nz)
  {
    if (what <= 2)
    {
      cout << head << "\n";
    }
    return;
  }

  first = 1;
  k = 0;
  z = "";
  *curcol = 0;

  while ( nz-- )
  {
    if ( !*z )
    {
      z = rdline(buf);
    }

    if ( first )
    {
      cout << head << "\n";
      first = 0;
    }

    while ( !(n = exindx(&z)) )
    {
      if ( k )
      {
        printf("    %-8.8s  %-8.8s  %.15s\n",
          curcol, rownm[0], rc1);
        k = 0;
      }

      if ( blanksubst )
      {
        if ( *z )
        {
          blankfix ( z );
        }
        else
        {
          z = head;
        }
      }

      strncpy ( curcol, z, 8 );

      if ( what == 1 )
      {
        namstore ( ++cn, z );
      }

      z = rdline ( buf );

    }

    if ( 4 <= what )
    {

      if ( 7 <= n )
      {
        sprintf(msgbuf, "bad bound type index = %d",n);
        scream(msgbuf, "");
      }

      if ( !*z )
      {
        z = rdline(buf);
      }
      namfetch ( (int)nrow + exindx(&z), rownm[0] );
      if ( 4 <= n-- )
      {
        printf(" %s %-8.8s  %.8s\n", bt[n], curcol, *rownm);
        continue;
      }
    }
    else namfetch ( n, rownm[k] );
    if ( !*z )
    {
      z = rdline(buf);
    }
    if (k) rc2 = exform(rcbuf2, &z);
    else rc1 = exform(rcbuf1, &z);
    if (what <= 3)
    {
      if ( just1 )
        printf("    %-8.8s  %-8.8s  %.15s\n",
          curcol, rownm[0], rc1);
      else
      {
        if (++k == 1) continue;
        printf(fmt2, curcol, rownm[0], rc1,
          rownm[1], rc2);
        k = 0;
      }
    }
    else
    {
      printf(" %s %-8.8s  %-8.8s  %.15s\n", bt[n], curcol,
        rownm[0], rc1);
    }
  }

  if ( k )
  {
    printf("    %-8.8s  %-8.8s  %.15s\n", curcol, *rownm, rc1);
  }

}
//****************************************************************************80

void early_eof ( void )

//****************************************************************************80
//
//  Purpose:
//
//    EARLY_EOF handles the situation when the input file ends prematurely.
//
{
  lastl = "";
  scream ( "Premature end of file", "" );
}
//****************************************************************************80

char *exform ( char *s0, char **Z )

//****************************************************************************80
//
//  Purpose:
//
//    EXFORM expands a compressed string.
//
{
  char *d;
  char db[32];
  int ex;
  int k;
  int nd;
  int nelim;
  char *s;
  char sbuf[32];
  long x, y = 0;
  register char *Tr = invtrtab;
  register char *z = *Z;

  d = db;
  k = tr(*z++);
//
//  Supersparse index.
//
  if ( k < 46 )
  {
    k = exindx ( Z );
    if ( kmax < k )
    {
      char msgbuf[64];
      sprintf(msgbuf, "index %u > kmax = %u in %%s", k, kmax);
      scream ( msgbuf, z-1 );
    }
    return ss + (k << 4);
  }

  s = sbuf;
  k -= 46;

  if ( 23 <= k )
  {
    *s++ = '-';
    k -= 23;
    nelim = 11;
  }
  else
  {
    nelim = 12;
  }
//
//  Integer floating point.
//
  if ( 11 <= k )
  {
    k -= 11;
    *d++ = '.';
    if ( k >= 6 )
    {
      x = k - 6;
    }
    else
    {
      x = k;
      for ( ; ; )
      {
        k = tr(*z++);
        x = x*46 + k;
        if ( k >= 46 )
        {
          x -= 46;
          break;
        }
      }
    }

    if (!x) *d++ = '0';
    else do
    {
      *d++ = '0' + x%10;
      x /= 10;
    } while(x);

    do *s++ = *--d; while(d > db);
    }
//
//  General floating point.
//
  else
  {
    ex = (int)tr(*z++) - 50;
    x = tr(*z++);
    while(--k >= 0)
    {
      if (x >= 100000000)
      {
        y = x;
        x = tr(*z++);
      }
      else x = x*92 + tr(*z++);
    }
    if (y)
    {
      while(x > 1) { *d++ = x%10 + '0'; x /= 10; }
      for(;; y /= 10)
      {
        *d++ = y%10 + '0';
        if (y < 10) break;
      }
    }
    else if (x) for(;; x /= 10)
    {
      *d++ = x%10 + '0';
      if (x < 10) break;
    }
    else *d++ = '0';
    nd = d - db + ex;
    if ( ex > 0 )
    {
      if (nd < nelim || ex < 3)
      {
        while(d > db) *s++ = *--d;
        do *s++ = '0'; while(--ex);
        *s++ = '.';
      }
      else goto Eout;
      }
    else if (nd >= 0)
    {
      while(--nd >= 0) *s++ = *--d;
      *s++ = '.';
      while(d > db) *s++ = *--d;
      }
    else if (ex > -nelim)
    {
      *s++ = '.';
      while(++nd <= 0) *s++ = '0';
      while(d > db) *s++ = *--d;
      }
    else {
Eout:
      ex += d - db - 1;
      if (ex == -10) ex = -9;
      else {
        if (ex > 9 && ex <= d - db + 8) {
          do { *s++ = *--d;
            } while (--ex > 9);
          }
        *s++ = *--d;
        }
      *s++ = '.';
      while(d > db) *s++ = *--d;
      *s++ = 'E';
      if (ex < 0) { *s++ = '-'; ex = -ex; }
      while(ex) { *d++ = '0' + ex%10; ex /= 10; }
      while(d > db) *s++ = *--d;
      }
    }
  *s = 0;
  k = s - sbuf;
  s = s0;
  while(k++ < 12) *s++ = ' ';
  strcpy(s, sbuf);
  *Z = z;

  return s0;
}
//****************************************************************************80

int exindx ( char **s )

//****************************************************************************80
//
//  Purpose:
//
//    EXINDX expands a supersparse index.
//
{
  register char *Tr = invtrtab;
  register char *z = *s;
  register int k;
  register int x;

  k = tr(*z++);

  if (k >= 46)
  {
    scream ( "exindx: Bad index in %s", z );
  }

  if ( k >= 23 )
  {
    x = k - 23;
  }
  else
  {
    x = k;
    for ( ; ; )
    {
      k = tr(*z++);
      x = x*46 + k;
      if ( 46 <= k )
      {
        x -= 46;
        break;
      }
    }
  }

  *s = z;

  return x;
}
//****************************************************************************80

void namfetch ( int i, char s[8] )

//****************************************************************************80
//
//  Purpose:
//
//    namfetch ???
//
//  The following routine extends the
//  size of problems that the small-memory MS-DOS version of emps can
//  handle.  If you have a compiler that makes "huge" pointers available
//  and can arrange for bigstore to be a huge pointer (one that can
//  address a region larger than 64 kilobytes), then you can use
//  suitably modified versions of the namfetch given below.
//  (If you are using the large memory model, then this only matters
//  for the larger problems, those for which the number of rows plus
//  the number of columns is more than 8191.)
//
{
  if ( i <= 0 || i > BSsize )
  {
    scream ( "Bad i to namfetch", "" );
  }
  strncpy ( s, bigstore + (i<<3), 8 );
}
//****************************************************************************80

void namstore ( int i, char s[8] )

//****************************************************************************80
//
//  Purpose:
//
//    namstore ???
//
//  The following routine extends the
//  size of problems that the small-memory MS-DOS version of emps can
//  handle.  If you have a compiler that makes "huge" pointers available
//  and can arrange for bigstore to be a huge pointer (one that can
//  address a region larger than 64 kilobytes), then you can use
//  suitably modified versions of the namstore given below.
//  (If you are using the large memory model, then this only matters
//  for the larger problems, those for which the number of rows plus
//  the number of columns is more than 8191.)
//
{
  if ( i <= 0 || i > BSsize )
  {
    scream ( "Bad i to namstore", "" );
  }

  strncpy ( bigstore + (i<<3), s, 8 );
}
//****************************************************************************80

FILE *newfile ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    newfile ???
{
  char **av;
  char *s1;
  FILE *f = 0;

  av = *xargv;

  if ( *av && ( s1 = *++av ) && *s1 != '-' )
  {
    fclose ( inf );
    f = fopen(s1, "r");
    if ( !f )
    {
      cerr << progname << ": can't open " << s1 << "\n";
      exit ( 1 );
    }
    inf = f;
    infile = s1;
    *xargv = av;
    nline = n;
  }
  return f;
}
//****************************************************************************80

void newofile ( char *buf )

//****************************************************************************80
//
//  Purpose:
//
//    newofile ???
{
  register unsigned char *s;
  register unsigned char *t;
  register int c;
  char namebuf[80];

  for ( s = (unsigned char *)buf + 4; *s <= ' '; s++)

    if ( !*s )
    {
      scream ( "Blank NAME line", "" );
    }

  t = (unsigned char *)namebuf;

  if (sflag == 2)
  {
    while((c = *s++) > ' ')
      *t++ = c >= 'A' && c <= 'Z' ? c + 'a' - 'A' : c;
  }
  else
  {
    while((c = *s++) > ' ')
      *t++ = c;
  }

  *t = 0;

  if ( !freopen(namebuf, "w", stdout) )
  {
    scream ( "can't open \"%s\"", namebuf );
  }

}
//****************************************************************************80

void process ( FILE *f, char *infile1 )

//****************************************************************************80
//
//  Purpose:
//
//    process ???
{
  char *b1;
  char buf[80];

  long ncol;
  long colmx;
  long nz, nrhs, rhsnz, nran, ranz, nbd, bdnz, ns;
  int i;
  char *s;
  char *ss0;
  char *z;

  infile = infile1;
  inf = f;
  nline = 0;
  canend = 0;
  ncs = 1;
  rdline(buf);
top:
  kmax = -1;
  ncs = 1;
//
//  NAME line
//
  while ( strncmp(buf,"NAME",4) )
    if ( !rdline(buf) )
      goto done;
  canend = 0;

  if ( sflag )
  {
    newofile ( buf );
  }

  cout << buf << "\n";
  ncs = 1;
//
//  problem statistics.
//
  rdline(buf);
  if ( sscanf(buf,"%ld %ld %ld %ld %ld %ld %ld %ld", &nrow, &ncol,
    &colmx, &nz, &nrhs, &rhsnz, &nran, &ranz) != 8 ||
    rdline(buf), sscanf(buf, "%ld %ld %ld", &nbd, &bdnz, &ns) != 3 )
  {
      scream ( "Bad statistics line:\n%s\n", buf );
  }

  ncs = 1;
  cn = nrow;
  i = cn + ncol;

  if (i != nrow + ncol)
  {
    scream ( "Problem too big", "" );
  }
//
//  read, expand number table
//
  BSsize = nrow + ncol;
  ss0 = (char *) malloc((unsigned)(ns<<4) + (BSsize<<3));
  bigstore = ss0 + (ns<<4) - 8;

  if ( !ss0 )
  {
    scream ( "malloc failure!", "" );
  }

  ss = ss0 - 16;
  z = "";
  for(s = ss0, i = ns; i--; s += 16)
  {
    if (!*z) z = rdline(buf);
    exform(s, &z);
  }

  kmax = ns;
//
//  read, print row names
//
  b1 = buf + 1;
  for ( i = 1; i <= nrow; i++)
  {
    rdline(buf);
    if ( i == 1 )
    {
      cout << "ROWS\n";
    }

    if ( blanksubst )
    {
      blankfix(b1);
    }

    cout << *buf << " " << b1 << "\n";

    namstore(i, b1);
  }
//
//  Read, print columns
//
  colout ( "COLUMNS", nz, 1 );
//
//  right-hand sides.
//
  colout ( "RHS", rhsnz, 2 );
//
//  Ranges.
//
  colout("RANGES", ranz, 3);
//
//  Bounds.
//
  colout("BOUNDS", bdnz, 4);
//
//  Final checksum line...
//
  if ( ncs > 1 )
  {
    checkline ( inf );
  }

  cout << "ENDATA\n";
//
//  See whether there's another LP in this file...
//
  free(ss0);
  canend = ncs = 1;

  if ( rdline( buf ) )
  {
    goto top;
  }

 done:
  fclose(inf);
}
//****************************************************************************80

char *rdline ( char s[77] )

//****************************************************************************80
//
//  Purpose:
//
//    rdline ???
{
  FILE *f = inf;

again:

  nline++;

  if (!fgets(s, 77, f))
  {
    if ( f = newfile(0) )
      goto again;
    if ( canend )
      return 0;
    early_eof();
  }

  checkchar(s);

  if ( ncs >= 72 )
  {
    checkline(f);
  }

  if ( *s == ':' )
  {
    if ( keepmyst )
    {
      cout << s+1 << "\n";
    }
    goto again;
  }

  return lastl = s;
}
//****************************************************************************80

void scream ( char *s, char *t )

//****************************************************************************80
//
//  Purpose:
//
//    SCREAM prints an error message and terminates execution.
//
{
  long c;

  cerr << progname << ": ";
  fprintf ( stderr, s, t);
  c = ftell(inf) - strlen(lastl);
  fprintf(stderr, ": line %ld (char %ld) of %s\n", nline, c, infile);

  exit ( 1 );
}
//****************************************************************************80

void usage ( char **o, int rc )

//****************************************************************************80
//
//  Purpose:
//
//    USAGE prints out usage information for the program.
//
{
  char *s;

  cerr << "Usage: " << progname << " [options] [file ...]\n";
  cerr << "Options:\n";

  while ( s = *o++ )
  {
    cerr << "  " << s << "\n";
  }

  exit ( rc );
}
