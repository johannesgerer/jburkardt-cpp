// This file is intended to write a bitmap image file corresponding to
// the output produced from the heatedplate.c or heatedplate.f90
// program files.
//
// The BmpImage class was written by Daniel LePage, 2007
// The main method was written by Aaron Bloomfield, 2008

# include <iostream>
# include <vector>
# include <fstream>
# include <stdio.h>

using namespace std;
using std::vector;

class Color {
 public:
  Color(double red, double green, double blue) 
    {
      r = (char)(red*255);
      g = (char)(green*255);
      b = (char)(blue*255);
    };
  char r, g, b;
};


class BmpImage {
 public:
  // size of this image
  const int width;
  const int height;

  // Create a new bitmap of size w x h, with a white background
  BmpImage(int w, int h);

  // Set the pixel at (x,y) to the color (r,g,b). r, g, and b should be
  // doubles between 0 and 1.
  void putpixel(int w, int h, double r, double g, double b);

  void writeToFile ( const char* filename );
  void writeToFile ( string filename );

 protected:
  enum endianness {BIG,LITTLE,UNKNOWN};
  static endianness endian;

  static bool isPlatformLittleEndian();

  void putpixel(int w, int h, Color c);
  vector<vector<Color > > data;

  void writeInt (int value, ofstream &f);
  void clear(Color c);

  void writeHeader(ofstream &f);
};


/** Public methods ***/

BmpImage::endianness BmpImage::endian = BmpImage::UNKNOWN;

BmpImage::BmpImage(int w, int h): width(w), height(h) {
  clear(Color(1,1,1));
}

void BmpImage::putpixel(int x, int y, double r, double g, double b) {
  putpixel(x,y,Color(r,g,b));
}

void BmpImage::writeToFile(string filename) {
  writeToFile(filename.c_str());
}

void BmpImage::writeToFile(const char* filename) {
  ofstream outfile;
  outfile.open (filename, ios::out | ios::binary);
  writeHeader(outfile);
  
  int w = width, h = height; //laziness
  int real_width = w*3+4-(w*3)%4;
  int total_data_bytes = real_width*h;

  int index = 0;
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      outfile << data[i][j].b;
      outfile << data[i][j].g;
      outfile << data[i][j].r;
    }
    //pad row with zeros
    for (int i = 0; i < (4-(w*3)%4)%4; i++)
      outfile << char(0);
  }
  outfile.close();
}


/** Private methods ***/

void BmpImage::putpixel(int x, int y, Color c) {
  data[x][y] = c;
}

void BmpImage::clear(Color background) {
  data.clear();
  for (int i = 0; i < width; i++) {
    vector<Color> c(height,background);
    data.push_back(c);
  }
}



bool BmpImage::isPlatformLittleEndian() {
  if (endian != UNKNOWN)
    return endian == LITTLE;

  union {
    char charword[4];
    unsigned int intword;
  } check;

  check.charword[0] = 1;    check.charword[1] = 2;
  check.charword[2] = 3;    check.charword[3] = 4;

  if (check.intword == 0x01020304)
    endian = LITTLE;
  else
    endian = BIG;
  return endian == LITTLE;
} 
  
void BmpImage::writeInt (int value, ofstream& f) {
  union {
    int intvalue;
    struct{
      char a, b, c, d;
    } c;
  } e;
  e.intvalue = value;
  if (!isPlatformLittleEndian()) {
    f<< e.c.a << e.c.b << e.c.c << e.c.d;
  } else {
    f<< e.c.d << e.c.c << e.c.b << e.c.a;
  }
}

//****************************************************************************80

void BmpImage::writeHeader ( ofstream &f ) 

//****************************************************************************80
{
  f << "BM"; // Magic number
  writeInt(0,f); // File Size
  writeInt(0,f); // Reserved
  writeInt(54,f); // Data offset (=length of all headers)
  writeInt(40,f); // Length of info header
  writeInt(width,f); // width
  writeInt(height,f); // height
  f << char(1); // Color planes
  f << char(0); // Color planes
  f << char(24) << char(0); // Bits per Pixel
  writeInt(0,f); // No compression
  writeInt(0,f); // image size
  writeInt(2835,f); // Horizontal resolution = 72 dpi
  writeInt(2835,f); // Vertical resolution = 72 dpi
  writeInt(0,f); // Number of colors (automatic)
  writeInt(0,f); // "Important" colors (all)

  return;
}
//****************************************************************************80

int main ( int argc, char** argv ) 

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GRID_TO_BMP.
//
//  Discussion:
//
//    GRID_TO_BMP creates a Microsoft BMP color image file that represents
//    scalar data read from a text file.
//
//    The text file should contain the values of a quantity on an M by N
//    grid, stored as an M by N array, which we will call "U".
//
//    The first two records of the text file should contain the values of
//    M and N, respectively.
//
//    There should follow M records, each of length N, containing, in order,
//    the rows of U.
//
//  Usage:
//
//    grid_to_bmp input_file output_file
//
//    where "input_file" is the text file to be read, and "output_file"
//    is the name of the BMP file to be created.  It is customary to use
//    a file extension of ".bmp" for the BMP file.
//
//  Modified:
//
//    22 July 2008
//
{
  double b;
  int col;
  float d;
  float d_max;
  float d_min;
  double f;
  double g;
  double h;
  int height;
  int hi;
  BmpImage *image;
  FILE *input_file;
  FILE *output_file;
  double p;
  double q;
  double r;
  int row;
  double s;
  double t;
  double v;
  int width;

  cout << "\n";
  cout << "GRID_TO_BMP:\n";
  cout << "  C++ version\n";

  if ( argc != 3 ) 
  {
    cerr << "\n";
    cerr << "GRID_TO_BMP\n";
    cerr << "  USAGE: " << argv[0] << " <input-filename> <output-filename>\n" << endl;
    return 1;
  }
//
//  Open the input file;
//  Read HEIGHT and WIDTH.
//  Determine the maximum and minimum values.
//  Close the file.
//
  input_file = fopen ( argv[1], "r" );

  if ( input_file == NULL ) 
  {
    cerr << "Error: unable to open the input file." << endl;
    return 2;
  }

  if ( fscanf (input_file, "%d\n", &height ) != 1 ) 
  {
    cerr << "Error: parse error reading input file (image size)" << endl;
    return 3;
  }
  if ( fscanf (input_file, "%d\n", &width ) != 1 ) 
  {
    cerr << "Error: parse error reading input file (image size)" << endl;
    return 3;
  }

  for ( row = 0; row < height; row++ ) 
  {
    for ( col = 0; col < width; col++ ) 
    {
      if ( fscanf ( input_file, "%f ", &d ) != 1 ) 
      {
        cerr << "Error: parse error reading in row " << row << ", col " << col << endl;
        return 4;
      }
      if ( row == 0 && col == 0 ) 
      {
        d_max = d;
        d_min = d;
      }
      else
      {
        if ( d < d_min )
        {
          d_min = d;
        }
        if ( d_max < d )
        {
          d_max = d;
        }
      }      
    }
  }
  fclose ( input_file );

  cout << "\n";
  cout << "  Minimum data value is " << d_min << "\n";
  cout << "  Maximum data value is " << d_max << "\n";
//
//  Make sure we can open the output file.
//
  output_file = fopen ( argv[2], "w" );

  if ( output_file == NULL ) 
  {
    cerr << "Error: unable to open output file." << endl;
    return 2;
  }
  fclose ( output_file );
//
//  Create the image.
//
  image = new BmpImage ( width, height );

  cout << "\n";
  cout << "  Creating image of size " << height << "x" << width << "\n";
//
//  Reopen the file.
//  Read the data, scaling each item, converting to a hue, and then
//  converting the hue to an (R,G,B) color.
//
  input_file = fopen ( argv[1], "r" );

  if ( fscanf (input_file, "%d\n", &height ) != 1 ) 
  {
    cerr << "Error: parse error reading input file (image size)" << endl;
    return 3;
  }
  if ( fscanf (input_file, "%d\n", &width ) != 1 ) 
  {
    cerr << "Error: parse error reading input file (image size)" << endl;
    return 3;
  }
  for ( row = 0; row < height; row++ ) 
  {
    for ( col = 0; col < width; col++ ) 
    {
      if ( fscanf ( input_file, "%f ", &d ) != 1 ) 
      {
        cerr << "Error: parse error reading in row " << row << ", col " << col << endl;
        return 4;
      }
//
//  Scale so maximum D has H = 0, minimum D has H = 1.
//
      h = 240.0 * ( d_max - d ) / ( d_max - d_min );
      s = 1.0;
      v = 1.0;
// 
//  Convert HSV color to RGB.
//
      hi = ( ( int ) ( h / 60.0 ) ) % 6;
      f = h / 60.0 - hi;
      p = v * ( 1.0 - s);
      q = v * ( 1.0 - f * s );
      t = v * ( 1.0 - ( 1.0 - f ) * s );

      switch ( hi ) 
      {
        case 0: r=v; g=t; b=p; break;
        case 1: r=q; g=v; b=p; break;
        case 2: r=p; g=v; b=t; break;
        case 3: r=p; g=q; b=v; break;
        case 4: r=t; g=p; b=v; break;
        case 5: r=v; g=p; b=q; break;
        default: cout << "Error: converting HSV to RGB" << endl; return 5; break;
      }
//
//  We implicitly invert the row index here.
//
      image->putpixel ( col, height-row-1, r, g, b );
    }
  }

  fclose ( input_file );
//
//  Write the image to the output file.
//
  image->writeToFile ( argv[2] );

  cout << "\n";
  cout << "  BMP graphics information written to \"" << argv[2] << "\".\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "GRID_TO_BMP:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
