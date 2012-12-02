namespace MPI 
{
  void Finalize (  );
  void Init ( int &argc, char **&argv );
  void Init ( );
  double Wtick ( );
  double Wtime ( );
  class Comm 
  {
    public:
      int Get_size (  ) const;
      int Get_rank ( ) const;
  };
  extern Comm COMM_WORLD;

  const int SOURCE = 1;
  const int TAG = 2;
  const int COUNT = 3;

  const int ANY_SOURCE = -1;
  const int ANY_TAG = -1;

  const int INT = 1;
  const int FLOAT = 2;
  const int DOUBLE = 3;
  const int BYTE = 4;

  const int SUM = 1;
  const int MAX = 2;
  const int MIN = 3;
  const int PRODUCT = 4;
};

