typedef struct 
{
  int owner;
  int count;
} omp_nest_lock_t;

typedef int omp_lock_t;

enum { UNLOCKED = -1, INIT, LOCKED };

enum { NOOWNER = -1, MASTER = 0 };

void omp_destroy_lock ( omp_lock_t *lock );
void omp_destroy_nest_lock ( omp_nest_lock_t *nlock );
int omp_get_dynamic ( );
int omp_get_max_threads ( );
int omp_get_nested ( );
int omp_get_num_procs ( );
int omp_get_num_threads ( );
int omp_get_thread_num ( );
double omp_get_wtick ( );
double omp_get_wtime ( );
int omp_in_parallel ( );
void omp_init_lock ( omp_lock_t *lock );
void omp_init_nest_lock ( omp_nest_lock_t *nlock );
void omp_set_dynamic ( int dynamic_threads );
void omp_set_lock ( omp_lock_t *lock );
void omp_set_nest_lock ( omp_nest_lock_t *nlock );
void omp_set_nested ( int nested );
void omp_set_num_threads ( int num_threads );
int omp_test_lock ( omp_lock_t *lock );
int omp_test_nest_lock ( omp_nest_lock_t *nlock );
void omp_unset_lock ( omp_lock_t *lock );
void omp_unset_nest_lock ( omp_nest_lock_t *nlock );


