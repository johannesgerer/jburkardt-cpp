uint32_t cong_seeded ( uint32_t &jcong );
uint32_t cong_value ( );
double cpu_time ( );
float efix ( );
uint32_t kiss_seeded ( uint32_t &jcong, uint32_t &jsr, uint32_t &w, uint32_t &z );
uint32_t kiss_value ( );
uint32_t mwc_seeded ( uint32_t &w, uint32_t &z );
uint32_t mwc_value ( );
float nfix ( );
void r4_exp_setup ( );
float r4_exp_value ( );
void r4_nor_setup ( );
float r4_nor_value ( );
float r4_uni_value ( );
uint32_t shr3_seeded ( uint32_t &jsr );
uint32_t shr3_value ( );
void timestamp ( );
void zigget ( uint32_t &jsr_value, uint32_t &jcong_value,
  uint32_t &w_value, uint32_t &z_value );
void zigset ( uint32_t jsr_value, uint32_t jcong_value,
  uint32_t w_value, uint32_t z_value );

