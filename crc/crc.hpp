//  crc.h

unsigned long crc ( unsigned char *buf, int len );
void make_crc_table ( void );
void print_crc_table ( void );
unsigned long update_crc_c ( unsigned long crc, unsigned char c );
unsigned long update_crc_s ( unsigned long crc, unsigned char *buf, int len );
void timestamp ( void );

