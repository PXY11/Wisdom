# include <cmath>
# include <cstdlib>
# include <ctime>
# include <stdint.h>
# include <iomanip>
# include <iostream>

//  Ziggurat algorithm of Marsaglia and Tsang
float r4_uni(uint32_t &jsr);
uint32_t shr3_seeded(uint32_t &jsr);
float r4_nor(uint32_t &jsr, uint32_t kn[128], float fn[128], float wn[128]);
void r4_nor_setup(uint32_t kn[128], float fn[128], float wn[128]);
