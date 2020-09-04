#ifndef __BLCA_DIGAMMA_H__
#define __BLCA_DIGAMMA_H__

#ifndef M_PIl
/** The constant Pi in high precision */
#define M_PIl 3.1415926535897932384626433832795029L
#endif
#ifndef M_GAMMAl
/** Euler's constant in high precision */
#define M_GAMMAl 0.5772156649015328606065120900824024L
#endif
#ifndef M_LN2l
/** the natural logarithm of 2 in high precision */
#define M_LN2l 0.6931471805599453094172321214581766L
#endif

#include <math.h>
#include "kncoe_table.h"

long double BLCA_digammal(long double x);

#endif