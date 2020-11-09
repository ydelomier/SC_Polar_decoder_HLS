/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_mtwister_
#define _sc_mtwister_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_mtwister)          	// module (class) declaration
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_out< float > s; 		// vers de decodeur

  SC_CTOR(sc_mtwister)
  {
    SC_CTHREAD(do_gen, clk.pos());
    reset_signal_is(reset, true);
  }

private:
    void do_gen();
};
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#define TINYMT32_MEXP 127
#define TINYMT32_SH0 1
#define TINYMT32_SH1 10
#define TINYMT32_SH8 8
#define TINYMT32_MASK ((unsigned int)(0x7fffffff))
#define TINYMT32_MUL (1.0f / 4294967296.0f)

struct TINYMT32_T {
    unsigned int status[4];
    unsigned int mat1;
    unsigned int mat2;
    unsigned int tmat;
};

typedef struct TINYMT32_T tinymt32_t;
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
inline void tinymt32_next_state(tinymt32_t * random) {
	unsigned int x;
	unsigned int y;

    y = random->status[3];
    x = (random->status[0] & TINYMT32_MASK) ^ random->status[1] ^ random->status[2];
    x ^= (x << TINYMT32_SH0);
    y ^= (y >> TINYMT32_SH0) ^ x;
    random->status[0] = random->status[1];
    random->status[1] = random->status[2];
    random->status[2] = x ^ (y << TINYMT32_SH1);
    random->status[3] = y;
    random->status[1] ^= -((unsigned int)(y & 1)) & random->mat1;
    random->status[2] ^= -((unsigned int)(y & 1)) & random->mat2;
//    printf("next = %X\n", random->status[0]);
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
inline unsigned int tinymt32_temper(tinymt32_t * random) {
	unsigned int t0, t1;
    t0  = random->status[3];
    t1  = random->status[0] ^ (random->status[2] >> TINYMT32_SH8);
    t0 ^= t1;
    t0 ^= -((unsigned int)(t1 & 1)) & random->tmat;
//    printf("integer = %d\n", t0);
    return t0;
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#define MIN_LOOP 8
#define PRE_LOOP 8

static int offset = 0;

void tinymt32_init(tinymt32_t* random) {
	random->status[0] = 'T' + offset;
	random->status[1] = 'I' + offset;
	random->status[2] = 'N' + offset;
	random->status[3] = 'Y' + offset;
	random->mat1      = 0x8f7011ee;
	random->mat1      = 0xfc78ff1f;
	random->tmat      = 0x3793fdff;
	offset += 1;

	for (int i = 0; i < PRE_LOOP; i++) {
		tinymt32_next_state(random);
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
void sc_mtwister::do_gen()
{
	tinymt32_t      random;
	tinymt32_init( &random );
	while( true )
	{
		tinymt32_next_state(&random);
		float f = tinymt32_temper(&random) * (float)TINYMT32_MUL;
//    printf("float = %f\n", f);
		s.write( f );
//    if( cpt++ == 8 ) exit( 0 );
	}
}

#endif
