/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_xorshift128plus_
#define _sc_xorshift128plus_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_xorshift128plus)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_out< float > s1;
  sc_fifo_out< float > s2;

  SC_CTOR(sc_xorshift128plus)
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
#define TINYMT32_MUL (1.0f / 4294967296.0f)

void sc_xorshift128plus::do_gen()
{
	sc_uint<64> t0 = 0x1234567809876543;
	sc_uint<64> t1 = 0x12121212FCBADEFF;
	wait();

	while( true )
	{
		#pragma HLS PIPELINE
		sc_uint<64> x = t0;
		const sc_uint<64> y = t1;
		t0 = y;
		x ^= x << 23;
		t1 = x ^ y ^ (x >> 17) ^ (y >> 26);

		sc_uint<64> w = t1 + y;

		sc_uint<32> rand1 = w.range(31,  0);
		sc_uint<32> rand2 = w.range(63, 32);

		float f1 = 1.0f - (float)rand1 * (float)TINYMT32_MUL;
		float f2 = 1.0f - (float)rand2 * (float)TINYMT32_MUL;

		s1.write( f1 );
		s2.write( f2 );
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
