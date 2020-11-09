/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_xorshift1024star_
#define _sc_xorshift1024star_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_xorshift1024star)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_out< float > s1;
  sc_fifo_out< float > s2;

  SC_CTOR(sc_xorshift1024star)
  {
    SC_CTHREAD(do_gen, clk.pos());
    reset_signal_is(reset, true);
  }

private:
    void do_gen();

};


#define TINYMT32_MUL (1.0f / 4294967296.0f)

void sc_xorshift1024star::do_gen()
{
	int p = 0;
	const sc_uint<64> valu = 1181783497276652981;
	sc_uint<64> s[16] = { 0x84242f96eca9c41d,
			0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
			0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
			0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
			0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
			0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
		};
	wait();

	while( true )
	{
		#pragma HLS PIPELINE
		const sc_uint<64> t0 = s[p];
		sc_uint<64> t1 = s[p = (p + 1) & 15];
		t1 ^= t1 << 31; // a
		s[p] = t1 ^ t0 ^ (t1 >> 11) ^ (t0 >> 30);

		const sc_uint<64> res = s[p] * valu;

		sc_uint<32> rand1 = res.range(31,  0);
		sc_uint<32> rand2 = res.range(63, 32);

		float f1 = 1.0f - (float)rand1 * (float)TINYMT32_MUL;
		float f2 = 1.0f - (float)rand2 * (float)TINYMT32_MUL;

		s1.write( f1 );
		s2.write( f2 );
	}
}

#undef TINYMT32_MUL

#endif
