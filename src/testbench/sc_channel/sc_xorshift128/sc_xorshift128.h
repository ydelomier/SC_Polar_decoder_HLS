/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 * 29/10/2016
 *   DANS SA FORME ACTUELLE LE MODULE SYNTHETISE SUR NEXYS 4 AVEC LES PARAMETRES
 *   SUIVANTS POUR LA LOOP 1 (CONTRAINTE = 10ns + 2ns uncertaincy)
 *   - LATENCY  = 12
 *   - INTERCAL = 1
 *   COUT DE L'ARCHITECTURE
 *   - 864 LUTS
 *    - 719 FFs
 *   -   3 DSPs
 *
**/

#ifndef _sc_xorshift128_
#define _sc_xorshift128_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_xorshift128)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_in < sc_uint<8> > seed;

  sc_fifo_out< float > s1;
  sc_fifo_out< float > s2;

	SC_CTOR(sc_xorshift128)
  {
    SC_CTHREAD(do_gen1, clk.pos());
    reset_signal_is(reset, true);

    SC_CTHREAD(do_gen2, clk.pos());
    reset_signal_is(reset, true);
  }

private:
  void do_gen1();
  void do_gen2();
};
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#define TINYMT32_MUL (1.0f / 4294967296.0f)

void sc_xorshift128::do_gen1()
{
	wait( );
	sc_uint<8> mask = seed.read();
	sc_uint<32> xMk = (mask, mask, mask, mask);
	sc_uint<32> x   = 0x12311178 & xMk;
	sc_uint<32> y   = 0x65498732 | xMk;
	sc_uint<32> z   = 0xFEDCAA01 ^ xMk;
	sc_uint<32> w   = 0xF489A179 + xMk;

  /*
    cout << x << endl;
    cout << y << endl;
    cout << z << endl;
    cout << w << endl;
  */

	while( true )
	{
		#pragma HLS PIPELINE
		sc_uint<32> t = x;
		t ^= t << 11;
		t ^= t >> 8;
		x  = y;
		y  = z;
		z  = w;
		w ^= w >> 19;
		w ^= t;
		float f = 1.0f - (float)w * (float)TINYMT32_MUL;
    //printf("%f\n", f);
		s1.write( f );
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
void sc_xorshift128::do_gen2()
{
	wait( );
	sc_uint<8> mask = seed.read();
	sc_uint<32> xMk = (mask, mask, mask, mask);
	sc_uint<32> x = 0x98765432 & xMk;
	sc_uint<32> y = 0x12345678 | xMk;
	sc_uint<32> z = 0xFCBADEFF ^ xMk;
	sc_uint<32> w = 0x12121212 + xMk;
/*
  cout << x << endl;
  cout << y << endl;
  cout << z << endl;
  cout << w << endl;
*/
	while( true )
	{
		#pragma HLS PIPELINE
		sc_uint<32> t = x;
		t ^= t << 11;
		t ^= t >> 8;
		x  = y;
		y  = z;
		z  = w;
		w ^= w >> 19;
		w ^= t;
		float f = 1.0f - (float)w * (float)TINYMT32_MUL;
		s2.write( f );
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
