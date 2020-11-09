/*
 *  sc_awgn.h
 *  sc_awgn
 *
 *  Created by Bertrand LE GAL on 06/11/2016.
 *  Copyright 2016 Bertrand LE GAL. All rights reserved.
 *
 */
#ifndef _sc_awgn_
#define _sc_awgn_

#include "systemc.h"
#include <cmath>

SC_MODULE(sc_awgn)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_in < float > e1;
  sc_fifo_in < float > e2;

  sc_fifo_out< float > s;

  SC_CTOR(sc_awgn)
  {
    SC_CTHREAD(do_action, clk.pos());
    reset_signal_is(reset,true);
  }

private:
  void do_action();
};
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#ifdef _SIMULATION_
	#define sincosf __sincosf
#endif
//
//
////////////////////////////////////////////////////////////////////////////////
//
//

float mult2(float value)
{
	unsigned int f1 = static_cast<unsigned int>(value);
	sc_uint<32>  f2 = f1;
	sc_uint< 8>  ex = f2.range(30, 23);
	sc_uint< 8>  eX = ex + 1;
	sc_uint<32>  f3 = (f2.test(32), eX, f2.range(22, 0));
	unsigned int s1 = f3.to_uint();
	float        f4 = static_cast<float>(s1);
	return f4;
}

void sc_awgn::do_action()
{
	const float _1PI = 3.14159265358979f;
	const float _2PI = 2.0f * _1PI;

	while( true ){
		#pragma HLS PIPELINE
		float r1 = e1.read();
		float y  = _2PI * e2.read();
		float x  = sqrtf (-2.0f * logf( r1 ) );
		float vsin, vcos;
#if SYNTHESIS
		sincosf( y, &vsin, &vcos);
#else
		vsin = sinf(y);
		vcos = cosf(y);
#endif
		float Ph = x * vsin;
		float Qu = x * vcos;
/*
    printf("%f\n", r1);
    printf("%f\n", y/_2PI);
    printf("%f\n", Ph);
    printf("%f\n", Qu);
    exit( 0 );
*/    
		s.write( Ph );
		s.write( Qu );
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
