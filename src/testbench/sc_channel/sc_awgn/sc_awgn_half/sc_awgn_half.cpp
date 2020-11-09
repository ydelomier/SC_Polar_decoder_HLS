/*
 *  sc_awgn.cpp
 *  sc_awgn
 *
 *  Created by Bertrand LE GAL on 06/11/2016.
 *  Copyright 2016 Bertrand LE GAL. All rights reserved.
 *
 */

#include "sc_awgn_half.h"

#include <cmath>

#ifdef _SIMULATION_
#define sincosf __sincosf
#endif

sc_lv<16> chalf;

// sc_uint< 1> csign;
// sc_int < 6> cexpo;
// sc_uint<10> cmant;

void sc_awgn_half::do_action()
{
	const half _1PI = 3.14159265358979f;
	const half _2PI = 2.0f * _1PI;

	while( true ){
#pragma HLS PIPELINE
		half r1 = e1.read();
		half y  = e2.read();
		half x  = sqrtf (-2.0f * logf( r1 ) );
		float vsin, vcos;
		sincosf(_2PI * y, &vsin, &vcos);
		half Ph = x * vsin;
		half Qu = x * vcos;
		s.write( Ph );
		s.write( Qu );
	}
}
