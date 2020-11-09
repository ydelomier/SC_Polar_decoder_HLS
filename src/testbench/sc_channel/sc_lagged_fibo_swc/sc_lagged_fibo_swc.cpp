/*
*  Adder.cpp
*  SystemC_SimpleAdder
*
*  Created by Le Gal on 07/05/07.
*  Copyright 2007 __MyCompanyName__. All rights reserved.
*
*/

#include "sc_lagged_fibo_swc.h"

/*
* 29/10/2016
*   DANS SA FORME ACTUELLE LE MODULE SYNTHETISE SUR NEXYS 4 AVEC LES PARAMETRES
*   SUIVANTS POUR LA LOOP 1 (CONTRAINTE = 10ns + 2ns uncertaincy)
*   - LATENCY  = X
*   - INTERCAL = X
*   COUT DE L'ARCHITECTURE
*   - X LUTS
*   - X FFs
*   - X DSPs
*
*/

#define TINYMT32_MUL (1.0f / 4294967296.0f)

void sc_lagged_fibo_swc::do_gen1()
{
//	sc_uint<32> x = 0x12345678;
//	sc_uint<32> y = 0X98765432;
//	sc_uint<32> z = 0x12121212;
//	sc_uint<32> w = 0xFCBADEFF;
	for(int i=0; i<1024; i++){
		sc_uint<64> pattern;
		pattern.range(31,  0) = rand();
		pattern.range(63, 32) = rand();
		lastElements[ i ]     = pattern;
	}

	sc_uint< 1> c = 0;
	wait( );

	int i = 0;
	while( true )
	{
//		#pragma HLS PIPELINE
// 273,607

		sc_uint<65> first  = lastElements[ (1024 + i - 273) & 1024 ];
		sc_uint<64> second = lastElements[ (1024 + i - 607) & 1024 ];
		sc_uint<65> result = first - second - c;
		c = result.test(65);
		lastElements[ i ] = result.range(63, 0);
		i = (i + 1) & 1024;

		float f1 = (float)c.range(31,  0) * (float)TINYMT32_MUL;
		float f2 = (float)c.range(63, 32) * (float)TINYMT32_MUL;
		s1.write( f1 );
		s2.write( f2 );
	}
}
