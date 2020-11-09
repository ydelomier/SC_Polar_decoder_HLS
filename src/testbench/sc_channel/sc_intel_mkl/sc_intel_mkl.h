/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_intel_mkl_
#define _sc_intel_mkl_
#ifdef _SIMULATION_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_intel_mkl)          	// module (class) declaration
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_fifo_out< float > s1; 		// vers de decodeur
    sc_fifo_out< float > s2; 		// vers de decodeur

    SC_CTOR(sc_intel_mkl)
    {
        SC_CTHREAD(do_gen, clk.pos());
        reset_signal_is(reset, true);
    }

private:
    void do_gen();
};

#include "mkl.h"
  #ifndef VSL_METHOD_SGAUSSIAN_BOXMULLER2
  #define VSL_METHOD_SGAUSSIAN_BOXMULLER2 1
#endif

void sc_intel_mkl::do_gen()
{
	const int nbData = 65536;
	float* noise = new float[2*nbData];

	VSLStreamStatePtr stream;

	int status = vslNewStream( &stream, VSL_BRNG_MT2203, rand() );
	if( status != VSL_STATUS_OK ){
		printf("(EE) Error during vslNewStream execution\n");
		exit( 0 );
	}

	while( true )
	{
		vsRngGaussian( VSL_METHOD_SGAUSSIAN_BOXMULLER2, stream, 2*nbData, noise, 0.0f, 1.0f );
		for(int i=0; i<nbData; i+=1)
		{
			s1.write( noise[2 * i + 0] );
			s2.write( noise[2 * i + 1] );
		}


	}
	vslDeleteStream( &stream );
	delete[] noise;
}

#endif
#endif
