/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_random_
#define _sc_random_

#include "systemc.h"

#include <iostream>
using namespace std;

SC_MODULE(sc_random)          	// module (class) declaration
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_fifo_out< float > s1; 		// vers de decodeur
    sc_fifo_out< float > s2; 		// vers de decodeur

    SC_CTOR(sc_random)
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
void sc_random::do_gen()
{
#ifdef _SIMULATION_
	while( true )
	{
		float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		s1.write( r1 );

		float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		s2.write( r2 );
	}
#endif
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
