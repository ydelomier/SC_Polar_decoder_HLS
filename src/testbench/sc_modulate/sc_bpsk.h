/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_bpsk_
#define _sc_bpsk_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#else
#include "../../module/config.h"
#endif


#include <iostream>
using namespace std;

SC_MODULE(sc_bpsk)          	// module (class) declaration
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_fifo_in < BIT   > e; 		// vers de decodeur
    sc_fifo_out< float > s; 		// vers de decodeur

    SC_CTOR(sc_bpsk)
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
void sc_bpsk::do_gen()           // process
{
	while( true )
	{
		#pragma HLS PIPELINE

		//float value = (e.read() == 1) ? +1.0f : -1.0f;
		float value = (e.read() == 1) ? -1.0f : +1.0f;
		s.write( value );
	}
}
//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
