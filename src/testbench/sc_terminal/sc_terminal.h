/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_terminal_
#define _sc_terminal_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#else
#include "../../module/config.h"
#endif


#include <iostream>

SC_MODULE(sc_terminal)          	// module (class) declaration
{
    sc_in <bool> clk;
    sc_in <bool> reset;

    sc_in< sc_uint<64> > err_bits;
  	sc_in< sc_uint<64> > err_frames;
  	sc_in< sc_uint<64> > proc_bits;
  	sc_in< sc_uint<64> > proc_frames;
  	sc_in< sc_uint<32> > mbps_value;
  	sc_in< sc_uint<32> > fer_value;

    SC_CTOR(sc_terminal)
    {
        SC_CTHREAD(do_gen, clk.pos());
        reset_signal_is(reset, true);
    }

private:

    void do_gen()           // process
    {
    	int start = 0;

    	sc_uint<32> last_frames = 0;
    	int counter = 0;

    	while( true )
    	{
            wait();
            if( last_frames != proc_frames.read() )
            {
    					unsigned long long errBER   = err_bits.read();
    					unsigned long long errFER   = err_frames.read();
    					unsigned long long nbBits   = proc_bits.read();
    					unsigned long long nbFrames = proc_frames.read();
    					unsigned long long cFERv    = fer_value.read();

    					float BER = (float)errBER/(float)nbBits;
    					float FER = (float)errFER/(float)nbFrames;
    					int end   = 0;
    					int eTime = 0;
    					printf("(DD) sc_terminal :: Run %6u | %10llu bits | %6llu BE(s) | %6llu frames | %6llu FE(s) | BER = %1.2e | FER = %1.2e (%8llu) | TIME = %d sec. | ", counter, nbBits, errBER, nbFrames, errFER, BER, FER, cFERv, (int)eTime );
    					cout << sc_time_stamp() << "      \r";
    					fflush( stdout );

    					last_frames = proc_frames.read();
            }
    	}
    }

};


#endif
