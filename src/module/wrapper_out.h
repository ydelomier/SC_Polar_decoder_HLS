#ifndef _wrapper_out_
#define _wrapper_out_

#include "systemc.h"

#include "polar_parameters.h"

#include "config.h"

#include "../../../shared/src/library.h"


SC_MODULE(wrapper_out)
{
	sc_in <bool> clk;
	sc_in <bool> reset;

    sc_fifo_in < TYPE_BITS > e;
    sc_fifo_out< BIT > s;

	SC_CTOR(wrapper_out){
		SC_CTHREAD(do_action, clk.pos());
		reset_signal_is(reset,true);
	}

	void do_action(){

		while(1){
			TYPE_BITS bit = e.read();

			for(COUNTER i = 0; i < PAR; i++){
				s.write( (BIT) bit[i] );
			}

		}
	}
};


#endif
