#ifndef _wrapper_in_
#define _wrapper_in_

#include "systemc.h"

#include "polar_parameters.h"

#include "config.h"

#include "../../../shared/src/library.h"


SC_MODULE(wrapper_in)
{
	sc_in <bool> clk;
	sc_in <bool> reset;

    sc_fifo_in < LLR > e;
    sc_fifo_out< TYPE_LLRS > s;

	SC_CTOR(wrapper_in){
		SC_CTHREAD(do_action, clk.pos());
		reset_signal_is(reset,true);
	}

	void do_action(){

		while(1){

			TYPE_LLRS llr_in = 0;

			for(COUNTER i = 0; i < PAR; i++){
				LLR value =  e.read();
				LLR value_adapt = Adapt_format<LLR_BITS> (value);

				//cout << "[WRAPPER] value = " << value << "  ( bin : "<< (sc_bv<LLR_BITS>) value << " ), adapt = " << (sc_bv<LLR_BITS>) value_adapt << endl;

				llr_in = (value_adapt,  llr_in.range( (LLR_BITS * PAR) - 1, LLR_BITS ) ) ;

			}

			s.write(llr_in);
		}
	}
};

#endif
