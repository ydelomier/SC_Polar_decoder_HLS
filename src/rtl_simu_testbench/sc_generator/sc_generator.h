#ifndef _sc_generator_
#define _sc_generator_

#include "systemc.h"
#include "math.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#include "polar_parameters.h"
#else
#include "../../module/config.h"
#include "../../module/polar_parameters.h"
#endif

#include <iostream>


SC_MODULE(sc_generator)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_out < BIT > FB;

  sc_fifo_out < TYPE_LLRS > e;

  SC_CTOR(sc_generator)
  {
	SC_CTHREAD(do_fb, clk.pos());
	reset_signal_is(reset, true);

	SC_CTHREAD(write_fifo, clk.pos());
	reset_signal_is(reset, true);
  }

private:

	void do_fb() // process
	{
		bool finish = false;
		wait();

		while(true){

			if(finish == false)
			{
				for(COUNTER i = 0; i < _NBITS; i++)
				{
					FB.write( (BIT) Frozen_Bits[i]);
				}
			}

			finish = true;
			wait();

		}
	}

	void write_fifo(){

		wait();
		while(true)
		{
			TYPE_LLRS value = (TYPE_LLRS) 21;
			e.write(value);
			wait();
		}
	}
};

#endif
