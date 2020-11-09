/*
*  Adder.h
*  SystemC_SimpleAdder
*
*  Created by Le Gal on 07/05/07.
*  Copyright 2007 __MyCompanyName__. All rights reserved.
*
*/
#ifndef _sc_mbps_compute_
#define _sc_mbps_compute_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#include "polar_parameters.h"
#else
#include "../../module/config.h"
#include "../../module/polar_parameters.h"
#endif



#include <iostream>

SC_MODULE(sc_mbps_compute)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_in < BIT > dec_frame;

  sc_out< sc_uint<32> > mbps;

  SC_CTOR(sc_mbps_compute)
  {
    SC_CTHREAD(do_gen1, clk.pos());
    reset_signal_is(reset, true);

    SC_CTHREAD(do_gen2, clk.pos());
    reset_signal_is(reset, true);
  }

private:

	sc_signal< sc_uint<64> > b_counter;

	void do_gen1()
	{
		//
		// GESTION DU RESET DU COMPOSANT
		//

		b_counter.write( 0 );
		wait( );

		//
		// COMPORTEMENT DU COMPOSANT
		//

		while( true )
		{
			wait();

			if( dec_frame.read() == 1 )
			{
				b_counter.write( b_counter.read() + _NBITS );
			}
		}
	}

	void do_gen2()
	{
		sc_uint<64> last_counter = 0;
		mbps.write     ( 0 );
		wait( );

		while( true )
		{
			sc_uint<32> c_counter = 99999999 - 6; // cf rapport de synthese HLS
			do
			{
				c_counter -= 1;
				wait();
			}while( c_counter != 0 );

			sc_uint<64> diff = (b_counter.read() - last_counter);
			sc_uint<32> bits = diff.range(31, 0) / 1000; // pour obtenir une conversion en kBPS
			last_counter     = b_counter.read();
			mbps.write( bits );
		}
	}

};

//
//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
