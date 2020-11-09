/*
*  Adder.h
*  SystemC_SimpleAdder
*
*  Created by Le Gal on 07/05/07.
*  Copyright 2007 __MyCompanyName__. All rights reserved.
*
*/
#ifndef _sc_xer_compute_
#define _sc_xer_compute_

#include "systemc.h"

#ifdef __RTL_SIMULATION__
#include "config.h"
#else
#include "../../module/config.h"
#endif

SC_MODULE(sc_xer_compute)
{
  sc_in < bool > clk;
  sc_in < bool > reset;

  sc_in < sc_uint<32> > decs;
  sc_in < sc_uint<32> > errs;

  sc_out< sc_uint<32> > xer;

  SC_CTOR(sc_xer_compute)
  {
    SC_CTHREAD(do_compute, clk.pos());
    reset_signal_is(reset, true);
  }

private:
  void do_compute()           // process
  {
  	xer.write( 0 );
  	wait( );

  	while( true )
  	{
  		wait( );
  		sc_uint<64> tmp = 100000000 * errs.read();
  		sc_uint<32> val = tmp / (decs.read() + 1);
  		xer.write( val );
  	}
  }
};

//
////////////////////////////////////////////////////////////////////////////////
//
//
#endif
