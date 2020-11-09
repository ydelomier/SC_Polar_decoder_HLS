/*
 *  sc_fake_awgn.h
 *  sc_fake_awgn
 *
 *  Created by Bertrand LE GAL on 06/11/2016.
 *  Copyright 2016 Bertrand LE GAL. All rights reserved.
 *
 *  USEFULL only for sc_intel_mkl channel simulation
 */
#ifndef _sc_fake_awgn_
#define _sc_fake_awgn_

#include "systemc.h"

SC_MODULE(sc_fake_awgn)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_in < float > e1;
  sc_fifo_in < float > e2;

  sc_fifo_out< float > s;

  SC_CTOR(sc_fake_awgn)
  {
    SC_CTHREAD(do_action, clk.pos());
    reset_signal_is(reset,true);
  }

private:
  void do_action();
};

void sc_fake_awgn::do_action()
{
	while( true ){
    s.write( e1.read() );
    s.write( e2.read() );
	}
}

#endif
