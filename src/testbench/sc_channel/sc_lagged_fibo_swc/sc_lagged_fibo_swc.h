/*
 *  Adder.h
 *  SystemC_SimpleAdder
 *
 *  Created by Le Gal on 07/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _sc_lagged_fibo_swc_
#define _sc_lagged_fibo_swc_

#include "systemc.h"

#include <iostream>
using namespace std;
SC_MODULE(sc_lagged_fibo_swc)
{
  sc_in <bool> clk;
  sc_in <bool> reset;

  sc_fifo_out< float > s1;
  sc_fifo_out< float > s2;

	SC_CTOR(sc_lagged_fibo_swc)
  {
    SC_CTHREAD(do_gen1, clk.pos());
    reset_signal_is(reset, true);

    SC_CTHREAD(do_gen2, clk.pos());
    reset_signal_is(reset, true);
  }

private:
  void do_gen1();
  void do_gen2();

  sc_uint<64> lastElements[1024];

};

#endif
