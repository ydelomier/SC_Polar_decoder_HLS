#ifndef _sc_top_module_
#define _sc_top_module_

#ifdef __RTL_SIMULATION__
  #include "my_module_rtl_wrapper.h"
  #define sc_fec my_module_rtl_wrapper
#else
  #include "../module/my_module.h"
  #define sc_fec my_module
#endif

#ifdef __RTL_SIMULATION__
	#include "sc_monitor.h"
	#include "sc_generator.h"
#else
	#include "sc_monitor/sc_monitor.h"
	#include "sc_generator/sc_generator.h"
#endif


SC_MODULE(sc_top_module)
{
    sc_in <bool> clk;
    sc_in <bool> reset;

   SC_CTOR(sc_top_module) :
	polar  ("polar" ),
	monitor ("monit"),
	generator("gene"),
	busy ("busy"),
	load("load"),
	store("store"),
	N_value("N_v"),
	Fct_ID ("F_ID"),
	f8 ("wrap_2_dec",  16),  // TYPE LLR
	f9 ("dec_2_wrap",  16),  // TYPE BITS
	FB ("fb", 1)
    {
	    generator.clk(clk);
        generator.reset(reset);
        generator.e(f8);
        generator.FB(FB);

        polar.clk   ( clk      );
        polar.reset ( reset    );
        polar.e     ( f8       );
        polar.s     ( f9       );
        polar.busy (busy);
        polar.load (load);
        polar.store (store);
        polar.N_value (N_value);
        polar.Fct_ID (Fct_ID);
        polar.FB    ( FB       );
        polar.cnt_i(cnt_i);
        polar.dgb_n(dgb_n);

        monitor.cnt_i(cnt_i);
        monitor.clk(clk);
        monitor.reset(reset);
        monitor.s(f9);
        monitor.busy (busy);
        monitor.load (load);
        monitor.store (store);
        monitor.N_value (N_value);
        monitor.Fct_ID (Fct_ID);

    }

private:

    sc_fec             polar;
    sc_monitor		   monitor;
    sc_generator	   generator;

    sc_fifo< TYPE_LLRS > f8;
    sc_fifo< TYPE_BITS > f9;
    sc_fifo< BIT > FB;
    sc_signal <bool> busy;
    sc_signal <bool> load;
    sc_signal <bool> store;
    sc_signal <int> N_value;
    sc_signal <COUNTER> cnt_i;
    sc_signal <COUNTER> dgb_n;
    sc_signal < sc_uint<8> > Fct_ID;

};

#endif
