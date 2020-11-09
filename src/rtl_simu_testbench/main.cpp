//#include "./sc_top_module/sc_dual_top_module.h"

#include "sc_top_module.h"
#include <iomanip>

#ifdef __RTL_SIMULATION__
#include "polar_parameters.h"
#else
#include "../module/polar_parameters.h"
#endif

#define _VSATN	-31
#define _VSATP	31
#define _BETA	4


int sc_main (int argc, char * argv []){

	sc_clock clock("clock", 10, SC_NS, 0.5);
	sc_signal< bool > reset ("8");

	sc_report_handler::set_actions (SC_ID_LOGIC_X_TO_BOOL_, SC_DO_NOTHING);

	sc_top_module top ( "top"  );
	top.clk      ( clock     );
	top.reset    ( reset     );

	reset = true;

	cout << "(II) Reset the system" << endl;

	sc_start( 20, SC_NS );

	reset  = false;

	cout << "(II) Running the system" << endl;

	sc_start(30000, SC_US);


  return 0;
}
