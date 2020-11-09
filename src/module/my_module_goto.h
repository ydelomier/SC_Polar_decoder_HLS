#ifndef _my_module_
#define _my_module_

#include "systemc.h"

#include "polar_parameters.h"

#include "config.h"

#include "../../../shared/src/library.h"

//#define DEBUG_LIGHT
//#define DEBUG

SC_MODULE(my_module)
{
	sc_in <bool> clk;
	sc_in <bool> reset;

	// Monitor
#ifdef _MONITORING_
	sc_out <bool> busy;
	sc_out <bool> load;
	sc_out <bool> store;
	sc_out <int> N_value;
	sc_out <COUNTER> dgb_n;
	sc_out <COUNTER> cnt_i;
	sc_out < sc_uint<8> > Fct_ID;  // F: 0x88, G: 0x48, H: 0x28, R: 0x18
	  	  	  	  	  	  	  	  // R0: 0x00, R1: 0x0F, REP: 0x02, SPC: 0x04
#endif
	//

    sc_fifo_in < BIT > FB;

    sc_fifo_in < TYPE_LLRS > e;
    sc_fifo_out< TYPE_BITS > s;

	SC_CTOR(my_module):
		enable ("en")
	{
		SC_CTHREAD(do_prunning, clk.pos());
		reset_signal_is(reset,true);

		SC_CTHREAD(do_action, clk.pos());
		reset_signal_is(reset,true);
	}

private:

	sc_signal < bool > enable;

	TYPE_LLRS llr_mem_a[N_DIVIDED];
	TYPE_LLRS llr_mem_b[N_DIVIDED];

	TYPE_BITS bit_mem_1[N_DIVIDED];
	TYPE_BITS bit_mem_2[N_DIVIDED];
	
	sc_bv<PAR> Bit_Frozen[N_DIVIDED];
	sc_uint<4> Node_Type[N_DIVIDED];  // 0x00 : R0, 0x0F : R1, 0x02 : REP, 0x03 : REP2, 0x04 : SPC, 0x05 : SPC2, 0x08 : ?

	void do_prunning(){

		bool finish = false;
		enable.write(false);
		wait();

		while(true){

			if(finish == false)
			{

				for( COUNTER i = 0; i < N_DIVIDED; i++)
				{
					BIT fb;
					BIT R0 = 0;
					BIT R1 = 1;

					BIT SPC_1st = 0;
					BIT SPC_R1 = 1;

					BIT SPC_2nd = 0;
					BIT SPC_2nd_R1 = 1;

					BIT REP_R0 = 0;
					BIT REP_last = 1;

					BIT REP_2Last_R0 = 0;
					BIT REP_2Last = 1;

					sc_bv<PAR> tab_FB;

					for(COUNTER k = 0; k < PAR; k++)
					{
						fb = FB.read();
						tab_FB = ( (sc_bv<1>)fb ,  tab_FB.range( PAR-1, 1 ) ) ;
						R0 = R0 | fb;
						R1 = R1 & fb;
						if(k == 0)
							SPC_1st = fb;
						else
							SPC_R1 = SPC_R1 & fb;

						if(k == 0)
							SPC_2nd = fb;
						else if (k == 1)
							SPC_2nd = SPC_2nd | fb;
						else
							SPC_2nd_R1 = SPC_2nd_R1 & fb;

						if ( k == (PAR-1) )
							REP_last = fb;
						else
							REP_R0 = REP_R0 | fb;

						if ( k == (PAR-2) )
							REP_2Last = fb;
						else if(k == (PAR-1) )
							REP_2Last = REP_2Last & fb;
						else
							REP_2Last_R0 = REP_2Last_R0 | fb;
					}
					SPC_1st = (~SPC_1st) & SPC_R1;
					SPC_2nd = (~SPC_2nd) & SPC_2nd_R1;
					REP_last = REP_last & (~REP_R0);
					REP_2Last = REP_2Last & (~REP_2Last_R0);

					Bit_Frozen[i] = tab_FB;

#if  PRUNING_LEVEL > 0
					if(R0 == 0)
						Node_Type[i] = NODE_R0;  // Node R0
	#if  ELAG_R1 == 1
					else if(R1 == 1)
						Node_Type[i] = NODE_R1;  // Node R1
	#endif
	#if  ELAG_REP == 1
					else if(REP_last == 1)
						Node_Type[i] = NODE_REP;  // Node REP
	#endif
	#if  ELAG_SPC == 1
					else if(SPC_1st == 1)
						Node_Type[i] = NODE_SPC;  // Node SPC
	#endif
	#if  ELAG_REP2 == 1
					else if(REP_2Last == 1)
						Node_Type[i] = NODE_REP2;  // Node REP2
	#endif
	#if  ELAG_SPC2 == 1
					else if(SPC_2nd == 1)
						Node_Type[i] = NODE_SPC2;  // Node SPC2
	#endif
					else
#endif
						Node_Type[i] = NODE_RN;  // Node ?

				}

				enable.write(true);
			}

			finish = true;
			wait();

		}

	}


	
/*********************************************************************************/
/*** 								DECODER PROCESS							   ***/
/*********************************************************************************/

	void do_action(){

		COUNTER ptr_FB = 0;
		sc_bv<PAR> Is_Frozen = 0;
		sc_uint<4> Node = 0;
		COUNTER N_REG = 0;
		COUNTER NB_ITER = 0;
		bool R_State_condition = 0;

		// Pointer llr_mem
		COUNTER adr_a = 0;
		COUNTER adr_b = 0;
		COUNTER adr_w_a = 0;
		COUNTER adr_w_b = 0;
		//Pointer bit_mem
		COUNTER adr_s = 0;
		COUNTER ps_adr_a = 0;
		COUNTER ps_adr_b = 0;
		COUNTER ps_adr = 0;

		COUNTER i = 0;

		TYPE_LLRS reg_result = 0;
		// stack
		sc_biguint<DEPTH_DIV * 2> stack = (sc_biguint<DEPTH_DIV * 2>) 0;
		sc_biguint< 2 > G_stack_value = 1;
		// Node type stack
		sc_biguint<DEPTH_DIV * 8> Node_type_stack = (sc_biguint<DEPTH_DIV * 8>) 0;

		/// ELAG NODE COMPUTATION ///
		sc_uint<4> Node_T = 0;
		sc_uint<4> l_R0 = 0x00;
		sc_uint<4> l_R1 = 0x0F;
		sc_uint<4> l_SPC_1st = 0x00;
		sc_uint<4> l_SPC_R1 = 0x0F;
		sc_uint<4> l_REP_R0 = 0x00;
		sc_uint<4> l_REP_last = 0x00;
		sc_uint<4> r_R0 = 0x00;
		sc_uint<4> r_R1 = 0x0F;
		sc_uint<4> r_SPC_1st = 0x00;
		sc_uint<4> r_SPC_R1 = 0x0F;
		sc_uint<4> r_REP_R0 = 0x00;
		sc_uint<4> r_REP_last = 0x00;
		sc_uint<4> left_Node = 0;
		sc_uint<4> right_Node = 0;

		//
		TYPE_LLRS la = 0;
		TYPE_LLRS lb = 0;
		TYPE_BITS sa = 0;
		TYPE_BITS sb = 0;
		TYPE_LLRS result = 0;

#ifdef _MONITORING_
		load.write(false);
		store.write(false);
		busy.write(false);
		N_value.write(0);
		Fct_ID.write(0);
#endif
		wait();
////////////// RESET  //////////////

		while(true){

			// Load Frame to memory
#ifdef _MONITORING_
			store.write(false);
			load.write(true);
#endif

			load_loop_a: for (i = 0; i < (N_DIVIDED >> 1); i++)
			{
			#pragma HLS PIPELINE
				llr_mem_a[i] = e.read();
#ifndef __SYNTHESIS__
				wait();
#endif
			}
			load_loop_b: for ( i = 0; i < (N_DIVIDED >> 1); i++)
			{
			#pragma HLS PIPELINE
				llr_mem_b[i] = e.read();
#ifndef __SYNTHESIS__
				wait();
#endif
			}

			////////////// INITIALIZATION  //////////////

			// Wait for Frozen Bit initialization
			while( enable.read() == false)
				wait();

			// Initialize pointer
			ptr_FB = 0;
			Is_Frozen = Bit_Frozen[ptr_FB];
			Node = Node_Type[ptr_FB];
			N_REG = N_DIVIDED;
			adr_a = 0;
			adr_b = 0;
			adr_w_a = (N_DIVIDED >> 1);
			adr_w_b = (N_DIVIDED >> 1);
			adr_s = 0;
			stack = 0;
			Node_type_stack = push_stack<8, DEPTH_DIV>( Node_type_stack , ((sc_biguint<4>)NODE_RN , (sc_biguint<4>)NODE_RN) );  // First 2 nodes are never pruned

			/**************************************************************/

				F_STATE:
				{
					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					NB_ITER = (N_REG >> 1);
					N_REG = (N_REG >> 1);
#ifdef _MONITORING_
					load.write(false);
					busy.write(true);
					N_value.write((int) N_REG );
					Fct_ID.write(0x88);
#endif

					stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management

#ifdef DEBUG_LIGHT
cout << "	+++++ [F_STATE]  N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", w_a : "<< adr_w_a << ", w_b : "<< adr_w_b << endl;
#endif

					/// ELAG NODE COMPUTATION ///
					l_R0 = 0x00;
					l_R1 = 0x0F;
					l_SPC_1st = 0x00;
					l_SPC_R1 = 0x0F;
					l_REP_R0 = 0x00;
					l_REP_last = 0x00;
					r_R0 = 0x00;
					r_R1 = 0x0F;
					r_SPC_1st = 0x00;
					r_SPC_R1 = 0x0F;
					r_REP_R0 = 0x00;
					r_REP_last = 0x00;
					/////////////////////////////

					f_loop: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE
						/*%%%%%%%%%%%%%%%%%%*/
						cnt_i.write(i );
						/*%%%%%%%%%%%%%%%%%%*/
						la = llr_mem_a[adr_a];
						lb = llr_mem_b[adr_b];

						result = PU_FUNCTION_F< PAR, LLR_BITS > (la, lb);

						reg_result = result; // Register for R_STATE

#ifdef DEBUG
cout << "	+++++ [F_STATE] " << endl;
cout << "					la =  "; SHOW_LLRS<PAR, LLR_BITS >(la); cout << endl;
cout << "					lb =  "; SHOW_LLRS<PAR, LLR_BITS >(lb); cout << endl;
cout << "					rs =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
#endif

						if(NB_ITER == 1){
							llr_mem_a[adr_w_a] = result;
						}
						else{
							if(i < (NB_ITER >> 1)){
								llr_mem_a[adr_w_a] = result;
								adr_w_a++;

								/// ELAG NODE COMPUTATION ///
								Node_T = Node_Type[ptr_FB + i];
								l_R0 = VECTOR_OR( l_R0 , Node_T ) ;
								l_R1 = VECTOR_AND( l_R1 , Node_T );
								if(i == 0)
									l_SPC_1st = Node_T;
								else
									l_SPC_R1 = VECTOR_AND( l_SPC_R1 , Node_T );

								if ( i == ((NB_ITER >> 1) - 1) )
									l_REP_last = Node_T;
								else
									l_REP_R0 = VECTOR_OR( l_REP_R0 , Node_T );
								/////////////////////////////
							}
							else{
								llr_mem_b[adr_w_b] = result;
								adr_w_b++;

								/// ELAG NODE COMPUTATION ///
								Node_T = Node_Type[ptr_FB + i];
								r_R0 = VECTOR_OR( r_R0 , Node_T ) ;
								r_R1 = VECTOR_AND( r_R1 , Node_T );
								if(i == (NB_ITER >> 1) )
									r_SPC_1st = Node_T;
								else
									r_SPC_R1 = VECTOR_AND( r_SPC_R1 , Node_T );

								if ( i == (NB_ITER - 1) )
									r_REP_last = Node_T;
								else
									r_REP_R0 = VECTOR_OR( r_REP_R0 , Node_T );
								/////////////////////////////
							}
						}

						adr_a++;
						adr_b++;

#ifndef __SYNTHESIS__
						wait();
#endif
					}

					/// ELAG NODE COMPUTATION ///
					if(l_R0 == NODE_R0)
						left_Node = NODE_R0;  // Node R0
					else if(l_R1 == NODE_R1)
						left_Node = NODE_R1;  // Node R1
					else if( ( l_REP_R0 == NODE_R0 ) && ( (sc_uint<3>)l_REP_last.range(3,1) == 0x01 ) )
						left_Node = l_REP_last;  // Node REP / REP2
					else if( ( l_SPC_R1 == NODE_R1 ) && ( (sc_uint<3>)l_SPC_1st.range(3,1) == 0x02 ) )
						left_Node = l_SPC_1st;  // Node SPC / SPC2
					else
						left_Node = NODE_RN;  // Node ?

					if(r_R0 == NODE_R0)
						right_Node = NODE_R0;  // Node R0
					else if(r_R1 == NODE_R1)
						right_Node = NODE_R1;  // Node R1
					else if( ( r_REP_R0 == NODE_R0 ) && ( (sc_uint<3>)r_REP_last.range(3,1) == 0x01 ) )
						right_Node = r_REP_last;  // Node REP / REP2
					else if( ( r_SPC_R1 == NODE_R1 ) && ( (sc_uint<3>)r_SPC_1st.range(3,1) == 0x02 ) )
						right_Node = r_SPC_1st;  // Node SPC / SPC2
					else
						right_Node = NODE_RN;  // Node ?

					Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)(left_Node,right_Node) ); // Node Stack management
					/////////////////////////////

					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					if(N_REG > 1){
#if  PRUNING_LEVEL == 2
						switch(left_Node){
#if ELAG_H0 == 1
							case NODE_R0 :{
								N_REG = (N_REG >> 1);
								stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management
								Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
								adr_s += N_REG;
								ptr_FB += N_REG;
								Is_Frozen = Bit_Frozen[ptr_FB];
								Node = Node_Type[ptr_FB];

								sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
								sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node
								ps_adr = adr_s - N_REG;
								G_stack_value = (sc_biguint<2>) 2;
								switch(right_Node){  // Node R0 and REP are not possible
			#if  ELAG_R1 == 1
									case NODE_R1 :
										goto G_R1_STATE; break;
			#endif
			#if  ELAG_SPC == 1
									case NODE_SPC :
										goto G_SPC_STATE; break;
			#endif
									default :
										goto G_STATE; break; //next state : G_STATE; Rate ?
								}
								break;
							}
#else
							case NODE_R0 :
								goto F_R0_STATE; break;
#endif
	#if ELAG_RARE == 1
		#if  ELAG_R1 == 1
							case NODE_R1 :
								goto F_R1_STATE; break;
		#endif
		#if  ELAG_SPC == 1
							case NODE_SPC :
								goto F_SPC_STATE; break;
		#endif
	#endif
	#if  ELAG_REP == 1
							case NODE_REP :
								goto F_REP_STATE; break;
	#endif
							default :
								goto F_STATE; break; //next state : F_STATE; Rate ?
						}
#else
						goto F_STATE; break; //next state : F_STATE; Rate ?
#endif
					}
					else{
						adr_a--;
						adr_b--;
						R_State_condition = true;
						goto R_STATE; //next state : R_STATE;
					}
				}

			/**************************************************************/

				R_STATE:
				{
					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/
					////
					sc_uint<3> type = (sc_uint<3>) Node.range(3,1);
					sc_uint<1> sel = (sc_uint<1>) Node[0];
#ifdef _MONITORING_
					N_value.write((int)N_REG);
					sc_uint<4> _id = Node;
					Fct_ID.write( ( (sc_uint<4>)(1), Node) );
#endif

					ptr_FB += 1;

					TYPE_BITS ps;

#ifdef DEBUG_LIGHT
cout << "	_____ [R_STATE]  N_value : " << N_REG << ", Node : " << Node << ", FB : "<< Is_Frozen << ", adr_s : " << adr_s <<endl;
#endif

#if PRUNING_LEVEL == 1
					switch( type ){    // Node = 0x00 : R0, 0x0F : R1, 0x02 : REP, 0x03 : REP2, 0x04 : SPC, 0x05 : SPC2, 0x08 : ?
									   // type = 0x00 : R0, 0x07 : R1, 0x01 : REP, 0x01 : REP2, 0x02 : SPC, 0x02 : SPC2, 0x04 : ?
						case 0x00 :
							ps = Spec_Node_R0<PAR, LLR_BITS> (reg_result); break;
	#if ELAG_R1 == 1
						case 0x07 :
							ps = Spec_Node_R1<PAR, LLR_BITS> (reg_result); break;
	#endif
	#if  ELAG_REP == 1
						case 0x01 :
		#if ELAG_REP2 == 0
							ps = Spec_REP_Node<PAR, LLR_BITS> (reg_result); break;
		#else
							ps = Spec_REP_REP2_Node<PAR, LLR_BITS> (reg_result, (bool)sel); break;
		#endif
	#endif
	#if  ELAG_SPC == 1
						case 0x02 :
		#if ELAG_SPC2 == 0
							ps = Spec_SPC_Node<PAR, LLR_BITS> (reg_result); break;
		#else
							ps = Spec_SPC_SPC2_Node<PAR, LLR_BITS> (reg_result, (bool)sel); break;
		#endif
	#endif
						default :
							ps = Spec_Polar_Decoder< PAR, LLR_BITS > (reg_result, Is_Frozen);break;
					}
#else
					ps = Spec_Polar_Decoder< PAR, LLR_BITS > (reg_result, Is_Frozen);
#endif

#ifdef DEBUG
cout << "	_____ [R_STATE] " << endl;
cout << "					la =  "; SHOW_LLRS<PAR, LLR_BITS >(reg_result); cout << endl;
cout << "					fb =  "; SHOW_BITS<PAR >(Is_Frozen); cout << endl;
cout << "					ps =  "; SHOW_BITS<PAR >(ps); cout << endl;
#endif

					bit_mem_1[adr_s] = ps;
					bit_mem_2[adr_s] = ps;
					adr_s++;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];

#ifndef __SYNTHESIS__
					wait();
#endif

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack

					////
					if(R_State_condition == true){
						ps_adr = adr_s - N_REG;
#if  PRUNING_LEVEL == 2
						switch(right_Node){
	#if ELAG_RARE == 1
							case NODE_R0 :
								G_stack_value = (sc_biguint<2>) 1;
								goto G_R0_STATE; break; // next state : Pruned G
		#if  ELAG_REP == 1
							case NODE_REP :
								G_stack_value = (sc_biguint<2>) 1;
								goto G_REP_STATE; break;
		#endif
	#endif
	#if  ELAG_R1 == 1
							case NODE_R1 :
								G_stack_value = (sc_biguint<2>) 1;
								goto G_R1_STATE; break;
	#endif
	#if  ELAG_SPC == 1
							case NODE_SPC :
								G_stack_value = (sc_biguint<2>) 1;
								goto G_SPC_STATE; break;
	#endif
							default :
								G_stack_value = (sc_biguint<2>) 1;
								goto G_STATE; break; //next state : G_STATE; Rate ?
						}
#else
						G_stack_value = (sc_biguint<2>) 1;
						goto G_STATE; break; //next state : G_STATE; Rate ?
#endif
					}
					else{
						ps_adr_a = adr_s - (N_REG << 1);
						ps_adr_b = adr_s - N_REG;
						if( condition == 1)
							goto H_STATE;
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
						else
							goto H0_STATE;
#endif
					}
				}

			/**************************************************************/

				G_STATE:
				{
					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					NB_ITER = N_REG;

#ifdef _MONITORING_
					//N_value.write( (int)(N_REG<<1) );
					N_value.write( (int) N_REG );
					Fct_ID.write(0x48);
#endif

					stack = write_stack< 2, DEPTH_DIV >( stack, G_stack_value); // Stack management

#ifdef DEBUG_LIGHT
cout << "	같같� [G_STATE]  N_value : " << (N_REG<<1) << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", ps_adr : " << ps_adr << ", w_a : "<< adr_w_a << ", w_b : "<< adr_w_b  << endl;
#endif

					/// ELAG NODE COMPUTATION ///
					l_R0 = 0x00;
					l_R1 = 0x0F;
					l_SPC_1st = 0x00;
					l_SPC_R1 = 0x0F;
					l_REP_R0 = 0x00;
					l_REP_last = 0x00;
					r_R0 = 0x00;
					r_R1 = 0x0F;
					r_SPC_1st = 0x00;
					r_SPC_R1 = 0x0F;
					r_REP_R0 = 0x00;
					r_REP_last = 0x00;
					/////////////////////////////

					g_loop: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE
						/*%%%%%%%%%%%%%%%%%%*/
						cnt_i.write(i );
						/*%%%%%%%%%%%%%%%%%%*/
						la = llr_mem_a[adr_a];
						lb = llr_mem_b[adr_b];
						if( G_stack_value == 2)
							sa = 0;
						else
							sa = bit_mem_1[ps_adr];

						result = PU_FUNCTION_G< PAR, LLR_BITS > (la, lb, sa);

						reg_result = result; // Register for R_STATE
#ifdef DEBUG
cout << "	같같� [G_STATE] " << endl;
cout << "					la =  "; SHOW_LLRS<PAR, LLR_BITS >(la); cout << endl;
cout << "					lb =  "; SHOW_LLRS<PAR, LLR_BITS >(lb); cout << endl;
cout << "					sa =  "; SHOW_BITS<PAR >(sa); cout << endl;
cout << "					rs =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout <<endl;
#endif

						if(NB_ITER == 1){
							llr_mem_b[adr_w_b] = result;
						}
						else{
							if(i < (NB_ITER >> 1)){
								llr_mem_a[adr_w_a] = result;
								adr_w_a++;

								/// ELAG NODE COMPUTATION ///
								Node_T = Node_Type[ptr_FB + i];
								l_R0 = VECTOR_OR( l_R0 , Node_T ) ;
								l_R1 = VECTOR_AND( l_R1 , Node_T );
								if(i == 0)
									l_SPC_1st = Node_T;
								else
									l_SPC_R1 = VECTOR_AND( l_SPC_R1 , Node_T );

								if ( i == ((NB_ITER >> 1) - 1) )
									l_REP_last = Node_T;
								else
									l_REP_R0 = VECTOR_OR( l_REP_R0 , Node_T );
								/////////////////////////////
							}
							else{
								llr_mem_b[adr_w_b] = result;
								adr_w_b++;

								/// ELAG NODE COMPUTATION ///
								Node_T = Node_Type[ptr_FB + i];
								r_R0 = VECTOR_OR( r_R0 , Node_T ) ;
								r_R1 = VECTOR_AND( r_R1 , Node_T );
								if(i == (NB_ITER >> 1) )
									r_SPC_1st = Node_T;
								else
									r_SPC_R1 = VECTOR_AND( r_SPC_R1 , Node_T );

								if ( i == (NB_ITER - 1) )
									r_REP_last = Node_T;
								else
									r_REP_R0 = VECTOR_OR( r_REP_R0 , Node_T );
								/////////////////////////////
							}
						}

						adr_a++;
						adr_b++;
						ps_adr++;

#ifndef __SYNTHESIS__
						wait();
#endif
					}

					/// ELAG NODE COMPUTATION ///
					if(l_R0 == NODE_R0)
						left_Node = NODE_R0;  // Node R0
					else if(l_R1 == NODE_R1)
						left_Node = NODE_R1;  // Node R1
					else if( ( l_REP_R0 == NODE_R0 ) && ( (sc_uint<3>)l_REP_last.range(3,1) == 0x01 ) )
						left_Node = l_REP_last;  // Node REP / REP2
					else if( ( l_SPC_R1 == NODE_R1 ) && ( (sc_uint<3>)l_SPC_1st.range(3,1) == 0x02 ) )
						left_Node = l_SPC_1st;  // Node SPC / SPC2
					else
						left_Node = NODE_RN;  // Node ?

					if(r_R0 == NODE_R0)
						right_Node = NODE_R0;  // Node R0
					else if(r_R1 == NODE_R1)
						right_Node = NODE_R1;  // Node R1
					else if( ( r_REP_R0 == NODE_R0 ) && ( (sc_uint<3>)r_REP_last.range(3,1) == 0x01 ) )
						right_Node = r_REP_last;  // Node REP / REP2
					else if( ( r_SPC_R1 == NODE_R1 ) && ( (sc_uint<3>)r_SPC_1st.range(3,1) == 0x02 ) )
						right_Node = r_SPC_1st;  // Node SPC / SPC2
					else
						right_Node = NODE_RN;  // Node ?

					Node_type_stack = write_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)(left_Node,right_Node) ); // Node Stack management
					/////////////////////////////

					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					////
					if(N_REG > 1){
#if  PRUNING_LEVEL == 2
						switch(left_Node){
#if ELAG_H0 == 1
							case NODE_R0 :{
								N_REG = (N_REG >> 1);
								stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management
								Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
								adr_s += N_REG;
								ptr_FB += N_REG;
								Is_Frozen = Bit_Frozen[ptr_FB];
								Node = Node_Type[ptr_FB];

								sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
								sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node
								ps_adr = adr_s - N_REG;
								G_stack_value = (sc_biguint<2>) 2;
								switch(right_Node){  // Node R0 and REP are not possible
			#if  ELAG_R1 == 1
									case NODE_R1 :
										goto G_R1_STATE; break;
			#endif
			#if  ELAG_SPC == 1
									case NODE_SPC :
										goto G_SPC_STATE; break;
			#endif
									default :
										goto G_STATE; break; //next state : G_STATE; Rate ?
								}
								break;
							}
#else
							case NODE_R0 :
								goto F_R0_STATE; break;
#endif
	#if ELAG_RARE == 1
		#if  ELAG_R1 == 1
							case NODE_R1 :
								goto F_R1_STATE; break;
		#endif
		#if  ELAG_SPC == 1
							case NODE_SPC :
								goto F_SPC_STATE; break;
		#endif
	#endif
	#if  ELAG_REP == 1
							case NODE_REP :
								goto F_REP_STATE; break;
	#endif
							default :
								goto F_STATE; break; //next state : F_STATE; Rate ?
						}
#else
						goto F_STATE; break; //next state : F_STATE; Rate ?
#endif
					}
					else{
						adr_a--;
						adr_b--;
						R_State_condition = false;
						goto R_STATE;
					}
				}

			/**************************************************************/

				H_STATE :
				{
					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

#ifdef _MONITORING_
					N_value.write( (int) N_REG );
					Fct_ID.write(0x28);
#endif

					NB_ITER = N_REG;
					N_REG = (N_REG << 1);

					stack = pop_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0);// Stack management

					Node_type_stack = pop_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>) 0x00 ); // Node Stack management

#ifdef DEBUG_LIGHT
cout << "	----- [H_STATE]  N_value : " << N_REG << ", iteration : "<< NB_ITER << ", ps_adr_a : " <<ps_adr_a << ", ps_adr_b : " <<ps_adr_b << endl;
#endif
					////
					h_loop: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=bit_mem_1 array inter RAW false
					#pragma HLS DEPENDENCE variable=bit_mem_2 array inter RAW false
					#pragma HLS PIPELINE
						/*%%%%%%%%%%%%%%%%%%*/
						cnt_i.write(i );
						/*%%%%%%%%%%%%%%%%%%*/
						TYPE_BITS sa = bit_mem_1[ps_adr_a];
						TYPE_BITS sb = bit_mem_2[ps_adr_b];

						TYPE_BITS a_xor_b = PU_FUNCTION_H< PAR > (sa, sb);

#ifdef DEBUG
cout << "	----- [H_STATE] " << endl;
cout << "					sa =  "; SHOW_BITS<PAR >(sa); cout << endl;
cout << "					sb =  "; SHOW_BITS<PAR >(sb); cout << endl;
cout << "					rs =  "; SHOW_BITS<PAR >(a_xor_b); cout << endl;
#endif

						bit_mem_1[ps_adr_a] = a_xor_b;
						bit_mem_2[ps_adr_a] = a_xor_b;

						ps_adr_a++;
						ps_adr_b++;

#ifndef __SYNTHESIS__
						wait();
#endif
					}
					adr_a 	-= N_REG;
					adr_b 	-= N_REG;
					adr_w_a -= NB_ITER;
					adr_w_b -= NB_ITER;

					////
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					if( condition == 1){
						ps_adr_a = adr_s - (N_REG << 1);
						ps_adr_b = adr_s - N_REG;
						goto H_STATE;
					}
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else if( condition == 2){
						ps_adr_a = adr_s - (N_REG << 1);
						ps_adr_b = adr_s - N_REG;
						goto H0_STATE;
					}
#endif
					else{
						if(ptr_FB == N_DIVIDED )
							goto END;
						else{
							ps_adr = adr_s - N_REG;
#if  PRUNING_LEVEL == 2
							switch(right_Node){
	#if ELAG_RARE == 1
								case NODE_R0 :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_R0_STATE; break; // next state : Pruned G
		#if  ELAG_REP == 1
								case NODE_REP :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_REP_STATE; break;
		#endif
	#endif
	#if  ELAG_R1 == 1
								case NODE_R1 :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_R1_STATE; break;
	#endif
	#if  ELAG_SPC == 1
								case NODE_SPC :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_SPC_STATE; break;
	#endif
								default :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_STATE; break; //next state : G_STATE; Rate ?
							}
#else
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
#endif
						}
					}
				}

#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
			/**************************************************************/
				H0_STATE :
				{
					NB_ITER = N_REG;
					N_REG = (N_REG << 1);

#ifdef _MONITORING_
					N_value.write( (int) N_REG );
					Fct_ID.write(0x20);
#endif

					stack = pop_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0);// Stack management

					Node_type_stack = pop_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>) 0x00 ); // Node Stack management

#ifdef DEBUG_LIGHT
cout << "	::::: [H0_STAT]  N_value : " << N_REG << ", iteration : "<< NB_ITER << ", ps_adr_a : " <<ps_adr_a << ", ps_adr_b : " <<ps_adr_b << endl;
#endif
					////
					h0_loop: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=bit_mem_1 array inter RAW false
					#pragma HLS DEPENDENCE variable=bit_mem_2 array inter RAW false
					#pragma HLS PIPELINE

						TYPE_BITS sb = bit_mem_2[ps_adr_b];

#ifdef DEBUG
cout << "	::::: [H0_STAT] " << endl;
cout << "					sb =  "; SHOW_BITS<PAR >(sb); cout << endl;
#endif

						bit_mem_1[ps_adr_a] = sb;
						bit_mem_2[ps_adr_a] = sb;

						ps_adr_a++;
						ps_adr_b++;

#ifndef __SYNTHESIS__
						wait();
#endif
					}
					adr_a 	-= N_REG;
					adr_b 	-= N_REG;
					adr_w_a -= NB_ITER;
					adr_w_b -= NB_ITER;

					////
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					if( condition == 1){
						ps_adr_a = adr_s - (N_REG << 1);
						ps_adr_b = adr_s - N_REG;
						goto H_STATE;
					}
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else if( condition == 2){
						ps_adr_a = adr_s - (N_REG << 1);
						ps_adr_b = adr_s - N_REG;
						goto H0_STATE;
					}
#endif
					else{
						if(ptr_FB == N_DIVIDED )
							goto END;
						else{
							ps_adr = adr_s - N_REG;
#if  PRUNING_LEVEL == 2
							switch(right_Node){
	#if ELAG_RARE == 1
								case NODE_R0 :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_R0_STATE; break; // next state : Pruned G
		#if  ELAG_REP == 1
								case NODE_REP :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_REP_STATE; break;
		#endif
	#endif
	#if  ELAG_R1 == 1
								case NODE_R1 :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_R1_STATE; break;
	#endif
	#if  ELAG_SPC == 1
								case NODE_SPC :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_SPC_STATE; break;
	#endif
								default :
									G_stack_value = (sc_biguint<2>) 1;
									goto G_STATE; break; //next state : G_STATE; Rate ?
							}
#else
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
#endif
						}
					}
				}
#endif
			/**************************************************************/

#if  PRUNING_LEVEL == 2

#if ELAG_H0 == 0
				F_R0_STATE :
				{

					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					NB_ITER = (N_REG >> 1);
					N_REG = (N_REG >> 1);

#ifdef _MONITORING_
					N_value.write((int)N_REG);
					Fct_ID.write(0x80);
#endif

					stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management

					Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	===== [F_ELAG_] ( R0 ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", adr_s : "<< adr_s   <<  endl;
#endif

					f_pruned_R0: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE
						/*%%%%%%%%%%%%%%%%%%*/
						cnt_i.write(i );
						/*%%%%%%%%%%%%%%%%%%*/
						TYPE_BITS ps = 0;
						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;

						adr_a++;
						adr_b++;
						adr_s++;
					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node R0 " << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					/*%%%%%%%%%%%%%%%%%%*/
					dgb_n.write(N_REG );
					/*%%%%%%%%%%%%%%%%%%*/

					ps_adr = adr_s - N_REG;
					switch(right_Node){
#if ELAG_RARE == 1
						case NODE_R0 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R0_STATE; break; // next state : Pruned G
	#if  ELAG_REP == 1
						case NODE_REP :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_REP_STATE; break;
	#endif
#endif
#if  ELAG_R1 == 1
						case NODE_R1 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R1_STATE; break;
#endif
#if  ELAG_SPC == 1
						case NODE_SPC :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_SPC_STATE; break;
#endif
						default :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
					}
				}
#endif

			/**************************************************************/

#if  ELAG_R1 == 1 && ELAG_RARE == 1
				F_R1_STATE :
				{

					NB_ITER = (N_REG >> 1);
					N_REG = (N_REG >> 1);
#ifdef _MONITORING_
					N_value.write((int)N_REG);
					Fct_ID.write(0x8F);
#endif
					stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management

					Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	===== [F_ELAG_] ( R1 ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", adr_s : "<< adr_s   <<  endl;
#endif

					f_pruned_R1: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE

						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];

						TYPE_LLRS result = PU_FUNCTION_F< PAR, LLR_BITS > (la, lb);

						TYPE_BITS ps = VECTOR_SIGN< PAR, LLR_BITS > (result);

						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;
						adr_a++;
						adr_b++;
						adr_s++;
					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node R1 " << endl;
cout << "					llr =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					ps_adr = adr_s - N_REG;
					switch(right_Node){
#if ELAG_RARE == 1
						case NODE_R0 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R0_STATE; break; // next state : Pruned G
	#if  ELAG_REP == 1
						case NODE_REP :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_REP_STATE; break;
	#endif
#endif
#if  ELAG_R1 == 1
						case NODE_R1 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R1_STATE; break;
#endif
#if  ELAG_SPC == 1
						case NODE_SPC :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_SPC_STATE; break;
#endif
						default :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
					}
				}
#endif

			/**************************************************************/
#if  ELAG_REP == 1
				F_REP_STATE :
				{
					NB_ITER = (N_REG >> 1);
					N_REG = (N_REG >> 1);
#ifdef _MONITORING_
					N_value.write((int)N_REG);
					Fct_ID.write(0x82);
#endif
					stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management

					Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	===== [F_ELAG_] (REP ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", adr_s : "<< adr_s   <<  endl;
#endif
					sc_bigint< LLR_BITS + LOG2_PAR + 1 > sum = 0;
					f_pruned_REP: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE

						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];

						TYPE_LLRS result = PU_FUNCTION_F< PAR, LLR_BITS > (la, lb);

						sum = ADD_TREE_FUNCTION< PAR, LLR_BITS >( result, sum );

						TYPE_BITS ps = 0 ;       // assumption : 0 is the value
						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;

						adr_a++;
						adr_b++;
						adr_s++;
					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node REP " << endl;
cout << "					llr =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					sc_biguint<1> sign = VECTOR_SIGN< 1, (LLR_BITS + LOG2_PAR + 1) > (sum);
					if (sign == 1)
					{
						adr_s-=NB_ITER;
						TYPE_BITS ps = VECTOR_INIT_1< PAR >() ;  // all 1 vector
						for( i = 0; i < NB_ITER ; i++){
							bit_mem_1[adr_s] = ps;
							bit_mem_2[adr_s] = ps;
							adr_s++;
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node REP " << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
						}
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					ps_adr = adr_s - N_REG;
					switch(right_Node){
#if ELAG_RARE == 1
						case NODE_R0 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R0_STATE; break; // next state : Pruned G
	#if  ELAG_REP == 1
						case NODE_REP :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_REP_STATE; break;
	#endif
#endif
#if  ELAG_R1 == 1
						case NODE_R1 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R1_STATE; break;
#endif
#if  ELAG_SPC == 1
						case NODE_SPC :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_SPC_STATE; break;
#endif
						default :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
					}
				}
#endif
			/**************************************************************/
#if  ELAG_SPC == 1 && ELAG_RARE == 1
				F_SPC_STATE :
				{
					NB_ITER = (N_REG >> 1);
					N_REG = (N_REG >> 1);

#ifdef _MONITORING_
					N_value.write((int)N_REG);
					Fct_ID.write(0x84);
#endif
					stack = push_stack< 2, DEPTH_DIV >( stack, (sc_biguint<2>) 0); // Stack management

					Node_type_stack = push_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	===== [F_ELAG_] (SPC ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", adr_s : "<< adr_s   <<  endl;
#endif

					sc_biguint< 1 > parity = 0;
					sc_biguint< LLR_BITS > old_min = 0xFF;
					TYPE_BITS old_mask;
					COUNTER adr_min = 0;
					TYPE_BITS sign_min;
					f_pruned_SPC: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE

						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];

						TYPE_LLRS result = PU_FUNCTION_F< PAR, LLR_BITS > (la, lb);

						TYPE_BITS sign = VECTOR_SIGN< PAR, LLR_BITS > (result);
						parity = PARITY_TREE_FUNCTION< PAR > (sign, parity);

						sc_bigint< LLR_BITS + PAR > min_mask = MIN_MASK_TREE_FCT< PAR, LLR_BITS> (result);
						sc_biguint< LLR_BITS > new_min = (sc_biguint< LLR_BITS >) min_mask.range( (LLR_BITS + PAR) - 1, PAR );
						TYPE_BITS new_mask = (TYPE_BITS) min_mask.range( PAR - 1, 0 );

						if ( new_min < old_min ){
							old_min = new_min;
							old_mask = new_mask;
							adr_min = i;           // Position of the min
							sign_min = sign;
						}

						bit_mem_1[adr_s] = sign;
						bit_mem_2[adr_s] = sign;
						adr_a++;
						adr_b++;
						adr_s++;
					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node SPC " << endl;
cout << "					llr    =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps     =  "; SHOW_BITS<PAR>(sign); cout << endl;
#endif
					}
					if( parity != 0 )
					{
						TYPE_BITS value = VECTOR_XOR< PAR >( sign_min, old_mask );  // Invert the min sign if parity check fail
						bit_mem_1[ adr_s - NB_ITER + adr_min ] = value;
						bit_mem_2[ adr_s - NB_ITER + adr_min ] = value;
#ifdef DEBUG
cout << "	===== [F_ELAG_] Node SPC " << endl;
cout << "					s_min  =  "; SHOW_BITS<PAR>(sign_min); cout << endl;
cout << "					o_mask =  "; SHOW_BITS<PAR>(old_mask); cout << endl;
cout << "					ps_min =  "; SHOW_BITS<PAR>(value); cout << endl;
cout << "					adr    =  "<< adr_min << endl;
#endif
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;

					sc_uint<8> Node_cond = read_stack< 8, DEPTH_DIV >( Node_type_stack, 2 ); // test 2nd element of stack
					sc_uint<4> right_Node = (sc_uint<4>) Node_cond.range(3,0);    // test next G Node

					ps_adr = adr_s - N_REG;
					switch(right_Node){
#if ELAG_RARE == 1
						case NODE_R0 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R0_STATE; break; // next state : Pruned G
	#if  ELAG_REP == 1
						case NODE_REP :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_REP_STATE; break;
	#endif
#endif
#if  ELAG_R1 == 1
						case NODE_R1 :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_R1_STATE; break;
#endif
#if  ELAG_SPC == 1
						case NODE_SPC :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_SPC_STATE; break;
#endif
						default :
							G_stack_value = (sc_biguint<2>) 1;
							goto G_STATE; break; //next state : G_STATE; Rate ?
					}
				}
#endif
			/**************************************************************/

#if ELAG_RARE == 1
				G_R0_STATE :
				{

					NB_ITER = N_REG;

#ifdef _MONITORING_
					N_value.write( (int)N_REG );
					Fct_ID.write(0x40);
#endif
					stack = write_stack< 2, DEPTH_DIV >( stack, G_stack_value); // Stack management

					Node_type_stack = write_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	##### [G_ELAG_] ( R0 ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", ps_adr : "<< ps_adr << ", adr_s : "<< adr_s  << endl;
#endif

					g_pruned_R0: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS PIPELINE
						TYPE_BITS ps = 0;
						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;

						adr_a++;
						adr_b++;
						ps_adr++;
						adr_s++;

					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	##### [G_ELAG_] Node R0 " << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;
					ps_adr_a = adr_s - (N_REG << 1);
					ps_adr_b = adr_s - N_REG;
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack
					if( condition == 1)
						goto H_STATE;
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else
						goto H0_STATE;
#endif
				}
#endif
			/**************************************************************/
#if  ELAG_R1 == 1
				G_R1_STATE :
				{
					NB_ITER = N_REG;

#ifdef _MONITORING_
					N_value.write( (int) N_REG );
					Fct_ID.write(0x4F);
#endif
					stack = write_stack< 2, DEPTH_DIV >( stack, G_stack_value); // Stack management

					Node_type_stack = write_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	##### [G_ELAG_] ( R1 ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", ps_adr : "<< ps_adr << ", adr_s : "<< adr_s  << endl;
#endif

					g_pruned_R1: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS DEPENDENCE variable=bit_mem_1 array inter RAW false
					#pragma HLS PIPELINE
						/*%%%%%%%%%%%%%%%%%%*/
						cnt_i.write(i );
						/*%%%%%%%%%%%%%%%%%%*/
						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];
						TYPE_BITS sa;
						if( G_stack_value == 2)
							sa = 0;
						else
							sa = bit_mem_1[ps_adr];

						TYPE_LLRS result = PU_FUNCTION_G< PAR, LLR_BITS > (la, lb, sa);

						TYPE_BITS ps = VECTOR_SIGN< PAR, LLR_BITS > (result);

						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;

						adr_a++;
						adr_b++;
						ps_adr++;
						adr_s++;

					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	##### [G_ELAG_] Node R1 " << endl;
cout << "					llr =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;
					ps_adr_a = adr_s - (N_REG << 1);
					ps_adr_b = adr_s - N_REG;
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack
					if( condition == 1)
						goto H_STATE;
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else
						goto H0_STATE;
#endif
				}
#endif
			/**************************************************************/
#if  ELAG_REP == 1 && ELAG_RARE == 1
				G_REP_STATE :
				{
					NB_ITER = N_REG;

#ifdef _MONITORING_
					N_value.write( (int)N_REG );
					Fct_ID.write(0x42);
#endif
					stack = write_stack< 2, DEPTH_DIV >( stack, G_stack_value); // Stack management

					Node_type_stack = write_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	##### [G_ELAG_] (REP ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", ps_adr : "<< ps_adr << ", adr_s : "<< adr_s  << endl;
#endif

					sc_bigint< LLR_BITS + LOG2_PAR + 1 > sum = 0;
					g_pruned_REP: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS DEPENDENCE variable=bit_mem_1 array inter RAW false
					#pragma HLS PIPELINE

						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];
						TYPE_BITS sa;
						if( G_stack_value == 2)
							sa = 0;
						else
							sa = bit_mem_1[ps_adr];

						TYPE_LLRS result = PU_FUNCTION_G< PAR, LLR_BITS > (la, lb, sa);

						sum = ADD_TREE_FUNCTION< PAR, LLR_BITS >( result, sum );

						TYPE_BITS ps = 0 ;       // assumption : 0 is the value
						bit_mem_1[adr_s] = ps;
						bit_mem_2[adr_s] = ps;

						adr_a++;
						adr_b++;
						ps_adr++;
						adr_s++;

					#ifndef __SYNTHESIS__
						wait();
					#endif
#ifdef DEBUG
cout << "	##### [G_ELAG_] Node REP " << endl;
cout << "					llr =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
					}
					sc_biguint<1> sign = VECTOR_SIGN< 1, (LLR_BITS + LOG2_PAR + 1) > (sum);
					if (sign == 1)
					{
						adr_s-=NB_ITER;
						TYPE_BITS ps = VECTOR_INIT_1< PAR >() ;  // all 1 vector
						for( i = 0; i < NB_ITER ; i++){
						#pragma HLS PIPELINE
							bit_mem_1[adr_s] = ps;
							bit_mem_2[adr_s] = ps;
							adr_s++;
#ifdef DEBUG
cout << "	##### [G_ELAG_] Node REP " << endl;
cout << "					ps =  "; SHOW_BITS<PAR>(ps); cout << endl;
#endif
						}
					}
					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;
					ps_adr_a = adr_s - (N_REG << 1);
					ps_adr_b = adr_s - N_REG;
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack
					if( condition == 1)
						goto H_STATE;
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else
						goto H0_STATE;
#endif
				}
#endif
			/**************************************************************/
#if  ELAG_SPC == 1
				G_SPC_STATE :
				{
					NB_ITER = N_REG;
#ifdef _MONITORING_
					N_value.write( (int)N_REG );
					Fct_ID.write(0x44);
#endif
					stack = write_stack< 2, DEPTH_DIV >( stack, G_stack_value); // Stack management

					Node_type_stack = write_stack< 8, DEPTH_DIV >( Node_type_stack, (sc_biguint<8>)0x00 ); // Node Stack management
					/////////////////////////////
#ifdef DEBUG_LIGHT
cout << "	##### [G_ELAG_] (SPC ) N_value : " << N_REG << ", iteration : "<< NB_ITER << ", adr_a : "<< adr_a << ", adr_b : "<< adr_b << ", ps_adr : "<< ps_adr << ", adr_s : "<< adr_s  << endl;
#endif

					sc_biguint< 1 > parity = 0;
					sc_biguint< LLR_BITS > old_min = 0xFFFF;
					TYPE_BITS old_mask;
					COUNTER adr_min = 0;
					TYPE_BITS sign_min;
					g_pruned_SPC: for( i = 0; i < NB_ITER ; i++)
					{
					#pragma HLS DEPENDENCE variable=llr_mem_b array inter RAW false
					#pragma HLS DEPENDENCE variable=llr_mem_a array inter RAW false
					#pragma HLS DEPENDENCE variable=bit_mem_1 array inter RAW false
					#pragma HLS PIPELINE

						TYPE_LLRS la = llr_mem_a[adr_a];
						TYPE_LLRS lb = llr_mem_b[adr_b];
						TYPE_BITS sa;
						if( G_stack_value == 2)
							sa = 0;
						else
							sa = bit_mem_1[ps_adr];

						TYPE_LLRS result = PU_FUNCTION_G< PAR, LLR_BITS > (la, lb, sa);

						TYPE_BITS sign = VECTOR_SIGN< PAR, LLR_BITS > (result);

						bit_mem_1[adr_s] = sign;
						bit_mem_2[adr_s] = sign;

						adr_a++;
						adr_b++;
						ps_adr++;
						adr_s++;

						parity = PARITY_TREE_FUNCTION< PAR > (sign, parity);

						sc_bigint< LLR_BITS + PAR > min_mask = MIN_MASK_TREE_FCT< PAR, LLR_BITS> (result);
						sc_biguint< LLR_BITS > new_min = (sc_biguint< LLR_BITS >) min_mask.range( (LLR_BITS + PAR) - 1, PAR );
						TYPE_BITS new_mask = (TYPE_BITS) min_mask.range( PAR - 1, 0 );

						if ( new_min < old_min ){
							old_min = new_min;
							old_mask = new_mask;
							adr_min = i;           // Position of the min
							sign_min = sign;
						}

					#ifndef __SYNTHESIS__
						wait();
					#endif

#ifdef DEBUG
cout << "	##### [G_ELAG_] Node SPC " << endl;
cout << "					la =  "; SHOW_LLRS<PAR, LLR_BITS >(la); cout << endl;
cout << "					lb =  "; SHOW_LLRS<PAR, LLR_BITS >(lb); cout << endl;
cout << "					sa =  "; SHOW_BITS<PAR >(sa); cout << endl;
cout << "					llr    =  "; SHOW_LLRS<PAR, LLR_BITS >(result); cout << endl;
cout << "					ps     =  "; SHOW_BITS<PAR>(sign); cout << endl;
#endif
					}

					if( parity != 0 )
					{
						TYPE_BITS value = VECTOR_XOR< PAR >( sign_min, old_mask );  // Invert the min sign if parity check fail
						bit_mem_1[ adr_s - NB_ITER + adr_min ] = value;
						bit_mem_2[ adr_s - NB_ITER + adr_min ] = value;
#ifdef DEBUG
cout << "	##### [G_ELAG_] Node SPC " << endl;
cout << "					s_min  =  "; SHOW_BITS<PAR>(sign_min); cout << endl;
cout << "					o_mask =  "; SHOW_BITS<PAR>(old_mask); cout << endl;
cout << "					ps_min =  "; SHOW_BITS<PAR>(value); cout << endl;
cout << "					adr    =  "<< adr_min << endl;
#endif
					}

					/////////////////////////////

					ptr_FB += NB_ITER;
					Is_Frozen = Bit_Frozen[ptr_FB];
					Node = Node_Type[ptr_FB];
					adr_a -= NB_ITER;
					adr_b -= NB_ITER;
					ps_adr_a = adr_s - (N_REG << 1);
					ps_adr_b = adr_s - N_REG;
					sc_uint<2> condition = read_stack< 2, DEPTH_DIV >( stack, 1 ); // test first element of stack
					if( condition == 1)
						goto H_STATE;
#if PRUNING_LEVEL == 2 && ELAG_H0 == 1
					else
						goto H0_STATE;
#endif
				}

#endif
			/**************************************************************/
#endif

			END:
			// Send decoded Frame
#ifdef DEBUG_LIGHT
cout << "(DEBUG) send decoded frame " << endl;
#endif

#ifdef _MONITORING_
			busy.write(false);
			store.write(true);
#endif
			send_loop: for ( i = 0; i < N_DIVIDED; i++)
			{
				#pragma HLS PIPELINE
				s.write(bit_mem_1[i]);
#ifndef __SYNTHESIS__
				wait();
#endif
			}

/////////// FINISH /////////////////

		}
	}
};


#endif
