//test_p.h

/*********************************************************************/
/* This source code file is part of PRIMoRDiA software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

#ifndef TEST_PRIMORDIA
#define TEST_PRIMORDIA

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "../include/common.h"

//==================================================================
/**
 * @class test_p
 * @author barden
 * @date 17/04/18
 * @file test_p.h
 * @brief Class to instantiate a object with a collection of memeber functions for PRIMoRDiA-libs unit testing.
 */
class test_p {
	public:			
		test_p();
		~test_p();
		void test_systems_support();
		void init_general_test();
		void unit_test();
		void test_reaction_analysis();
		void test_traj_analysis();
};

#endif
//================================================================================
//END OF FILE
//================================================================================