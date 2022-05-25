//test_p.cpp

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

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <omp.h>


#include "../include/test_p.h"
#include "../include/autoprimordia.h"

//#include <gtest/gtest.h>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::move;

/************************************************************************************/
test_p::test_p(){
}
/************************************************************************************/
test_p::~test_p(){}
/************************************************************************************/
void test_p::init_general_test(){
	//AutoPrimordia run_test("/home/igorchem/Documents/tutorials_1.2/Tutorial_2/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/Documents/tutorials_1.2/Tutorial_3/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/Documents/tutorials_1.2/Tutorial_5/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/Documents/tutorials_1.2/Tutorial_1/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/Documents/tutoriais_primordia/tutorials_1_2/Tutorial_1/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/CCDIR/primordia_FeDFO/primordia.input");
	AutoPrimordia run_test("/home/igorchem/CCDIR/FURG/4qsu_protein_results/primordia.input");
	//AutoPrimordia run_test("/home/igorchem/CCDIR/calcs_turco/primordia.input");
	run_test.init();
}
/************************************************************************************/
void test_p::unit_test(){
	
}
/************************************************************************************/
void test_p::test_reaction_analysis(){
	AutoPrimordia run_test1d("/home/igorchem/CCDIR/primordia_FeDFO/primordia.input");
	//AutoPrimordia run_test1d("/media/igorchem/CCDIR/ccdir/primordia_testes/deHalo_big335_PM7/primordia.input");
	run_test1d.init();
	//AutoPrimordia run_test2d("/home/igorchem/Documents/PM6_16x16/primordia.input");
	//run_test2d.init();
}

/************************************************************************************/
void test_p::test_traj_analysis(){
	AutoPrimordia Run("/home/igorchem/CCDIR/cmp2-2/primordia.input");
	Run.init();
}
//================================================================================
//END OF FILE
//================================================================================