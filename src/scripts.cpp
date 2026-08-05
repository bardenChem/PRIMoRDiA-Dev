//scripts.cpp

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
#include <cstring>
#include <string>
#include <fstream>

#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/scripts.h"
#include "../include/local_rd.h"
#include "../include/Icube.h"
#include "../include/ReactionAnalysis.h"
#include "../include/reaction_coord.h"
#include "../include/pos_traj.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::to_string;

/************************************************************************/
scripts::scripts()		:
	file_name("noname")	,
	s_type("notype")	{
		
	script_file.open("default.R");

}
/************************************************************************/
scripts::scripts(	string Nm	,
					string _type 	):
	file_name(Nm)					,
	s_type(_type)					{
	
	string fname = file_name;
	fname += "_";
	fname += s_type;
	if ( s_type == "pymols" || s_type == "pymols_pdb" ){
		fname += ".pym";
	}else{
		fname += ".R";
	}
	script_file.open( fname.c_str() );
	
	if( s_type == "pymols" || s_type == "pymols_pdb" ){
		script_file << "preset.publication(selection='all')\n"
					<< "set sphere_scale, 0.2\n"
					<< "set bg_rgb, white \n"
					<< "set stick_radius, 0.18\n";
	}else{
		script_file << "#!/usr/bin/env Rscript\n";
	}	
}
/**********************************************************/
scripts::scripts(const scripts& rhs):
	file_name(rhs.file_name)     	,
	s_type(rhs.s_type)				{
	
}
/***********************************************************/
scripts& scripts::operator=(const scripts& rhs){
	if ( this != &rhs){
		file_name	= rhs.file_name;
		s_type		= rhs.s_type;
	}
	return *this;
}
/**************************************************************************/
scripts::~scripts(){
	script_file.close();
}
/*************************************************************************/
void scripts::write_r_dos(vector<double>& energies){
	
	string Name = file_name;
	Name 		+= ".DOS";
	std::ofstream dos_file( Name.c_str() );
	dos_file << "Energy\n";
	for( unsigned i=0; i<energies.size(); i++ ) { dos_file << energies[i] << endl; }
	
	m_log->input_message("Outputing Density of States information to file.");
	
	script_file	 << "require(ggplot2) \n"
				 << "dos = read.table( '" << Name << "',header=T)\n"
				 << "attach(dos) \n"
				 << "p <-ggplot(dos, aes( x=Energy) )+\n"
				 << "\tgeom_density(fill='blue',bw=1) \n"
				 << "png('"<< Name << ".png',width = 5, height = 3.5, units ='"  << "in', res = 400)\n"
				 << "p\ndev.off()";
}
/****************************************************************************/
void scripts::write_pymol_cube(local_rd& lrdVol){
	
	if ( lrdVol.lrds[5].voxelN > 0 ){
		string pdb_name = file_name + ".pdb";
		std::string typestr;
		std::string typestr2;
		if ( lrdVol.FD ){ 
			typestr = lrdVol.name+"FD_";
			typestr2 = get_file_name( typestr.c_str() );
		}
		else{ 
			typestr = lrdVol.name+"FOA_";
			typestr2 = get_file_name( typestr.c_str() );
		}
		
		double mean,max,sum,min = 0.0;
		
		//isos used in the fukui functions and MOs
		double iso1 = 0.0005;
		double iso2 = 0.005;
		double iso3 = -0.0005;
		double iso4 = -0.005;
		
		//iso for LH LCP
		lrdVol.lrds[9].get_cube_stats(mean,min,max);
		double iso5 = mean*2;
		double iso6 = iso5*30;
		
		//iso for LH Vee
		lrdVol.lrds[10].get_cube_stats(mean,min,max);
		double iso7 = mean;
		double iso8 = iso7*10;
		
		//iso for potential Fukui
		lrdVol.lrds[11].get_cube_stats(mean,min,max);
		double iso9 = mean;
		double iso10 = iso9*10;
		
		
		//iso for LH int 
		lrdVol.lrds[14].get_cube_stats(mean,min,max);
		double iso11 = mean;
		double iso12 = iso11*10;
		
		
		lrdVol.lrds[15].get_cube_stats(mean,min,max);
		double iso13 = min*0.1;
		double iso14 = min*0.3;
		double iso15 = max*0.1;
		double iso16 = max*0.3;
		
		//iso for right potential Fukui
		lrdVol.lrds[12].get_cube_stats(mean,min,max);
		double iso17 = mean*4;
		double iso18 = iso17*10;
		
		//iso for zero potential Fukui
		lrdVol.lrds[13].get_cube_stats(mean,min,max);
		double iso19 = mean*2;
		double iso20 = iso19*15;
		
		//iso for 
		lrdVol.lrds[16].get_cube_stats(mean,min,max);
		double iso21 = mean;
		double iso22 = iso21*10;
		
		//iso for 
		lrdVol.lrds[17].get_cube_stats(mean,min,max);
		double iso23 = 0.0002;
		double iso24 = iso23*5;
		
		//iso for multiphilicity
		lrdVol.lrds[18].get_cube_stats(mean,min,max);
		double iso25 = -max*0.3;
		double iso26 = -max*0.6;
		double iso27 = max*0.3;
		double iso28 = max*0.6;
		
		//iso for MEP
		lrdVol.lrds[19].get_cube_stats(mean,min,max);
		double iso29 = 0.05;
		double iso30 = 0.5;
		double iso31 = -0.5;
		double iso32 = -0.05;
		
		//iso for 
		lrdVol.lrds[20].get_cube_stats(mean,min,max);
		double iso33 = max*0.01;
		double iso34 = max*0.2;
		
		
		script_file	<< "load "	<< pdb_name									<< " \n"
					<< "load "	<< typestr << lrdVol.lrds[0].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[0].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[1].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[1].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[5].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[6].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[7].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[8].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[8].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[9].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[10].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[11].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[14].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[15].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[15].name			<< "_ph2.cube\n"
					<< "volume " << typestr2 << lrdVol.lrds[0].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[0].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[0].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[0].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[1].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[1].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[1].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[1].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[5].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[5].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[6].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[6].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[7].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[7].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[8].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[8].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[8].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[8].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[9].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[9].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[10].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[10].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[11].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[11].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[14].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[14].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[15].name		<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[15].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[15].name		<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[15].name	<<	"_ph2"	<< " \n" 
					<< "volume_color " << typestr2 << lrdVol.lrds[0].name	<< "_ph1_volume, "	<< iso1		<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[0].name	<< "_ph2_volume, "	<< iso4		<< " red 0.04 "		<< iso3		<< " orange 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[1].name	<< "_ph1_volume, "	<< iso1 	<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[1].name	<< "_ph2_volume, "	<< iso4		<< " red 0.04 "		<< iso3		<< " orange 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[5].name	<< "_volume, "		<< iso1		<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[6].name	<< "_volume, "		<< iso1		<< " pink 0.02 "	<< iso2		<< " red  0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[7].name	<< "_volume, "		<< iso1		<< " yellow 0.02 "	<< iso2		<< " green 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[8].name	<< "_ph1_volume, "	<< iso1		<< " pink 0.03 "	<< iso2		<< " red 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[8].name	<< "_ph2_volume, "	<< iso4		<< " blue 0.05 "	<< iso3		<< " cyan 0.02 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[9].name	<< "_volume, "		<< iso5		<< " limon 0.02 "	<< iso6		<< " purpleblue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[10].name	<< "_volume, "		<< iso7		<< " brightorange 0.001 "	<< iso8	<< " purple 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[11].name	<< "_volume, "		<< iso9		<< " aquamarine 0.002 "		<< iso10 << " purple 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[14].name	<< "_volume, "		<< iso11	<< " salmon 0.02 "			<< iso12	<< " slate 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[15].name	<< "_ph1_volume, "	<< iso14	<< " skyblue  0.04 "		<< iso13	<<" greencyan 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[15].name	<< "_ph2_volume, "	<< iso15	<< " yelloworange 0.02 "	<< iso16	<< " warmpink 0.05 \n";
		
		if ( extra_RD ){
			script_file	<< "load "	<< typestr << lrdVol.lrds[12].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[13].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[16].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[17].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[18].name			<< "_ph1.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[18].name			<< "_ph2.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[19].name			<< "_ph1.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[19].name			<< "_ph2.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[20].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[21].name			<< ".cube\n"
						<< "volume " << typestr2 << lrdVol.lrds[12].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[12].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[13].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[13].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[16].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[16].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[17].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[17].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[18].name		<< "_ph1_volume, "	<< typestr2	<< lrdVol.lrds[18].name	<< "_ph1"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[18].name		<< "_ph2_volume, "	<< typestr2	<< lrdVol.lrds[18].name	<< "_ph2"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[19].name		<< "_ph1_volume, "	<< typestr2	<< lrdVol.lrds[19].name	<< "_ph1"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[19].name		<< "_ph2_volume, "	<< typestr2	<< lrdVol.lrds[19].name	<< "_ph2"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[20].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[20].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[21].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[21].name	<< " \n" 
						<< "volume_color " << typestr2 << lrdVol.lrds[12].name	<< "_volume, "		<< iso17	<< " cyan 0.0 "				<< iso18	<< " blue 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[13].name	<< "_volume, "		<< iso19	<< " yellow 0.0 "			<< iso20	<< " orange 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[16].name	<< "_volume, "		<< iso21 	<< " cyan 0.0 "				<< iso22	<< " blue 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[17].name	<< "_volume, "		<< iso23	<< " yellow 0.02 "			<< iso24	<< " orange 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[18].name	<< "_ph1_volume, "	<< iso26	<< " blue 0.04 "			<< iso25	<< " cyan 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[18].name	<< "_ph2_volume, "	<< iso27	<< " red 0.05 "				<< iso28	<< " pink  0.02 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[19].name	<< "_ph1_volume, "	<< iso29	<< " pink 0.0 "				<< iso30	<< " red  0.04 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[19].name	<< "_ph2_volume, "	<< iso31	<< " blue 0.04 "			<< iso32 << " aquamarine  0.0 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[20].name	<< "_volume, "		<< iso33	<< " greencyan 0.02 "		<< iso34	<< " green  0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[21].name	<< "_volume, "		<< iso7		<< " brightorange 0.001 "	<< iso8 	<< " purple 0.05 \n";
		}
	}
	
	
}
/**************************************************************************/
void scripts::write_pymol_pdb(){
	
	
	//std::ofstream script_f;
	string fname = get_file_name( file_name.c_str() );
	fname += ".pym";
	//script_f.open( fname.c_str() );

	string fname2 = fname.substr( 0,fname.size()-4 );

	string load_pdb_basic = "load "+ fname2 + "_PDB_RD/" + fname2;
	
	script_file	<< load_pdb_basic  << "_nucleophilicity.pdb\n"	//0
				<< load_pdb_basic  << "_electrophilicity.pdb\n"	//1
				<< load_pdb_basic  << "_radicality.pdb\n"		//2
				<< load_pdb_basic  << "_netphilicity.pdb\n"		//3
				<< load_pdb_basic  << "_hardness_Vee.pdb\n"		//4
				<< load_pdb_basic  << "_hardness_lcp.pdb\n"		//5
				<< load_pdb_basic  << "_fukui_pot_left.pdb\n"	//6
				<< load_pdb_basic  << "_fukui_pot_right.pdb\n"	//7
				<< load_pdb_basic  << "_fukui_pot_zero.pdb\n"	//8
				<< load_pdb_basic  << "_softness_dual.pdb\n"	//9
				<< load_pdb_basic  << "_hyper_softness.pdb\n"	//10
				<< load_pdb_basic  << "_fukushima.pdb\n"		//11
				<< load_pdb_basic  << "_mep.pdb\n"				//12
				<< load_pdb_basic  << "_hardness_TFD.pdb\n"		//13
				<< load_pdb_basic  << "_softness_avg.pdb\n"		//14
				<< load_pdb_basic  << "_hardness_int.pdb\n"		//15
				<< load_pdb_basic  << "_multiphilicity.pdb\n"	//16
				<< load_pdb_basic  << "_charge.pdb\n"			//17
				<< load_pdb_basic  << "_electron_density.pdb\n"	//18
				<< "spectrum b, blue_white_red, minimum=-0.5, maximum=0.5\n"
				<< "spectrum b, white_yellow_orange_red_black, minimum=0.1, maximum=0.5\n"
				<< "spectrum b, white_cyan_blue, minimum=0, maximum=0.1\n"
				<< "spectrum b, white_magenta_purple_purpleblue_black, minimum=1, maximum=5\n"
				<< "spectrum b, white_pink_red, minimum=0, maximum=0.1\n";
}
/***************************************************************************/
void scripts::write_r_heatmap(vector< vector<double> > rd_numerical	, 
							vector<string> rds						,
							vector<string> residues)				{
								
	string sname = file_name;
	sname+="residuesRD_4_R";
	std::ofstream data_txt(sname.c_str() );
	
	data_txt << "residue ";
	for(unsigned i=0;i<rds.size();i++){
		data_txt << rds[i] << " ";
	}
	data_txt << "\n";
	
	for(unsigned i=0;i<residues.size();i++){
		data_txt << residues[i] << " ";
		for(unsigned j=0;j<rds.size();j++){
			data_txt << rd_numerical[i][j] << " "; 
		}
		data_txt << "\n";
	}
	
	script_file << "library(pheatmap)\n"
				<< "data1 <-read.table('" << sname << "',header=T)\n"
				<< "attach(data1)\n"
				<< "ras_s <-scale(RAS,scale=T,center=F)\n"
				<< "hard_s <-scale(Hardness,scale=T,center=F)\n"
				<< "d = ras_s + hard_s\n"
				<< "d_vec <-c(d)\n"
				<< "data2 <-data.frame(data1,d_vec)\n"
				<< "detach(data1)\n"
				<< "attach(data2)\n"
				<< "data3 <-subset(data2,d_vec>1.5*mean(d_vec))\n"
				<< "detach(data2)\n"
				<< "attach(data3)\n"
				<< "pro1 <-data.matrix(data3[2:7])\n"
				<< "pro1 <-scale(pro1,scale=T,center=F)\n"
				<< "rownames(pro1) <-res\n"
				<< "detach(data3)\n"
				<< "pheatmap(pro1,color = colorRampPalette(c('navy','white','red'))(100),border_color=NA,cluster_cols=F,fontsize=8,main='\n',filenam='"
				<< sname << "heatmap.png'";
}
/*********************************************************************/
void scripts::write_r_residuos_barplot(){
	script_file << "require(ggpubr)\n"
				<< "data_stat <-read.table('residues_data_stat',header=T) \n"
				<< "avg <-subset(data_stat,type=='AVG')\n"
				<< "sd <-subset(data_stat,type=='SD')\n"
				<< "#-------------------------------------\n" 
				<< "p1 <-ggplot( avg, aes(x=res, y=Nucleophilicity) ) +\n"
				<< "geom_bar(stat='identity', color='black',fill='blue', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Nucleophilicity')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Nucleophilicity-sd$Nucleophilicity,ymax=Nucleophilicity+sd$Nucleophilicity), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p2 <-ggplot(avg, aes(x=res, y=Electrophilicity)) +\n" 
				<< "\tgeom_bar(stat='identity', color='black',fill='red', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Electrophilicity')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Electrophilicity-sd$Electrophilicity,ymax=Electrophilicity+sd$Electrophilicity), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p3 <-ggplot(avg, aes(x=res, y=Netphilicity)) +\n"
				<< "\tgeom_bar(stat='identity', color='black',fill='lightgreen', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Netphilicity')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Netphilicity-sd$Netphilicity,ymax=Netphilicity+sd$Netphilicity), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p4 <-ggplot(avg, aes(x=res, y=Hardness_Vee)) +\n"
				<< "\tgeom_bar(stat='identity', color='black',fill='orange', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Hardness Vee')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Hardness_Vee-sd$Hardness_Vee,ymax=Hardness_Vee+sd$Hardness_Vee), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p5 <-ggplot(avg, aes(x=res, y=Hardness_LCP)) +\n"
				<< "\tgeom_bar(stat='identity', color='black',fill='aquamarine3', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Hardness LCP')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Hardness_LCP-sd$Hardness_LCP,ymax=Hardness_LCP+sd$Hardness_LCP), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p6 <-ggplot(avg, aes(x=res, y=Fukui_pot_left)) +\n"
				<< "\tgeom_bar(stat='identity', color='black',fill='blueviolet', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Fukui Potential Left')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Fukui_pot_left-sd$Fukui_pot_left,ymax=Fukui_pot_left+sd$Fukui_pot_left), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p7 <-ggplot(avg, aes(x=res, y=Electron_Density)) +\n"
				<< "\tgeom_bar(stat='identity', color='black',fill='springgreen3', position=position_dodge() ) +\n"
				<< "\ttheme_minimal()+\n"
				<< "\tylab('Electron Density')+\n"
				<< "\txlab('Residues')+\n"
				<< "\tgeom_errorbar(aes(ymin=Electron_Density-sd$Electron_Density,ymax=Electron_Density+sd$Electron_Density), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "png('Nph_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p1\n"
				<< "dev.off()\n"
				<< "png('Eph_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p2\n"
				<< "dev.off()\n"
				<< "png('NET_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p3\n"
				<< "dev.off()\n"
				<< "png('Hardness_Vee_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p4\n"
				<< "dev.off()\n"
				<< "png('Hardness_LCP_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p5\n"
				<< "dev.off()\n"
				<< "png('Fukui_pot_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p6\n"
				<< "dev.off()\n"
				<< "png('Electron_dens_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p7\n"
				<< "dev.off()\n"
				<< "#=====================================================\n"
				<< "dat <-read.table('residues_data_frames',header=T)\n"
				<< "#-------------------------------------\n"	
				<< "pm1<-ggplot(dat, aes(x = frame, y = Electrophilicity) )+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Electrophiliicty')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm2<-ggplot(dat, aes(x = frame, y = Nucleophilicity ) )+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Nucleophilicity')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm3<-ggplot(dat, aes(x = frame, y =Netphilicity ))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Netphilicity')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm4<-ggplot(dat, aes(x = frame, y =Hardness_Vee ))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Hardness Vee')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm5<-ggplot(dat, aes(x = frame, y =Hardness_LCP ))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Hardness LCP')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm6<-ggplot(dat, aes(x = frame, y =Fukui_pot_left))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Fukui Potential')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm7<-ggplot(dat, aes(x = frame, y =hyper_softness ))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('hyper_softness')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm8<-ggplot(dat, aes(x = frame, y =Electron_Density ))+\n"
				<< "\ttheme_minimal()+\n"
				<< "\tgeom_point(aes(color=res))+\n"
				<< "\tgeom_smooth(aes(color=res),method='loess')+\n"
				<< "\tylab('Electron Density')+\n"
				<< "\txlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "wdht = 5\n"
				<< "png('nucleophilicity_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm1\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('electrophilicity_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm2\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"				
				<< "png('net_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm3\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('hard_Vee_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm4\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('hard_LCP_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm5\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('Fukui_pot_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm6\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"				
				<< "png('Hsoftness_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm7\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('ED_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm8\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n";			
}
/*********************************************************************/
void scripts::write_r_reaction_analysis(traj_rd& path_rd			,
										ReactionAnalysis& r_info	,
										string& nameb				){
											
	
	vector<string> gl_names		= {"HOF","ECP","Hardness","Softness","Electrophilicity","Energy"};
	vector<string> gl_legends	= {"HOF \\n(kCal/mol)",
									"ECP (eV)",
									"Hardness \\n(eV)",
									"Softness \\n(eV)",
									"Electrophilicity \\n(eV)",
									"Energy \\n(eV)"};
	
	string delim = "#====================================================\n";
	string delim2 = "#---------------------------------------------------\n";
	
	script_file		<< delim
					<< "require(ggpubr)\n"
					<< "jet.colors <-colorRampPalette(c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow','#FF7F00','red','#7F0000'))\n\n"
					<< delim
					<< "atom_lrd='" << nameb <<  "'\n"
					<< "rc1_name='" << r_info.RCs[0].rc_label << "'\n"
					<< "df1 <-read.table(atom_lrd,header=T)\n";
					if ( r_info.nrcs == 2 ){
						script_file << "rc2_name='" << r_info.RCs[1].rc_label << "'\n";
					}
	script_file << delim;

	string pr_obj;
	
	unsigned int a = 0;
	unsigned int c = 0;	
	
	if ( r_info.ndim == 1 ){
		for( unsigned i=0; i<gl_names.size(); i++){
			script_file << "grc1_" << to_string(i) << " <-ggplot(df1,aes(x=RC1,y="
						<< gl_names[i] << ") )+\n "
						<< "\tgeom_point()+ \n"
						<< "\tgeom_line()+ \n"
						<< "\ttheme_minimal()+ \n"
						<< "\tylab('" << gl_legends[i] << "') + \n"
						<< "\txlab(rc1_name)\n\n"
						<< delim2;
		}
	
		script_file << "png('global_rc1.png',width=6.5,height=5,units='in',res=1000)\n";
		script_file << "ggarrange(grc1_0,grc1_1,grc1_2,grc1_3,grc1_4,grc1_5,ncol=3,nrow=2)\n";
		script_file << "dev.off()\n";
		script_file << delim;
		
		for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
			pr_obj = "la_";
			script_file<< pr_obj
						<< to_string(i) << " <-ggplot(df1,aes(x=RC1,y="
						<< path_rd.rds_labels[i] << ") )+\n "
						<< "\tgeom_point()+ \n"
						<< "\tgeom_line()+ \n"
						<< "\ttheme_minimal()+ \n"
						<< "\tylab('" << path_rd.atoms_labels[i] << "') + \n"
						<< "\txlab(rc1_name)  \n"
						<< delim2;
			c++;
			if( c % 13 == 0 ){
				script_file << "png('" <<pr_obj << to_string(a++) << ".png',width=8,height=7,units='in',res=1000)\n";
				script_file << "ggarrange(" 
							<< pr_obj << to_string(c-13) << ","
							<< pr_obj << to_string(c-12) << ","
							<< pr_obj << to_string(c-11) << ","
							<< pr_obj << to_string(c-10) << ","
							<< pr_obj << to_string(c-9) << ","
							<< pr_obj << to_string(c-8) << ","
							<< pr_obj << to_string(c-7) << ","
							<< pr_obj << to_string(c-6) << ","
							<< pr_obj << to_string(c-5) << ","
							<< "ncol=3,nrow=3)\n";
				script_file << "dev.off()\n";
				script_file << delim;
			}		
		}	
		c = 0;
		if ( r_info.nrcs == 2 ){
			c = 0;
			a = 0;
			for( unsigned i=0;i<gl_names.size(); i++){
				script_file << " grc2_" << to_string(i) << " <-ggplot(df1,aes(x=RC2,y="
							<< gl_names[i] << ") )+\n "
							<< "\tgeom_point()+ \n"
							<< "\tgeom_line()+ \n"
							<< "\ttheme_minimal()+ \n"
							<< "\tylab('" << gl_legends[i] << "') + \n"
							<< "\txlab(rc2_name) \n"
							<< delim2;
			}
			script_file << "png('global_rc2.png',width=6.5,height=5,units='in',res=1000)\n"
						<< "ggarrange(grc2_0,grc2_1,grc2_2,grc2_3,grc2_4,grc2_5,ncol=3,nrow=2)\n"
						<< "dev.off()\n"
						<< delim;
		
			for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
				pr_obj = "la_";
				script_file << pr_obj
							<< to_string(i) << " <-ggplot(df1,aes(x=RC2,y="
							<< path_rd.rds_labels[i] << ") )+\n "
							<< "\tgeom_point()+ \n"
							<< "\tgeom_line()+ \n"
							<< "\ttheme_minimal()+ \n"
							<< "\tylab('" << path_rd.atoms_labels[i] << "') + \n"
							<< "\txlab(rc2_name)  \n"
							<< delim2;
				c++;
				if( c % 13 == 0 ){
					script_file	<< "png('" <<pr_obj << "rc2_" << to_string(a++) <<".png',width=8,height=7,units='in',res=1000)\n";
					script_file	<< "ggarrange(" 
								<< pr_obj << to_string(c-13) << ","
								<< pr_obj << to_string(c-12) << ","
								<< pr_obj << to_string(c-11) << ","
								<< pr_obj << to_string(c-10) << ","
								<< pr_obj << to_string(c-9) << ","
								<< pr_obj << to_string(c-8) << ","
								<< pr_obj << to_string(c-7) << ","
								<< pr_obj << to_string(c-6) << ","
								<< pr_obj << to_string(c-5) << ","					 
								<< "ncol=3,nrow=3)\n";					 
					script_file<< "dev.off()\n";		
					script_file<< delim;
				}		
			}	
			c = 0;
		}			
	}else if ( r_info.ndim == 2 ){
		for( unsigned i=0; i<gl_names.size(); i++){
			script_file << "grc1_" << to_string(i) << " <-ggplot(df1,aes(x=rc1,y=rc2,z="
						<< gl_names[i] << ") )+\n "
						<< "\tstat_contour(geom='polygon', aes(fill = ..level..))+ \n" 
						<< "\tgeom_tile(aes(fill =" << gl_names[i] << " ))+ \n"
						<< "\tstat_contour(bins = 6,color='black')+ \n"
						<< "\tscale_fill_gradientn(colours=jet.colors(70))+ \n"
						<< "\tylab(rc1_name) + \n"
						<< "\txlab(rc2_name) + \n"
						<< "\tguides(fill = guide_colorbar('" << gl_legends[i] <<"'))\n"
						<< delim2;
		}
	
		script_file << "png('global_2d.png',width=8,height=8,units='in',res=1000)\n";
		script_file << "ggarrange(grc1_0,grc1_1,grc1_2,grc1_3,grc1_4,grc1_5,ncol=2,nrow=3)\n";
		script_file << "dev.off()\n";
		script_file << delim;
		
		for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
				pr_obj = "la_";
				script_file << pr_obj
							<< to_string(i) << " <-ggplot(df1,aes(x=rc1,y=rc2,z="
							<< path_rd.rds_labels[i] << ") )+\n "
							<< "\tstat_contour(geom='polygon', aes(fill = ..level..))+ \n" 
							<< "\tstat_contour(bins = 6,color='black')+ \n"
							<< "\tscale_fill_gradientn(colours=jet.colors(7))+\n"
							<< "\tgeom_tile(aes(fill =" << path_rd.rds_labels[i] << "))+\n"
							<< "\tylab(rc2_name) + \n"
							<< "\txlab(rc1_name)  +\n"
							<< "\tguides(fill = guide_colorbar(title ='" << path_rd.atoms_labels[i] << "'))\n"
							<< delim2;
				c++;
			if( c % 13 == 0 ){
				script_file	<< "png('" <<pr_obj << "rc2_" << to_string(a++) <<".png',width=11,height=7,units='in',res=1000)\n";
				script_file	<< "ggarrange(" 
							<< pr_obj << to_string(c-13) << ","
							<< pr_obj << to_string(c-12) << ","
							<< pr_obj << to_string(c-11) << ","
							<< pr_obj << to_string(c-10) << ","
							<< pr_obj << to_string(c-9) << ","
							<< pr_obj << to_string(c-8) << ","
							<< pr_obj << to_string(c-7) << ","
							<< pr_obj << to_string(c-6) << ","
							<< pr_obj << to_string(c-5) << ","
							<< "ncol=3,nrow=3)\n";
					script_file<< "dev.off()\n";
					script_file<< delim;
			}		
		}
	}
}
/***************************************************************************
 * Writes an R script for PCA and KRLS analyses of the final
 * protein-ligand descriptor matrix.
 *
 * The input matrix must already contain the summed descriptor differences
 * calculated by PRIMoRDiA.
 *
 * Expected metadata columns:
 *
 *     frame : complex identifier
 *     Type  : categorical activity/class label
 *     DG    : experimental binding free energy
 *
 * All remaining numerical columns are treated as molecular descriptors.
 ***************************************************************************/
void scripts::write_r_complex_analysis(
    const std::string& matrix_file,
    const std::string& output_prefix
){
    m_log->input_message(
        "Writing PCA and KRLS protein-ligand analysis R script."
    );

    /*
     * Configuration transferred from C++ to the generated R script.
     *
     * The names can be edited directly in the generated R file if the
     * matrix header changes later.
     */
    script_file
        << "\n"
        << "#============================================================\n"
        << "# PRIMoRDiA protein-ligand statistical analysis\n"
        << "# PCA and Kernel Regularized Least Squares\n"
        << "#============================================================\n\n"
        << "input_file <- \"" << matrix_file << "\"\n"
        << "output_prefix <- \"" << output_prefix << "\"\n\n";

    script_file << R"PRIMORDIA_R(


#============================================================
# Analysis configuration
#============================================================

id_column <- "frame"
pdb_column <- "pdb_code"
class_column <- "Type"
response_column <- "DG"

residue_column <- "res"
residue_type_column <- "res_typ"

# Number of principal components calculated and displayed.
pca_components <- 5

# Number of descriptors selected automatically for the second PCA.
pca2_number_variables <- 7

# KRLS regularization parameter.
krls_lambda <- 0.1

manual_pca2_variables <- character(0)

#============================================================
# DG classification limits
#
# DG is expressed in kcal/mol.
# More negative values indicate stronger binding.
#
# DG <= -10             Very Strong
# -10 < DG <= -8        Strong
# -8  < DG <= -6        Medium
# DG > -6               Weak
#============================================================

very_strong_limit <- -10.0
strong_limit <- -8.0
medium_limit <- -6.0

#============================================================
# Packages
#============================================================

library(ggplot2)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(dplyr)
library(caret)
library(KRLS)

#============================================================
# Read the residue descriptor-difference matrix
#============================================================

residue_data <- read.table(
    input_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

#============================================================
# Extract the PDB code from the frame name
#
# Example:
#
# complex_5NXG_residues -> 5NXG
#============================================================

residue_data[[pdb_column]] <- toupper(
    sub(
        "^complex_([[:alnum:]]{4})_residues.*$",
        "\\1",
        residue_data[[id_column]]
    )
)

#============================================================
# Identify the 19 numerical descriptor columns
#============================================================

metadata_columns <- c(
    id_column,
    pdb_column,
    residue_column,
    residue_type_column
)

descriptor_names <- setdiff(
    names(residue_data),
    metadata_columns
)

descriptor_names <- descriptor_names[
    vapply(
        residue_data[
            ,
            descriptor_names,
            drop = FALSE
        ],
        is.numeric,
        logical(1)
    )
]

#============================================================
# Sum residue descriptor differences for each complex
#
# The input contains one line per selected residue.
# The resulting table contains one line per complex.
#============================================================

complex_descriptor_data <- aggregate(
    residue_data[
        ,
        descriptor_names,
        drop = FALSE
    ],
    by = list(
        frame = residue_data[[id_column]],
        pdb_code = residue_data[[pdb_column]]
    ),
    FUN = sum,
    na.rm = TRUE
)

# Save the matrix produced after summing the selected residues.

write.table(
    complex_descriptor_data,
    file = paste0(
        output_prefix,
        "_summed_descriptor_matrix.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#============================================================
# Read experimental DG
#
# File structure:
#
# 5NXG  -11.700
# 5NXI   -9.100
#
# Lines beginning with # are ignored.
#============================================================

experimental_data <- read.table(
    experimental_file,
    header = FALSE,
    comment.char = "#",
    stringsAsFactors = FALSE,
    col.names = c(
        pdb_column,
        response_column
    )
)

experimental_data[[pdb_column]] <- toupper(
    trimws(
        experimental_data[[pdb_column]]
    )
)

experimental_data[[response_column]] <- as.numeric(
    experimental_data[[response_column]]
)

#============================================================
# Join calculated descriptors and experimental DG
#============================================================

original_pdb_order <- complex_descriptor_data[[pdb_column]]

analysis_data <- merge(
    complex_descriptor_data,
    experimental_data,
    by = pdb_column,
    all.x = TRUE,
    sort = FALSE
)

# merge() may change the row order. Restore the original matrix order.

analysis_data <- analysis_data[
    match(
        original_pdb_order,
        analysis_data[[pdb_column]]
    ),
    ,
    drop = FALSE
]

rownames(analysis_data) <- NULL

#============================================================
# Create the affinity category
#============================================================

analysis_data[[class_column]] <- cut(
    analysis_data[[response_column]],
    breaks = c(
        -Inf,
        very_strong_limit,
        strong_limit,
        medium_limit,
        Inf
    ),
    labels = c(
        "Very Strong",
        "Strong",
        "Medium",
        "Weak"
    ),
    right = TRUE,
    ordered_result = TRUE
)

# Remove categories that are not present in this particular family.

analysis_data[[class_column]] <- droplevels(
    analysis_data[[class_column]]
)

#============================================================
# Save descriptors, DG and classification
#============================================================

write.table(
    analysis_data,
    file = paste0(
        output_prefix,
        "_analysis_matrix.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#============================================================
# Show the association in the terminal
#============================================================

print(
    analysis_data[
        ,
        c(
            id_column,
            pdb_column,
            response_column,
            class_column
        )
    ]
)

print(
    table(
        analysis_data[[class_column]],
        useNA = "ifany"
    )
)

#============================================================
# Prepare descriptor data for PCA and KRLS
#
# DG must not be included among the predictors.
#============================================================

numeric_columns <- names(analysis_data)[
    vapply(
        analysis_data,
        is.numeric,
        logical(1)
    )
]

descriptor_names <- setdiff(
    numeric_columns,
    response_column
)

descriptor_data <- analysis_data[
    ,
    descriptor_names,
    drop = FALSE
]

# Remove descriptors with zero variance.

near_zero_variance <- nearZeroVar(
    descriptor_data,
    saveMetrics = TRUE
)

selected_descriptor_names <- rownames(
    near_zero_variance
)[
    !near_zero_variance$zeroVar
]

descriptor_data <- descriptor_data[
    ,
    selected_descriptor_names,
    drop = FALSE
]

descriptor_names <- colnames(
    descriptor_data
)

analysis_class <- analysis_data[[class_column]]
experimental_dg <- analysis_data[[response_column]]


descriptor_data <- analysis_data[
    ,
    descriptor_names,
    drop = FALSE
]

# Remove constant or nearly constant variables because they cannot be
# meaningfully standardized for PCA or KRLS.

near_zero_variance <- nearZeroVar(
    descriptor_data,
    saveMetrics = TRUE
)

selected_descriptor_names <- rownames(
    near_zero_variance
)[
    !near_zero_variance$zeroVar
]

descriptor_data <- descriptor_data[
    ,
    selected_descriptor_names,
    drop = FALSE
]

descriptor_names <- colnames(descriptor_data)

# Save the descriptors effectively included in the analysis.

write.table(
    data.frame(
        descriptor = descriptor_names
    ),
    file = paste0(
        output_prefix,
        "_descriptors_used.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#============================================================
# 1. Descriptor distributions
#============================================================

create_descriptor_visualization <- function(
    data,
    descriptor,
    class_name
) {

    strip_plot <- ggplot(
        data,
        aes(
            x = .data[[class_name]],
            y = .data[[descriptor]],
            color = .data[[class_name]]
        )
    ) +
        geom_jitter(
            width = 0.15,
            height = 0,
            alpha = 0.70,
            size = 2
        ) +
        stat_summary(
            fun = mean,
            geom = "point",
            color = "black",
            size = 3,
            shape = 18
        ) +
        theme_minimal() +
        labs(
            x = "Type",
            y = descriptor,
            color = "Type"
        )

    histogram_plot <- ggplot(
        data,
        aes(
            x = .data[[descriptor]]
        )
    ) +
        geom_histogram(
            bins = 10,
            color = "black",
            fill = "grey75"
        ) +
        geom_vline(
            xintercept = mean(
                data[[descriptor]],
                na.rm = TRUE
            ),
            linetype = "dashed"
        ) +
        theme_minimal() +
        labs(
            x = descriptor,
            y = "Count"
        )

    ggarrange(
        strip_plot,
        histogram_plot,
        ncol = 2,
        nrow = 1
    )
}

for (descriptor in descriptor_names) {

    descriptor_plot <- create_descriptor_visualization(
        analysis_data,
        descriptor,
        class_column
    )

    ggsave(
        filename = paste0(
            output_prefix,
            "_",
            descriptor,
            "_distribution.png"
        ),
        plot = descriptor_plot,
        width = 8,
        height = 4,
        dpi = 600
    )
}

#============================================================
# 2. Descriptor correlation matrix
#============================================================

descriptor_correlation <- cor(
    descriptor_data,
    use = "pairwise.complete.obs"
)

png(
    paste0(
        output_prefix,
        "_descriptor_correlations.png"
    ),
    width = 2400,
    height = 2200,
    units = "px",
    res = 300
)

corrplot(
    descriptor_correlation,
    method = "color",
    type = "upper",
    order = "hclust",
    tl.cex = 0.65,
    tl.col = "black",
    diag = FALSE
)

dev.off()

write.table(
    descriptor_correlation,
    file = paste0(
        output_prefix,
        "_descriptor_correlations.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

#============================================================
# 3. First PCA: all nonconstant numerical descriptors
#============================================================

number_pca_components <- min(
    pca_components,
    ncol(descriptor_data),
    nrow(descriptor_data) - 1
)

pca_result <- PCA(
    descriptor_data,
    scale.unit = TRUE,
    graph = FALSE,
    ncp = number_pca_components
)

#------------------------------------------------------------
# 3.1 Eigenvalues and explained variance
#------------------------------------------------------------

pca_eigenvalues <- get_eigenvalue(
    pca_result
)

write.table(
    pca_eigenvalues,
    file = paste0(
        output_prefix,
        "_pca_eigenvalues.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

png(
    paste0(
        output_prefix,
        "_pca_explained_variance.png"
    ),
    width = 1600,
    height = 1200,
    units = "px",
    res = 300
)

print(
    fviz_eig(
        pca_result,
        addlabels = TRUE
    )
)

dev.off()

#------------------------------------------------------------
# 3.2 Variable coordinates, contributions and cos2
#------------------------------------------------------------

pca_variables <- get_pca_var(
    pca_result
)

write.table(
    pca_variables$coord,
    file = paste0(
        output_prefix,
        "_pca_variable_coordinates.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

write.table(
    pca_variables$contrib,
    file = paste0(
        output_prefix,
        "_pca_variable_contributions.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

write.table(
    pca_variables$cos2,
    file = paste0(
        output_prefix,
        "_pca_variable_cos2.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

png(
    paste0(
        output_prefix,
        "_pca_variables_cos2.png"
    ),
    width = 1800,
    height = 1500,
    units = "px",
    res = 300
)

print(
    fviz_pca_var(
        pca_result,
        col.var = "cos2",
        gradient.cols = c(
            "#00AFBB",
            "#E7B800",
            "#FC4E07"
        ),
        repel = TRUE
    )
)

dev.off()

png(
    paste0(
        output_prefix,
        "_pca_contribution_matrix.png"
    ),
    width = 1400,
    height = 1800,
    units = "px",
    res = 300
)

corrplot(
    pca_variables$contrib,
    is.corr = FALSE
)

dev.off()

#------------------------------------------------------------
# 3.3 PCA scores
#------------------------------------------------------------

pca_scores <- data.frame(
    frame = analysis_data[[id_column]],
    Type = analysis_data[[class_column]],
    pca_result$ind$coord,
    check.names = FALSE
)

write.table(
    pca_scores,
    file = paste0(
        output_prefix,
        "_pca_scores.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#------------------------------------------------------------
# 3.4 PCA biplots
#------------------------------------------------------------

png(
    paste0(
        output_prefix,
        "_pca_biplot_PC1_PC2.png"
    ),
    width = 1900,
    height = 1400,
    units = "px",
    res = 300
)

print(
    fviz_pca_biplot(
        pca_result,
        axes = c(1, 2),
        col.ind = analysis_data[[class_column]],
        palette = "jco",
        addEllipses = TRUE,
        label = "var",
        col.var = "black",
        repel = TRUE,
        legend.title = "Activity"
    )
)

dev.off()

if (number_pca_components >= 3) {

    png(
        paste0(
            output_prefix,
            "_pca_biplot_PC1_PC3.png"
        ),
        width = 1900,
        height = 1400,
        units = "px",
        res = 300
    )

    print(
        fviz_pca_biplot(
            pca_result,
            axes = c(1, 3),
            col.ind = analysis_data[[class_column]],
            palette = "jco",
            addEllipses = TRUE,
            label = "var",
            col.var = "black",
            repel = TRUE,
            legend.title = "Activity"
        )
    )

    dev.off()
}

#============================================================
# 4. Second PCA
#
# The original script selected columns by fixed numerical positions.
# Here, variables can be provided manually or selected automatically
# from the largest summed contributions to the first three PCs.
#============================================================

if (length(manual_pca2_variables) > 0) {

    pca2_descriptor_names <- manual_pca2_variables

} else {

    number_selection_axes <- min(
        3,
        ncol(pca_variables$contrib)
    )

    contribution_scores <- rowSums(
        pca_variables$contrib[
            ,
            seq_len(number_selection_axes),
            drop = FALSE
        ]
    )

    ordered_contributions <- sort(
        contribution_scores,
        decreasing = TRUE
    )

    number_selected_variables <- min(
        pca2_number_variables,
        length(ordered_contributions)
    )

    pca2_descriptor_names <- names(
        ordered_contributions
    )[
        seq_len(number_selected_variables)
    ]
}

pca2_data <- descriptor_data[
    ,
    pca2_descriptor_names,
    drop = FALSE
]

write.table(
    data.frame(
        descriptor = pca2_descriptor_names
    ),
    file = paste0(
        output_prefix,
        "_pca2_selected_descriptors.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

number_pca2_components <- min(
    pca_components,
    ncol(pca2_data),
    nrow(pca2_data) - 1
)

pca2_result <- PCA(
    pca2_data,
    scale.unit = TRUE,
    graph = FALSE,
    ncp = number_pca2_components
)

pca2_eigenvalues <- get_eigenvalue(
    pca2_result
)

write.table(
    pca2_eigenvalues,
    file = paste0(
        output_prefix,
        "_pca2_eigenvalues.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

png(
    paste0(
        output_prefix,
        "_pca2_explained_variance.png"
    ),
    width = 1600,
    height = 1200,
    units = "px",
    res = 300
)

print(
    fviz_eig(
        pca2_result,
        addlabels = TRUE
    )
)

dev.off()

pca2_variables <- get_pca_var(
    pca2_result
)

write.table(
    pca2_variables$contrib,
    file = paste0(
        output_prefix,
        "_pca2_variable_contributions.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

png(
    paste0(
        output_prefix,
        "_pca2_contribution_matrix.png"
    ),
    width = 1400,
    height = 1800,
    units = "px",
    res = 300
)

corrplot(
    pca2_variables$contrib,
    is.corr = FALSE
)

dev.off()

pca2_scores <- data.frame(
    frame = analysis_data[[id_column]],
    Type = analysis_data[[class_column]],
    pca2_result$ind$coord,
    check.names = FALSE
)

write.table(
    pca2_scores,
    file = paste0(
        output_prefix,
        "_pca2_scores.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

png(
    paste0(
        output_prefix,
        "_pca2_biplot_PC1_PC2.png"
    ),
    width = 1900,
    height = 1400,
    units = "px",
    res = 300
)

print(
    fviz_pca_biplot(
        pca2_result,
        axes = c(1, 2),
        col.ind = analysis_data[[class_column]],
        palette = "jco",
        addEllipses = TRUE,
        label = "var",
        col.var = "black",
        repel = TRUE,
        legend.title = "Activity"
    )
)

dev.off()

if (number_pca2_components >= 3) {

    png(
        paste0(
            output_prefix,
            "_pca2_biplot_PC1_PC3.png"
        ),
        width = 1900,
        height = 1400,
        units = "px",
        res = 300
    )

    print(
        fviz_pca_biplot(
            pca2_result,
            axes = c(1, 3),
            col.ind = analysis_data[[class_column]],
            palette = "jco",
            addEllipses = TRUE,
            label = "var",
            col.var = "black",
            repel = TRUE,
            legend.title = "Activity"
        )
    )

    dev.off()
}

#============================================================
# 5. KRLS regression
#
# The selected PCA2 descriptors are used as KRLS predictors.
# They are centered and scaled before fitting.
#============================================================

krls_predictor_data <- pca2_data

krls_preprocess <- preProcess(
    krls_predictor_data,
    method = c(
        "center",
        "scale"
    )
)

krls_scaled_data <- predict(
    krls_preprocess,
    krls_predictor_data
)

experimental_dg <- analysis_data[
    [
        response_column
    ]
]

krls_model <- krls(
    X = as.matrix(krls_scaled_data),
    y = experimental_dg,
    lambda = krls_lambda
)

capture.output(
    summary(krls_model),
    file = paste0(
        output_prefix,
        "_krls_summary.txt"
    )
)

png(
    paste0(
        output_prefix,
        "_krls_diagnostic.png"
    ),
    width = 1800,
    height = 1400,
    units = "px",
    res = 300
)

plot(krls_model)

dev.off()

fitted_dg <- as.numeric(
    krls_model$fitted
)

krls_predictions <- data.frame(
    frame = analysis_data[[id_column]],
    Type = analysis_data[[class_column]],
    Experimental_DG = experimental_dg,
    Fitted_DG = fitted_dg,
    Residual = experimental_dg - fitted_dg
)

write.table(
    krls_predictions,
    file = paste0(
        output_prefix,
        "_krls_fitted_values.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#------------------------------------------------------------
# 5.1 KRLS fit metrics
#
# These metrics describe the fit to the same data used to train the
# model. They are not yet external-validation metrics.
#------------------------------------------------------------

residual_sum_squares <- sum(
    (
        experimental_dg -
        fitted_dg
    )^2
)

total_sum_squares <- sum(
    (
        experimental_dg -
        mean(experimental_dg)
    )^2
)

r_squared_fit <- 1.0 - (
    residual_sum_squares /
    total_sum_squares
)

rmse_fit <- sqrt(
    mean(
        (
            experimental_dg -
            fitted_dg
        )^2
    )
)

mae_fit <- mean(
    abs(
        experimental_dg -
        fitted_dg
    )
)

spearman_fit <- cor(
    experimental_dg,
    fitted_dg,
    method = "spearman"
)

krls_metrics <- data.frame(
    metric = c(
        "R2_fit",
        "RMSE_fit",
        "MAE_fit",
        "Spearman_fit",
        "Lambda",
        "Number_samples",
        "Number_predictors"
    ),
    value = c(
        r_squared_fit,
        rmse_fit,
        mae_fit,
        spearman_fit,
        krls_lambda,
        nrow(krls_scaled_data),
        ncol(krls_scaled_data)
    )
)

write.table(
    krls_metrics,
    file = paste0(
        output_prefix,
        "_krls_metrics.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

#------------------------------------------------------------
# 5.2 Experimental versus fitted plot
#------------------------------------------------------------

krls_plot <- ggplot(
    krls_predictions,
    aes(
        x = Experimental_DG,
        y = Fitted_DG,
        color = Type
    )
) +
    geom_point(
        size = 3,
        alpha = 0.80
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "black",
        linewidth = 0.8
    ) +
    annotate(
        "text",
        x = Inf,
        y = -Inf,
        hjust = 1.1,
        vjust = -1.2,
        label = paste0(
            "R² fit = ",
            round(
                r_squared_fit,
                4
            ),
            "\nRMSE = ",
            round(
                rmse_fit,
                4
            ),
            "\nMAE = ",
            round(
                mae_fit,
                4
            )
        )
    ) +
    theme_minimal() +
    labs(
        title = "Experimental versus fitted values — KRLS",
        x = "Experimental binding free energy",
        y = "Fitted binding free energy",
        color = "Activity"
    )

ggsave(
    filename = paste0(
        output_prefix,
        "_krls_fitted.png"
    ),
    plot = krls_plot,
    width = 6,
    height = 5,
    dpi = 600
)

ggsave(
    filename = paste0(
        output_prefix,
        "_krls_fitted.pdf"
    ),
    plot = krls_plot,
    width = 6,
    height = 5
)

#============================================================
# 6. Analysis summary
#============================================================

cat("\n")
cat("PRIMoRDiA PCA and KRLS analysis completed.\n")
cat("Input matrix:", input_file, "\n")
cat("Number of complexes:", nrow(analysis_data), "\n")
cat("Number of descriptors in PCA1:", ncol(descriptor_data), "\n")
cat("Number of descriptors in PCA2/KRLS:", ncol(pca2_data), "\n")
cat("KRLS lambda:", krls_lambda, "\n")
cat("KRLS fitted R2:", r_squared_fit, "\n")
cat("KRLS fitted RMSE:", rmse_fit, "\n")
cat("KRLS fitted MAE:", mae_fit, "\n")
cat("\n")

)PRIMORDIA_R";
}


/********************************************************************************/

//================================================================================
//END OF FILE
//================================================================================

