/*
 * ntuple_content.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */


#include "../../../DeepTau/DeepTauID/interface/ntuple_content.h"

#include <stdexcept>


bool ntuple_content::useoffsets=true;
bool ntuple_content::debug=false;

ntuple_content::~ntuple_content(){


}



const double * ntuple_content::rhoInfo()const{
	return rhoInfo_;
}
