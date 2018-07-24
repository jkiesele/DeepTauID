/*
 * ntuple_global.cc
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */


#include "../interface/ntuple_global.h"
#include "../interface/ntuple_content.h"

void ntuple_global::clear(){
	event=0;
	run=0;
	lumi=0;
	numPileUp=0;
	rho=0;
}

void ntuple_global::initBranches(TTree* t){

	ADDBRANCH(t, event);
	ADDBRANCH(t, run);
	ADDBRANCH(t, lumi);
	ADDBRANCH(t, numPileUp);
	ADDBRANCH(t, rho);
}
void ntuple_global::readEvent(const edm::Event& iEvent){

	event=iEvent.id().event();
	lumi=iEvent.id().luminosityBlock();
	run=iEvent.id().run();
	numPileUp=vertices()->size();
	if(rhoInfo())
		rho=rhoInfo()[0];
}

