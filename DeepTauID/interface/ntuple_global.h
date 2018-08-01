/*
 * ntuple_global.h
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_

#include "ntuple_content.h"

class ntuple_global: public ntuple_content{
public:

	ntuple_global();

    void getInput(const edm::ParameterSet& iConfig){};
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup){};
    //use either of these functions

    bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){return true;}

    void clear();


private:
	int event;
	int run;
	int lumi;
	float numPileUp;
	float rho;

};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_GLOBAL_H_ */
