/*
 * ntuple_JetInfo.h
 *
 *  Created on: 20 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_

/**
 *
 * Contains global jet information
 *
 */
class ntuple_JetInfo: public ntuple_content{
public:
    ntuple_JetInfo():ntuple_content(){}


    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);
    //use either of these functions

    bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau);


private:




};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_JETINFO_H_ */
