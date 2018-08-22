
#include "../interface/ntuple_JetInfo.h"



ntuple_JetInfo::ntuple_JetInfo():ntuple_content(){
	edm::ParameterSet cfgPFJetIdAlgo;
	cfgPFJetIdAlgo.addParameter<std::string>("version", "FIRSTDATA");
	cfgPFJetIdAlgo.addParameter<std::string>("quality", "LOOSE");
	loosePFJetIdAlgo_ = new PFJetIDSelectionFunctor(cfgPFJetIdAlgo);
	clear();
}

ntuple_JetInfo::~ntuple_JetInfo(){
	delete loosePFJetIdAlgo_;
}

void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig){
}

void ntuple_JetInfo::initBranches(TTree* t){

    ADDBRANCH(t,recJet_pt);
    ADDBRANCH(t,recJet_eta);
    ADDBRANCH(t,recJet_phi);
    ADDBRANCH(t,recJet_mass);
    ADDBRANCH(t,recJetLooseId);

    ADDBRANCH(t,genJet_pt);
    ADDBRANCH(t,genJet_eta);
    ADDBRANCH(t,genJet_phi);
    ADDBRANCH(t,genJet_mass);

    ADDBRANCH(t,hasGenMatch);

    ADDBRANCH(t,hasElecMatch05);
    ADDBRANCH(t,hasMuonMatch05);
    ADDBRANCH(t,hasTauMatch05);
}
void ntuple_JetInfo::clear(){

    recJet_pt=0;
    recJet_eta=0;
    recJet_phi=0;
    recJet_mass=0;
    genJet_pt=0;
    genJet_eta=0;
    genJet_phi=0;
    genJet_mass=0;
    recJetLooseId=0;
    hasGenMatch=0;


    hasElecMatch05=0;
    hasMuonMatch05=0;
    hasTauMatch05=0;
}

bool ntuple_JetInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau,
		const std::vector<const reco::GenParticle*> * gen){
	if(!recJet){
		clear();
		return true;
	}
	if(gen){
		for(auto genParticle: *gen){
			if(!(genParticle->isLastCopy() && genParticle->statusFlags().fromHardProcess()))continue;
			if(reco::deltaR2(genParticle->p4(),recJet->p4()) > 0.5*0.5) continue;

			if(TMath::Abs(genParticle->pdgId()) == 11) hasElecMatch05++;
			if(TMath::Abs(genParticle->pdgId()) == 13) hasMuonMatch05++;
			if(TMath::Abs(genParticle->pdgId()) == 15) hasTauMatch05++;

		}
	}
	return fillBranches(recTau,recJet,genTau);
}
bool ntuple_JetInfo::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){
	if(!recJet){
		clear();
		return true;
	}
	recJet_pt=recJet->pt();
	recJet_eta=recJet->eta();
	recJet_phi=recJet->phi();
	recJet_mass=recJet->mass();

	recJetLooseId = ( (*loosePFJetIdAlgo_)(*recJet) ) ? 1 : 0;

	hasGenMatch=0;
	if(recJet->genJet()){
		hasGenMatch=1;
	    genJet_pt=recJet->genJet()->pt();
	    genJet_eta=recJet->genJet()->eta();
	    genJet_phi=recJet->genJet()->phi();
	    genJet_mass=recJet->genJet()->mass();

	}

	return true;
}
