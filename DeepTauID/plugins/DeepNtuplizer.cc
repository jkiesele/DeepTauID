// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include <TRandom.h>

//CMSSW includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"

#include <algorithm>

//trash?

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include "../../DeepTauID/interface/ntuple_content.h"
#include "../../DeepTauID/interface/ntuple_genInfo.h"
#include "../../DeepTauID/interface/ntuple_JetInfo.h"
#include "../../DeepTauID/interface/ntuple_pfCands.h"
#include "../../DeepTauID/interface/ntuple_recTau.h"
#include "../../DeepTauID/interface/genDecayHelper.h"



#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

struct MagneticField;


class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit DeepNtuplizer(const edm::ParameterSet&);
    ~DeepNtuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //can return 0
    pat::Jet* findMatchingJet(const reco::GenParticle* tau, const edm::Handle<edm::View<pat::Jet> >& jets)const;
    pat::Tau* findMatchingRecTau(const reco::GenParticle* tau, const edm::Handle<edm::View<pat::Tau> >& taus)const;
    pat::Tau* findMatchingRecTau(const reco::Jet* jet, const edm::Handle<edm::View<pat::Tau> >& taus)const;


    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
    edm::EDGetTokenT<edm::View<pat::Tau> >      tauToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >      genToken_;
    edm::EDGetTokenT<double> rhoToken_;

    std::string t_qgtagger;

    edm::Service<TFileService> fs;
    TTree *tree_;


    ntuple_content * addModule(ntuple_content *m){
        modules_.push_back(m);
        return m;
    }
    std::vector<ntuple_content* > modules_;

};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
                                    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
									tauToken_(consumes<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("taus"))),
                                    genToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genparticles"))),
                                    rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"))),
									tree_(0)
{
    /*
     *  Initialise the modules here
     *  Everything else does not need to be changed if
     *  modules don't interact.
     */

	// prunedGenParticles


    ntuple_genInfo* genInfo=new ntuple_genInfo();
    // optional configuration of the module
    // ...
    // add module
    addModule(genInfo);


    ntuple_recTau* recTau=new ntuple_recTau();
    // optional configuration of the module
    // ...
    // add module
    addModule(recTau);


    ntuple_JetInfo* JetInfo=new ntuple_JetInfo();
    // optional configuration of the module
    // ...
    // add module
    addModule(JetInfo);


    ntuple_pfCands* pfCands=new ntuple_pfCands();
    // optional configuration of the module
    // ...
    // add module
    addModule(pfCands);


    /*
     *
     * Modules initialized
     *
     * parse the input parameters (if any)
     */
    for(auto& m: modules_){
        m->getInput(iConfig);
        m->setIsRead(false);
    }
}


DeepNtuplizer::~DeepNtuplizer()
{
    return;
    for(auto& m:modules_)
        delete m;
}

pat::Jet* DeepNtuplizer::findMatchingJet(const reco::GenParticle* tau, const edm::Handle<edm::View<pat::Jet> >& jets)const{
	pat::Jet* ret=0;
	double mindr=0.3;
	for(size_t t=0;t<jets->size();t++){
		auto jet=jets->at(t);
		if(deltaR(jet.p4(),tau->p4())<mindr)
			ret=&jet;
	}
	return ret;
}
pat::Tau* DeepNtuplizer::findMatchingRecTau(const reco::GenParticle* gentau, const edm::Handle<edm::View<pat::Tau> >& taus)const{
	pat::Tau* ret=0;
	double mindr=0.3;
	for(size_t t=0;t<taus->size();t++){
		auto rectau=taus->at(t);
		if(deltaR(gentau->p4(),rectau.p4())<mindr)
			ret=&rectau;
	}
	return ret;
}
pat::Tau* DeepNtuplizer::findMatchingRecTau(const reco::Jet* jet, const edm::Handle<edm::View<pat::Tau> >& taus)const{
	pat::Tau* ret=0;
	double mindr=0.3;
	for(size_t t=0;t<taus->size();t++){
		auto tau=taus->at(t);
		if(deltaR(tau.p4(),jet->p4())<mindr)
			ret=&tau;
	}
	return ret;
}


// ------------ method called for each event  ------------
void
DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //global info


    edm::Handle<double> rhoInfo;
    iEvent.getByToken(rhoToken_,rhoInfo);

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<edm::View<pat::Tau> > taus;
    iEvent.getByToken(tauToken_, taus);

    edm::Handle<edm::View<reco::GenParticle> > gens;
    iEvent.getByToken(genToken_, gens);

    for(auto& m:modules_){
        m->readSetup(iSetup);
        m->readEvent(iEvent);
    }


    // loop over gen taus and associated jets
    // save all associated jets
    std::vector<pat::Jet *> tauJets;

    genDecayHelper dec_helper;
    for(size_t i=0;i<gens->size();i++){
    	auto p = gens->at(i);

    	if(!dec_helper.selectTau(p))continue;
    	genDecayHelper::decayModes dec=dec_helper.getGenTauDecayMode(&p);
    	if(dec == genDecayHelper::electron || dec == genDecayHelper::muon) continue;

    	pat::Jet * jet = findMatchingJet(&p,jets);
    	pat::Tau * rectau = findMatchingRecTau(&p,taus);

    	tauJets.push_back(jet);

    	bool writetau=true;
    	for(auto& m:modules_){
    		if(! m->fillBranches(rectau,jet,&p)){
    			writetau=false;
    		}
    	}
    	if(writetau)
    		tree_->Fill();

    }
    for(size_t j=0;j<jets->size();j++){
    	auto jet = jets->at(j);
    	if(std::find(tauJets.begin(),tauJets.end(),&jet) != tauJets.end()) continue;

    	//these have no gen Tau // maybe a fake reco tau
    	reco::GenParticle * notau=0;
    	pat::Tau * rectau = findMatchingRecTau(&jet,taus);

    	bool writejet=true;
    	for(auto& m:modules_){
    		if(! m->fillBranches(rectau,&jet,notau)){
    			writejet=false;
    		}
    	}
    	if(writejet)
    		tree_->Fill();

    }



}


// ------------ method called once each job just before starting event loop  ------------
void
DeepNtuplizer::beginJob()
{
    if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
    tree_=(fs->make<TTree>("tree" ,"tree" ));

    for(auto& m:modules_)
        m->initBranches(tree_);

}

// ------------ method called once each job just after ending the event loop  ------------
void
DeepNtuplizer::endJob()
{


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
