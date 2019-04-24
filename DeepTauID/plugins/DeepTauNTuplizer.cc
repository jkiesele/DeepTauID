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
#include "../../DeepTauID/interface/ntuple_global.h"
#include "../../DeepTauID/interface/ntuple_JetInfo.h"
#include "../../DeepTauID/interface/ntuple_pfCands.h"
#include "../../DeepTauID/interface/ntuple_recTau.h"
#include "../../DeepTauID/interface/genDecayHelper.h"

#include "TRandom3.h"


#if defined( __GXX_EXPERIMENTAL_CXX0X__)
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#endif

struct MagneticField;


class DeepTauNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit DeepTauNTuplizer(const edm::ParameterSet&);
    ~DeepTauNTuplizer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //can return 0
    const pat::Jet* findMatchingJet(const reco::GenParticle* tau,    const std::vector<const pat::Jet *> & jets, bool requireGen)const;
    const pat::Tau* findMatchingRecTau(const reco::GenParticle* tau, const std::vector<const pat::Tau *> & taus)const;
    const pat::Tau* findMatchingRecTau(const reco::Jet* jet,         const std::vector<const pat::Tau *> & taus)const;



    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<edm::View<pat::Jet> >      jetToken_;
    edm::EDGetTokenT<edm::View<pat::Tau> >      tauToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >      genToken_;
    edm::EDGetTokenT<double> rhoToken_;

    TRandom3 rand;

    double gentau_minpt_;
    double gentau_maxeta_;
    double rectau_minpt_;
    double rectau_maxeta_;
    double jet_minpt_;
	double jet_maxeta_;
	double noTauNoGen_reduction_;
	double promptLepton_reduction_;
	double backgroundjet_reduction_;
	double oversample_;
	bool onlyrectaus_,onlytaus_;

    edm::Service<TFileService> fs;
    TTree *tree_;

    genDecayHelper dec_helper;

    ntuple_content * addModule(ntuple_content *m){
        modules_.push_back(m);
        return m;
    }
    std::vector<ntuple_content* > modules_;

};

DeepTauNTuplizer::DeepTauNTuplizer(const edm::ParameterSet& iConfig):
				vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
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

    gentau_minpt_=18;
    gentau_maxeta_=3;
    rectau_minpt_=20;
    rectau_maxeta_=3;
    jet_minpt_=20;
	jet_maxeta_=3;
	noTauNoGen_reduction_=0.0; //reduce non-tau no gen (pileup) contribution to 0%
	promptLepton_reduction_=0.6; //reduce to 0% muon/electron contamination for tests
	backgroundjet_reduction_=0.6;
	onlyrectaus_=iConfig.getParameter<bool>("onlyRecTaus");
	onlytaus_=iConfig.getParameter<bool>("onlyTaus");
	oversample_=iConfig.getParameter<int>("overSample");

	ntuple_global * globals = new ntuple_global();
	addModule(globals);

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


DeepTauNTuplizer::~DeepTauNTuplizer()
{
    return;
    for(auto& m:modules_)
        delete m;
}

const pat::Jet* DeepTauNTuplizer::findMatchingJet(const reco::GenParticle* tau, const std::vector<const pat::Jet *> & jets, bool requireGen)const{
	const pat::Jet* ret=0;
	double mindr=0.3;
	reco::Candidate::LorentzVector vismomentum=tau->p4();
	if(dec_helper.isPromptTau(*tau)){
		vismomentum=dec_helper.getVisMomentum(tau);
	}
	for(size_t t=0;t<jets.size();t++){
		auto jet=jets.at(t);
		if(requireGen && !jet->genJet()) continue;
		double deltar=deltaR(jet->p4(),vismomentum);
		if(deltar<mindr){
			mindr=deltar;
			ret=jet;
		}
	}
	return ret;
}
const pat::Tau* DeepTauNTuplizer::findMatchingRecTau(const reco::GenParticle* gentau, const std::vector<const pat::Tau *> & taus)const{
	const pat::Tau* ret=0;
	double mindr=0.3;
	reco::Candidate::LorentzVector vismomentum=gentau->p4();
	if(dec_helper.isPromptTau(*gentau)){
		vismomentum=dec_helper.getVisMomentum(gentau);
	}
	for(size_t t=0;t<taus.size();t++){
		auto rectau=taus.at(t);
		double deltar=deltaR(vismomentum,rectau->p4());
		if(deltar<mindr){
			mindr=deltar;
			ret=rectau;
		}
	}
	return ret;
}
const pat::Tau* DeepTauNTuplizer::findMatchingRecTau(const reco::Jet* jet, const std::vector<const pat::Tau *> & taus)const{
	const pat::Tau* ret=0;
	double mindr=0.15;
	for(size_t t=0;t<taus.size();t++){
		auto tau=taus.at(t);
		double deltar=deltaR(tau->p4(),jet->p4());
		if(deltar<mindr){
			mindr=deltar;
			ret=tau;
		}
	}
	return ret;
}


// ------------ method called for each event  ------------
void
DeepTauNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


    edm::Handle<edm::View<reco::GenParticle> > h_gens;
    iEvent.getByToken(genToken_, h_gens);

    if(! h_gens.product()) return;
    std::vector<const reco::GenParticle *> gens;
    for(size_t i=0;i<h_gens->size();i++)
    	gens.push_back(&h_gens->at(i));

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return;


    edm::Handle<double> rhoInfo;
    iEvent.getByToken(rhoToken_,rhoInfo);

    edm::Handle<edm::View<pat::Jet> > h_jets;
    iEvent.getByToken(jetToken_, h_jets);
    std::vector<const pat::Jet *> jets;
    if(h_jets.product()){
    	for(size_t i=0;i<h_jets->size();i++)
    		jets.push_back(&h_jets->at(i));
    }

    edm::Handle<edm::View<pat::Tau> > h_taus;
    iEvent.getByToken(tauToken_, h_taus);
    std::vector<const pat::Tau *> taus;
    if(h_taus.product()){
    	for(size_t i=0;i<h_taus->size();i++)
    		taus.push_back(&h_taus->at(i));
    }

    if(gens.size())
    	rand.SetSeed((int) (1000*gens.at(0)->pt()));


    for(auto& m:modules_){
    	m->clear();
        m->setPrimaryVertices(vertices.product());
        m->readSetup(iSetup);
        m->readEvent(iEvent);
        m->setRhoInfo(rhoInfo.product());
    }


    // loop over gen taus and associated jets
    // save all associated jets
    std::vector<const pat::Jet *> tauJets;
    std::vector<const pat::Tau *> recTaus;
    std::vector<const reco::GenParticle *> leptons; //for overlap check

    for(auto p: gens){

    	if(dec_helper.isLooseGenLepton(*p)){
    		leptons.push_back(p);
    	}
    	// this also allows electrons and muons and tau decays to electrons and muons - later flagged as !isTau
    	if(!dec_helper.isPromptLepton(*p))continue;


    	const pat::Jet * jet = findMatchingJet(p,jets,true);
		tauJets.push_back(jet);

    	const pat::Tau * rectau = findMatchingRecTau(p,taus);

    	//if(jet && (jet->pt()<jet_minpt_ || fabs(jet->eta())>jet_maxeta_)) jet=0;
    	recTaus.push_back(rectau);

    	if(dec_helper.isPromptTau(*p)){
    		reco::Candidate::LorentzVector vismomentum=dec_helper.getVisMomentum(p);
    		if(( vismomentum.pt() < gentau_minpt_ || fabs(vismomentum.eta())> gentau_maxeta_)) continue;
    	}

    	if(! dec_helper.isPromptTau(*p) && rand.Uniform(0.,1.) > promptLepton_reduction_ ) continue;

    	bool writetau=true;
    	if(onlyrectaus_ && (!rectau || rectau->pt()<rectau_minpt_ || fabs(rectau->eta()) > rectau_maxeta_)) continue;

    	for(auto& m:modules_){
    		if(! m->fillBranches(rectau,jet,p,&gens)){ //MAKE TAU GEN CUTS HERE
    			writetau=false;
    		}
    	}
    	if(writetau){
    		for(int i=0;i<oversample_;i++)
    			tree_->Fill();
    	}

    }

    if(onlytaus_) return;

    for(auto& m:modules_){
    	m->clear();
    }

    std::vector<const pat::Tau *> unmatched_recTaus;
    for(const auto t:taus){
    	if(std::find(recTaus.begin(), recTaus.end(), t) != recTaus.end())continue;
    	unmatched_recTaus.push_back(t);
    }

    for(auto jet: jets){
    	if(std::find(tauJets.begin(),tauJets.end(),jet) != tauJets.end()) continue;


    	//remove some PU jets
    	if(! jet->genJet() && rand.Uniform(0.,1.)>noTauNoGen_reduction_)continue;

    	//these have no gen Tau // maybe a fake reco tau
    	reco::GenParticle * nolepton=0;
    	const pat::Tau * rectau = findMatchingRecTau(jet,unmatched_recTaus);


    	if(onlyrectaus_ && (!rectau || rectau->pt()<rectau_minpt_ || fabs(rectau->eta()) > rectau_maxeta_)) continue;
    	if(!onlyrectaus_ && (jet->pt()<jet_minpt_ || fabs(jet->eta())>jet_maxeta_)) continue;

    	// remove any jet that could overlap with a tau or lepton (dR=0.5)
    	// avoiding ambiguous cases
    	bool usejet=true;
    	for(auto lep: leptons){
    		if(reco::deltaR2(lep->p4(),jet->p4()) < 0.5*0.5){
    			usejet=false;
    			break;
    		}
    	}
    	if(!usejet) continue;

    	if(rand.Uniform(0.,1.)>backgroundjet_reduction_)continue;

    	// the implicit cuts that CAN be implemented in fillBranches are NOT used here
    	bool writejet=true;
    	if(onlyrectaus_ && ! rectau) writejet=false;
    	for(auto& m:modules_){
    		if(! m->fillBranches(rectau,jet,nolepton,&gens)){
    			writejet=false;
    		}
    	}
    	if(writejet)
    		for(int i=0;i<oversample_;i++)
    			tree_->Fill();

    }



}


// ------------ method called once each job just before starting event loop  ------------
void
DeepTauNTuplizer::beginJob()
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
DeepTauNTuplizer::endJob()
{


}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DeepTauNTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepTauNTuplizer);
