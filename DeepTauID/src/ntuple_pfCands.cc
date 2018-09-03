/*
 * ntuple_pfCands.cc
 *
 *  Created on: 24 Jul 2018
 *      Author: jkiesele
 */


#include "../interface/ntuple_pfCands.h"
#include "../interface/ntuple_content.h"
#include "../interface/ntuple_pfCands.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "../interface/sorting_modules.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"

class TrackInfoBuilder{
public:
	TrackInfoBuilder(edm::ESHandle<TransientTrackBuilder> & build):
		builder(build),
		trackMomentum_(0),
		trackEta_(0),
		trackEtaRel_(0),
		trackPtRel_(0),
		trackPPar_(0),
		trackDeltaR_(0),
		trackPtRatio_(0),
		trackPParRatio_(0),
		trackSip2dVal_(0),
		trackSip2dSig_(0),
		trackSip3dVal_(0),
		trackSip3dSig_(0),

		trackJetDistVal_(0),
		trackJetDistSig_(0)
{


}

	void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_ ,const math::XYZVector&  jetDir, GlobalVector refjetdirection, const reco::Vertex & pv){
		const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();

		reco::TransientTrack transientTrack;
		transientTrack=builder->build(PseudoTrack);
		Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
		Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
		Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
		math::XYZVector trackMom = PseudoTrack.momentum();
		double trackMag = std::sqrt(trackMom.Mag2());
		TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());
		TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());


		trackMomentum_=std::sqrt(trackMom.Mag2());
		trackEta_= trackMom.Eta();
		trackEtaRel_=reco::btau::etaRel(jetDir, trackMom);
		trackPtRel_=trackMom3.Perp(jetDir3);
		trackPPar_=jetDir.Dot(trackMom);
		trackDeltaR_=reco::deltaR(trackMom, jetDir);
		trackPtRatio_=trackMom3.Perp(jetDir3) / trackMag;
		trackPParRatio_=jetDir.Dot(trackMom) / trackMag;
		trackSip2dVal_=(meas_ip2d.value());

		trackSip2dSig_=(meas_ip2d.significance());
		trackSip3dVal_=(meas_ip3d.value());


		trackSip3dSig_=meas_ip3d.significance();
		trackJetDistVal_= jetdist.value();
		trackJetDistSig_= jetdist.significance();

	}

	const float& getTrackDeltaR() const {return trackDeltaR_;}
	const float& getTrackEta() const {return trackEta_;}
	const float& getTrackEtaRel() const {return trackEtaRel_;}
	const float& getTrackJetDistSig() const {return trackJetDistSig_;}
	const float& getTrackJetDistVal() const {return trackJetDistVal_;}
	const float& getTrackMomentum() const {return trackMomentum_;}
	const float& getTrackPPar() const {return trackPPar_;}
	const float& getTrackPParRatio() const {return trackPParRatio_;}
	const float& getTrackPtRatio() const {return trackPtRatio_;}
	const float& getTrackPtRel() const {return trackPtRel_;}
	const float& getTrackSip2dSig() const {return trackSip2dSig_;}
	const float& getTrackSip2dVal() const {return trackSip2dVal_;}
	const float& getTrackSip3dSig() const {return trackSip3dSig_;}
	const float& getTrackSip3dVal() const {return trackSip3dVal_;}

private:

	edm::ESHandle<TransientTrackBuilder>& builder;

	float trackMomentum_;
	float trackEta_;
	float trackEtaRel_;
	float trackPtRel_;
	float trackPPar_;
	float trackDeltaR_;
	float trackPtRatio_;
	float trackPParRatio_;
	float trackSip2dVal_;
	float trackSip2dSig_;
	float trackSip3dVal_;
	float trackSip3dSig_;

	float trackJetDistVal_;
	float trackJetDistSig_;

};




void ntuple_pfCands::readSetup(const edm::EventSetup& iSetup){

	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}

void ntuple_pfCands::getInput(const edm::ParameterSet& iConfig){

}


float ntuple_pfCands::catchInfsAndBound_track(const pat::PackedCandidate* cand, const float& in,const float& replace_value,
        const float& lowerbound, const float& upperbound,const float offset)const{
	if(cand->hasTrackDetails())
		return catchInfsAndBound(in,replace_value,lowerbound,upperbound,offset);
	else
		return replace_value;
}


void ntuple_pfCands::initBranches(TTree* t){
	ADDBRANCH(t,nCpfcan);
	ADDBRANCH(t,Cpfcan_pt);
	ADDBRANCH(t,Cpfcan_eta);
	ADDBRANCH(t,Cpfcan_phi);
	ADDBRANCH(t,Cpfcan_ptrel);
	ADDBRANCH(t,Cpfcan_erel);
	ADDBRANCH(t,Cpfcan_phirel);
	ADDBRANCH(t,Cpfcan_etarel);
	ADDBRANCH(t,Cpfcan_deltaR);
	ADDBRANCH(t,Cpfcan_puppiw);
	ADDBRANCH(t,Cpfcan_VTX_ass);

	ADDBRANCH(t,Cpfcan_fromPV);

	// empty ADDBRANCH(t,Cpfcan_vertexChi2);
	// empty ADDBRANCH(t,Cpfcan_vertexNdof);
	// empty ADDBRANCH(t,Cpfcan_vertexNormalizedChi2);

	ADDBRANCH(t,Cpfcan_vertex_rho);
	ADDBRANCH(t,Cpfcan_vertex_phirel);
	ADDBRANCH(t,Cpfcan_vertex_etarel);
	// empty ADDBRANCH(t,Cpfcan_vertexRef_mass);


	ADDBRANCH(t,Cpfcan_dz);
	ADDBRANCH(t,Cpfcan_dxy);

	ADDBRANCH(t,Cpfcan_dxyerrinv);
	ADDBRANCH(t,Cpfcan_dxysig);


	ADDBRANCH(t,Cpfcan_BtagPf_trackMomentum);
	ADDBRANCH(t,Cpfcan_BtagPf_trackEta);
	ADDBRANCH(t,Cpfcan_BtagPf_trackEtaRel);
	ADDBRANCH(t,Cpfcan_BtagPf_trackPtRel);
	ADDBRANCH(t,Cpfcan_BtagPf_trackPPar);
	ADDBRANCH(t,Cpfcan_BtagPf_trackDeltaR);
	ADDBRANCH(t,Cpfcan_BtagPf_trackPtRatio);
	ADDBRANCH(t,Cpfcan_BtagPf_trackPParRatio);
	ADDBRANCH(t,Cpfcan_BtagPf_trackSip3dVal);
	ADDBRANCH(t,Cpfcan_BtagPf_trackSip3dSig);
	ADDBRANCH(t,Cpfcan_BtagPf_trackSip2dVal);
	ADDBRANCH(t,Cpfcan_BtagPf_trackSip2dSig);


	ADDBRANCH(t,Cpfcan_BtagPf_trackJetDistVal);

	ADDBRANCH(t,Cpfcan_isMu);
	ADDBRANCH(t,Cpfcan_isEl);
	ADDBRANCH(t,Cpfcan_pdgID);

	ADDBRANCH(t,Cpfcan_charge);


	ADDBRANCH(t,Cpfcan_lostInnerHits);
	ADDBRANCH(t,Cpfcan_numberOfPixelHits);

	ADDBRANCH(t,Cpfcan_chi2);
	ADDBRANCH(t,Cpfcan_quality);


	ADDBRANCH(t,nNpfcand);
	ADDBRANCH(t,Npfcan_pt);
	ADDBRANCH(t,Npfcan_eta);
	ADDBRANCH(t,Npfcan_phi);
	ADDBRANCH(t,Npfcan_ptrel);
	ADDBRANCH(t,Npfcan_erel);
	ADDBRANCH(t,Npfcan_puppiw);
	ADDBRANCH(t,Npfcan_phirel);
	ADDBRANCH(t,Npfcan_etarel);
	ADDBRANCH(t,Npfcan_deltaR);
	ADDBRANCH(t,Npfcan_isGamma);
	ADDBRANCH(t,Npfcan_HadFrac);
	ADDBRANCH(t,Npfcan_pdgID);

}

void ntuple_pfCands::clear(){
	nCpfcan=0;
	Cpfcan_pt.clear();
	Cpfcan_eta.clear();
	Cpfcan_phi.clear();
	Cpfcan_ptrel.clear();
	Cpfcan_erel.clear();
	Cpfcan_phirel.clear();
	Cpfcan_etarel.clear();
	Cpfcan_deltaR.clear();
	Cpfcan_puppiw.clear();
	Cpfcan_VTX_ass.clear();

	Cpfcan_fromPV.clear();

	Cpfcan_vertexChi2.clear();
	Cpfcan_vertexNdof.clear();
	Cpfcan_vertexNormalizedChi2.clear();
	Cpfcan_vertex_rho.clear();
	Cpfcan_vertex_phirel.clear();
	Cpfcan_vertex_etarel.clear();
	Cpfcan_vertexRef_mass.clear();


	Cpfcan_dz.clear();
	Cpfcan_dxy.clear();

	Cpfcan_dxyerrinv.clear();
	Cpfcan_dxysig.clear();


	Cpfcan_BtagPf_trackMomentum.clear();
	Cpfcan_BtagPf_trackEta.clear();
	Cpfcan_BtagPf_trackEtaRel.clear();
	Cpfcan_BtagPf_trackPtRel.clear();
	Cpfcan_BtagPf_trackPPar.clear();
	Cpfcan_BtagPf_trackDeltaR.clear();
	Cpfcan_BtagPf_trackPtRatio.clear();
	Cpfcan_BtagPf_trackPParRatio.clear();
	Cpfcan_BtagPf_trackSip3dVal.clear();
	Cpfcan_BtagPf_trackSip3dSig.clear();
	Cpfcan_BtagPf_trackSip2dVal.clear();
	Cpfcan_BtagPf_trackSip2dSig.clear();


	Cpfcan_BtagPf_trackJetDistVal.clear();

	Cpfcan_isMu.clear();
	Cpfcan_isEl.clear();
	Cpfcan_pdgID.clear();
	Cpfcan_charge.clear();


	Cpfcan_lostInnerHits.clear();
	Cpfcan_numberOfPixelHits.clear();

	Cpfcan_chi2.clear();
	Cpfcan_quality.clear();


	nNpfcand=0;
	Npfcan_pt.clear();
	Npfcan_eta.clear();
	Npfcan_phi.clear();
	Npfcan_ptrel.clear();
	Npfcan_erel.clear();
	Npfcan_puppiw.clear();
	Npfcan_phirel.clear();
	Npfcan_etarel.clear();
	Npfcan_deltaR.clear();
	Npfcan_isGamma.clear();
	Npfcan_HadFrac.clear();
	Npfcan_drminsv.clear();
	Npfcan_pdgID.clear();
}

bool ntuple_pfCands::fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau){
	clear();
	if(!recJet)
		return true;

	const pat::Jet & jet=*recJet;

	reco::Candidate::LorentzVector refVector=jet.p4();
	if(recTau)
		refVector=recTau->p4();

	float etasign = 1.;
	if (refVector.eta()<0) etasign =-1.;


	math::XYZVector jetDir = refVector.Vect().Unit();
	GlobalVector jetRefTrackDir(refVector.px(),refVector.py(),refVector.pz());


	TrackInfoBuilder trackinfo(builder);

	for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
		const pat::PackedCandidate* PackedCandidate_ = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
		if(!PackedCandidate_)continue;

		if(PackedCandidate_->charge()!=0 ){

			Cpfcan_pt.push_back(PackedCandidate_->pt());
			Cpfcan_eta.push_back(PackedCandidate_->eta());
			Cpfcan_phi.push_back(PackedCandidate_->phi());
			Cpfcan_ptrel.push_back(catchInfsAndBound(PackedCandidate_->pt()/refVector.pt(),0,-1,0,-1));
			Cpfcan_erel.push_back(catchInfsAndBound(PackedCandidate_->energy()/refVector.pt(),0,-1,0,-1));
			Cpfcan_phirel.push_back(catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),refVector.phi())),0,-2,0,-0.5));
			Cpfcan_etarel.push_back(catchInfsAndBound(fabs(PackedCandidate_->eta()-refVector.eta()),0,-2,0,-0.5));
			Cpfcan_deltaR.push_back(catchInfsAndBound(reco::deltaR(*PackedCandidate_,refVector),0,-0.6,0,-0.6));

			Cpfcan_VTX_ass.push_back(PackedCandidate_->pvAssociationQuality());

			Cpfcan_fromPV.push_back(PackedCandidate_->fromPV());

			Cpfcan_puppiw.push_back(PackedCandidate_->puppiWeight());

			Cpfcan_isMu.push_back(abs(PackedCandidate_->pdgId())==13);
			Cpfcan_isEl.push_back(abs(PackedCandidate_->pdgId())==11);
			Cpfcan_pdgID.push_back(PackedCandidate_->pdgId());

			Cpfcan_charge.push_back(PackedCandidate_->charge());

			Cpfcan_lostInnerHits.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->lostInnerHits(),2,0,10));
			Cpfcan_numberOfPixelHits.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->numberOfPixelHits(),-1,0,10));

			Cpfcan_vertexChi2.push_back(PackedCandidate_->vertexChi2());
			Cpfcan_vertexNdof.push_back(PackedCandidate_->vertexNdof());
			//divided
			Cpfcan_vertexNormalizedChi2.push_back(PackedCandidate_->vertexNormalizedChi2());
			Cpfcan_vertex_rho.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->vertex().rho(),0,-1,50));
			Cpfcan_vertex_phirel.push_back(catchInfsAndBound_track(PackedCandidate_,reco::deltaPhi(PackedCandidate_->vertex().phi(),refVector.phi()),0,0,10));
			Cpfcan_vertex_etarel.push_back(catchInfsAndBound_track(PackedCandidate_,etasign*(PackedCandidate_->vertex().eta()-refVector.eta()),0,-5,5));
			Cpfcan_vertexRef_mass.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->vertexRef()->p4().M(),0,0,10));


			if(PackedCandidate_->hasTrackDetails()){// && PackedCandidate_->pt()<0.75){
				Cpfcan_dxy.push_back(catchInfsAndBound_track(PackedCandidate_, fabs(PackedCandidate_->dxy()),0,-50,50));
				Cpfcan_dxyerrinv.push_back(catchInfsAndBound_track(PackedCandidate_,1/PackedCandidate_->dxyError(),0,-1, 10000.));
				Cpfcan_dxysig.push_back(catchInfsAndBound_track(PackedCandidate_,fabs(PackedCandidate_->dxy()/PackedCandidate_->dxyError()),0.,-2000,2000));
				Cpfcan_dz.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->dz(),0,-100,100));

				trackinfo.buildTrackInfo(PackedCandidate_,jetDir,jetRefTrackDir,vertices()->at(0));

				Cpfcan_BtagPf_trackMomentum.push_back(catchInfsAndBound(trackinfo.getTrackMomentum(),0,0 ,1000));
				Cpfcan_BtagPf_trackEta.push_back(catchInfsAndBound(trackinfo.getTrackEta()   ,  0,-5,5));
				Cpfcan_BtagPf_trackEtaRel.push_back(     catchInfsAndBound(trackinfo.getTrackEtaRel(),  0,-5,15));
				Cpfcan_BtagPf_trackPtRel.push_back(      catchInfsAndBound(trackinfo.getTrackPtRel(),   0,-1,4));
				Cpfcan_BtagPf_trackPPar.push_back(       catchInfsAndBound(trackinfo.getTrackPPar(),    0,-1e5,1e5 ));
				Cpfcan_BtagPf_trackDeltaR.push_back(     catchInfsAndBound(trackinfo.getTrackDeltaR(),  0,-5,5 ));
				Cpfcan_BtagPf_trackPtRatio.push_back(    catchInfsAndBound(trackinfo.getTrackPtRatio(), 0,-1,10 ));
				Cpfcan_BtagPf_trackPParRatio.push_back(  catchInfsAndBound(trackinfo.getTrackPParRatio(),0,-10,100));
				Cpfcan_BtagPf_trackSip3dVal.push_back(   catchInfsAndBound(trackinfo.getTrackSip3dVal(), 0, -1,1e5 ));
				Cpfcan_BtagPf_trackSip3dSig.push_back(   catchInfsAndBound(trackinfo.getTrackSip3dSig(), 0, -1,4e4 ));
				Cpfcan_BtagPf_trackSip2dVal.push_back(   catchInfsAndBound(trackinfo.getTrackSip2dVal(), 0, -1,70 ));
				Cpfcan_BtagPf_trackSip2dSig.push_back(   catchInfsAndBound(trackinfo.getTrackSip2dSig(), 0, -1,4e4 ));

				Cpfcan_BtagPf_trackJetDistVal.push_back( catchInfsAndBound(trackinfo.getTrackJetDistVal(),0,-20,1 ));

				Cpfcan_chi2.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->pseudoTrack().normalizedChi2(),300,-1,300));
				Cpfcan_quality.push_back(catchInfsAndBound_track(PackedCandidate_,PackedCandidate_->pseudoTrack().qualityMask(),0,0,100));
			}
			else{

				Cpfcan_dxy.push_back(0);
				Cpfcan_dxyerrinv.push_back(0);
				Cpfcan_dxysig.push_back(0);
				Cpfcan_dz.push_back(0);


				Cpfcan_BtagPf_trackMomentum.push_back(  0);
				Cpfcan_BtagPf_trackEta.push_back(       0);
				Cpfcan_BtagPf_trackEtaRel.push_back(    0);
				Cpfcan_BtagPf_trackPtRel.push_back(     0);
				Cpfcan_BtagPf_trackPPar.push_back(      0);
				Cpfcan_BtagPf_trackDeltaR.push_back(    0);
				Cpfcan_BtagPf_trackPtRatio.push_back(   0);
				Cpfcan_BtagPf_trackPParRatio.push_back( 0);
				Cpfcan_BtagPf_trackSip3dVal.push_back(  0);
				Cpfcan_BtagPf_trackSip3dSig.push_back(  0);
				Cpfcan_BtagPf_trackSip2dVal.push_back(  0);
				Cpfcan_BtagPf_trackSip2dSig.push_back(  0);
				Cpfcan_BtagPf_trackJetDistVal.push_back(0);

				Cpfcan_chi2.push_back(0);
				Cpfcan_quality.push_back(0);
			}

		}
		else{// neutral candidates
			Npfcan_pt.push_back(PackedCandidate_->pt());
			Npfcan_eta.push_back(PackedCandidate_->eta());
			Npfcan_phi.push_back(PackedCandidate_->phi());
			Npfcan_ptrel.push_back(catchInfsAndBound(PackedCandidate_->pt()/refVector.pt(),0,-1,0,-1));
			Npfcan_erel.push_back(catchInfsAndBound(PackedCandidate_->energy()/refVector.pt(),0,-1,0,-1));
			Npfcan_puppiw.push_back(PackedCandidate_->puppiWeight());
			Npfcan_phirel.push_back(catchInfsAndBound(fabs(reco::deltaPhi(PackedCandidate_->phi(),refVector.phi())),0,-2,0,-0.5));
			Npfcan_etarel.push_back(catchInfsAndBound(fabs(PackedCandidate_->eta()-refVector.eta()),0,-2,0,-0.5));
			Npfcan_deltaR.push_back(catchInfsAndBound(reco::deltaR(*PackedCandidate_,refVector),0,-0.6,0,-0.6));
			Npfcan_isGamma.push_back(fabs(PackedCandidate_->pdgId())==22);
			Npfcan_HadFrac.push_back(PackedCandidate_->hcalFraction());
			Npfcan_pdgID.push_back(PackedCandidate_->pdgId());

		}

	} // end loop over jet.numberOfDaughters()


	nCpfcan=Cpfcan_pt.size();
	nNpfcand=Npfcan_pt.size();

	return true; //for making cuts
}







