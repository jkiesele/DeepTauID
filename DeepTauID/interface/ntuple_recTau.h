/*
 * ntuple_recTau.h
 *
 *  Created on: 20 Jul 2018
 *      Author: jkiesele
 */

#ifndef DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_RECTAU_H_
#define DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_RECTAU_H_


#include "ntuple_content.h"
/**
 * Contains gen information for the reco tau object (hi-level variables)
 */
class ntuple_recTau: public ntuple_content{
public:
	ntuple_recTau();


    void getInput(const edm::ParameterSet& iConfig){}
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent){}
    void readSetup(const edm::EventSetup& iSetup){}
    //use either of these functions

    bool fillBranches(const pat::Tau* recTau, const pat::Jet* recJet, const reco::GenParticle* genTau);

    void clear();

private:

    float pt_weighted_dx(const pat::Tau& tau, int mode = 0, int var = 0, int decaymode = -1)const;
    float pt_weighted_dr_signal(const pat::Tau& tau, int dm)const;
    float pt_weighted_deta_strip(const pat::Tau& tau, int dm)const;
    float pt_weighted_dphi_strip(const pat::Tau& tau, int dm)const;
    float pt_weighted_dr_iso(const pat::Tau& tau, int dm)const;

    unsigned int n_photons_total(const pat::Tau& tau)const;
    bool dynamic_isInside(float photon_pt, float deta, float dphi)const;


    int isRecTau;
    int recTauDecayMode_i;
    float recTauDecayMode;
    float recTau_pt;
    float recTau_eta;
    float recTau_phi;
    float recTau_M;
    float recTauVtxZ;
    float recImpactParamPCA_x;
    float recImpactParamPCA_y;
    float recImpactParamPCA_z;
    float recImpactParam;
    float recImpactParamSign;
    //float recImpactParamPCA3D_x;
    //float recImpactParamPCA3D_y;
    //float recImpactParamPCA3D_z;
    float recImpactParam3D;
    float recImpactParamSign3D;
    //float recImpactParamZ;
    //float recImpactParamSignZ;
    //float recImpactParamTk2;
    //float recImpactParamSignTk2;
    //float recImpactParam3DTk2;
    //float recImpactParamSignZTk2;
    //float recImpactParamTk3;
    //float recImpactParamSignTk3;
    //float recImpactParam3DTk3;
    //float recImpactParamSign3DTk3;
    //float recImpactParamZTk3;
    //float recImpactParamSignZTk3;
    //float recDecayLengthTk1;
    //float recDecayLengthSignTk1;
    //float recDecayLengthTk2;
    //float recDecayLengthSignTk2;
    //float recDecayLengthTk3;
    //float recDecayLengthSignTk3;
    //float recDecayDist2D;
    //float recDecayDistSign2D;
    float recChi2DiffEvtVertex;
    //
    float hasRecDecayVertex;
    //float recDecayVertex_x; //not in miniAOD
    //float recDecayVertex_y;
    //float recDecayVertex_z;
    //missing recDecayVertexCov
    //float recDecayVertexChi2;
    float recDecayDist_x;
    float recDecayDist_y;
    float recDecayDist_z;
    //missing recDecayDistCov
    float recDecayDistSign;
    //float recEvtVertex_x;
    //float recEvtVertex_y;
    //float recEvtVertex_z;
    //missing recEvtVertexCov
    float recTauPtWeightedDetaStrip;
    float recTauPtWeightedDphiStrip;
    float recTauPtWeightedDrSignal;
    float recTauPtWeightedDrIsolation;
    float recTauNphoton;
    float recTauEratio;
    float recTauLeadingTrackChi2;
    float recTauNphotonSignal;
    float recTauNphotonIso;
    float recTauNphotonTotal;

    //FIXME TBI
//    float recTauGJangleMeasured;
//    float recTauGJangleDiff;

    float recTauLeadPFCand_pt;
    float recTauLeadPFCand_eta;
    float recTauLeadPFCand_phi;
    float recTauLeadPFCand_M;

    float recTauLeadPFChargedHadrCand_pt;
    float recTauLeadPFChargedHadrCand_eta;
    float recTauLeadPFChargedHadrCand_phi;
    float recTauLeadPFChargedHadrCand_M;

    //don't consider pi0 given that it is easy to reconstruct from gammas in the DNN

    float chargedIsoPtSum;
    float neutralIsoPtSum;
    float puCorrPtSum;
    float neutralIsoPtSumWeight;
    float footprintCorrection;
    float photonPtSumOutsideSignalCone;
    float chargedIsoPtSumdR03;
    float neutralIsoPtSumdR03;
    float neutralIsoPtSumWeightdR03;
    float footprintCorrectiondR03;
    float photonPtSumOutsideSignalConedR03;



};



#endif /* DEEPTAU_DEEPTAUID_INTERFACE_NTUPLE_RECTAU_H_ */
