//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 31 23:06:44 2018 by ROOT version 5.34/30
// from TTree ModFixedStripTaus/ModFixedStripTaus
// found on file: GGHTT_PU200.root
//////////////////////////////////////////////////////////

#ifndef signal_pu0_h
#define signal_pu0_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxtauPt = 1;
   const Int_t kMaxtauEta = 1;
   const Int_t kMaxgenTauPt = 1;
   const Int_t kMaxgenTauEta = 1;
   const Int_t kMaxgenTauMatch = 1;
   const Int_t kMaxnvtx = 1;
   const Int_t kMaxvtxX = 1;
   const Int_t kMaxvtxY = 1;
   const Int_t kMaxvtxZ = 1;
   const Int_t kMaxvtxIndex = 1;
   const Int_t kMaxdm = 1;
   const Int_t kMaxgoodReco = 1;
   const Int_t kMaxtauMass = 1;

class signal_pu0 : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   Double_t        tauPt;
   Double_t        tauEta;
   Double_t        genTauPt;
   Double_t        genTauEta;
   Int_t           genTauMatch;
   Int_t           nvtx;
   Double_t        vtxX;
   Double_t        vtxY;
   Double_t        vtxZ;
   Int_t           vtxIndex;
   Double_t        dz_tt;
   Double_t        dxy_tt;
   Int_t           Nparticles;
   Float_t         DR_IsoCandTau[31];   //[Nparticles]
   Float_t         cand_pt[31];   //[Nparticles]
   Float_t         cand_eta[31];   //[Nparticles]
   Float_t         cand_phi[31];   //[Nparticles]
   Int_t           cand_q[31];   //[Nparticles]
   Int_t           cand_numHits[31];   //[Nparticles]
   Float_t         cand_dz_new[31];   //[Nparticles]
   Float_t         cand_dxy_new[31];   //[Nparticles]
   Int_t           Ngammas;
   Float_t         DR_gamma[75];   //[Ngammas]
   Float_t         gamma_pt[75];   //[Ngammas]
   Float_t         gamma1_pt[75];   //[Ngammas]
   Float_t         gamma2_pt[75];   //[Ngammas]
   Float_t         gamma_eta[75];   //[Ngammas]
   Float_t         gamma_phi[75];   //[Ngammas]
   Int_t           dm;
   Int_t           goodReco;
   Double_t        tauMass;
   Bool_t          taupfTausDiscriminationByDecayModeFinding;
   Bool_t          taupfTausDiscriminationByDecayModeFindingNewDMs;
   Float_t         tauByIsolationMVArun2v1DBnewDMwLTraw;
   Float_t         tauByIsolationMVArun2v1DBoldDMwLTraw;
   Float_t         tauByIsolationMVArun2v1PWnewDMwLTraw;
   Float_t         tauByIsolationMVArun2v1PWoldDMwLTraw;
   Float_t         tauChargedIsoPtSum;
   Float_t         tau_iso_pv;
   Float_t         tau_iso_dz0015;
   Float_t         tau_iso_dz0015_dz;
   Float_t         tau_iso_dz0015_dxy;
   Float_t         tau_iso_dz01_dz;
   Float_t         tau_iso_dz01_dxy;
   Float_t         tau_iso_dz01;
   Float_t         tauNeutralIsoPtSum;
   Bool_t          tauByLooseCombinedIsolationDeltaBetaCorr3Hits;
   Bool_t          tauByMediumCombinedIsolationDeltaBetaCorr3Hits;
   Bool_t          tauByTightCombinedIsolationDeltaBetaCorr3Hits;
   Float_t         tauCombinedIsolationDeltaBetaCorrRaw3Hits;
   Float_t         tauPuCorrPtSum;
   Float_t         taufootprintCorrection;
   Float_t         tauphotonPtSumOutsideSignalCone;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_tauPt_;   //!
   TBranch        *b_tauEta_;   //!
   TBranch        *b_genTauPt_;   //!
   TBranch        *b_genTauEta_;   //!
   TBranch        *b_genTauMatch_;   //!
   TBranch        *b_nvtx_;   //!
   TBranch        *b_vtxX_;   //!
   TBranch        *b_vtxY_;   //!
   TBranch        *b_vtxZ_;   //!
   TBranch        *b_vtxIndex_;   //!
   TBranch        *b_dz_tt;   //!
   TBranch        *b_dxy_tt;   //!
   TBranch        *b_Nparticles;   //!
   TBranch        *b_DR_IsoCandTau;   //!
   TBranch        *b_cand_pt;   //!
   TBranch        *b_cand_eta;   //!
   TBranch        *b_cand_phi;   //!
   TBranch        *b_cand_q;   //!
   TBranch        *b_cand_numHits;   //!
   TBranch        *b_cand_dz_new;   //!
   TBranch        *b_cand_dxy_new;   //!
   TBranch        *b_Ngammas;   //!
   TBranch        *b_DR_gamma;   //!
   TBranch        *b_gamma_pt;   //!
   TBranch        *b_gamma1_pt;   //!
   TBranch        *b_gamma2_pt;   //!
   TBranch        *b_gamma_eta;   //!
   TBranch        *b_gamma_phi;   //!
   TBranch        *b_dm_;   //!
   TBranch        *b_goodReco_;   //!
   TBranch        *b_tauMass_;   //!
   TBranch        *b_taupfTausDiscriminationByDecayModeFinding;   //!
   TBranch        *b_taupfTausDiscriminationByDecayModeFindingNewDMs;   //!
   TBranch        *b_tauByIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_tauByIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_tauByIsolationMVArun2v1PWnewDMwLTraw;   //!
   TBranch        *b_tauByIsolationMVArun2v1PWoldDMwLTraw;   //!
   TBranch        *b_tauChargedIsoPtSum;   //!
   TBranch        *b_tau_iso_pv;   //!
   TBranch        *b_tau_iso_dz0015;   //!
   TBranch        *b_tau_iso_dz0015_dz;   //!
   TBranch        *b_tau_iso_dz0015_dxy;   //!
   TBranch        *b_tau_iso_dz01_dz;   //!
   TBranch        *b_tau_iso_dz01_dxy;   //!
   TBranch        *b_tau_iso_dz01;   //!
   TBranch        *b_tauNeutralIsoPtSum;   //!
   TBranch        *b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauByTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tauCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_tauPuCorrPtSum;   //!
   TBranch        *b_taufootprintCorrection;   //!
   TBranch        *b_tauphotonPtSumOutsideSignalCone;   //!

   signal_pu0(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~signal_pu0() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(signal_pu0,0);
};

#endif

#ifdef signal_pu0_cxx
void signal_pu0::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("tauPt", &tauPt, &b_tauPt_);
   fChain->SetBranchAddress("tauEta", &tauEta, &b_tauEta_);
   fChain->SetBranchAddress("genTauPt", &genTauPt, &b_genTauPt_);
   fChain->SetBranchAddress("genTauEta", &genTauEta, &b_genTauEta_);
   fChain->SetBranchAddress("genTauMatch", &genTauMatch, &b_genTauMatch_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx_);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX_);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY_);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ_);
   fChain->SetBranchAddress("vtxIndex", &vtxIndex, &b_vtxIndex_);
   fChain->SetBranchAddress("dz_tt", &dz_tt, &b_dz_tt);
   fChain->SetBranchAddress("dxy_tt", &dxy_tt, &b_dxy_tt);
   fChain->SetBranchAddress("Nparticles", &Nparticles, &b_Nparticles);
   fChain->SetBranchAddress("DR_IsoCandTau", DR_IsoCandTau, &b_DR_IsoCandTau);
   fChain->SetBranchAddress("cand_pt", cand_pt, &b_cand_pt);
   fChain->SetBranchAddress("cand_eta", cand_eta, &b_cand_eta);
   fChain->SetBranchAddress("cand_phi", cand_phi, &b_cand_phi);
   fChain->SetBranchAddress("cand_q", cand_q, &b_cand_q);
   fChain->SetBranchAddress("cand_numHits", cand_numHits, &b_cand_numHits);
   fChain->SetBranchAddress("cand_dz_new", cand_dz_new, &b_cand_dz_new);
   fChain->SetBranchAddress("cand_dxy_new", cand_dxy_new, &b_cand_dxy_new);
   fChain->SetBranchAddress("Ngammas", &Ngammas, &b_Ngammas);
   fChain->SetBranchAddress("DR_gamma", DR_gamma, &b_DR_gamma);
   fChain->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
   fChain->SetBranchAddress("gamma1_pt", gamma1_pt, &b_gamma1_pt);
   fChain->SetBranchAddress("gamma2_pt", gamma2_pt, &b_gamma2_pt);
   fChain->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
   fChain->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
   fChain->SetBranchAddress("dm", &dm, &b_dm_);
   fChain->SetBranchAddress("goodReco", &goodReco, &b_goodReco_);
   fChain->SetBranchAddress("tauMass", &tauMass, &b_tauMass_);
   fChain->SetBranchAddress("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding, &b_taupfTausDiscriminationByDecayModeFinding);
   fChain->SetBranchAddress("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs, &b_taupfTausDiscriminationByDecayModeFindingNewDMs);
   fChain->SetBranchAddress("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw, &b_tauByIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw, &b_tauByIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw, &b_tauByIsolationMVArun2v1PWnewDMwLTraw);
   fChain->SetBranchAddress("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw, &b_tauByIsolationMVArun2v1PWoldDMwLTraw);
   fChain->SetBranchAddress("tauChargedIsoPtSum", &tauChargedIsoPtSum, &b_tauChargedIsoPtSum);
   fChain->SetBranchAddress("tau_iso_pv", &tau_iso_pv, &b_tau_iso_pv);
   fChain->SetBranchAddress("tau_iso_dz0015", &tau_iso_dz0015, &b_tau_iso_dz0015);
   fChain->SetBranchAddress("tau_iso_dz0015_dz", &tau_iso_dz0015_dz, &b_tau_iso_dz0015_dz);
   fChain->SetBranchAddress("tau_iso_dz0015_dxy", &tau_iso_dz0015_dxy, &b_tau_iso_dz0015_dxy);
   fChain->SetBranchAddress("tau_iso_dz01_dz", &tau_iso_dz01_dz, &b_tau_iso_dz01_dz);
   fChain->SetBranchAddress("tau_iso_dz01_dxy", &tau_iso_dz01_dxy, &b_tau_iso_dz01_dxy);
   fChain->SetBranchAddress("tau_iso_dz01", &tau_iso_dz01, &b_tau_iso_dz01);
   fChain->SetBranchAddress("tauNeutralIsoPtSum", &tauNeutralIsoPtSum, &b_tauNeutralIsoPtSum);
   fChain->SetBranchAddress("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits, &b_tauByLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits, &b_tauByMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits, &b_tauByTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits, &b_tauCombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("tauPuCorrPtSum", &tauPuCorrPtSum, &b_tauPuCorrPtSum);
   fChain->SetBranchAddress("taufootprintCorrection", &taufootprintCorrection, &b_taufootprintCorrection);
   fChain->SetBranchAddress("tauphotonPtSumOutsideSignalCone", &tauphotonPtSumOutsideSignalCone, &b_tauphotonPtSumOutsideSignalCone);
}

Bool_t signal_pu0::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef signal_pu0_cxx
