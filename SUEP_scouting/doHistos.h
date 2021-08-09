//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep  9 20:50:52 2020 by ROOT version 6.12/06
// from TTree tree/tree
// found on file: input/scoutingQCD500to700.root
//////////////////////////////////////////////////////////

#ifndef doHistos_h
#define doHistos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include "SUEP_scouting/PhysicsObjects.h"
#include "SUEP_scouting/PlotHelper.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

// Header file for the classes stored in the TTree if any.
#include "vector"

PlotHelper plotter("");//set up the plotter
TCanvas *c1 = new TCanvas("c1","c1",800,800);

// helpful global variables
float evt_wght = 1.0;
float pu_wght = 1.0;
bool pass_DST_HT410_PFScouting = 0;
bool pass_DST_HT450_PFScouting = 0;
int ntracks=0;
int ntracks_09=0;
int ntracks_08=0;
int ntracks_07=0;
int ntracks_2=0;
int njets=0;
float ht=0;
float lead_jet_pt=0;
std::vector<Jet> jets;
std::vector<Track> tracks; 
std::vector<Vertex> vertices; 

class doHistos {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          lumSec;
   UInt_t          run;
   UChar_t         trig;
   vector<bool>    *l1Result;
   UInt_t          n_ele;
   std::vector<float>   *Electron_pt;
   std::vector<float>   *Electron_eta;
   std::vector<float>   *Electron_phi;
   std::vector<float>   *Electron_charge;
   std::vector<float>   *Electron_m;
   std::vector<float>   *Electron_tkiso;
   std::vector<float>   *Electron_HoE;
   std::vector<float>   *Electron_sigmaietaieta;
   std::vector<float>   *Electron_dphiin;
   std::vector<float>   *Electron_detain;
   std::vector<float>   *Electron_mHits;
   std::vector<float>   *Electron_ooEMOop;
   UInt_t          n_pho;
   std::vector<float>   *Photon_pt;
   std::vector<float>   *Photon_eta;
   std::vector<float>   *Photon_phi;
   std::vector<float>   *Photon_m;
   std::vector<float>   *Photon_hcaliso;
   std::vector<float>   *Photon_ecaliso;
   std::vector<float>   *Photon_HoE;
   std::vector<float>   *Photon_sigmaietaieta;
   UInt_t          n_pvs;
   std::vector<float>   *Vertex_x;
   std::vector<float>   *Vertex_y;
   std::vector<float>   *Vertex_z;
   std::vector<float>   *Vertex_tracksSize;
   std::vector<float>   *Vertex_chi2;
   std::vector<float>   *Vertex_ndof;
   std::vector<float>   *Vertex_isValidVtx;
   UInt_t          n_pfcand;
   std::vector<float>   *PFcand_pt;
   std::vector<float>   *PFcand_eta;
   std::vector<float>   *PFcand_phi;
   std::vector<float>   *PFcand_m;
   std::vector<float>   *PFcand_pdgid;
   std::vector<float>   *PFcand_vertex;
   UInt_t          n_mu;
   std::vector<float>   *Muon_pt;
   std::vector<float>   *Muon_eta;
   std::vector<float>   *Muon_phi;
   std::vector<float>   *Muon_m;
   std::vector<float>   *Muon_ecaliso;
   std::vector<float>   *Muon_hcaliso;
   std::vector<float>   *Muon_trkiso;
   std::vector<float>   *Muon_chi2;
   std::vector<float>   *Muon_ndof;
   std::vector<float>   *Muon_charge;
   std::vector<float>   *Muon_dxy;
   std::vector<float>   *Muon_dz;
   std::vector<float>   *Muon_nvalidmuon_hits;
   std::vector<float>   *Muon_validpixelhits;
   std::vector<float>   *Muon_nmatchedstations;
   std::vector<float>   *Muon_type;
   std::vector<float>   *Muon_nvalidstriphits;
   std::vector<float>   *Muon_trkqoverp;
   std::vector<float>   *Muon_trklambda;
   std::vector<float>   *Muon_trkpt;
   std::vector<float>   *Muon_trkphi;
   std::vector<float>   *Muon_trketa;
   std::vector<float>   *Muon_trkqoverperror;
   std::vector<float>   *Muon_trklambdaerror;
   std::vector<float>   *Muon_trkpterror;
   std::vector<float>   *Muon_trkphierror;
   std::vector<float>   *Muon_trketaerror;
   std::vector<float>   *Muon_trkdszerror;
   std::vector<float>   *Muon_trkdsz;
   UInt_t          n_jet;
   std::vector<float>   *Jet_pt;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_phi;
   std::vector<float>   *Jet_m;
   std::vector<float>   *Jet_area;
   std::vector<float>   *Jet_chargedHadronEnergy;
   std::vector<float>   *Jet_neutralHadronEnergy;
   std::vector<float>   *Jet_photonEnergy;
   std::vector<float>   *Jet_electronEnergy;
   std::vector<float>   *Jet_muonEnergy;
   std::vector<float>   *Jet_HFHadronEnergy;
   std::vector<float>   *Jet_HFEMEnergy;
   std::vector<float>   *Jet_HOEnergy;
   std::vector<float>   *Jet_chargedHadronMultiplicity;
   std::vector<float>   *Jet_neutralHadronMultiplicity;
   std::vector<float>   *Jet_photonMultiplicity;
   std::vector<float>   *Jet_electronMultiplicity;
   std::vector<float>   *Jet_muonMultiplicity;
   std::vector<float>   *Jet_HFHadronMultiplicity;
   std::vector<float>   *Jet_HFEMMultiplicity;
   std::vector<float>   *Jet_csv;
   std::vector<float>   *Jet_mvaDiscriminator;
   std::vector<vector<short> > *Jet_constituents;
   std::vector<float>   *FatJet_area;
   std::vector<float>   *FatJet_eta;
   std::vector<float>   *FatJet_n2b1;
   std::vector<float>   *FatJet_n3b1;
   std::vector<float>   *FatJet_phi;
   std::vector<float>   *FatJet_pt;
   std::vector<float>   *FatJet_tau1;
   std::vector<float>   *FatJet_tau2;
   std::vector<float>   *FatJet_tau3;
   std::vector<float>   *FatJet_tau4;
   std::vector<float>   *FatJet_mass;
   std::vector<float>   *FatJet_msoftdrop;
   std::vector<float>   *FatJet_mtrim;

   // List of branches
   TBranch        *b_lumSec;   //!
   TBranch        *b_run;   //!
   TBranch        *b_trig;   //!
   TBranch        *b_l1Result;   //!
   TBranch        *b_n_ele;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_m;   //!
   TBranch        *b_Electron_tkiso;   //!
   TBranch        *b_Electron_HoE;   //!
   TBranch        *b_Electron_sigmaietaieta;   //!
   TBranch        *b_Electron_dphiin;   //!
   TBranch        *b_Electron_detain;   //!
   TBranch        *b_Electron_mHits;   //!
   TBranch        *b_Electron_ooEMOop;   //!
   TBranch        *b_n_pho;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_m;   //!
   TBranch        *b_Photon_hcaliso;   //!
   TBranch        *b_Photon_ecaliso;   //!
   TBranch        *b_Photon_HoE;   //!
   TBranch        *b_Photon_sigmaietaieta;   //!
   TBranch        *b_n_pvs; //!
   TBranch        *b_Vertex_x; //!
   TBranch        *b_Vertex_y; //!
   TBranch        *b_Vertex_z; //!
   TBranch        *b_Vertex_tracksSize; //!
   TBranch        *b_Vertex_chi2; //!
   TBranch        *b_Vertex_ndof; //!
   TBranch        *b_Vertex_isValidVtx; //!
   TBranch        *b_n_pfcand;   //!
   TBranch        *b_PFcand_pt;   //!
   TBranch        *b_PFcand_eta;   //!
   TBranch        *b_PFcand_phi;   //!
   TBranch        *b_PFcand_m;   //!
   TBranch        *b_PFcand_pdgid;   //!
   TBranch        *b_PFcand_vertex;   //!
   TBranch        *b_n_mu;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_m;   //!
   TBranch        *b_Muon_ecaliso;   //!
   TBranch        *b_Muon_hcaliso;   //!
   TBranch        *b_Muon_trkiso;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_ndof;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_nvalidmuon_hits;   //!
   TBranch        *b_Muon_validpixelhits;   //!
   TBranch        *b_Muon_nmatchedstations;   //!
   TBranch        *b_Muon_type;   //!
   TBranch        *b_Muon_nvalidstriphits;   //!
   TBranch        *b_Muon_trkqoverp;   //!
   TBranch        *b_Muon_trklambda;   //!
   TBranch        *b_Muon_trkpt;   //!
   TBranch        *b_Muon_trkphi;   //!
   TBranch        *b_Muon_trketa;   //!
   TBranch        *b_Muon_trkqoverperror;   //!
   TBranch        *b_Muon_trklambdaerror;   //!
   TBranch        *b_Muon_trkpterror;   //!
   TBranch        *b_Muon_trkphierror;   //!
   TBranch        *b_Muon_trketaerror;   //!
   TBranch        *b_Muon_trkdszerror;   //!
   TBranch        *b_Muon_trkdsz;   //!
   TBranch        *b_n_jet;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_m;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_chargedHadronEnergy;   //!
   TBranch        *b_Jet_neutralHadronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_muonEnergy;   //!
   TBranch        *b_Jet_HFHadronEnergy;   //!
   TBranch        *b_Jet_HFEMEnergy;   //!
   TBranch        *b_Jet_HOEnergy;   //!
   TBranch        *b_Jet_chargedHadronMultiplicity;   //!
   TBranch        *b_Jet_neutralHadronMultiplicity;   //!
   TBranch        *b_Jet_photonMultiplicity;   //!
   TBranch        *b_Jet_electronMultiplicity;   //!
   TBranch        *b_Jet_muonMultiplicity;   //!
   TBranch        *b_Jet_HFHadronMultiplicity;   //!
   TBranch        *b_Jet_HFEMMultiplicity;   //!
   TBranch        *b_Jet_csv;   //!
   TBranch        *b_Jet_mvaDiscriminator;   //!
   TBranch        *b_Jet_constituents;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_mtrim;   //!

   doHistos(TTree *tree=0, bool isMC=0);
   virtual ~doHistos();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string s_sample,bool isMC);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef doHistos_cxx
doHistos::doHistos(TTree *tree, bool isMC) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("input/scoutingQCD500to700.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("input/scoutingQCD500to700.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("input/scoutingQCD500to700.root:/mmtree");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

doHistos::~doHistos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t doHistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t doHistos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void doHistos::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   l1Result = 0;
   Electron_pt = 0;
   Electron_eta = 0;
   Electron_phi = 0;
   Electron_charge = 0;
   Electron_m = 0;
   Electron_tkiso = 0;
   Electron_HoE = 0;
   Electron_sigmaietaieta = 0;
   Electron_dphiin = 0;
   Electron_detain = 0;
   Electron_mHits = 0;
   Electron_ooEMOop = 0;
   Photon_pt = 0;
   Photon_eta = 0;
   Photon_phi = 0;
   Photon_m = 0;
   Photon_hcaliso = 0;
   Photon_ecaliso = 0;
   Photon_HoE = 0;
   Photon_sigmaietaieta = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   Vertex_tracksSize = 0;
   Vertex_chi2 = 0;
   Vertex_ndof = 0;
   Vertex_isValidVtx = 0;
   PFcand_pt = 0;
   PFcand_eta = 0;
   PFcand_phi = 0;
   PFcand_m = 0;
   PFcand_pdgid = 0;
   PFcand_vertex = 0;
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_m = 0;
   Muon_ecaliso = 0;
   Muon_hcaliso = 0;
   Muon_trkiso = 0;
   Muon_chi2 = 0;
   Muon_ndof = 0;
   Muon_charge = 0;
   Muon_dxy = 0;
   Muon_dz = 0;
   Muon_nvalidmuon_hits = 0;
   Muon_validpixelhits = 0;
   Muon_nmatchedstations = 0;
   Muon_type = 0;
   Muon_nvalidstriphits = 0;
   Muon_trkqoverp = 0;
   Muon_trklambda = 0;
   Muon_trkpt = 0;
   Muon_trkphi = 0;
   Muon_trketa = 0;
   Muon_trkqoverperror = 0;
   Muon_trklambdaerror = 0;
   Muon_trkpterror = 0;
   Muon_trkphierror = 0;
   Muon_trketaerror = 0;
   Muon_trkdszerror = 0;
   Muon_trkdsz = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_m = 0;
   Jet_area = 0;
   Jet_chargedHadronEnergy = 0;
   Jet_neutralHadronEnergy = 0;
   Jet_photonEnergy = 0;
   Jet_electronEnergy = 0;
   Jet_muonEnergy = 0;
   Jet_HFHadronEnergy = 0;
   Jet_HFEMEnergy = 0;
   Jet_HOEnergy = 0;
   Jet_chargedHadronMultiplicity = 0;
   Jet_neutralHadronMultiplicity = 0;
   Jet_photonMultiplicity = 0;
   Jet_electronMultiplicity = 0;
   Jet_muonMultiplicity = 0;
   Jet_HFHadronMultiplicity = 0;
   Jet_HFEMMultiplicity = 0;
   Jet_csv = 0;
   Jet_mvaDiscriminator = 0;
   Jet_constituents = 0;
   FatJet_area = 0;
   FatJet_eta = 0;
   FatJet_n2b1 = 0;
   FatJet_n3b1 = 0;
   FatJet_phi = 0;
   FatJet_pt = 0;
   FatJet_tau1 = 0;
   FatJet_tau2 = 0;
   FatJet_tau3 = 0;
   FatJet_tau4 = 0;
   FatJet_mass = 0;
   FatJet_msoftdrop = 0;
   FatJet_mtrim = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("lumSec", &lumSec, &b_lumSec);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("trig", &trig, &b_trig);
   fChain->SetBranchAddress("l1Result", &l1Result, &b_l1Result);
   fChain->SetBranchAddress("n_ele", &n_ele, &b_n_ele);
   fChain->SetBranchAddress("Electron_pt", &Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_eta", &Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_phi", &Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_charge", &Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_m", &Electron_m, &b_Electron_m);
   fChain->SetBranchAddress("Electron_tkiso", &Electron_tkiso, &b_Electron_tkiso);
   fChain->SetBranchAddress("Electron_HoE", &Electron_HoE, &b_Electron_HoE);
   fChain->SetBranchAddress("Electron_sigmaietaieta", &Electron_sigmaietaieta, &b_Electron_sigmaietaieta);
   fChain->SetBranchAddress("Electron_dphiin", &Electron_dphiin, &b_Electron_dphiin);
   fChain->SetBranchAddress("Electron_detain", &Electron_detain, &b_Electron_detain);
   fChain->SetBranchAddress("Electron_mHits", &Electron_mHits, &b_Electron_mHits);
   fChain->SetBranchAddress("Electron_ooEMOop", &Electron_ooEMOop, &b_Electron_ooEMOop);
   fChain->SetBranchAddress("n_pho", &n_pho, &b_n_pho);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_m", &Photon_m, &b_Photon_m);
   fChain->SetBranchAddress("Photon_hcaliso", &Photon_hcaliso, &b_Photon_hcaliso);
   fChain->SetBranchAddress("Photon_ecaliso", &Photon_ecaliso, &b_Photon_ecaliso);
   fChain->SetBranchAddress("Photon_HoE", &Photon_HoE, &b_Photon_HoE);
   fChain->SetBranchAddress("Photon_sigmaietaieta", &Photon_sigmaietaieta, &b_Photon_sigmaietaieta);
   fChain->SetBranchAddress("n_pvs", &n_pvs, &b_n_pvs);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex_tracksSize", &Vertex_tracksSize, &b_Vertex_tracksSize);
   fChain->SetBranchAddress("Vertex_chi2", &Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_ndof", &Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_isValidVtx", &Vertex_isValidVtx, &b_Vertex_isValidVtx);
   fChain->SetBranchAddress("n_pfcand", &n_pfcand, &b_n_pfcand);
   fChain->SetBranchAddress("PFcand_pt", &PFcand_pt, &b_PFcand_pt);
   fChain->SetBranchAddress("PFcand_eta", &PFcand_eta, &b_PFcand_eta);
   fChain->SetBranchAddress("PFcand_phi", &PFcand_phi, &b_PFcand_phi);
   fChain->SetBranchAddress("PFcand_m", &PFcand_m, &b_PFcand_m);
   fChain->SetBranchAddress("PFcand_pdgid", &PFcand_pdgid, &b_PFcand_pdgid);
   fChain->SetBranchAddress("PFcand_vertex", &PFcand_vertex, &b_PFcand_vertex);
   fChain->SetBranchAddress("n_mu", &n_mu, &b_n_mu);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_m", &Muon_m, &b_Muon_m);
   fChain->SetBranchAddress("Muon_ecaliso", &Muon_ecaliso, &b_Muon_ecaliso);
   fChain->SetBranchAddress("Muon_hcaliso", &Muon_hcaliso, &b_Muon_hcaliso);
   fChain->SetBranchAddress("Muon_trkiso", &Muon_trkiso, &b_Muon_trkiso);
   fChain->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   fChain->SetBranchAddress("Muon_ndof", &Muon_ndof, &b_Muon_ndof);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_nvalidmuon_hits", &Muon_nvalidmuon_hits, &b_Muon_nvalidmuon_hits);
   fChain->SetBranchAddress("Muon_validpixelhits", &Muon_validpixelhits, &b_Muon_validpixelhits);
   fChain->SetBranchAddress("Muon_nmatchedstations", &Muon_nmatchedstations, &b_Muon_nmatchedstations);
   fChain->SetBranchAddress("Muon_type", &Muon_type, &b_Muon_type);
   fChain->SetBranchAddress("Muon_nvalidstriphits", &Muon_nvalidstriphits, &b_Muon_nvalidstriphits);
   fChain->SetBranchAddress("Muon_trkqoverp", &Muon_trkqoverp, &b_Muon_trkqoverp);
   fChain->SetBranchAddress("Muon_trklambda", &Muon_trklambda, &b_Muon_trklambda);
   fChain->SetBranchAddress("Muon_trkpt", &Muon_trkpt, &b_Muon_trkpt);
   fChain->SetBranchAddress("Muon_trkphi", &Muon_trkphi, &b_Muon_trkphi);
   fChain->SetBranchAddress("Muon_trketa", &Muon_trketa, &b_Muon_trketa);
   fChain->SetBranchAddress("Muon_trkqoverperror", &Muon_trkqoverperror, &b_Muon_trkqoverperror);
   fChain->SetBranchAddress("Muon_trklambdaerror", &Muon_trklambdaerror, &b_Muon_trklambdaerror);
   fChain->SetBranchAddress("Muon_trkpterror", &Muon_trkpterror, &b_Muon_trkpterror);
   fChain->SetBranchAddress("Muon_trkphierror", &Muon_trkphierror, &b_Muon_trkphierror);
   fChain->SetBranchAddress("Muon_trketaerror", &Muon_trketaerror, &b_Muon_trketaerror);
   fChain->SetBranchAddress("Muon_trkdszerror", &Muon_trkdszerror, &b_Muon_trkdszerror);
   fChain->SetBranchAddress("Muon_trkdsz", &Muon_trkdsz, &b_Muon_trkdsz);
   fChain->SetBranchAddress("n_jet", &n_jet, &b_n_jet);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_m", &Jet_m, &b_Jet_m);
   fChain->SetBranchAddress("Jet_area", &Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_chargedHadronEnergy", &Jet_chargedHadronEnergy, &b_Jet_chargedHadronEnergy);
   fChain->SetBranchAddress("Jet_neutralHadronEnergy", &Jet_neutralHadronEnergy, &b_Jet_neutralHadronEnergy);
   fChain->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
   fChain->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
   fChain->SetBranchAddress("Jet_muonEnergy", &Jet_muonEnergy, &b_Jet_muonEnergy);
   fChain->SetBranchAddress("Jet_HFHadronEnergy", &Jet_HFHadronEnergy, &b_Jet_HFHadronEnergy);
   fChain->SetBranchAddress("Jet_HFEMEnergy", &Jet_HFEMEnergy, &b_Jet_HFEMEnergy);
   fChain->SetBranchAddress("Jet_HOEnergy", &Jet_HOEnergy, &b_Jet_HOEnergy);
   fChain->SetBranchAddress("Jet_chargedHadronMultiplicity", &Jet_chargedHadronMultiplicity, &b_Jet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jet_neutralHadronMultiplicity", &Jet_neutralHadronMultiplicity, &b_Jet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jet_photonMultiplicity", &Jet_photonMultiplicity, &b_Jet_photonMultiplicity);
   fChain->SetBranchAddress("Jet_electronMultiplicity", &Jet_electronMultiplicity, &b_Jet_electronMultiplicity);
   fChain->SetBranchAddress("Jet_muonMultiplicity", &Jet_muonMultiplicity, &b_Jet_muonMultiplicity);
   fChain->SetBranchAddress("Jet_HFHadronMultiplicity", &Jet_HFHadronMultiplicity, &b_Jet_HFHadronMultiplicity);
   fChain->SetBranchAddress("Jet_HFEMMultiplicity", &Jet_HFEMMultiplicity, &b_Jet_HFEMMultiplicity);
   fChain->SetBranchAddress("Jet_csv", &Jet_csv, &b_Jet_csv);
   fChain->SetBranchAddress("Jet_mvaDiscriminator", &Jet_mvaDiscriminator, &b_Jet_mvaDiscriminator);
   fChain->SetBranchAddress("Jet_constituents", &Jet_constituents, &b_Jet_constituents);
   fChain->SetBranchAddress("FatJet_area", &FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_eta", &FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_n2b1", &FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", &FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_phi", &FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", &FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_tau1", &FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", &FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", &FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", &FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_mass", &FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", &FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_mtrim", &FatJet_mtrim, &b_FatJet_mtrim);
   Notify();
}

Bool_t doHistos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void doHistos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t doHistos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef doHistos_cxx
