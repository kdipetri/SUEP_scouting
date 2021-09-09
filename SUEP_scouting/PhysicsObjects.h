#ifndef physicsObjects_h
#define physicsObjects_h

// Structures for physics objects that don't require doHistos class 
// Useful for packing up n-tuple variables for a single muon
// Can loop through muons more easily to figure 

struct Vertex {
  unsigned int index;
  TLorentzVector p4;
  int tracksSize;
  float chi2;
  float ndof;
  float isValidVtx;
};

struct Jet {
  unsigned int index;
  TLorentzVector p4;
  bool id;
  float NHF; 
  float NEMF; 
  float CHF; 
  float MUF; 
  float CEMF; 
  float NumConst; 
  float NumNeutralParticles; 
  float CHM;  
};

struct FatJet {
  unsigned int index;
  TLorentzVector p4;
  int nconstituents;
};

struct Track {
  unsigned int index;
  TLorentzVector p4;
  int pdgId;
  bool fromPV;
  bool highPurity;
  bool isPFcand;
  bool isSuep;
};

struct SUEP_particle {
  unsigned int index;
  TLorentzVector p4;
  int pdgId;
  bool pt1GeV;
  bool eta2p5;
  bool isReco;
  bool charge;
};

struct SUEP_Jet {// jet of all truth particle pt > 1 GeV eta, charge from suep
  TLorentzVector p4;
  int nTruthTracks=0;
};

#endif