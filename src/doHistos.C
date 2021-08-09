#define doHistos_cxx
#include "SUEP_scouting/doHistos.h"
#include "src/kinematics.C"
#include "src/eventShapes.C"
#include "src/jetStudies.C"
#include "src/pileupRW.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


int get_charge(int pdgId){
  if (abs(pdgId) == 11) return 1; // electron
  else if (abs(pdgId) == 13) return 1; // muon
  else if (abs(pdgId) == 211) return 1; // pion
  return 0;
  // 130 = KLong - neutral hadron 
  // 22 = photon 
  // workbook also says
  // 1 = HF hadron, where HF means forward calo
  // 2 = HF em particle, where HF means forward calo
}
void doHistos::Loop(std::string s_sample,bool isMC)
{

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	std::cout << " Analyzing Tree : " << s_sample << " with " << nentries << " entries" << std::endl;

   // disable all branches
   fChain->SetBranchStatus("*",0);
   // only activate branches that are used (faster)
   std::vector<string> branches_used{
      "PFcand_pt",
      "PFcand_eta",
      "PFcand_phi",
      "PFcand_m",
      "PFcand_pdgid",
      "PFcand_vertex",
      "Jet_pt",
      "Jet_eta",
      "Jet_phi",
      "Jet_m",
      "Vertex_x",
      "Vertex_y",
      "Vertex_z",
      "Vertex_tracksSize", 
      "Vertex_chi2", 
      "Vertex_ndof", 
      "Vertex_isValidVtx",
      "trig", 
      "n_pvs",
      "Jet_chargedHadronEnergy", 
      "Jet_neutralHadronEnergy", 
      "Jet_photonEnergy", 
      "Jet_electronEnergy", 
      "Jet_muonEnergy", 
      "Jet_HFHadronEnergy", 
      "Jet_HFEMEnergy", 
      "Jet_HOEnergy", 
      "Jet_chargedHadronMultiplicity", 
      "Jet_neutralHadronMultiplicity", 
      "Jet_photonMultiplicity", 
      "Jet_electronMultiplicity", 
      "Jet_muonMultiplicity", 
      "Jet_HFHadronMultiplicity", 
      "Jet_HFEMMultiplicity", 
   };
   for(const auto& branch : branches_used) fChain->SetBranchStatus(branch.c_str(),1);

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		if (ientry%1000==0) std::cout << "Processed " << ientry << " events!" << std::endl;



    // Trigger info is an 8-bit word...
    // 
    //pass_DST_HT410_PFScouting = "{0:08b}".format(trig)[1]
    //pass_DST_HT450_PFScouting = "{0:08b}".format(trig)[0]
    //int trigger = byte((int)trig);
    //std::cout << "trig "  << trigger << std::endl;
    //std::cout << "pass_DST_HT410_PFScouting" << pass_DST_HT410_PFScouting << std::endl;
    //std::cout << "pass_DST_HT450_PFScouting" << pass_DST_HT450_PFScouting << std::endl;

    // * 
    // Packup jets
    // * 
    njets = 0;
    ht = 0;
    TLorentzVector jet_p4;
    jets.clear();
    for(unsigned int i=0; i < Jet_pt->size(); i++)
    {
      // basic cuts
      if ( Jet_pt->at(i) < 30 ) continue;
      if ( abs(Jet_eta->at(i) ) > 2.4 ) continue;

      // plot jet phi to double check ...

      
      // look into loose ID

      Jet jet;

      jet_p4.SetPtEtaPhiM(Jet_pt->at(i), Jet_eta->at(i), Jet_phi->at(i), Jet_m->at(i)); // assume pion mass
      jet.p4 = jet_p4;

      // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2018
      // For jet ID - confirmed energies sum up to jet_p4.E()
      jet.NHF  = Jet_neutralHadronEnergy->at(i) / jet_p4.E();
      jet.NEMF = Jet_photonEnergy->at(i)  / jet_p4.E();
      jet.CHF  = Jet_chargedHadronEnergy->at(i) / jet_p4.E();
      jet.MUF  = Jet_muonEnergy->at(i)  / jet_p4.E();
      jet.CEMF = Jet_electronEnergy->at(i) / jet_p4.E();
      jet.NumConst = Jet_chargedHadronMultiplicity->at(i) + Jet_neutralHadronMultiplicity->at(i) + Jet_photonMultiplicity->at(i) + Jet_electronMultiplicity->at(i) + Jet_muonMultiplicity->at(i) + Jet_HFHadronMultiplicity->at(i) + Jet_HFEMMultiplicity->at(i) ;
      jet.NumNeutralParticles = Jet_neutralHadronMultiplicity->at(i) + Jet_photonMultiplicity->at(i) ;  // not including HF? not used?
      jet.CHM      = Jet_chargedHadronMultiplicity->at(i) + Jet_electronMultiplicity->at(i) + Jet_muonMultiplicity->at(i) ; 

      jet.id = ( jet.CEMF<0.8 && jet.CHM>0 && jet.CHF>0 && jet.NumConst>1 && jet.NEMF<0.9 && jet.MUF <0.8 && jet.NHF < 0.9 );

      //std::cout << i << "pass ID " << jet.id << std::endl;

      if (jet.id) {
        jets.push_back(jet);
        njets+=1;
        ht += Jet_pt->at(i);        
      }
      
    }

    

		// * 
		// Packup final state particles
		// * 
		ntracks=0;
		TLorentzVector trk_p4;
		tracks.clear();
		for (unsigned int i = 0; i <PFcand_pt->size(); i++)
		{
	
			if (PFcand_vertex->at(i) != 0) continue;     
			if (abs(PFcand_eta->at(i)) > 2.5) continue; 
			if (PFcand_pt->at(i) < 1) continue;
			if (get_charge(PFcand_pdgid->at(i)) == 0) continue;
      //std::cout << i << std::endl;
      //std::cout << "eta " << PFcand_eta->at(i) << std::endl;
      //std::cout << "phi " << PFcand_phi->at(i) << std::endl;
      //std::cout << "pt " << PFcand_pt->at(i) << std::endl;
      //std::cout << "mass " << PFcand_m->at(i) << std::endl;
      //std::cout << "pdgid " << PFcand_pdgid->at(i) << std::endl;
      //std::cout << "vertex " << PFcand_vertex->at(i) << std::endl; 

			Track track;

			trk_p4.SetPtEtaPhiM(PFcand_pt->at(i), PFcand_eta->at(i), PFcand_phi->at(i), 0.13957); // assume pion mass
			track.p4 = trk_p4; 
			track.pdgId = PFcand_pdgid->at(i); 
		  	//track.isSuep = false;
			
			tracks.push_back(track);
			ntracks+=1;
      //plotter.Plot1D(Form("%s_all_tracks_pt"  ,s_sample.c_str()),";track p_{T} [GeV]", track.p4.Pt(), 100,0,10);
      //plotter.Plot1D(Form("%s_all_tracks_ptL" ,s_sample.c_str()),";track p_{T} [GeV]", track.p4.Pt(), 100,0,100);
		}


    vertices.clear();
    TLorentzVector vertex_p4;
    for (unsigned int i = 0; i <Vertex_x->size(); i++)
    {

      Vertex vertex;

      vertex_p4.SetXYZM(Vertex_x->at(i), Vertex_y->at(i), Vertex_z->at(i), 0); // assume pion mass
      vertex
      .p4 = vertex_p4; 

      vertex.tracksSize = Vertex_tracksSize->at(i);
      vertex.chi2       = Vertex_chi2->at(i);
      vertex.ndof       = Vertex_ndof->at(i);
      vertex.isValidVtx = Vertex_isValidVtx->at(i);
      vertices.push_back(vertex);


    }

    // Get Weight info...
    pu_wght = (isMC) ? get_PU( vertices.size() ) : 1.0;
    evt_wght = pu_wght;  
    ///std::cout << pu_wght << " " << evt_wght << " " << evt_wght/pu_wght << std::endl;


    // * 
    // Do Plotting  
    // * 
    basic_kinematics(s_sample,"all"); // no trigger requirement

    // Basic Requirements
    if (vertices.size() == 0) continue;

    if (ht < 500) continue; // may need to up this to be fully in the turn on?

    basic_kinematics(s_sample,"scouting");
    //plotEventShapes(s_sample,"scouting", tracks);
    fatjet_plots(s_sample,"scouting", tracks, jentry, 1.5);

    //if (ht < 1200) continue;
    //basic_kinematics(s_sample,"offline");
    //plotEventShapes(s_sample,"scouting",tracks);



   }// end loop over events
}

int main(int argc, char* argv[]){

    // defaults 
    std::string tree_name = "mmtree/tree";
    std::string file_name = "input/scoutingData18.root";
    std::string sample_name = "data18";
    std::string output_name = "data18";

    // Pick file
    // ./doHistos sample_name
    std::cout << argc << std::endl;
    if (argc > 2){
        sample_name = argv[1];
        output_name = argv[1];
        file_name = argv[2]; 

      //std::cout << "file_name " << file_name << std::endl;
      //std::cout << "sample_name " << sample_name << std::endl;
      //std::cout << "output_name " << output_name << std::endl;
    }

    std::cout << "Starting SUEP Studies!" << std::endl;

    gROOT->SetBatch();
    gStyle->SetOptStat(0);
    PlotHelper::setPlotStyle();

    TFile *file = TFile::Open(file_name.c_str());
    TTree *tree = (TTree*)file->Get(tree_name.c_str());

	// Figure out if is MC
    bool isMC = 1;
    if (sample_name=="data18") isMC = 0;

    // Do analysis
    doHistos analysis(tree,isMC);
    analysis.Loop(sample_name,isMC);

    // Save histograms here
    TFile *output_file;
    output_file = TFile::Open(Form("output/%s.root",output_name.c_str()),"RECREATE");
    c1->SetTickx(true);
    c1->SetTicky(true);
    plotter.DrawAll1D(c1);
    plotter.DrawAll2D(c1, "colz");

}
