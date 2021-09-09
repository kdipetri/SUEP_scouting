#define jetStudies_cxx

using namespace fastjet;
using namespace contrib;

std::string evtntrackregion(int ntracks)
{
	if      (ntracks < 10) return "0to10";
	else if (ntracks < 30) return "10to30";
	else if (ntracks < 50) return "30to50";
	else if (ntracks < 70) return "50to70";
	else return "70toinf";
}
std::string jetntrackregion(int ntracks)
{
	if      (ntracks < 5 ) return "0to5";
	else if (ntracks < 10) return "5to10";
	else if (ntracks < 20) return "10to20";
	else if (ntracks < 30) return "20to30";
	else if (ntracks < 50) return "30to50";
	else return "50toinf";
}
std::string rhoregion(float rho)
{
	if      (rho < 0.5) return "0to0p5";
	else if (rho < 1.0) return "0p5to1";
	else if (rho < 2.0) return "1to2";
	else if (rho < 5.0) return "2to5";
	else if (rho < 10.0) return "5to10";
	else return "inf";
}
std::string rho0region(float rho)
{
	if (rho < 12.0) return "low";
	else return "high";
}
std::string ptregion(float pt)
{
	if      (pt < 175) return "150to175";
	else if (pt < 200) return "175to200";
	else if (pt < 250) return "200to250";
	else return "inf";
}
double delta_r(PseudoJet jet1, PseudoJet jet2)
{
    double dphi = abs(jet1.phi() - jet2.phi());
    if (dphi > pi) {dphi = twopi - dphi;}
    double deta = jet1.eta() - jet2.eta();
    return (dphi*dphi + deta*deta);
}
/// This function needs the jet under consideration, a radius, and deltaR
/// defines a ring with thickness deltaR around the jet axis,
/// 	inner radius ra = r-deltaR/2
/// 	outer radius rb = r+deltaR/2
/// The distance of a given particle i from the jet axis can be expressed in terms of the azimuthal angle ϕi and the pseudorapidity ηi as
/// In the current analysis we used a binning of δr=0.05 and restricted ourselves to the fiducial region r < 0.6 to avoid edge effects
// from https://arxiv.org/pdf/2008.08500.pdf
vector<float> jet_rhos(PseudoJet jet, float dR, int nSteps) 
{
	// define radii & initialize
	vector<float> Rs;
	vector<float> pTsums;
	float R = dR/2.0;
	for (int step=0; step<nSteps; step++){
		Rs.push_back(R);
		pTsums.push_back(0);
		R += dR;
	}

    // get constituents
    std::vector<PseudoJet> constituents = jet.constituents();
    if (constituents.size() == 0 ) return pTsums;
    
    // scalar sum of contituent pT if inside our ring
    for (PseudoJet constituent : constituents ) 
    {
    	float dist = delta_r(jet,constituent);
    	for (int i=0; i<nSteps; i++){
      		if ( dist > ( Rs[i] + dR/2.0 ) ) continue; 
      		if ( dist < ( Rs[i] - dR/2.0 ) ) continue; 
      		pTsums[i] += constituent.pt();
      		break;
    	}

    }

    // compute rho
    vector<float> rhos;
    for (float pTsum : pTsums)
    {	
    	rhos.push_back(pTsum / jet.pt() / dR);
    	//std::cout << Rs[i] << std::endl;
    	//std::cout << pTsums[i]  << std::endl;
    	//float rho = pTsum / jet.pt() / dR;
    }
    
    
    return rhos;
}
// mean interparticle distance
float meanDR(PseudoJet jet){
	float mean_dR = 0;
	float n_pairs = 0;
	for (auto track1 : jet.constituents()){
		for (auto track2 : jet.constituents()){
			if (track1 == track2) continue;
			float dist = delta_r(track1,track2);
			mean_dR += dist;
			n_pairs += 1;

		}
	}
	mean_dR = mean_dR/n_pairs;
	return mean_dR;
}
//
// mean interparticle distance
float meanMinDR(PseudoJet jet){
	float mean_min_dR = 0;
	float n_pairs = 0;
	for (auto track1 : jet.constituents()){
		float min_dR = 10;
		for (auto track2 : jet.constituents()){
			if (track1 == track2) continue;

			float dR = delta_r(track1,track2);
			if (dR < min_dR) min_dR = dR; 
		}
		mean_min_dR += min_dR;
		n_pairs += 1;
	}
	mean_min_dR = mean_min_dR/n_pairs;
	return mean_min_dR;
}
//
float nsubjettiness(PseudoJet jet, int n, float R) 
{// https://githubmemory.com/repo/pkomiske/Nsubjettiness

	double beta = 1.0; // can try beta = 2

	Nsubjettiness nSub(n, 
		OnePass_WTA_KT_Axes(), 
		NormalizedMeasure(beta, R));

    
    return nSub(jet);
}
float nsubjettinessratio(PseudoJet jet, int n1, int n2, float R) 
{// https://githubmemory.com/repo/pkomiske/Nsubjettiness

	double beta = 1.0; // can try beta = 2

	NsubjettinessRatio nSubRatio(n1,n2,
                          OnePass_WTA_KT_Axes(),
                          NormalizedMeasure(beta, R));
    
    return nSubRatio(jet);
}
double width(PseudoJet jet)
{
	float num =0;
	float den =0;
    // get constituents
    std::vector<PseudoJet> constituents = jet.constituents();
    if (constituents.size() == 0 ) return -1;
  
    for (PseudoJet constituent: constituents) 
    {
  
  	    num += delta_r(constituent,jet)*constituent.pt();
  	    den += constituent.pt();
  
    }
    return num/den; 
}

void makeAKDisplay(PseudoJet suep_jet, const vector<fastjet::PseudoJet> jets, std::string sample, Long64_t ievent){

	std::vector<TGraph*> graphs;
	std::vector<TGraph*> jet_graphs; 
	
	//std::vector<TGraph2D> graphs;
	int k=0;
	for (unsigned i = 0; i < jets.size(); i++){
		
		// draw all constituents
		TGraph *graph = new TGraph();
		graph->SetTitle(Form("g_%s_%lli_%i",sample.c_str(), ievent, i));
		//TGraph2D *graph = new TGraph2D();

		std::vector<PseudoJet> constituents = jets[i].constituents();
		for (unsigned j = 0; j < constituents.size(); j++) {
      		graph->SetPoint(j,constituents[j].eta(),constituents[j].phi_std());
   		}

   		//graph->SetMarkerSize(color.at(i));
   		if (jets[i] == suep_jet)  graph->SetMarkerColor(kBlack);
   		else graph->SetMarkerColor(color.at(i));
   		graph->SetMarkerStyle(21);
		graphs.push_back(graph);

		// jet graphs, only draw high pT jets
		if (jets[i].pt() < 100 ) continue;

		TGraph *jet_graph = new TGraph();
		jet_graph->SetTitle(Form("g_jets_%s_%lli_%i",sample.c_str(), ievent, k));
		jet_graph->SetPoint(k, jets[i].eta(), jets[i].phi_std());
		jet_graph->SetMarkerStyle(24);
		jet_graph->SetMarkerSize(40);
		if (jets[i] == suep_jet)  jet_graph->SetMarkerColor(kBlack);
   		else jet_graph->SetMarkerColor(color.at(i));
   		jet_graphs.push_back(jet_graph);
   		k++;
	}
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle(";#eta; #phi");
	for (auto gr: graphs){
	//for (TGraph2D gr: graphs){
		mg->Add(gr);		
	}
	for (auto gr: jet_graphs){
		mg->Add(gr);
	}
        

	c1->cd();
	mg->Draw("ap");
	mg->GetYaxis()->SetLimits(-3.5,3.5);          
	mg->GetXaxis()->SetLimits(-3.5,3.5);// along  
	c1->Update();
	c1->Print(Form("plots/rhoStudy/eventDisplay/%s_%lli.png", sample.c_str(), ievent));
}

vector<FatJet> fatjet_plots(std::string sample, std::string sel, std::vector<Track> tracks, Long64_t ievent, float R=0.8){

	// Get the particles ready
	std::vector<PseudoJet> particles;
	// an event with  particles: px py pz E
	for (auto track : tracks){
		particles.push_back( PseudoJet( track.p4.Px(), track.p4.Py(), track.p4.Pz(), track.p4.E()) );
	}
	
	// choose a jet definition
	int cone = R*10;
	JetDefinition jet_def(antikt_algorithm, R);
	double ghost_maxrap = 2.5; // e.g. if particles go up to y=5
	//AreaDefinition area_def(passive_area, GhostedAreaSpec(ghost_maxrap));

	// run the clustering, extract the jets
	ClusterSequence cs(particles, jet_def);
	//ClusterSequenceArea cs(particles, jet_def, area_def);
	std::vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
	std::vector<FatJet> fat_jets = {}; // for using outside of this function

	// *
	// Make some nice plots...
	// *
	int n_fatjets=0;
	float dR = 6.0;
	PseudoJet suep_jet;
	int max_nconsit=0;

	///if (ievent < 200) makeAKDisplay(jets, sample, cone, ievent);
	FatJet fat_jet;
	for (unsigned i = 0; i < jets.size(); i++) {

		// basic cleaning cuts 
		if ( jets[i].pt() < 30 ) continue; // need some min jet cut
		// some min # tracks cut ? > 1 

		n_fatjets+=1;

		//For saving & plotting outside this function
		fat_jet.p4.SetPtEtaPhiM( jets[i].pt(), jets[i].eta(), jets[i].phi_std(), jets[i].m());
		fat_jet.nconstituents = jets[i].constituents().size();
		fat_jets.push_back(fat_jet);

		// the basics...
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_pt" , sample.c_str(),sel.c_str(),cone),";jet pt" , jets[i].pt()      , 100, 0, 1000 );
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_eta", sample.c_str(),sel.c_str(),cone),";jet eta", jets[i].eta()     , 100, -3.5, 3.5 );
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_phi", sample.c_str(),sel.c_str(),cone),";jet phi", jets[i].phi_std() , 100, -3.5, 3.5 );
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_m"  , sample.c_str(),sel.c_str(),cone),";jet m"  , jets[i].m() 		, 100, 0, 2000 );
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_moverpt"  , sample.c_str(),sel.c_str(),cone),";jet mass/pt"  , jets[i].m()/jets[i].pt() , 100, 0, 5 );

		// add plots of other variables here
		TLorentzVector jet_p4;
		jet_p4.SetPtEtaPhiM(jets[i].pt(),jets[i].eta(),jets[i].phi_std(),jets[i].m());
		//float dR_tmp= jet_p4.DeltaR(scalar);
		//plotter.Plot1D(Form( "%s_%s_jetsAK%i_dRscalar"  , sample.c_str(),sel.c_str(),cone),";jet,scalar dR"  , dR_tmp		, 50, 0, 6.0 );

		// Event displays
		//plotter.Plot2D(Form("%s_evt%lli_event_display_jetsAK%i",sample_name.c_str(),ievent,cone),";eta;phi;pt", jets[i].eta(), jets[i].phi_std(), 100, 3.5,3.5,100,3.5,3.5 , jets[i].perp());

		// Constituent based plots
		vector<PseudoJet> constituents = jets[i].constituents();
		for (unsigned j = 0; j < constituents.size(); j++) {
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_constit_pt", sample.c_str(),sel.c_str(),cone),";constit pt", constituents[j].pt(), 100, 0, 100 );

		}
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_nconstit", sample.c_str(),sel.c_str(),cone),";n constit.", constituents.size(), 100, 0, 500 );
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_width"   , sample.c_str(),sel.c_str(),cone),";jet width" , width(jets[i]), 100, 0, 3.0 );
		//plotter.Plot1D(Form( "%s_%s_jetsAK%i_rho"     , sample.c_str(),sel.c_str(),cone),";jet rho"   , rho(jets[i]), 100, 0, 10.0 );
		
		if (constituents.size() > max_nconsit) {
			max_nconsit = constituents.size();
			suep_jet = jets[i];
		}
	}
	plotter.Plot1D(Form("%s_%s_jetsAK%i_njets", sample.c_str(),sel.c_str(),cone),";njets", n_fatjets, 20, -0.5, 19.5 );

	// * 
	// Now make suep jet plots
	// * 
	if (max_nconsit > 0){// then we found a suep jet
		// the basics...
		plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_pt" , sample.c_str(),sel.c_str(),cone),";jet pt" , suep_jet.pt()      , 100, 0, 2000 );
		
		if ( suep_jet.pt() > 150 ) {

			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_eta", sample.c_str(),sel.c_str(),cone),";jet eta", suep_jet.eta()     , 100, -3.5, 3.5 );
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_phi", sample.c_str(),sel.c_str(),cone),";jet phi", suep_jet.phi_std() , 100, -3.5, 3.5 );
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_m"  , sample.c_str(),sel.c_str(),cone),";jet m"  , suep_jet.m() 		, 100, 0, 2000 );
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_moverpt" , sample.c_str(),sel.c_str(),cone),";jet mass/pt" , suep_jet.m()/suep_jet.pt(), 100, 0, 5.0 );
			
			// add plots of other variables here
			TLorentzVector suep_p4;
			suep_p4.SetPtEtaPhiM(suep_jet.pt(),suep_jet.eta(),suep_jet.phi_std(),suep_jet.m());
			//float dR_tmp= suep_p4.DeltaR(scalar);
			//plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_dRscalar"  , sample.c_str(),sel.c_str(),cone),";jet,scalar dR"  , dR_tmp		, 50, 0, 6.0 );
		
			// Constituent based plots
			vector<PseudoJet> constituents = suep_jet.constituents();
            for (PseudoJet constituent : constituents)
            {
  	            plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_constit_pt", sample.c_str(),sel.c_str(),cone),";constit pt", constituent.pt(), 100, 0, 100 );

            }
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nconstit", sample.c_str(),sel.c_str(),cone),";n constit.", constituents.size(), 100, 0, 500 );
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_width"   , sample.c_str(),sel.c_str(),cone),";jet width" , width(suep_jet), 100, 0, 3.0 );		
            
            //
            // jet momentum profiles	
            //
			vector<float> rhos = jet_rhos(suep_jet,0.05,15);
			std::string evtntracksel = evtntrackregion(tracks.size());
			std::string jetntracksel = jetntrackregion(constituents.size());
			std::string jetptsel     = ptregion(suep_jet.pt());
			std::string rhosel;    
			std::string rho0sel = rho0region(rhos[0]);
            for (int i=0; i<rhos.size() ; i++){
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i" , sample.c_str(),sel.c_str(),cone,i),Form(";jet rho%i",i) , rhos[i], 100, 0, 20.0 );		

                // 2D plots 
			    plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_rho_dR0p05_%i" , sample.c_str(),sel.c_str(),cone,i),Form(";jet rho%i;ntrk"  ,i) , rhos[i], tracks.size()      , 100, 0, 20.0, 100, 0, 100.0 );		
			    plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_jetntrk_v_rho_dR0p05_%i" , sample.c_str(),sel.c_str(),cone,i),Form(";jet rho%i;nconst",i) , rhos[i], constituents.size(), 100, 0, 20.0, 100, 0, 100.0 );		
				
				// slice ntracks in rho regions
				rhosel = rhoregion(rhos[i]);
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_evtntrk" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str()),";ntrk"   , tracks.size()      , 100, 0, 100.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_jetntrk" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str()),";nconst" , constituents.size(), 100, 0, 100.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_jet_pt"  , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str()),";jet pt" , suep_jet.pt() , 100, 0, 500.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_jet_eta" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str()),";jet eta", suep_jet.eta(), 100, -3.0, 3.0 );	

			    // slice rho in ntrack regions
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_evtntrk%s" , sample.c_str(),sel.c_str(),cone,i,evtntracksel.c_str()),Form(";jet rho%i",i) , rhos[i], 200, 0, 20.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_jetntrk%s" , sample.c_str(),sel.c_str(),cone,i,jetntracksel.c_str()),Form(";jet rho%i",i) , rhos[i], 200, 0, 20.0 );	

			    // 
			    // rho 0 region study
			    // 
			    // split into dijet/middle/single jet
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho0%s" , sample.c_str(),sel.c_str(),cone,i, rho0sel.c_str()),Form(";jet rho%i",i) , rhos[i], 100, 0, 20.0 );		
			    
			    plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_rho_dR0p05_%i_rho0%sL" , sample.c_str(),sel.c_str(),cone,i, rho0sel.c_str()),Form(";jet rho%i;ntrk"  ,i) , rhos[i], tracks.size()      , 200, 0, 20.0, 250, 0, 250.0 );		
			    plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_rho_dR0p05_%i_rho0%s" , sample.c_str(),sel.c_str(),cone,i, rho0sel.c_str()),Form(";jet rho%i;ntrk"  ,i) , rhos[i], tracks.size()      , 200, 0, 20.0, 100, 0, 100.0 );		
			    plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_jetntrk_v_rho_dR0p05_%i_rho0%s" , sample.c_str(),sel.c_str(),cone,i, rho0sel.c_str()),Form(";jet rho%i;nconst",i) , rhos[i], constituents.size(), 200, 0, 20.0, 100, 0, 100.0 );		

			    // slice ntracks in rho & rho0 regions
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_rho0%s_evtntrk" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str(),rho0sel.c_str()),";ntrk"   , tracks.size()      , 100, 0, 100.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_rho0%s_jetntrk" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str(),rho0sel.c_str()),";nconst" , constituents.size(), 100, 0, 100.0 );		
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_rho0%s_jet_pt"  , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str(),rho0sel.c_str()),";jet pt" , suep_jet.pt() , 100, 0, 500.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho%s_rho0%s_jet_eta" , sample.c_str(),sel.c_str(),cone,i,rhosel.c_str(),rho0sel.c_str()),";jet eta", suep_jet.eta(), 100, -3.0, 3.0 );				    
			    // slice rho in ntrack regions
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho0%s_evtntrk%s" , sample.c_str(),sel.c_str(),cone,i,rho0sel.c_str(),evtntracksel.c_str()),Form(";jet rho%i",i) , rhos[i], 200, 0, 20.0 );	
			    plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_rho_dR0p05_%i_rho0%s_jetntrk%s" , sample.c_str(),sel.c_str(),cone,i,rho0sel.c_str(),jetntracksel.c_str()),Form(";jet rho%i",i) , rhos[i], 200, 0, 20.0 );	            
			    // 
			    if (ievent < 400 && i==0) makeAKDisplay(suep_jet, jets, Form("%s_%s_jetsAK%i_suep_rho0%s",sample.c_str(),sel.c_str(),cone,rho0sel.c_str()), ievent);
            }				
	
            //
            // Nsubjettiness
            //
            float tau21 = nsubjettinessratio(suep_jet,2,1, R);
            float tau32 = nsubjettinessratio(suep_jet,3,2, R);
            float tau31 = nsubjettinessratio(suep_jet,3,1, R);

           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub1"   , sample.c_str(),sel.c_str(),cone),";jet tau_{1}"  , nsubjettiness(suep_jet,1, R)   , 100, 0, 1.0 );		
           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub2"   , sample.c_str(),sel.c_str(),cone),";jet tau_{2}"  , nsubjettiness(suep_jet,2, R)   , 100, 0, 1.0 );		
           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub3"   , sample.c_str(),sel.c_str(),cone),";jet tau_{3}"  , nsubjettiness(suep_jet,3, R)   , 100, 0, 1.0 );		
           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub32"  , sample.c_str(),sel.c_str(),cone),";jet tau_{32}" , tau32, 100, 0, 1.5 );		
           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub31"  , sample.c_str(),sel.c_str(),cone),";jet tau_{32}" , tau31, 100, 0, 1.5 );		
           	plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub21"  , sample.c_str(),sel.c_str(),cone),";jet tau_{21}" , tau21, 100, 0, 1.5 );		
           	
           	// compare tau21 and ntracks
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_nsub21L", sample.c_str(),sel.c_str(),cone),";jet tau_{21};ntrk"   , tau21, tracks.size()      , 100, 0, 1.0, 250, 0, 250.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_nsub21" , sample.c_str(),sel.c_str(),cone),";jet tau_{21};ntrk"   , tau21, tracks.size()      , 100, 0, 1.0, 100, 0, 100.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_jetntrk_v_nsub21" , sample.c_str(),sel.c_str(),cone),";jet tau_{21};nconst" , tau21, constituents.size(), 100, 0, 1.0, 100, 0, 100.0 );
			
			// slice tau in ntrack regions 
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub21_evtntrk%s" , sample.c_str(),sel.c_str(),cone,evtntracksel.c_str()), ";jet tau_{21}", tau21, 100, 0, 1.0 );	
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub21_jetntrk%s" , sample.c_str(),sel.c_str(),cone,jetntracksel.c_str()), ";jet tau_{21}", tau21, 100, 0, 1.0 );	

            // compare tau and rho0/1
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_tau21_v_rho0_dR0p05"   , sample.c_str(),sel.c_str(),cone),Form(";jet rho0;jet tau_{21}" ) , rhos[0], tau21   , 100, 0, 20.0, 100, 0, 1.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_tau21_v_rho1_dR0p05"   , sample.c_str(),sel.c_str(),cone),Form(";jet rho1;jet tau_{21}" ) , rhos[1], tau21   , 100, 0, 20.0, 100, 0, 1.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_tau32_v_rho0_dR0p05"   , sample.c_str(),sel.c_str(),cone),Form(";jet rho0;jet tau_{32}" ) , rhos[0], tau32   , 100, 0, 20.0, 100, 0, 1.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_tau32_v_rho1_dR0p05"   , sample.c_str(),sel.c_str(),cone),Form(";jet rho1;jet tau_{32}" ) , rhos[1], tau32   , 100, 0, 20.0, 100, 0, 1.0 );		

            // slice tau in ntrack/rho0 regions
            plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub21_evtntrk%s_rho0%s" , sample.c_str(),sel.c_str(),cone,evtntracksel.c_str(), rho0sel.c_str() ), ";jet tau_{21}", tau21, 100, 0, 1.0 );	
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_nsub21_rho0%s" 		   , sample.c_str(),sel.c_str(),cone,rho0sel.c_str()					   ), ";jet tau_{21}", tau21, 100, 0, 1.0 );	

            // compare tau21 and ntracks for different rho0 regions
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_nsub21L_rho0%s", sample.c_str(),sel.c_str(),cone,rho0sel.c_str()),";jet tau_{21};ntrk"   , tau21, tracks.size()      , 100, 0, 1.0, 250, 0, 250.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_v_nsub21_rho0%s" , sample.c_str(),sel.c_str(),cone,rho0sel.c_str()),";jet tau_{21};ntrk"   , tau21, tracks.size()      , 100, 0, 1.0, 100, 0, 100.0 );		
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_jetntrk_v_nsub21_rho0%s" , sample.c_str(),sel.c_str(),cone,rho0sel.c_str()),";jet tau_{21};nconst" , tau21, constituents.size(), 100, 0, 1.0, 100, 0, 100.0 );

			// interparticle distances
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_meanDR"    , sample.c_str(),sel.c_str(),cone ), ";jet mean dR"    , meanDR(suep_jet)   , 100, 0, 4.0 );	
			plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_meanMinDR" , sample.c_str(),sel.c_str(),cone ), ";jet mean min dR", meanMinDR(suep_jet), 100, 0, 1.0 );	
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_meanDR"    , sample.c_str(),sel.c_str(),cone ), ";jet mean dR;ntrk"    , meanDR(suep_jet)   , tracks.size() , 100, 0, 4.0, 250, 0, 250 );	
			plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_evtntrk_meanMinDR" , sample.c_str(),sel.c_str(),cone ), ";jet mean min dR;ntrk", meanMinDR(suep_jet), tracks.size() , 100, 0, 1.0, 250, 0, 250 );	
			//if (cone == 15){
			//
			// COMPARISON WITH TRUTH SUEP "jet"
			//
			//plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_dRtruth", sample.c_str(),sel.c_str(),cone),";dR", truth_suep_jet.p4.DeltaR(suep_p4) , 100, 0, 5.0 );
			//plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_truth_ratio_mass"    , sample.c_str(),sel.c_str(),cone),";reco/truth m "		, suep_jet.m()		 /truth_suep_jet.p4.M() 		, 100, 0, 2 );
			//plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_truth_ratio_pt"      , sample.c_str(),sel.c_str(),cone),";reco/truth pt"		, suep_jet.pt()		 /truth_suep_jet.p4.Pt() 		, 100, 0, 2 );
			//plotter.Plot1D(Form( "%s_%s_jetsAK%i_suep_truth_ratio_nconstit", sample.c_str(),sel.c_str(),cone),";reco/truth nconst"	, (float)constituents.size()/(float)truth_suep_jet.nTruthTracks 	, 100, 0, 2 );
	
			//plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_truth_mass"    , sample.c_str(),sel.c_str(),cone),";truth m;reco m "			, truth_suep_jet.p4.M() 		, suep_jet.m()	, 100, 0, 2000, 100, 0, 2000 );
			//plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_truth_pt"      , sample.c_str(),sel.c_str(),cone),";truth pt;reco pt"			, truth_suep_jet.p4.Pt() 		, suep_jet.pt()	, 100, 0, 2000, 100, 0, 2000 );
			//plotter.Plot2D(Form( "%s_%s_jetsAK%i_suep_truth_nconstit", sample.c_str(),sel.c_str(),cone),";truth nconst;reco nconst"	, truth_suep_jet.nTruthTracks 	, constituents.size()	, 100, 0,  500, 100, 0,  500 );				
			//}



		}



	}


	return fat_jets;
}
