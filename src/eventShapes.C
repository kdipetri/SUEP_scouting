#define eventShapes_cxx

// 
void makeEventDisplay(std::vector<Track> tracks, std::string sample, std::string sel, Long64_t ievent, double sphericity){

  std::vector<TGraph*> graphs;
  
  int i=0;
  for (Track track: tracks){
    
    // draw all constituents 
    TGraph *graph = new TGraph();
    graph->SetTitle(Form("g_%s_%s_%lli_%i",sample.c_str(), sel.c_str(),ievent, i));
    graph->SetPoint(0,track.p4.Eta(),track.p4.Phi());
      
    graph->SetMarkerSize( std::log( track.p4.Pt() + 1) ); // test
    graph->SetMarkerColor(kBlue+1);
    graph->SetMarkerStyle(21);
    graphs.push_back(graph);

    i++;
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";#eta; #phi");
  for (auto gr: graphs){
    mg->Add(gr);    
  }   

  c1->cd();
  mg->Draw("ap");
  mg->GetYaxis()->SetLimits(-3.5,3.5);          
  mg->GetXaxis()->SetLimits(-3.5,3.5);// along  
  c1->Update();

  TLatex l;
  l.SetTextSize(0.03);
  l.DrawLatex(0.1,0.85,Form("sphericity %.2f",sphericity) );

  c1->Print(Form("plots/boostingABCD/eventDisplay/%s_%s_%lli.png", sample.c_str(), sel.c_str(), ievent));
}

/// the return value is 1 for spherical events and 0 for events linear in r-phi. This function
/// needs the number of steps to determine how fine the granularity of the algorithm in phi
/// should be
double isotropy(std::vector<Track> tracks, const unsigned int& numberOfSteps=1000) 
{
  //std::cout << "isotropy " << tracks.size() << std::endl;
  const double deltaPhi = 2 * TMath::Pi() / numberOfSteps;
  double phi = 0, eIn = -1., eOut = -1.;
  for (unsigned int i = 0; i < numberOfSteps; ++i) {
    phi += deltaPhi;
    double sum = 0;
    double cosphi = TMath::Cos(phi);
    double sinphi = TMath::Sin(phi);
    for (const auto& track : tracks) {
      // sum over inner product of unit vectors and momenta
      sum += TMath::Abs(cosphi * track.p4.X() + sinphi * track.p4.Y());
    }
    if (eOut < 0. || sum < eOut)
      eOut = sum;
    if (eIn < 0. || sum > eIn)
      eIn = sum;
  }
  return 1.0 - (eIn - eOut) / eIn;
}

/// the return value is 1 for spherical and 0 linear events in r-phi. This function needs the
/// number of steps to determine how fine the granularity of the algorithm in phi should be
double circularity(std::vector<Track> tracks, const unsigned int& numberOfSteps=1000) 
{
  //std::cout << "circularity " << tracks.size() << std::endl;
  const double deltaPhi = 2 * TMath::Pi() / numberOfSteps;
  double circularity = -1, phi = 0, area = 0;
  for (const auto& track : tracks) {
    area += TMath::Sqrt(track.p4.X() * track.p4.X() + track.p4.Y() * track.p4.Y());
  }
  for (unsigned int i = 0; i < numberOfSteps; ++i) {
    phi += deltaPhi;
    double sum = 0, tmp = 0.;
    double cosphi = TMath::Cos(phi);
    double sinphi = TMath::Sin(phi);
    for (const auto& track : tracks) {
      sum += TMath::Abs(cosphi * track.p4.X() + sinphi * track.p4.Y());
    }
    tmp = TMath::Pi() / 2 * sum / area;
    if (circularity < 0 || tmp < circularity) {
      circularity = tmp;
    }
  }
  return circularity;
}

/// helper function to fill the 3 dimensional momentum tensor from the inputVectors where needed
/// also fill the 3 dimensional vectors of eigen-values and eigen-vectors;
/// the largest (smallest) eigen-value is stored at index position 0 (2)
std::vector<double> compTensorsAndVectors(std::vector<Track> tracks) {
  //if (tensors_computed_)
  //  return;
  //
  //if (inputVectors_.size() < 2) {
  //  tensors_computed_ = true;
  //  return;
  //}
  // 
  //std::cout << "tensors " << tracks.size() << std::endl;
  double r_=2.0;
  std::vector<double> eigenValues_ = std::vector<double>(3, 0);
  std::vector<double> eigenValuesNoNorm_ = std::vector<double>(3, 0);

  TMatrixD eigenVectors_;
  TVectorD eigenValuesTmp_, eigenValuesNoNormTmp_;

  TMatrixDSym momentumTensor(3);
  momentumTensor.Zero();

  if (tracks.size() == 0) return eigenValues_;

  // fill momentumTensor from inputVectors
  double norm = 0.;
  TVector3 vec; 
  for (const auto& track : tracks) {

  	// get 3 vector
  	vec = track.p4.Vect();
    double p2 = (vec).Dot(vec);
    double pR = (r_ == 2.) ? p2 : TMath::Power(p2, 0.5 * r_);
    //std::cout << " pR " << p2 <<  std::endl;
    norm += pR;
    double pRminus2 = (r_ == 2.) ? 1. : TMath::Power(p2, 0.5 * r_ - 1.);
    momentumTensor(0, 0) += pRminus2 * vec.X() * vec.X();
    momentumTensor(0, 1) += pRminus2 * vec.X() * vec.Y();
    momentumTensor(0, 2) += pRminus2 * vec.X() * vec.Z();
    momentumTensor(1, 0) += pRminus2 * vec.Y() * vec.X();
    momentumTensor(1, 1) += pRminus2 * vec.Y() * vec.Y();
    momentumTensor(1, 2) += pRminus2 * vec.Y() * vec.Z();
    momentumTensor(2, 0) += pRminus2 * vec.Z() * vec.X();
    momentumTensor(2, 1) += pRminus2 * vec.Z() * vec.Y();
    momentumTensor(2, 2) += pRminus2 * vec.Z() * vec.Z();
  }

  if (momentumTensor.IsSymmetric() && (momentumTensor.NonZeros() != 0)) {
    momentumTensor.EigenVectors(eigenValuesNoNormTmp_);
  }
  eigenValuesNoNorm_[0] = eigenValuesNoNormTmp_(0);
  eigenValuesNoNorm_[1] = eigenValuesNoNormTmp_(1);
  eigenValuesNoNorm_[2] = eigenValuesNoNormTmp_(2);

  // momentumTensor normalized to determinant 1
  momentumTensor *= (1. / norm);
  //std::cout << "norm " << norm << std::endl;
  //std::cout << " mt(0, 0)" <<  momentumTensor(0, 0) << std::endl;
  //std::cout << " mt(0, 1)" <<  momentumTensor(0, 1) << std::endl;
  //std::cout << " mt(0, 2)" <<  momentumTensor(0, 2) << std::endl;
  //std::cout << " mt(1, 0)" <<  momentumTensor(1, 0) << std::endl;
  //std::cout << " mt(1, 1)" <<  momentumTensor(1, 1) << std::endl;
  //std::cout << " mt(1, 2)" <<  momentumTensor(1, 2) << std::endl;
  //std::cout << " mt(2, 0)" <<  momentumTensor(2, 0) << std::endl;
  //std::cout << " mt(2, 1)" <<  momentumTensor(2, 1) << std::endl;
  //std::cout << " mt(2, 2)" <<  momentumTensor(2, 2) << std::endl;

  // now get eigens
  if (momentumTensor.IsSymmetric() && (momentumTensor.NonZeros() != 0)) {
    momentumTensor.EigenVectors(eigenValuesTmp_);
  }
  eigenValues_[0] = eigenValuesTmp_(0);
  eigenValues_[1] = eigenValuesTmp_(1);
  eigenValues_[2] = eigenValuesTmp_(2);
  //std::cout << eigenValues_[0] << " " << eigenValues_[1] << " " << eigenValues_[2] << std::endl;
  return eigenValues_;
}

/// 1.5*(q1+q2) where q0>=q1>=q2>=0 are the eigenvalues of the momentum tensor sum{p_j[a]*p_j[b]}/sum{p_j**2}
/// normalized to 1. Return values are 1 for spherical, 3/4 for plane and 0 for linear events
double sphericity(std::vector<double> eigenValues_) {
  return 1.5 * (eigenValues_[1] + eigenValues_[2]);
}

/// 1.5*q2 where q0>=q1>=q2>=0 are the eigenvalues of the momentum tensor sum{p_j[a]*p_j[b]}/sum{p_j**2}
/// normalized to 1. Return values are 0.5 for spherical and 0 for plane and linear events
double aplanarity(std::vector<double> eigenValues_) {
  return 1.5 * eigenValues_[2];
}

/// 3.*(q0*q1+q0*q2+q1*q2) where q0>=q1>=q2>=0 are the eigenvalues of the momentum tensor sum{p_j[a]*p_j[b]}/sum{p_j**2}
/// normalized to 1. Return value is between 0 and 1
/// and measures the 3-jet structure of the event (C vanishes for a "perfect" 2-jet event)
double C(std::vector<double> eigenValues_) {
  return 3. *
         (eigenValues_[0] * eigenValues_[1] + eigenValues_[0] * eigenValues_[2] + eigenValues_[1] * eigenValues_[2]);
}

/// 27.*(q0*q1*q2) where q0>=q1>=q2>=0 are the eigenvalues of the momemtum tensor sum{p_j[a]*p_j[b]}/sum{p_j**2}
/// normalized to 1. Return value is between 0 and 1
/// and measures the 4-jet structure of the event (D vanishes for a planar event)
double D(std::vector<double> eigenValues_) {
  return 27. * eigenValues_[0] * eigenValues_[1] * eigenValues_[2];
}

void getFWmoment(){
	double fwmom_maxl_ = 10.0;
  	std::vector<double> fwmom_ = std::vector<double>(fwmom_maxl_, 0.);
}

void plotEventShapes(std::string sample_name, std::string sel, std::vector<Track> tracks, Long64_t ievent)
{
  if ( tracks.size() == 0 ) return;

	double iso = isotropy(tracks);
	double cir = circularity(tracks);
	plotter.Plot1D(Form("%s_%s_evtshape_isotropy"    ,sample_name.c_str(), sel.c_str()),";isotropy"   , iso, 100,0,1, evt_wght);
	plotter.Plot1D(Form("%s_%s_evtshape_circularity" ,sample_name.c_str(), sel.c_str()),";circularity", cir, 100,0,1, evt_wght);

	std::vector<double> eigenValues = compTensorsAndVectors(tracks);
	double sphere = sphericity(eigenValues);
	double aplan  = aplanarity(eigenValues);
	double evt_c  = C(eigenValues);
	double evt_d  = D(eigenValues);

	plotter.Plot1D(Form("%s_%s_evtshape_sphericity"  ,sample_name.c_str(), sel.c_str()),";sphericity" , sphere, 100,0,1, evt_wght);
	plotter.Plot1D(Form("%s_%s_evtshape_aplanarity"  ,sample_name.c_str(), sel.c_str()),";aplanarity"	, aplan , 100,0,1, evt_wght);
	plotter.Plot1D(Form("%s_%s_evtshape_c"    		   ,sample_name.c_str(), sel.c_str()),";c"   		    , evt_c , 100,0,1, evt_wght);
	plotter.Plot1D(Form("%s_%s_evtshape_d" 			     ,sample_name.c_str(), sel.c_str()),";d"			    , evt_d , 100,0,1, evt_wght);

  // useful for bkg estimation
  plotter.Plot2D(Form("%s_%s_evtshape_ht_v_circularity"     ,sample_name.c_str(), sel.c_str()),";circularity;ht [GeV]", cir,ht      , 50,0,1, 50,0,3000    , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_njets_v_circularity"  ,sample_name.c_str(), sel.c_str()),";circularity;njets"   , cir,njets   , 50,0,1, 20,-0.5, 19.5, evt_wght);
  
  plotter.Plot2D(Form("%s_%s_evtshape_ht_v_sphericity"     ,sample_name.c_str(), sel.c_str()),";sphericity;ht [GeV]", sphere,ht      , 50,0,1, 50,0,3000    , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_njets_v_sphericity"  ,sample_name.c_str(), sel.c_str()),";sphericity;njets"   , sphere,njets   , 50,0,1, 20,-0.5, 19.5, evt_wght);

  // useful for bkg estimation
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_isotropy"   ,sample_name.c_str(), sel.c_str()),";isotropy;ntracks"    , iso   ,ntracks , 50,0,1, 50,0,200     , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_circularity",sample_name.c_str(), sel.c_str()),";circularity;ntracks" , cir   ,ntracks , 50,0,1, 50,0,200     , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_sphericity" ,sample_name.c_str(), sel.c_str()),";sphericity;ntracks"  , sphere,ntracks , 50,0,1, 50,0,200     , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_aplanarity" ,sample_name.c_str(), sel.c_str()),";aplanarity;ntracks" , aplan ,ntracks , 50,0,1, 50,0,200     , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_c"          ,sample_name.c_str(), sel.c_str()),";c;ntracks"          , evt_c ,ntracks , 50,0,1, 50,0,200     , evt_wght);
  plotter.Plot2D(Form("%s_%s_evtshape_ntracks_v_d"          ,sample_name.c_str(), sel.c_str()),";d;ntracks"          , evt_d ,ntracks , 50,0,1, 50,0,200     , evt_wght);

  if (ievent < 100) makeEventDisplay(tracks, sample_name,  sel, ievent, sphere);
}




