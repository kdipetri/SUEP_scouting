#define kinematics_cxx




//void scalar_plots(std::string sample, std::string sel){
//	  plotter.Plot1D(Form("%s_%s_scalar_pt"  ,sample.c_str(), sel.c_str()),";scalar pT"  , scalar.Pt() , 100,0,1000);
//    plotter.Plot1D(Form("%s_%s_scalar_eta" ,sample.c_str(), sel.c_str()),";scalar eta" , scalar.Eta(), 100,-3.5,3.5);
//    plotter.Plot1D(Form("%s_%s_scalar_phi" ,sample.c_str(), sel.c_str()),";scalar phi" , scalar.Phi(), 100,-3.5,3.5);
//    plotter.Plot1D(Form("%s_%s_scalar_m"   ,sample.c_str(), sel.c_str()),";scalar mass", scalar.M()  , 100,0,1000);
//}
void vertex_plots(std::string sample, std::string sel){
    plotter.Plot1D(Form("%s_%s_nvtx_after"         ,sample.c_str(),sel.c_str()),";n_{PV}", vertices.size()                  , 50,0,50  , evt_wght);
    plotter.Plot1D(Form("%s_%s_nvtx_before"        ,sample.c_str(),sel.c_str()),";n_{PV}", vertices.size()                  , 50,0,50  , 1.0);
    plotter.Plot1D(Form("%s_%s_vertex0_ntracks"    ,sample.c_str(),sel.c_str()),";n_{tracks}", vertices.at(0).tracksSize    , 100,0,100, evt_wght);
    for (auto vertex : vertices){
      plotter.Plot1D(Form("%s_%s_vertex_ntracks"     ,sample.c_str(),sel.c_str()),";n_{tracks};", vertex.tracksSize         , 100,0,100, evt_wght);
    }
    
    //plotter.Plot1D(Form("%s_%s_nchpfs_07" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_07, 100,0,1000);
    //plotter.Plot1D(Form("%s_%s_nchpfs_08" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_08, 100,0,1000);
    //plotter.Plot1D(Form("%s_%s_nchpfs_09" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_09, 100,0,1000);      
    //plotter.Plot1D(Form("%s_%s_nchpfs_2"  ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_2 , 100,0,1000);  

    for (auto track : tracks){
      plotter.Plot1D(Form("%s_%s_track_pt"  ,sample.c_str(),sel.c_str()) ,";track p_{T} [GeV]", track.p4.Pt()   , 100, 0, 10    , evt_wght);
      plotter.Plot1D(Form("%s_%s_track_ptL" ,sample.c_str(),sel.c_str()) ,";track p_{T} [GeV]", track.p4.Pt()   , 100, 0, 100   , evt_wght);
      plotter.Plot1D(Form("%s_%s_track_eta" ,sample.c_str(),sel.c_str()) ,";track eta"        , track.p4.Eta()  , 100, -3.5, 3.5, evt_wght);
      plotter.Plot1D(Form("%s_%s_track_phi" ,sample.c_str(),sel.c_str()) ,";track phi"        , track.p4.Phi()  , 100, -3.5, 3.5, evt_wght);
    }
}
void track_plots(std::string sample, std::string sel){
    plotter.Plot1D(Form("%s_%s_ntracksS"    ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks   , 100,0,100 , evt_wght);
    plotter.Plot1D(Form("%s_%s_ntracks"     ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks   , 100,0,500 , evt_wght);
    plotter.Plot1D(Form("%s_%s_ntracksL"    ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks   , 100,0,1000, evt_wght);
    //plotter.Plot1D(Form("%s_%s_nchpfs_07" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_07, 100,0,1000);
    //plotter.Plot1D(Form("%s_%s_nchpfs_08" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_08, 100,0,1000);
    //plotter.Plot1D(Form("%s_%s_nchpfs_09" ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_09, 100,0,1000);      
    //plotter.Plot1D(Form("%s_%s_nchpfs_2"  ,sample.c_str(),sel.c_str()),";n_{tracks}", ntracks_2 , 100,0,1000);  

    for (auto track : tracks){
      plotter.Plot1D(Form("%s_%s_track_pt"  ,sample.c_str(),sel.c_str()) ,";track p_{T} [GeV]", track.p4.Pt()   , 100, 0, 10 , evt_wght);
      plotter.Plot1D(Form("%s_%s_track_ptL" ,sample.c_str(),sel.c_str()) ,";track p_{T} [GeV]", track.p4.Pt()   , 100, 0, 100, evt_wght);
      plotter.Plot1D(Form("%s_%s_track_eta" ,sample.c_str(),sel.c_str()) ,";track eta"        , track.p4.Eta()  , 100, -3.5, 3.5, evt_wght);
      plotter.Plot1D(Form("%s_%s_track_phi" ,sample.c_str(),sel.c_str()) ,";track phi"        , track.p4.Phi()  , 100, -3.5, 3.5, evt_wght);
    }
}
void jet_plots(std::string sample, std::string sel){
      plotter.Plot1D(Form("%s_%s_ht"        ,sample.c_str(),sel.c_str()),";H_{T} [GeV]"    , ht         , 30,0,1000 , evt_wght);
      plotter.Plot1D(Form("%s_%s_htL"       ,sample.c_str(),sel.c_str()),";H_{T} [GeV]"    , ht         , 30,0,3000 , evt_wght);
      plotter.Plot1D(Form("%s_%s_njets"     ,sample.c_str(),sel.c_str()),";n_{jets}"       , njets      , 11,-0.5,10.5 , evt_wght);
      plotter.Plot1D(Form("%s_%s_njetsL"    ,sample.c_str(),sel.c_str()),";n_{jets}"       , njets      , 21,-0.5,20.5 , evt_wght);
      //plotter.Plot1D(Form("%s_%s_leadjetpt" ,sample.c_str(),sel.c_str()),";jet1 pT [GeV]" , lead_jet_pt, 20,0, 1000);
      for (auto jet : jets){
        plotter.Plot1D(Form("%s_%s_jet_pt"  ,sample.c_str(),sel.c_str()) ,";jet p_{T} [GeV]", jet.p4.Pt()   , 100, 0, 100 , evt_wght);
        plotter.Plot1D(Form("%s_%s_jet_ptL" ,sample.c_str(),sel.c_str()) ,";jet p_{T} [GeV]", jet.p4.Pt()   , 100, 0, 1000, evt_wght);
        plotter.Plot1D(Form("%s_%s_jet_eta" ,sample.c_str(),sel.c_str()) ,";jet eta"        , jet.p4.Eta()  , 100, -3.5, 3.5, evt_wght);
        plotter.Plot1D(Form("%s_%s_jet_phi" ,sample.c_str(),sel.c_str()) ,";jet phi"        , jet.p4.Phi()  , 100, -3.5, 3.5, evt_wght);
      }
}
void combined_plots(std::string sample, std::string sel){
	plotter.Plot2D(Form("%s_%s_ntracks_v_ht"      , sample.c_str(),sel.c_str()), ";H_{T} [GeV];n_{tracks}"     , ht   , ntracks, 30,0,3000   , 50,0,500, evt_wght);
  plotter.Plot2D(Form("%s_%s_ntracks_v_njets"   , sample.c_str(),sel.c_str()), ";n_{jets};n_{tracks}"        , njets, ntracks, 20,-0.5,19.5, 50,0,500, evt_wght);
  plotter.Plot2D(Form("%s_%s_ntracks_v_pvtracks", sample.c_str(),sel.c_str()), ";n_{tracks};PV n_{tracks}"   , ntracks, vertices.at(0).tracksSize, 50,0,100, 50, 0 ,100, evt_wght);
  plotter.Plot2D(Form("%s_%s_ntracks_v_npvs", sample.c_str(),sel.c_str()), ";n_{PV};n_{tracks}"   , vertices.size(), ntracks, 50,0,50, 50, 0 ,200, evt_wght);
  plotter.Plot2D(Form("%s_%s_njets_v_npvs", sample.c_str(),sel.c_str()), ";n_{PV};n_{jets}"   , vertices.size(), njets, 50,0,50, 20,-0.5,19.5, evt_wght);
}
void basic_kinematics(std::string sample, std::string sel){//, std::vector<Jet> jets, std::vector<Track> tracks, std::vector){
  jet_plots(sample,sel);
  vertex_plots(sample,sel);
	if (sel!="all") track_plots(sample,sel);
	//scalar_plots(sample,sel);
	combined_plots(sample,sel);
}
