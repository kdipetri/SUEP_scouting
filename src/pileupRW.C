#define pileupRW_cxx


TFile *file_prw   = TFile::Open("plots/safe/pu_reweighting.root");
TH1F *h_prw        = (TH1F*)file_prw->Get("ratio_scouting_nvtx");

float get_PU(int nvtx){

	float pu_weight  = h_prw->GetBinContent( h_prw->FindBin( nvtx ) );

	return pu_weight;
}