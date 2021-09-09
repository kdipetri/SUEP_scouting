import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?

# 
# Make 2D profiles of event shapes v ntracks
# 

def clean(hist):
    name = hist.GetName()
    if "circularity" in name : hist.GetXaxis().SetTitle("circularity")
    elif "isotropy" in name : hist.GetXaxis().SetTitle("isotropy")
    elif "_c" in name : hist.GetXaxis().SetTitle("C")
    elif "_d" in name : hist.GetXaxis().SetTitle("D")
    return

def profile2D(sample,dist):
    filename="{}_{}".format(sample,dist)
    c = ROOT.TCanvas(filename,"",800,700)
    c.SetLeftMargin(0.18)
    c.SetRightMargin(0.2)

    hist = get1D(sample,dist)
    clean(hist)
    hist.GetXaxis().SetRangeUser(0,30)
    hist.GetZaxis().SetTitle("Events [AU]")
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetNdivisions(505)
    hist.GetZaxis().SetTitleOffset(1.1)
    hist.Draw("COLZ")
    profile = hist.ProfileX("profx_"+filename)
    profile.Draw("same")
    
    ROOT.gPad.Modified()
    
    c.SetLogz(1)
    c.Print("plots/boostingABCD/{}_log.png".format(filename))
    c.SetLogz(0)
    c.Print("plots/boostingABCD/{}_lin.png".format(filename))

    return

profile2D("data18","scouting_evtshape_ntracks_v_circularity")
profile2D("data18","scouting_evtshape_ntracks_v_sphericity")

profile2D("data18","scouting_boosted_evtshape_ntracks_v_circularity")
profile2D("data18","scouting_boosted_evtshape_ntracks_v_sphericity")
profile2D("data18","scouting_boosted_evtshape_ntracks_v_isotropy")
profile2D("data18","scouting_boosted_evtshape_ntracks_v_c")
profile2D("data18","scouting_boosted_evtshape_ntracks_v_d")
profile2D("data18","scouting_boosted_evtshape_ntracks_v_aplanarity")

profile2D("QCD","scouting_boosted_evtshape_ntracks_v_circularity")
profile2D("QCD","scouting_boosted_evtshape_ntracks_v_sphericity")
profile2D("QCD","scouting_boosted_evtshape_ntracks_v_isotropy")
profile2D("QCD","scouting_boosted_evtshape_ntracks_v_c")
profile2D("QCD","scouting_boosted_evtshape_ntracks_v_d")
profile2D("QCD","scouting_boosted_evtshape_ntracks_v_aplanarity")

profile2D("QCD","scouting_evtshape_ntracks_v_circularity")
profile2D("QCD","scouting_evtshape_ntracks_v_sphericity")
profile2D("QCD","scouting_evtshape_ntracks_v_isotropy")
profile2D("QCD","scouting_evtshape_ntracks_v_c")
profile2D("QCD","scouting_evtshape_ntracks_v_d")
profile2D("QCD","scouting_evtshape_ntracks_v_aplanarity")

profile2D("QCD","scouting_jetsAK15_suep_evtntrk_meanDR")
profile2D("QCD","scouting_jetsAK15_suep_evtntrk_meanMinDR")

profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_circularity")
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_sphericity")
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_isotropy")
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_c")
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_d")
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_boosted_evtshape_ntracks_v_aplanarity")

