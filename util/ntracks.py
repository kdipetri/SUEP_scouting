import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?

def profile2D(sample,dist):
    filename="{}_{}".format(sample,dist)
    c = ROOT.TCanvas(filename,"",800,700)
    c.SetLeftMargin(0.18)
    c.SetRightMargin(0.2)

    hist = get1D(sample,dist)
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
    c.Print("plots/ntracks/{}_log.png".format(filename))
    c.SetLogz(0)
    c.Print("plots/ntracks/{}_lin.png".format(filename))

    return

dist = "scouting_ntracks_v_npvs"
profile2D("data18",dist)
profile2D("QCD",dist)
profile2D("mMed-200_mDark-2_temp-2_decay-darkPho",dist)
profile2D("mMed-300_mDark-2_temp-2_decay-darkPho",dist)
profile2D("mMed-400_mDark-2_temp-2_decay-darkPho",dist)
