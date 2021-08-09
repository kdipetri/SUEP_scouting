import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

def clean1D(hist):
    # Clean
    #adjust(hist)
    hist.SetLineWidth(2)
    hist.GetYaxis().SetNdivisions(505)
    hist.GetXaxis().SetNdivisions(505)
    hist.SetDirectory(0)
    hist.Scale(1.0/hist.Integral(0,-1))
    return hist

def get1D(sample,dist):

    # Get hist
    if "QCD" in sample: hist = getQCD(dist)
    else : 
        filename = "output/{}.root".format(sample)
        f = ROOT.TFile.Open(filename)
        histname = "{}_{}".format(sample,dist)
        #print(histname)
        hist = f.Get(histname)
    
    # if hist exists
    if hist : 
        clean1D(hist)
        return hist
    else : return 0

def qcd_xs(sample):

    if "QCD_HT200to300"   in sample : return 1559000
    if "QCD_HT300to500"   in sample : return 311900
    if "QCD_HT500to700"   in sample : return 29070
    if "QCD_HT700to1000"  in sample : return 5962
    if "QCD_HT1000to1500" in sample : return 1207 
    if "QCD_HT1500to2000" in sample : return 119.9 
    if "QCD_HT2000toInf"  in sample : return 25.24 

def getQCD(dist):

    # Get hist
    samples =[]
    samples.append("QCD_HT200to300")# do slicing later
    samples.append("QCD_HT300to500")# do slicing later
    samples.append("QCD_HT500to700")# do slicing later
    samples.append("QCD_HT700to1000")# do slicing later
    samples.append("QCD_HT1000to1500")# do slicing later
    samples.append("QCD_HT1500to2000")# do slicing later
    samples.append("QCD_HT2000toInf")# do slicing later
    
    hists = []
    for sample in samples: 
        f = ROOT.TFile.Open("output/{}.root".format(sample))
        if not f: continue
        h = f.Get("{}_{}".format(sample,dist))
        if not h:
            print("Missing hist: {}_{}".format(sample,dist)) 
            continue
        #scale to xs * lumi * totalweights
        h.Scale(qcd_xs(sample))
        h.SetDirectory(0)
        hists.append(h)

    hist_final = hists[0].Clone("QCD_"+dist)
    for i,hist in enumerate(hists):
        if i>0: hist_final.Add(hist)

    clean1D(hist_final)

    return hist_final



h1 = get1D("QCD","scouting_nvtx")
h2 = get1D("data18","scouting_nvtx")
ratio = h2.Clone("ratio_scouting_nvtx")
ratio.Divide(h1)

c = ROOT.TCanvas("ht","",800,800)

#dy = 0.05*2.0
#leg = ROOT.TLegend(0.18,0.86-dy,0.86,0.86)
#leg.SetTextSize(0.035)
#leg.SetBorderSize(0)


ratio.Draw("histe")

    #hist.SetLineColor(colors[i])
    #if "QCD" in labels[i]: 
    #    hist.SetLineColor(ROOT.kGray+1) 
    #    hist.SetFillStyle(1001)
    #    hist.SetFillColorAlpha(ROOT.kGray+1,0.8)
    #if "data" in labels[i]:
    #    hist.SetLineColor(ROOT.kBlack)
    #    hist.SetMarkerColor(ROOT.kBlack)
    #    hist.SetMarkerStyle(20)
    #    hist.SetMarkerSize(1.0)

#c.SetLogy(1)
#hists[0].GetYaxis().SetRangeUser(0.0001,ymax*100)
#c.Print("plots/{}_log.png".format(filename))
ratio.GetYaxis().SetRangeUser(0,1.5)
c.SetLogy(0)
c.Print("plots/pu_reweighting.png")

# save output
fout = ROOT.TFile("plots/safe/pu_reweighting.root","RECREATE")
fout.cd()
ratio.Write()
fout.Close()


