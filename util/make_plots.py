import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?

def adjust(hist):
    name = hist.GetName()
    if "nvtx" in name: hist.GetXaxis().SetTitle("n_{PV}")
    if "njetsL" in name: hist.GetXaxis().SetRangeUser(-0.5,12.5)
    if "track_pt" in name: hist.GetXaxis().SetRangeUser(0,5)
    #if "ntracks" in name: hist.Rebin()
    return

def clean1D(hist):
    # Clean
    adjust(hist)
    hist.SetLineWidth(2)
    hist.GetYaxis().SetNdivisions(505)
    hist.GetXaxis().SetNdivisions(505)
    hist.SetDirectory(0)
    hist.Scale(1.0/hist.Integral(0,-1))
    return hist


def label(sample):
    #return "(m_{S},m_{#phi},T)=(%i,%i,%i), %s"%(mMed,mDark,temp,decay_label(decay))
    if "QCD" in sample: return "QCD"
    elif "data18" in sample: return "data18"
    elif "200" in sample: return  "m_{S}=200, m_{#phi}=2, T=2"
    elif "300" in sample: return  "m_{S}=300, m_{#phi}=2, T=2"
    elif "400" in sample: return  "m_{S}=400, m_{#phi}=2, T=2"
    else : return ""

def doROC(histname):
    if "nchpfs" in histname: return 1
    if "evtshape" in histname: return 1
    if "suep" in histname and "dRtruth" not in histname: return 1
    else: return 0 

    
def compare1D(hists,labels,filename):
    c = ROOT.TCanvas(filename,"",800,800)

    dy = 0.05*len(hists)
    leg = ROOT.TLegend(0.5,0.86-dy,0.86,0.86)
    leg.SetTextSize(0.035)
    leg.SetBorderSize(0)

    ymax = 0
    for i,hist in enumerate(hists):
        #print(i,hist) 
        hist.SetLineColor(colors[i+1])
        if "QCD" in labels[i]: 
            hist.SetLineColor(ROOT.kGray+1) 
            hist.SetFillStyle(1001)
            hist.SetFillColorAlpha(ROOT.kGray+1,0.8)
        if "mMed-300" in labels[i]: 
            hist.SetLineColor(ROOT.kBlue+1) 
        if "mMed-200" in labels[i]: 
            hist.SetLineColor(ROOT.kRed+1) 
        if "data" in labels[i]:
            hist.SetLineColor(ROOT.kBlack)
            hist.SetMarkerColor(ROOT.kBlack)
            hist.SetMarkerStyle(20)
            hist.SetMarkerSize(1.0)
        if i==0: 
            if "data" in labels[i] : hist.Draw("e")
            else : hist.Draw("hist")
        else :
            if "data" in labels[i] : hist.Draw("e same")
            else : hist.Draw("hist same")

        if hist.GetMaximum() > ymax: ymax=hist.GetMaximum()

        if "data" in labels[i]: leg.AddEntry(hist,labels[i],"pl")
        elif "QCD" in labels[i] :  leg.AddEntry(hist,labels[i],"f")
        else : leg.AddEntry(hist,labels[i],"l")
        

    leg.Draw()
    
    c.SetLogy(1)
    hists[0].GetYaxis().SetRangeUser(0.0001,ymax*100)
    if "ntracks" in filename: hists[0].GetYaxis().SetRangeUser(0.00001,ymax*100)
    c.Print("plots/{}_log.png".format(filename))
    hists[0].GetYaxis().SetRangeUser(0,ymax*1.8)
    c.SetLogy(0)
    c.Print("plots/{}_lin.png".format(filename))

def ratio1D(histMC,histData,filename):

    c = ROOT.TCanvas(filename,"",600,800)

    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "", 0, 0.4, 1, 1.0)
    pad1.SetBottomMargin(0.05)  # joins upper and lower plot
    pad1.SetLeftMargin(0.2)
    #pad1.SetGridx()
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.4)
    pad2.SetTopMargin(0.05)  # joins upper and lower plot
    pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.2)
    pad2.SetGridy()
    pad2.Draw()

    # top pad
    pad1.cd()

    dy = 0.2
    leg = ROOT.TLegend(0.25,0.86-dy,0.86,0.86)
    leg.SetTextSize(0.05)
    leg.SetBorderSize(0)

    ymax = max(histMC.GetMaximum(),histData.GetMaximum())

    histMC.SetLineColor(ROOT.kGray+1) 
    histMC.SetFillStyle(1001)
    histMC.SetFillColorAlpha(ROOT.kGray+1,0.8)
    histMC.Draw("hist")

    histData.SetLineColor(ROOT.kBlack)
    histData.SetMarkerColor(ROOT.kBlack)
    histData.SetMarkerStyle(20)
    histData.SetMarkerSize(1.0)
    histData.Draw("e same")

    # asthetics
    histMC.GetYaxis().SetLabelSize(0.07)
    histMC.GetXaxis().SetLabelSize(0)

    leg.AddEntry(histData,"data","pl")
    leg.AddEntry(histMC,"MC","f")
    leg.Draw()
    
    # bottom pad
    pad2.cd()
    histRatio = histData.Clone("ratio"+histData.GetName().strip("data"))
    histRatio.Divide(histMC)
    histRatio.Clone()
    histRatio.Draw("e")

    # asthetics
    histRatio.GetYaxis().SetRangeUser(0,2)
    histRatio.GetYaxis().SetNdivisions(505)
    histRatio.GetYaxis().SetTitleOffset(0.6)
    histRatio.GetXaxis().SetTitleOffset(0.85)
    histRatio.GetYaxis().SetTitleSize(0.12)
    histRatio.GetXaxis().SetTitleSize(0.12)
    histRatio.GetYaxis().SetLabelSize(0.1)
    histRatio.GetXaxis().SetLabelSize(0.1)
    histRatio.GetYaxis().SetTitle("Data/MC")

    pad1.cd()

    c.Print("plots/{}_lin.png".format(filename))

    ROOT.gPad.SetLogy(1)
    histMC.GetYaxis().SetRangeUser(0.0001,ymax*100)
    c.Print("plots/{}_log.png".format(filename))    

    ROOT.gPad.SetLogy(0)

def compareDataMC(dist):
    samples = []
    samples.append("QCD")
    samples.append("data18")
    samples.append("mMed-200_mDark-2_temp-2_decay-darkPho")
    samples.append("mMed-300_mDark-2_temp-2_decay-darkPho")

    hists=[]
    labels=[]
    for sample in samples:

        hist = get1D(sample,dist)
        if hist : 
            clean1D(hist)
            if dist == "scouting_ntracks": 
                hist.GetXaxis().SetRangeUser(0,200)
                if sample == "data18" :continue

            hists.append(hist)
            labels.append(label(sample))

    compare1D(hists,labels,"compareDataMC/compare1D_{}".format(dist))

    ratio1D(hists[0],hists[1],"compareDataMC/ratio1D_{}".format(dist))
    #if histname=="h_pf_ntracks": 
    #if doROC(histname)  : makeROC(hists,labels,"roc_curve/temp{}_mDark{}_decay_{}_{}".format(temp,mDark,decay,histname))



sels = []
#sels.append("all")
sels.append("scouting")

dists=[]
dists.append("njets") 
dists.append("njets") 
dists.append("njetsL") 
dists.append("ht") 
dists.append("htL") 
dists.append("ntracksS")
dists.append("ntracks")
dists.append("ntracksL")

dists.append("jet_ptL") 
dists.append("jet_pt")
dists.append("jet_eta")
dists.append("jet_phi")

dists.append("track_ptL") 
dists.append("track_pt")
dists.append("track_eta")
dists.append("track_phi")

#dists.append("evtshape_aplanarity")   
#dists.append("evtshape_c")    
#dists.append("evtshape_circularity")  
#dists.append("evtshape_d")    
#dists.append("evtshape_isotropy") 
#dists.append("evtshape_sphericity")



dists.append("vertex_ntracks")
dists.append("vertex0_ntracks")
dists.append("nvtx_before")
dists.append("nvtx_after")

dists.append("jetsAK15_suep_eta")
dists.append("jetsAK15_suep_phi")
dists.append("jetsAK15_suep_m")
dists.append("jetsAK15_suep_moverpt")
dists.append("jetsAK15_suep_pt")
dists.append("jetsAK15_suep_width")
dists.append("jetsAK15_suep_nconstit")

for i in range(15):
    dists.append("jetsAK15_suep_rho_dR0p05_%i"%i)
#
#dists.append("trig_suepjet_min_dR") 
#dists.append("trig_suepjet_pt") 
#dists.append("trig_suepjet_scalar_dPhi")    
#dists.append("trig_suepjet_scalar_dR")
#dists.append("trig_suepjet_isolation")
#dists.append("trig_suepjet_multiplicity")



for sel in sels: 
    for dist in dists:
        name=sel+"_"+dist
        if sel=="all" and "track" in dist: continue
        compareDataMC(name)
        #compareMass(2,2,"darkPho",name)
        #compareDecay(750,2,2,name)

compareDataMC("offline_htL")

