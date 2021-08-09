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
    c.Print("plots/tauStudy/{}_log.png".format(filename))
    c.SetLogz(0)
    c.Print("plots/tauStudy/{}_lin.png".format(filename))

    return
    
def compare1D(hists,labels,filename,opt=""):
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
            elif opt=="err" : hist.Draw("hist e")
            else : hist.Draw("hist")
        else :
            if "data" in labels[i] : hist.Draw("e same")
            elif opt=="err" : hist.Draw("hist e same")
            else : hist.Draw("hist same")

        if hist.GetMaximum() > ymax: ymax=hist.GetMaximum()

        if "data" in labels[i]: leg.AddEntry(hist,labels[i],"pl")
        elif opt=="err" : leg.AddEntry(hist,labels[i],"pl")
        elif "QCD" in labels[i] :  leg.AddEntry(hist,labels[i],"f")
        else : leg.AddEntry(hist,labels[i],"l")
        

    leg.Draw()
    
    c.SetLogy(1)
    hists[0].GetYaxis().SetRangeUser(0.0001,ymax*100)
    c.Print("plots/{}_log.png".format(filename))
    hists[0].GetYaxis().SetRangeUser(0,ymax*1.8)
    c.SetLogy(0)
    c.Print("plots/{}_lin.png".format(filename))

def compareDataMC(dist):
    samples = []
    #samples.append("QCD")
    samples.append("data18")
    samples.append("mMed-200_mDark-2_temp-2_decay-darkPho")
    samples.append("mMed-300_mDark-2_temp-2_decay-darkPho")
    samples.append("mMed-400_mDark-2_temp-2_decay-darkPho")

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

    compare1D(hists,labels,"tauStudy/compare1D_{}".format(dist))

#
# Make nsubjettiness signal and background shape comparison plots
#
taus=[]
taus.append("nsub1")
taus.append("nsub2")
taus.append("nsub3")
taus.append("nsub21")
taus.append("nsub32")
for tau in taus:
    compareDataMC("scouting_jetsAK15_suep_"+tau)

# 
# Make 2D profiles of tau21 v ntracks
# 

profile2D("data18","scouting_jetsAK15_suep_evtntrk_v_nsub21")
profile2D("data18","scouting_jetsAK15_suep_jetntrk_v_nsub21")
    
profile2D("mMed-200_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_jetntrk_v_nsub21")
profile2D("mMed-200_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_evtntrk_v_nsub21")

#
## now compare tau in ntrack slices for data

def compareTauVNtrack(dist="jetsAK15_suep_nsub21",kind="evtntrk", sample="data18",sel="scouting"):

    ntrks = []

    if "jet" in kind: 
        #ntrks.append(kind+"0to5");
        ntrks.append(kind+"5to10");
        ntrks.append(kind+"10to20");
        ntrks.append(kind+"20to30");
        ntrks.append(kind+"30to50");
        ntrks.append(kind+"50toinf");
    else :     
        #ntrks.append("ntrks0to10");
        ntrks.append(kind+"10to30");
        ntrks.append(kind+"30to50");
        ntrks.append(kind+"50to70");
        ntrks.append(kind+"70toinf");

    hists=[]
    labels=[]
    for ntrk in ntrks:
        name = sel+"_"+dist+"_"+ntrk
        print(name)
        hist = get1D(sample,name)
        print(hist)
        if hist : 
            hist.Rebin(10)
            clean1D(hist)
            hists.append(hist)
            labels.append(ntrk)

    compare1D(hists,labels,"tauStudy/compare_correlation_{}_{}".format(dist,kind),opt="err")


dists = compareTauVNtrack()



