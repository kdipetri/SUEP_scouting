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
    c.Print("plots/rhoStudy/{}_log.png".format(filename))
    c.SetLogz(0)
    c.Print("plots/rhoStudy/{}_lin.png".format(filename))

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

    compare1D(hists,labels,"rhoStudy/compare1D_{}".format(dist))


# 
# Make rho signal & background shape comparison plots
#
r0bins = ["rho0low", "rho0high"]
dists=[]
for i in range(15):
    dists.append("scouting_jetsAK15_suep_rho_dR0p05_%i"%i)
#     
# 
# for dist in dists: 
#     compareDataMC(dist)
#     for r0bin in r0bins: 
#         compareDataMC(dist+"_"+r0bin)

#
# Make rho profiles for data
# 

def rhoProfiles(dists,sample="data18",dR=0.05):
    
    # initialize arrays
    x, y = array('d'), array('d')
    xerr, yerr = array('d'), array('d')

    # fill arrays
    r = dR/2.0
    for dist in dists:
        hist = get1D(sample,dist)
        if hist : 
            step = int(dist.split("_")[-1])
            x.append( r )
            xerr.append( dR/2.0 )
            y.append( hist.GetMean() ) 
            yerr.append( hist.GetRMS() )
        r += dR 

    # make graph
    g_rho_v_r = ROOT.TGraphErrors(len(x),x,y,xerr,yerr)
    g_rho_v_r.SetName("{}_mean_rho_v_r".format(sample)) 
    g_rho_v_r.SetTitle("")
    g_rho_v_r.GetXaxis().SetTitle("r")
    g_rho_v_r.GetYaxis().SetTitle("#rho(r)")


    # draw
    c = ROOT.TCanvas(sample+"mean_rho_v_r","",600,600)
    g_rho_v_r.Draw("ACP")
    c.Print("plots/rhoStudy/{}_mean_rho_v_r.png".format(sample))

    return g_rho_v_r

samples = []
samples.append("data18")
samples.append("QCD")
samples.append("mMed-200_mDark-2_temp-2_decay-darkPho")
samples.append("mMed-300_mDark-2_temp-2_decay-darkPho")
samples.append("mMed-400_mDark-2_temp-2_decay-darkPho")

dy = 0.05*len(samples)
leg = ROOT.TLegend(0.5,0.86-dy,0.86,0.86)
leg.SetTextSize(0.035)
leg.SetBorderSize(0)

multi = ROOT.TMultiGraph()
for i,sample in enumerate(samples):

    graph = rhoProfiles(dists,sample)
    graph.SetLineColor(ROOT.kBlack if i==0 else colors[i+1])

    leg.AddEntry(graph,label(sample))
    multi.Add(graph) 

c = ROOT.TCanvas("compare_mean_rho_v_r","",600,600)

multi.Draw("ACP")
leg.Draw()
c.Print("plots/rhoStudy/compare_mean_rho_v_r.png")


#
# now compare Rho V ntrack for data
# 
def compareRhoVNtrack(dist="jetsAK15_suep_rho_dR0p05_1",kind="evtntrk", sample="data18",sel="scouting"):

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

    compare1D(hists,labels,"rhoStudy/compare_correlation_{}_{}".format(dist,kind),opt="err")

def compareNtrackVRho(dist="jetsAK15_suep_rho_dR0p05_1",kind="evtntrk",sample="data18",sel="scouting"):

    rhos = []
    rhos.append("rho0to0p5")
    rhos.append("rho0p5to1")
    rhos.append("rho1to2")
    rhos.append("rho2to5")
    rhos.append("rho5to10")
    #rhos.append("rhoinf")

    hists=[]
    labels=[]
    for rho in rhos:
        name = sel+"_"+dist+"_"+rho+"_"+kind
        print(name)
        hist = get1D(sample,name)
        print(hist)
        if hist : 
            hist.Rebin(5)
            clean1D(hist)
            hists.append(hist)
            labels.append(rho)

    compare1D(hists,labels,"rhoStudy/compare_correlation_{}_{}".format(dist,kind),opt="err")

dists = []
dists.append("jetntrk")
dists.append("evtntrk")
#dists.append("jet_pt")
#dists.append("jet_eta")
for dist in dists:
    compareNtrackVRho("jetsAK15_suep_rho_dR0p05_1",dist)
    compareNtrackVRho("jetsAK15_suep_rho_dR0p05_2",dist)
    compareNtrackVRho("jetsAK15_suep_rho_dR0p05_1","rho0low_"+dist)
    compareNtrackVRho("jetsAK15_suep_rho_dR0p05_2","rho0low_"+dist)

compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_1","jetntrk")
compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_2","jetntrk")

compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_1","evtntrk")
compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_2","evtntrk")
compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_1","rho0low_evtntrk")
compareRhoVNtrack("jetsAK15_suep_rho_dR0p05_2","rho0low_evtntrk")

r0bins = ["rho0low", "rho0high"]
for rbin in range(3):

    profile2D("data18","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i"%rbin)
    profile2D("data18","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i"%rbin)
    profile2D("data18","scouting_jetsAK15_suep_jetntrk_v_rho_dR0p05_%i"%rbin)
    #profile2D("data18","scouting_jetsAK15_suep_jetpt_v_rho_dR0p05_%i"%rbin)

    for r0bin in r0bins: 
        profile2D("mMed-300_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%sL"%(rbin,r0bin))


        profile2D("QCD","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%sL"%(rbin,r0bin))
        profile2D("QCD","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))
        profile2D("QCD","scouting_jetsAK15_suep_jetntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))

        profile2D("data18","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%sL"%(rbin,r0bin))
        profile2D("data18","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))
        profile2D("data18","scouting_jetsAK15_suep_jetntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))
        #profile2D("data18","scouting_jetsAK15_suep_jetpt_v_rho_dR0p05_%i_%s"%(rbin,r0bin))
    
    #profile2D("mMed-200_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i"%rbin)
    #profile2D("mMed-200_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_jetntrk_v_rho_dR0p05_%i"%rbin)
    #profile2D("mMed-200_mDark-2_temp-2_decay-darkPho","scouting_jetsAK15_suep_jetpt_v_rho_dR0p05_%i"%rbin)


