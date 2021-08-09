import ROOT

hists1D = {}
hists2D = {}
plotdir = "plots/hists"

doPng=True
doPdf=False
doC=False

one   = ROOT.TColor(2001,143/255.,45 /255.,86/255.,"darkPurple")#quinacridone magenta
two   = ROOT.TColor(2002,119/255.,104/255.,174/255.,"blahBlue")#blue-violet
three = ROOT.TColor(2003,239/255.,71 /255.,111/255.,"pinkRed")#paradise pink
four  = ROOT.TColor(2004,247/255.,178/255.,103/255.,"orange")#orange
five  = ROOT.TColor(2005,42 /255.,157/255.,143/255.,"PersianGreen")# persian green
six   = ROOT.TColor(2006,38 /255.,70 /255.,83 /255.,"Charcol")# charcol
seven = ROOT.TColor(2007,116/255.,165/255.,127/255.,"Green")#forest green
eight = ROOT.TColor(2008,233/255.,196/255.,106/255.,"Maize")# maize
nine  = ROOT.TColor(2009,8/255.,103/255.,136/255.,"RussianViolet")#russian violet 
ten   = ROOT.TColor(2010,231/255.,111/255.,81 /255.,"TerraCotta")# terra cotta
colors = [] #[2001,2002,2003,2004,2005,2006,2007,2008,2009,2010]
colors.append(2003)#paradise
colors.append(2004)#orange
colors.append(2005)#persian green
colors.append(2002)#blue-violet
colors.append(2001)#quinacridone magenta
colors.append(2010)#terra cotta
colors.append(2008)#maize
colors.append(2007)#forest green
colors.append(2009)#bluesapphire
colors.append(2006)#charcol

def setStyle():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetLabelFont(42,"xyz")
    ROOT.gStyle.SetLabelSize(0.05,"xyz")
    ROOT.gStyle.SetTitleFont(42,"xyz")
    ROOT.gStyle.SetTitleFont(42,"t")
    ROOT.gStyle.SetTitleSize(0.06,"xyz")
    ROOT.gStyle.SetTitleSize(0.06,"t")

    ROOT.gStyle.SetPadBottomMargin(0.14)
    ROOT.gStyle.SetPadLeftMargin(0.14)

    ROOT.gStyle.SetPadGridX(0)
    ROOT.gStyle.SetPadGridY(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetTitleOffset(1,'y')
    ROOT.gStyle.SetLegendTextSize(0.04)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridColor(14)

    ROOT.gStyle.SetMarkerSize(1.0) #large markers
    ROOT.gStyle.SetHistLineWidth(2) # bold lines
    ROOT.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes

    #one = ROOT.TColor(2001,0.906,0.153,0.094)
    #two = ROOT.TColor(2002,0.906,0.533,0.094)
    #three = ROOT.TColor(2003,0.086,0.404,0.576)
    #four =ROOT.TColor(2004,0.071,0.694,0.18)
    #five =ROOT.TColor(2005,0.388,0.098,0.608)
    #six=ROOT.TColor(2006,0.906,0.878,0.094)
    #seven=ROOT.TColor(2007,0.99,0.677,0.614)
    #colors = [1,2001,2002,2003,2004]
    return       


def plot1D(name, title, x, nbinsx, xmin, xmax, weight=1.):
    if name in hists1D:
        # fill
        hists1D[name].Fill(x,weight)
    else : 
        # create and fill
        hist = ROOT.TH1F(name, title, nbinsx, xmin, xmax)
        hist.SetDirectory(0)
        hist.Fill(x,weight)
        hists1D[name] = hist
    return

def plot2D(name, title, x, y, nbinsx, xmin, xmax, nbinsy, ymin, ymax, weight=1.):
    if name in hists2D:
        # fill
        hists2D[name].Fill(x,y,weight)
    else : 
        # create and fill
        hist = ROOT.TH2F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
        hist.SetDirectory(0)
        hist.Fill(x,weight)
        hists2D[name] = hist

    return

def draw1D(c1,h, drawopt="hist"):  
    c1.cd() 
    c1.Clear()
    h.Draw(drawopt)
    if doPng: c1.Print("{}/{}.png".format(plotdir,h.GetName()))
    if doPdf: c1.Print("{}/{}.pdf".format(plotdir,h.GetName()))
    if doC  : c1.Print("{}/{}.C".format(plotdir,h.GetName()))
    h.Write()
    return 

def draw2D(c2, h, drawopt="COLZ"):  
    c2.cd() 
    c2.Clear()
    c2.SetTopMargin(0.05)
    c2.SetLeftMargin(0.2)
    c2.SetBottomMargin(0.2)
    c2.SetRightMargin(0.2);
    h.Draw(drawopt)

    if doPng: c2.Print("{}/{}.png".format(plotdir,h.GetName()))
    if doPdf: c2.Print("{}/{}.pdf".format(plotdir,h.GetName()))
    if doC  : c2.Print("{}/{}.C".format(plotdir,h.GetName()))
    h.Write()
    return 

def drawAll1D(c1,drawopt="hist"):
    for n, h in hists1D.items():
        draw1D(c1,h, drawopt);
    return

def drawAll2D(c2,drawopt="hist"):
    for n, h in hists2D.items():
        draw2D(c2,h, drawopt);
    return 


def get1D(sample,dist):

    # Get hist
    if "QCD" in sample: hist = getQCD(dist)
    else : 
        filename = "output/{}.root".format(sample)
        f = ROOT.TFile.Open(filename)
        histname = "{}_{}".format(sample,dist)
        hist = f.Get(histname)
    
    # if hist exists
    if hist : 
        hist.SetDirectory(0)
        return hist
    else : return 0

def qcd_xs(sample):

    if "QCD_HT200to300"   in sample : return 1559000 / 100000 * 135000
    if "QCD_HT300to500"   in sample : return 347700 / 100000 * 135000 #311900 - was missing k-factor
    if "QCD_HT500to700"   in sample : return 32100 / 100000 * 135000  #29070 - was missing k-factor
    if "QCD_HT700to1000"  in sample : return 6831 / 100000 * 135000  #5962 - was missing k-factor
    if "QCD_HT1000to1500" in sample : return 1207 / 100000 * 135000  
    if "QCD_HT1500to2000" in sample : return 119.9 / 7242 * 135000  # * 100000./7242. # missing some events
    if "QCD_HT2000toInf"  in sample : return 25.24 / 45606 * 135000  # * 100000./45606. # missing some events

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

    if len(hists)==0: return 0 
    hist_final = hists[0].Clone("QCD_"+dist)
    for i,hist in enumerate(hists):
        if i>0: hist_final.Add(hist)

    #if "scouting_htL" in dist:
    #    print("KARRI TEST") 
    #    print("INT {}".format( hist_final.Integral() ) ) 

    #if "offline_htL" in dist:
    #    print("KARRI TEST") 
    #    print("INT {}".format( hist_final.Integral() ) ) 
    #hist_final.Scale(40048/hist_final.Integral()) #40048 events in data passing hT cut...
    #clean1D(hist_final)
    hist_final.SetDirectory(0)

    return hist_final
