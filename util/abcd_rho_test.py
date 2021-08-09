import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?

#
# now compare Rho V ntrack for data
# 

def ABCDtest(rhocut = 0.9, lowntrk = "ntrk30to50", highntrk = "ntrk70toinf", dist="jetsAK15_suep_rho_dR0p05_1"):
    
    sample="data18"
    sel="scouting"

    hist_lowntrk   = get1D(sample,sel+"_"+lowntrk+"_"+dist)
    hist_highntrk  = get1D(sample,sel+"_"+highntrk+"_"+dist)

    a = hist_lowntrk.Integral(0,hist_lowntrk.FindBin( rhocut ))
    b = hist_lowntrk.Integral(hist_lowntrk.FindBin( rhocut ),-1)

    c = hist_highntrk.Integral(0,hist_highntrk.FindBin( rhocut ))
    d = hist_highntrk.Integral(hist_highntrk.FindBin( rhocut ),-1)

    dpred = round(b/a*c,2)
    print(a,b,c,dpred,d)
    return

#rhocuts = [0.75,1.0,1.5,2.0]
#for rhocut in rhocuts:
#    print(rhocut)
#    ABCDtest(rhocut, "ntrk30to50","ntrk70toinf", "jetsAK15_suep_rho_dR0p05_1")
#    ABCDtest(rhocut, "ntrk30to50","ntrk50to70" , "jetsAK15_suep_rho_dR0p05_1")
#    ABCDtest(rhocut, "ntrk50to70","ntrk70toinf", "jetsAK15_suep_rho_dR0p05_1")
#
#print("")
#print("")
#print("")
#for rhocut in rhocuts:
#    print(rhocut)
#    ABCDtest(rhocut, "ntrk30to50","ntrk70toinf", "jetsAK15_suep_rho_dR0p05_1")
#    ABCDtest(rhocut, "ntrk30to50","ntrk50to70" , "jetsAK15_suep_rho_dR0p05_1")
#    ABCDtest(rhocut, "ntrk50to70","ntrk70toinf", "jetsAK15_suep_rho_dR0p05_1")

#def ABCDtest2(ntrkcut, lowrho , highrho , dist="jetsAK15_suep_rho_dR0p05_1", kind="evtntrk"):
#
#    #"0to0p5";
#    #"0p5to1";
#    #"1to2";
#    #"2to5";
#    #"5to10";
#    
#    sample="data18"
#    sel="scouting"
#    #print(sel+"_"+lowrho +"_"+dist+"_"+kind)
#
#    hist_lowrho   = get1D(sample,sel+"_"+dist+"_"+lowrho +"_"+kind)
#    hist_highrho  = get1D(sample,sel+"_"+dist+"_"+highrho+"_"+kind)
#
#    a = hist_lowrho.Integral(0,hist_lowrho.FindBin( ntrkcut ))
#    c = hist_lowrho.Integral(hist_lowrho.FindBin( ntrkcut ),-1)
#
#    b = hist_highrho.Integral(0,hist_lowrho.FindBin( ntrkcut ))
#    d = hist_highrho.Integral(hist_lowrho.FindBin( ntrkcut ),-1)
#
#    dpred = round(b/a*c,2)
#    print(a,b,c,dpred,d)
#    return
#
#r0bins = ["rho0low", "rho0med", "rho0high"]
#ntrkcuts = [70]
#for r0bin in r0bins:
#    for ntrkcut in ntrkcuts:
#        print(r0bin, ntrkcut)
#        ABCDtest2(ntrkcut, "rho0to0p5_%s"%r0bin, "rho5to10_%s"%r0bin, "jetsAK15_suep_rho_dR0p05_1", "evtntrk")
#        ABCDtest2(ntrkcut, "rho0p5to1_%s"%r0bin, "rho5to10_%s"%r0bin, "jetsAK15_suep_rho_dR0p05_1", "evtntrk")
#        ABCDtest2(ntrkcut, "rho0to0p5_%s"%r0bin, "rho2to5_%s"%r0bin , "jetsAK15_suep_rho_dR0p05_1", "evtntrk")
#        ABCDtest2(ntrkcut, "rho0p5to1_%s"%r0bin, "rho2to5_%s"%r0bin , "jetsAK15_suep_rho_dR0p05_1", "evtntrk")

def ABCDtest3(ntrklow, ntrkcut, ntrkhigh, rholow, rhocut, rhohigh , dist="jetsAK15_suep_rho_dR0p05_1", kind="evtntrk"):

    
    sample="data18"
    hist  = get1D(sample,dist)
    #print(dist)

    xlow  = hist.GetXaxis().FindBin(rholow)
    xcut  = hist.GetXaxis().FindBin(rhocut)
    xhigh = hist.GetXaxis().FindBin(rhohigh)
    ylow  = hist.GetYaxis().FindBin(ntrklow)
    ycut  = hist.GetYaxis().FindBin(ntrkcut)
    yhigh = hist.GetYaxis().FindBin(ntrkhigh)
    #print (xlow,xcut,xhigh,ylow,ycut,yhigh)

    a = hist.Integral( xlow, xcut , ylow, ycut ) # lwo ntrk low rho
    b = hist.Integral( xcut, xhigh, ylow, ycut ) # lwo ntrk high rho

    c = hist.Integral( xlow, xcut , ycut, -1 ) # high ntrk, low rho
    d = hist.Integral( xcut, xhigh, ycut, -1 ) # high ntrk, high rho

    dpred = round(b/a*c,2)
    print(a,b,c,dpred,d)
    return

r0bin = "rho0low"
rbin=2
ntrkcuts = [70,75,80,85,90]
rhocuts = [0.7,1.0,1.5]
for rhocut in rhocuts:
    #print(rhocut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest3(0, ntrkcut,500, 0.,rhocut,7.0,"scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))

print("")
print("")
r0bin = "rho0low"
rbin=0
ntrkcuts = [70,75,80,85,90]
rhocuts = [0.7,1.0,1.5]
for rhocut in rhocuts:
    #print(rhocut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest3(0, ntrkcut,500, 0.,rhocut,10.0,"scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))
