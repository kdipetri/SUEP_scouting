import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?



def ABCDtest(ntrklow, ntrkcut, ntrkhigh, taulow, taucut, tauhigh , dist="jetsAK15_suep_rho_dR0p05_1", kind="evtntrk"):

    
    sample="data18"
    hist  = get1D(sample,dist)
    #print(dist)

    xlow  = hist.GetXaxis().FindBin(taulow)
    xcut  = hist.GetXaxis().FindBin(taucut)
    xhigh = hist.GetXaxis().FindBin(tauhigh)
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

#
# Event ntracks
#
ntrkcuts = [70,75,80,85,90]
taucuts = [0.55,0.60,0.65,0.7]
for taucut in taucuts:
    print(taucut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest(10, ntrkcut,500, 0.0,taucut,1.0,"scouting_jetsAK15_suep_evtntrk_v_nsub21")


#
# Jet ntracks
#
#ntrkcuts = [40,45,50]
#taucuts = [0.6,0.65,0.7]
#for taucut in taucuts:
#    print(taucut, "vary ntrk")
#    for ntrkcut in ntrkcuts: 
#        ABCDtest(0, ntrkcut,500, 0.1,taucut,1.0,"scouting_jetsAK15_suep_jetntrk_v_nsub21")
#