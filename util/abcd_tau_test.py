import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?



def ABCDtest(ntrklow, ntrkcut, ntrkhigh, taulow, taucut, tauhigh , dist="jetsAK15_suep_rho_dR0p05_1", kind="evtntrk"):

    #
    # background
    #
    sample="QCD"
    hist  = get1D(sample,dist)

    xlow  = hist.GetXaxis().FindBin(taulow)
    xcut  = hist.GetXaxis().FindBin(taucut)
    xhigh = hist.GetXaxis().FindBin(tauhigh)
    ylow  = hist.GetYaxis().FindBin(ntrklow)
    ycut  = hist.GetYaxis().FindBin(ntrkcut)
    yhigh = hist.GetYaxis().FindBin(ntrkhigh)
    #print (xlow,xcut,xhigh,ylow,ycut,yhigh)

    a = round(hist.Integral( xlow, xcut , ylow, ycut ),1) # lwo ntrk low tau 
    b = round(hist.Integral( xcut, xhigh, ylow, ycut ),1) # lwo ntrk high tau 
    c = round(hist.Integral( xlow, xcut , ycut, -1   ),1) # high ntrk, low tau 
    d = round(hist.Integral( xcut, xhigh, ycut, -1   ),1) # high ntrk, high tau 

    dpred = round(b/a*c,2)

    # 
    # signal
    sample="mMed-300_mDark-2_temp-2_decay-darkPho"
    histSig  = get1D(sample,dist)
    dSig = round(histSig.Integral( xcut, xhigh, ycut, -1   ),1) # high ntrk, high tau 

    signif = round(dSig / (dpred + dSig + (0.5*dpred)**2 )**0.5 ,1) 
    # 
    # Output
    #

    print(a,b,c,dpred,d,dSig, signif)

    return

#
# Event ntracks
#
ntrkcuts = [90,100,110,120,130]
taucuts = [0.65,0.7,0.75]
for taucut in taucuts:
    #print(taucut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest(10, ntrkcut,500, 0.0,taucut,1.0,"scouting_jetsAK15_suep_evtntrk_v_nsub21L")


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
