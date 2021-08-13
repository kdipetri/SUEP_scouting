import ROOT
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?


def ABCDtest3(ntrklow, ntrkcut, ntrkhigh, rholow, rhocut, rhohigh , dist="jetsAK15_suep_rho_dR0p05_1", kind="evtntrk"):

    # background
    sample="QCD"
    hist  = get1D(sample,dist)

    xlow  = hist.GetXaxis().FindBin(rholow)
    xcut  = hist.GetXaxis().FindBin(rhocut)
    xhigh = hist.GetXaxis().FindBin(rhohigh)
    ylow  = hist.GetYaxis().FindBin(ntrklow)
    ycut  = hist.GetYaxis().FindBin(ntrkcut)
    yhigh = hist.GetYaxis().FindBin(ntrkhigh)
    #print (xlow,xcut,xhigh,ylow,ycut,yhigh)

    a = round ( hist.Integral( xlow, xcut , ylow, ycut ) , 1) # lwo ntrk low rho
    b = round ( hist.Integral( xcut, xhigh, ylow, ycut ) , 1) # lwo ntrk high rho
    c = round ( hist.Integral( xlow, xcut , ycut, -1 ) , 1) # high ntrk, low rho
    d = round ( hist.Integral( xcut, xhigh, ycut, -1 ) , 1) # high ntrk, high rho

    dpred = round(b/a*c,2)


    # signal
    sample="mMed-300_mDark-2_temp-2_decay-darkPho"
    histSig  = get1D(sample,dist)
    dSig = round(histSig.Integral( xcut, xhigh, ycut, -1   ),1) # high ntrk, high tau 

    signif = round(dSig / (dpred + dSig + (0.5*dpred)**2 )**0.5 ,1) 

    # output
    print(a,b,c,dpred,d,dSig, signif)

    return

#r0bin = "rho0low"
#rbin=2
#ntrkcuts = [70,75,80,85,90]
#rhocuts = [0.7,1.0,1.5]
#for rhocut in rhocuts:
#    #print(rhocut, "vary ntrk")
#    for ntrkcut in ntrkcuts: 
#        ABCDtest3(0, ntrkcut,500, 0.,rhocut,7.0,"scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%s"%(rbin,r0bin))

print("")
print("")
r0bin = "rho0low"
rbin=0
ntrkcuts = [90, 100, 110, 120, 130]
rhocuts  = [0.7,1.0,1.5]
for rhocut in rhocuts:
    #print(rhocut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest3(0, ntrkcut,500, 0.1,rhocut,10.0,"scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%sL"%(rbin,r0bin))



