import ROOT
import ctypes
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()

lumi=0.1 #fb?


def getStatUnc(a,b,c,d,erra,errb,errc,errd):

    part1 = erra.value/a
    part2 = errb.value/b
    part3 = errc.value/c
    err = (part1**2+part2**2+part3**2)**0.5
    #print(part1,part2,part3,err)


    return b/a*c*err


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

    erra=ctypes.c_double(0.)
    errb=ctypes.c_double(0.)
    errc=ctypes.c_double(0.)
    errd=ctypes.c_double(0.)
    a = round(hist.IntegralAndError( xlow, xcut , ylow, ycut ,erra),1) # lwo ntrk low tau 
    b = round(hist.IntegralAndError( xcut, xhigh, ylow, ycut ,errb),1) # lwo ntrk high tau 
    c = round(hist.IntegralAndError( xlow, xcut , ycut, -1   ,errc),1) # high ntrk, low tau 
    d = round(hist.IntegralAndError( xcut, xhigh, ycut, -1   ,errd),1) # high ntrk, high tau 

    dpred = round(b/a*c,2)

    # 
    # signal
    sample="mMed-300_mDark-2_temp-2_decay-darkPho"
    histSig  = get1D(sample,dist)
    dSig = round(histSig.Integral( xcut, xhigh, ycut, -1   ),1) # high ntrk, high tau 

    

    stat_unc = round( getStatUnc(a,b,c,d,erra,errb,errc,errd),2)

    errdobs = errd.value
    if errdobs==0: errdobs = 1.0 

    pull = round( abs(dpred-d)/(stat_unc + errdobs ),2 ) 


    nonclosure = round( max(  max( abs(dpred-d)-stat_unc ,0) -errdobs ,0)/((dpred+d)/2.), 2 )

    #signif = round(dSig / (dpred + dSig + ( dpred**2 ) )**0.5 ,1) 
    signif = round(dSig / (dpred + dSig + (nonclosure*dpred)**2 +stat_unc**2 )**0.5 ,1) 
    # 
    # Output
    #

    print(a,",",b,",",c,",",dpred,"pm",stat_unc,",",d,"pm",round(errdobs,1),",",dSig,",",pull,",",nonclosure,",",signif)
    #print(nonclosure, stat_unc, signif)

    return

#
# Event ntracks
#
dists=[]
dists.append("sphericity")
dists.append("circularity")
dists.append("c")
dists.append("d")

ntrkcuts = [90,100,110,120,130
]
evtshapecuts = [0.5,0.6,0.7,0.8]

for dist in dists: 
    print(dist)

    for evtshapecut in evtshapecuts:
        print(evtshapecut, "vary ntrk")
        for ntrkcut in ntrkcuts: 
            ABCDtest(10, ntrkcut,500, 0.05,evtshapecut,1.0,"scouting_boosted_evtshape_ntracks_v_{}".format(dist))


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
