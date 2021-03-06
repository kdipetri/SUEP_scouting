import ROOT
import ctypes
from array import array
from plothelper import *
ROOT.gROOT.SetBatch(ROOT.kTRUE)

setStyle()



def getStatUnc(a,b,c,d,erra,errb,errc,errd):

    part1 = erra.value/a
    part2 = errb.value/b
    part3 = errc.value/c
    err = (part1**2+part2**2+part3**2)**0.5
    #print(part1,part2,part3,err)


    return b/a*c*err

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

    erra=ctypes.c_double(0.)
    errb=ctypes.c_double(0.)
    errc=ctypes.c_double(0.)
    errd=ctypes.c_double(0.)
    a = round ( hist.IntegralAndError( xlow, xcut , ylow, ycut ,erra) , 1) # lwo ntrk low rho
    b = round ( hist.IntegralAndError( xcut, xhigh, ylow, ycut ,errb) , 1) # lwo ntrk high rho
    c = round ( hist.IntegralAndError( xlow, xcut , ycut, -1 ,errc) , 1) # high ntrk, low rho
    d = round ( hist.IntegralAndError( xcut, xhigh, ycut, -1 ,errd
    ) , 1) # high ntrk, high rho

    dpred = round(b/a*c,2)

    # signal
    sample="mMed-300_mDark-2_temp-2_decay-darkPho"
    histSig  = get1D(sample,dist)
    dSig = round(histSig.Integral( xcut, xhigh, ycut, -1   ),1) # high ntrk, high tau 

    # metrics
    stat_unc = round( getStatUnc(a,b,c,d,erra,errb,errc,errd),2)
    errdobs = errd.value
    if errdobs==0: errdobs = 1.0 

    pull = round( abs(dpred-d)/(stat_unc + errdobs ),2 ) 

    nonclosure = round( max(  max( abs(dpred-d)-stat_unc ,0) -errdobs ,0)/((dpred+d)/2.), 2 )

    #signif = round(dSig / (dpred + dSig + ( dpred**2 ) )**0.5 ,1) 
    signif = round(dSig / (dpred + dSig + (nonclosure*dpred)**2 +stat_unc**2 )**0.5 ,1) 
    
    # Output
    print(a,",",b,",",c,",",dpred,"pm",stat_unc,",",d,"pm",round(errdobs,1),",",dSig,",",pull,",",nonclosure,",",signif)
    #print(nonclosure, stat_unc, signif)



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
    print(rhocut, "vary ntrk")
    for ntrkcut in ntrkcuts: 
        ABCDtest3(0, ntrkcut,500, 0.1,rhocut,10.0,"scouting_jetsAK15_suep_evtntrk_v_rho_dR0p05_%i_%sL"%(rbin,r0bin))



