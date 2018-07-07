import numpy as np
from scipy import optimize
from numpy import mean

def Ftest(C1,C2,a,b):
    """return True if the null Hypothesis apply. ie, if the C2
       fit isnt statistically better than the C1 fit"""
    alpha=alphadet(a,b)
    if (C1/C2-1)<alpha*a/b:
        return True
    else:
        return False

def Ffunc(x,a,b):
    from scipy import special as SP
    a=float(a)
    b=float(b)
    j=a*x/(a*x+b)
    F=SP.betainc(a/2,b/2,j)
    return F
    
def alphadet(a,b):
    acc=0.99
    C=0
    up=0.0
    while C<acc:
        up+=1.0
        C=Ffunc(up,a,b)
    dwn=up-1.0
    mid=(dwn+up)/2.0
    while(abs(C-acc)>0.00001):
        mid=(dwn+up)/2.0
        C=Ffunc(mid,a,b)
        if C<acc:
            dwn=mid
        else:
            up=mid
        #print C
    return mid


def Fourier(time,Y,sigma=None,order=5,P=1):
    """Assumes that the LC is phased, if not, a Period different from 1 is expected"""
    global func
    global func_a
    
    if(order<0):order=0
    cmd="func_a= lambda t, A: A[0]"
    for i in range(order):
        cmd+="+A[{0}]*np.sin(2*np.pi*{1}*t/{3}+A[{2}])".format(2*i+1,i+1,2*i+2,P)
    exec cmd in globals()
    A1=""
    A2=""
    for i in range(order):
        A1+=", A{0}, A{1}".format(2*i+1,2*i+2)
        A2+="+A{0}*np.sin(2*np.pi*{1}*t/{3}+A{2})".format(2*i+1,i+1,2*i+2,P)
    cmd="func= lambda t, A0{0}: A0{1}".format(A1,A2)
    exec cmd in globals()
    if sigma==None:
        p1, pcov=optimize.curve_fit(func, time, Y)
    else:
        p1, pcov=optimize.curve_fit(func, time, Y,sigma=sigma)
    
    
    if sigma==None:
        A=np.power(Y-func_a(time,p1),2)
    else:
        A=np.power(Y-func_a(time,p1),2)/np.power(sigma,2)
    Chi2=sum(A)/(len(time)-(2*i+1)-1)

    return func_a, p1, Chi2

def sFourier(time,Y,sigma=None,Print=False,sOrder=1,P=1,maxorder=6):
    #sOrder = Starting order
    if sOrder>maxorder or sOrder<0:
        sOrder=1
    Chi2=[]
    
    funca, p1a, C=Fourier(time, Y,sigma,order=sOrder,P=P)
    sOrder+=1
    Chi2.append(C)
    while sOrder<=maxorder:
        func, p1, C=Fourier(time, Y,sigma,order=sOrder,P=P)
        Chi2.append(C)
        a=2
        b=len(time)-(sOrder*2+1)
        if Ftest(Chi2[-2],Chi2[-1],a,b):
            sOrder+=1
            if sOrder>maxorder: return funca, p1a, Chi2[-2]
            func, p1, C=Fourier(time, Y,sigma,order=sOrder,P=P)
            Chi2.append(C)
            a=4
            b=len(time)-(sOrder*2+1)
            if Ftest(Chi2[-3],Chi2[-1],a,b):
                return funca,p1a,Chi2[-3]
        if Print==True: print Chi2
        funca=func
        p1a=p1
        sOrder+=1
    return func, p1, Chi2[-1]

def cutsFourier(time,Y,sigma=None,order=None,sOrder=1,LCcut=False,P=1,Gfill=False,ret=False):
    
    #if P=1, phased curve is assumed
    
    
    if order==None:
        func, p1, Chi2=sFourier(time, Y,sigma,sOrder=sOrder,P=P)
    else:
        func, p1, Chi2=Fourier(time,Y,sigma,order=order,P=P)
    T=np.linspace(0,1,100)
    Max=max(func(T,p1))
    Min=min(func(T,p1))

    l=len(time)

    amp=6*np.std(Y)
    Mean=(Max+Min)/2
    diff=abs(Y-Mean)
    time=time[diff<amp/2]
    Y=Y[diff<amp/2]
    if sigma!=None:
        sigma=sigma[diff<amp/2]
        
    
    if l!=len(time):
        if order==None:
            func, p1, Chi2=sFourier(time,Y,sigma,sOrder=sOrder,P=P)
        else:
            func, p1, Chi2=Fourier(time,Y,sigma,order=order,P=P)

    l=len(time)
    if LCcut and (sigma!= None):
        print sigma
        merr=mean(sigma)
        print merr
        diff2=abs(Y-func(time,p1))
        time=time[diff2<merr*3]
        Y=Y[diff2<merr*3]
        sigma=sigma[diff2<merr*3]
        if l!=len(time):
            if order==None:
                func, p1, Chi2=sFourier(time,Y,sigma,sOrder=sOrder,P=P)
            else:
                func, p1, Chi2=Fourier(time,Y,sigma,order=order,P=P)
    l=len(time)
    if Gfill:
        if P!=1: print "Gap Filling with non phased LightCurves isn't available yet."
        else:
            T=np.array(time)
            T.sort()
            L=len(time)
            for i in range(L-1):
                l=T[i+1]-T[i]
                if l>0.1:
                    fp=(T[i]+T[i+1])/2
                    a=i+2
                    b=i+3
                    C=0
                    if i+2==L:
                        a=0
                        b=1
                        C=2                        
                    elif i+3==L:
                        b=0
                        C=1
                    ep=(T[i+1]+T[a]+T[b]+C)/3
                    C=0
                    if i-1<0:
                        C=2
                    elif i-2<0:
                        C=1
                    sp=(T[i]+T[i-1]+T[i-2]-C)/3
                    #val=(func(T[i+1],p1)+func(T[i],p1))/2
                    val=(func(ep,p1)+func(sp,p1))/2
                    time=np.append(time,fp)
                    Y=np.append(Y,val)
                    if sigma!=None:
                        merr= mean(sigma)
                        sigma=np.append(sigma,merr/2)
                    
                    #print "Gap found: {0}, {1}, {2}".format(l,fp,val)
                    #print "         : {0}, {1}".format(sp,ep)
            l=1+T[0]-T[-1]
            if l>0.1:
                fp=(T[0]+1+T[-1])/2
                if fp>1:fp-=1
                sp=(T[-1]+T[-2]+T[-3])/3
                ep=(T[0]+T[1]+T[2])/3
                #val=(func(T[0],p1)+func(T[-1],p1))/2
                val=(func(ep,p1)+func(sp,p1))/2
                time=np.append(time,fp)
                Y=np.append(Y,val)    
                sigma=np.append(sigma,merr/2)
                
                #print "Gap found: {0}, {1}, {2}".format(l,fp,val)
        if l!=len(time):
            if order==None:
                func, p1, Chi2=sFourier(time,Y,sigma,sOrder=sOrder,P=P)
            else:
                func, p1, Chi2=Fourier(time,Y,sigma,order=order,P=P)

    if ret:
        return func,p1,Chi2,{"phased":time,"Var2":Y,"Var3":sigma}
    else:
        return func, p1, Chi2

def phase(time,period):
    Phased = (time%period)/period
    return Phased
