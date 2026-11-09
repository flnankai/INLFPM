"""Finite-sample tail calculations for the archived panel generator."""
import numpy as np
from scipy.optimize import nnls, brentq
from scipy.special import ndtr, stdtr
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from panel_math import spatial_median_fit, pqr_scale, in_fit, FitError

def sidak(p,n):
    return -np.expm1(n*np.log1p(-p)) if p<1 else 1.

class Reference:
    def __init__(self,y,pr):
        self.pr = pr
        self.z = pr.project(y)
        self.e = self.z - np.outer(pr.h, pr.h@self.z/pr.h2)
        self.M = np.eye(len(y))-pr.q@pr.q.T
        self.MX = self.M-np.outer(pr.h,pr.h)/pr.h2
        self.AX = self.MX**2
        self.chol = cho_factor(self.AX)

    def nuisance(self,d):
        e = self.e / d
        s = np.mean(e**2,axis=1)
        lam = cho_solve(self.chol,s)
        T,N=e.shape
        tol=64*np.finfo(float).eps*T
        if lam.min() < -tol*max(1.,abs(lam).max()):
            lam=nnls(self.AX,s,maxiter=30*T)[0]
        else:
            lam=np.maximum(lam,0)
        C=(self.MX*(lam/np.mean(lam)))@self.MX
        ev,U=np.linalg.eigh(C)
        keep=ev>tol*max(1.,abs(ev).max())
        r=keep.sum()
        if r<=2: raise FitError('insufficient whitening rank')
        W=(U[:,keep].T@e)/np.sqrt(ev[keep,None])
        gram=W@W.T
        g1=np.trace(gram)
        g2=np.sum(gram**2)
        gamma=max(1/N,(g2/g1**2-1/r)/(1+1/r-2/r**2))
        return lam,gamma

    def spatial(self,z):
        T,N=z.shape
        fit=spatial_median_fit(z)
        d,r=fit['d'],fit['r']
        a=1/r
        eta=self.pr.eta
        B=1-2*eta*a.mean()*r.mean()+eta*np.mean(r**2)*a.mean()**2
        numerator=T*N*a.mean()**2*np.max((fit['theta']/d)**2)
        raw=numerator/B
        previous=numerator/(1-2*eta*a.mean()*r.mean()+eta*np.mean(r**2)*np.mean(a**2))
        lam,gamma=self.nuisance(d)
        cj=1-gamma*(1-eta*a.mean()*r.mean())/B
        if cj<=0: raise FitError('nonpositive SS curvature')
        a=a/a.mean()
        v=self.M@a
        C=self.M-v[None,:]/a.sum()
        profile=np.sqrt(lam/lam.mean())
        C=(C*profile)*a[:,None]
        v=v*profile
        Q=np.outer(v,v)-(cj**2*raw*B)*(C.T@C)
        ev=np.linalg.eigvalsh(Q)
        tol=64*np.finfo(float).eps*T*max(1,abs(ev).max())
        if ev[-1]<=tol: return 0.,previous
        ratio=-ev[ev < -tol]/ev[-1]
        def integrand(x):
            return np.exp(-.5*np.sum(np.log1p(ratio/np.sin(x)**2)))
        marginal=2/np.pi*quad(integrand,0,np.pi/2,epsabs=1e-15,epsrel=1e-10)[0]
        return sidak(marginal,N),previous

    def inverse(self):
        d=pqr_scale(self.e)
        fit=in_fit(self.z,self.pr,d)
        lam,gamma=self.nuisance(d)
        r=fit['r']
        v=self.M@(self.pr.h/r**2)
        T,N=self.z.shape
        nu=self.pr.df
        k=(nu*N/T)*np.sum(lam*v**2)/np.sum((r*v)**2)
        kappa=pqr_kappa(lam)
        if kappa<2-1e-6: raise FitError('PQR moment out of range')
        rho=1+(kappa-2)/nu
        m=4+2/(rho*(nu-2)/(nu-4)-1)
        scale=(nu/(nu-2))/(m/(m-2))
        cj=1-2*gamma
        if cj<=0: raise FitError('nonpositive IN curvature')
        return fit['alpha'],np.sqrt(nu*N)*fit['A']/(d*np.sqrt(fit['Bc'])),cj/np.sqrt(k*scale),m

def pqr_kappa(lam):
    T=len(lam)
    lm=lam/lam.mean()
    pair=np.sqrt(lm[:,None]+lm[None,:])
    ps=pair[np.triu_indices(T,1)]
    total=len(ps)
    ps=ps[ps>0]
    nz=total-len(ps)
    if nz/total>=.5: raise FitError('zero PQR median')
    def cdf(q):
        return (nz+np.sum(2*ndtr(q/ps)-1))/total-.5
    q=brentq(cdf,0,ps.max()*.69,xtol=1e-13)
    u=2/np.sqrt(2*np.pi)*np.sum(np.exp(-.5*(q/ps)**2)/ps)/total
    P=np.ones_like(pair)
    keep=pair>0
    P[keep]=2*ndtr(q/pair[keep])-1
    np.fill_diagonal(P,0)
    sd=np.sqrt(lm[lm>0])
    pj=P.sum(axis=1)[lm>0]/(T-1)
    zeros=T-len(sd)
    coeff=16/(T*(q*u)**2)
    def fn(y,below):
        x=np.exp(y)
        a=ndtr((x+q)/sd)-ndtr((x-q)/sd)
        c=(a.sum()+zeros*below-a)/(T-1)-pj
        z=x/sd
        return np.sum(z*np.exp(-.5*z*z)*c*c)*np.sqrt(2/np.pi)
    kwargs=dict(epsabs=min(1e-10,2e-9/coeff/10),epsrel=2e-11,limit=300)
    l=quad(fn,np.log(sd.min())-48,np.log(q),args=(True,),**kwargs)[0]
    h=quad(fn,np.log(q),np.log(12*sd.max()),args=(False,),**kwargs)[0]
    return coeff*(l+h)
