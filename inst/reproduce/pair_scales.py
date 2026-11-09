"""Original Stat inverse-norm pair kernel with explicit design correction.

Zhao, Chen and Zi (2022), Stat 11 e490, pp. 3--4.  Each pair uses its own
diagonal sign-scale fit from the full unrestricted OLS residual array with
the two rows removed.  Both directions and both inverse radial weights use
restricted residuals Z=M_F Y and the SAME pair diagonal.

The published variance denominator T(T-1) is replaced, explicitly, by
H^2, H=sum(h^2), for noncentered factor designs. The literal
published-denominator statistic is also returned for an audit comparison.
No bootstrap centering, SS-SUM trace estimate or full-sample diagonal is
substituted for the displayed Stat construction.
"""
from __future__ import annotations
import numpy as np
from scipy.special import ndtr


class StatFitError(RuntimeError):
    pass


def full_diagonal(E,tol=1e-7,max_iter=500):
    T,N=E.shape;sq=E*E;D=np.mean(sq,axis=0)
    if np.any(D<=0): raise StatFitError("nonpositive unrestricted residual scale")
    D/=np.mean(D)
    for it in range(1,max_iter+1):
        r2=np.sum(sq/D,axis=1)
        if np.any(r2<=0): raise StatFitError("zero unrestricted residual radius")
        new=N*np.mean(sq/r2[:,None],axis=0)
        equation=np.max(np.abs(new/D-1))
        new/=np.mean(new);change=np.max(np.abs(np.log(new/D)))
        D=new
        if max(equation,change)<=tol:
            return D
    raise StatFitError("unrestricted full-diagonal warm start did not converge")


def pair_diagonals(E,tol=1e-7,max_iter=500,block_size=2048):
    """Exact fixed points computed simultaneously by matrix products.

    The number of retained observations is T-2.  Mean(D_ij)=1 fixes the
    otherwise unidentified common radial scale, identically for every pair.
    Full-sample D is only an initializer; every pair is iterated and checked.
    """
    E=np.asarray(E,float);T,N=E.shape;s,t=np.triu_indices(T,1)
    sq=E*E;initial=full_diagonal(E,tol,max_iter)
    Dall=np.empty((len(s),N));iterations=[];equations=[]
    for start in range(0,len(s),block_size):
        stop=min(start+block_size,len(s));m=stop-start
        Di=np.broadcast_to(initial,(m,N)).copy();rows=np.arange(m)
        left,right=s[start:stop],t[start:stop]
        for iteration in range(1,max_iter+1):
            radius2=(1/Di)@sq.T
            if np.any(radius2<=0) or not np.all(np.isfinite(radius2)):
                raise StatFitError("zero/nonfinite leave-two-out scale radius")
            inverse=1/radius2;inverse[rows,left]=0;inverse[rows,right]=0
            updated=(N/(T-2))*(inverse@sq)
            if np.any(updated<=0) or not np.all(np.isfinite(updated)):
                raise StatFitError("nonpositive leave-two-out diagonal update")
            equation=float(np.max(np.abs(updated/Di-1)))
            updated/=np.mean(updated,axis=1)[:,None]
            change=float(np.max(np.abs(np.log(updated/Di))))
            Di=updated
            if max(equation,change)<=tol:
                # Recheck at the returned, updated iterate.
                radius2=(1/Di)@sq.T
                inverse=1/radius2;inverse[rows,left]=0;inverse[rows,right]=0
                check=(N/(T-2))*(inverse@sq)/Di
                final_equation=float(np.max(np.abs(check-1)))
                if final_equation<=tol:break
        else:
            raise StatFitError("leave-two-out diagonal fixed point did not converge")
        Dall[start:stop]=Di;iterations.append(iteration);equations.append(final_equation)
    return {"D":Dall,"s":s,"t":t,"iterations":max(iterations),
            "equation":max(equations),"minimum_diagonal":float(Dall.min())}


def signal_pvalues(z, pr, settings, cache):
    T,N=z.shape
    left,right=cache["s"],cache["t"];W=1/cache["D"]
    zi=z[left];zj=z[right]
    rleft=np.sum(zi*zi*W,axis=1);rright=np.sum(zj*zj*W,axis=1)
    inner=np.sum(zi*zj*W,axis=1)
    smax=max(x[1] for x in settings)
    if smax:
        Pi=np.cumsum(zi[:,:smax]*W[:,:smax],axis=1)
        Pj=np.cumsum(zj[:,:smax]*W[:,:smax],axis=1)
        Pw=np.cumsum(W[:,:smax],axis=1)
    hs,ht=pr.h[left],pr.h[right]
    H=pr.h2;Dh=H*H
    if Dh<=0:raise RuntimeError("nonpositive design variance denominator")
    p_values=[]
    for design,s,delta in settings:
        if s:
            a=np.sqrt(delta/s)
            ai=Pi[:,s-1]*a;aj=Pj[:,s-1]*a;aa=Pw[:,s-1]*a*a
            ri=rleft+2*hs*ai+hs*hs*aa
            rj=rright+2*ht*aj+ht*ht*aa
            dot=inner+ht*ai+hs*aj+hs*ht*aa
        else:
            ri,rj,dot=rleft,rright,inner
        if np.any(ri<=0)|np.any(rj<=0):raise RuntimeError("zero restricted radius")
        q=dot/(ri*rj);weighted=hs*ht*q
        total=float(np.sum(weighted));square=float(np.sum(weighted*weighted))
        if square<=0:raise RuntimeError("nonpositive pair variance")
        Q=2*N/H*total;variance=4*N*N/Dh*square
        p_values.append(float(ndtr(-Q/np.sqrt(variance))))
    return p_values
