"""Panel generation and numerical fits used in the manuscript experiment."""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.special import ndtr, ndtri
from scipy.signal import lfilter

class FitError(RuntimeError):
    pass


@dataclass
class Projection:
    factors: np.ndarray
    q: np.ndarray
    h: np.ndarray
    h2: float
    eta: float
    df: int

    def project(self, y):
        return y - self.q @ (self.q.T @ y)


def prepare_projection(factors):
    factors = np.asarray(factors, float)
    if factors.ndim == 1:
        factors = factors[:, None]
    T, p = factors.shape
    if p and np.linalg.matrix_rank(factors) != p:
        raise FitError("rank-deficient factor matrix")
    q = np.linalg.qr(factors, mode="reduced")[0]
    h = np.ones(T) - q @ (q.T @ np.ones(T))
    h2 = float(h @ h)
    if h2 <= 1 or T-p-1 <= 4:
        raise FitError("unidentified intercept or insufficient residual degrees of freedom")
    return Projection(factors, q, h, h2, 1-h2/T, T-p-1)


def gumbel_p(maximum, N):
    centered = float(maximum - 2*np.log(N) + np.log(np.log(N)))
    log_rate = -.5*(centered + np.log(np.pi))
    p = 1.0 if log_rate > 6.6 else float(-np.expm1(-np.exp(log_rate)))
    return p, centered


def cauchy_p(pvalues, truncated=False):
    p = np.asarray(pvalues, float)
    if not np.all(np.isfinite(p)) or np.any((p<0)|(p>1)):
        return np.nan, np.nan
    if np.any(p == 0):
        return 0.0, np.inf
    if not truncated and np.any(p == 1):
        return 1.0, -np.inf
    active = p < .5 if truncated else np.ones(p.size, bool)
    terms = np.zeros(p.size)
    small = (p < 1e-8) & active
    terms[small] = 1/(np.pi*p[small])
    regular = active & ~small
    terms[regular] = np.tan(np.pi*(.5-p[regular]))
    C = float(np.mean(terms))
    if C > 1:
        tail = float(np.arctan(1/C)/np.pi)
    elif C < -1:
        tail = float(1+np.arctan(1/C)/np.pi)
    else:
        tail = float(.5-np.arctan(C)/np.pi)
    return tail, C


def ols_projected(z, pr):
    alpha = pr.h @ z / pr.h2
    residual = z - pr.h[:, None]*alpha
    sse = np.einsum("tn,tn->n", residual, residual)
    if np.any(sse <= 0):
        raise FitError("nonpositive OLS variance")
    sd = np.sqrt(sse/pr.df)
    t = np.sqrt(pr.h2)*alpha/sd
    scaled = residual/np.sqrt(sse)
    corr = scaled.T @ scaled
    N = z.shape[1]
    rho = corr[np.triu_indices(N, 1)]
    threshold = ndtri(1-.1/(2*N))
    retained = np.abs(np.sqrt(pr.df)*rho) > threshold
    a2 = 1+2*np.sum(rho[retained]**2)/N
    center = pr.df/(pr.df-2)
    scale = center*np.sqrt(2*(pr.df-1)/(pr.df-4)*a2)
    S = float(np.sum(t*t-center)/np.sqrt(N)/scale)
    p_sum = float(ndtr(-S))
    p_max, M = gumbel_p(np.max(t*t), N)
    p_cc, C = cauchy_p([p_sum, p_max])
    return {"alpha":alpha, "residual":residual, "sd":sd, "t":t,
            "S":S, "M":M, "C":C, "p":[p_sum,p_max,p_cc]}


def diagonal_sign_fit(z, tol=1e-7, max_iter=500):
    """Restricted, zero-location diagonal sign fit: mean diagonal is one."""
    T, N = z.shape
    sq = z*z
    D = np.mean(sq, axis=0)
    if np.any(D <= 0):
        raise FitError("degenerate initial sign scale")
    D /= np.mean(D)
    for iteration in range(1, max_iter+1):
        radius2 = np.sum(sq/D, axis=1)
        if np.any(radius2 <= 0) or not np.all(np.isfinite(radius2)):
            raise FitError("zero restricted sign radius")
        updated = N*np.mean(sq/radius2[:, None], axis=0)
        moments = updated/D
        updated /= np.mean(updated)
        eq = float(np.max(np.abs(moments-1)))
        update = float(np.max(np.abs(np.log(updated/D))))
        D = updated
        if max(eq,update) <= tol:
            break
    else:
        raise FitError("restricted sign scale did not converge")
    e = z/np.sqrt(D)
    r = np.linalg.norm(e, axis=1)
    U = e/r[:, None]
    final_eq = float(np.max(np.abs(N*np.mean(U*U,axis=0)-1)))
    if final_eq > tol:
        raise FitError("restricted sign equation residual exceeds tolerance")
    return {"D":D,"U":U,"r":r,"iterations":iteration,"equation":final_eq}


def sign_q(U, h):
    wsum = h @ U
    return float(U.shape[1]/(h@h)*(wsum@wsum - np.sum(h*h*np.sum(U*U,axis=1))))


def split_trace(y, pr, D):
    """Exact chronological split leave-two-out, including both ordered pairs.

    The three groups exploit the algebra of an even T. No observations,
    products, or scale fits are omitted. An independent literal implementation
    is used in validation.
    """
    T, N = y.shape
    if T % 2:
        return split_trace_literal(y, pr, D)
    f, h, h2 = pr.factors, pr.h, pr.h2
    s_all,t_all = np.triu_indices(T,1)
    H = T//2
    group_specs = [
        (t_all < H, np.arange(H+1), np.arange(H+1,T), (0,1), ()),
        ((s_all < H)&(t_all >= H), np.arange(H), np.arange(H,T), (0,), (1,)),
        (s_all >= H, np.arange(H-1), np.arange(H-1,T), (), (0,1)),
    ]
    total = 0.0
    invd = 1/np.sqrt(D)
    for mask, rows1, rows2, del1, del2 in group_specs:
        s,t = s_all[mask],t_all[mask]
        if not len(s):
            continue
        ff = (f[s],f[t]); yy = (y[s],y[t])

        def half_residuals(rows, deleted):
            A = np.broadcast_to(f[rows].T@f[rows], (len(s),f.shape[1],f.shape[1])).copy()
            B = f[rows].T@y[rows]
            for k in deleted:
                A -= np.einsum("mi,mj->mij",ff[k],ff[k])
            invA = np.linalg.inv(A)
            out=[]
            for target in (0,1):
                coeff = np.einsum("mi,mij->mj",ff[target],invA)
                e = yy[target] - coeff@B
                for k in deleted:
                    e += np.einsum("mi,mi->m",coeff,ff[k])[:,None]*yy[k]
                e *= invd
                radius = np.linalg.norm(e,axis=1)
                if np.any(radius <= 0):
                    raise FitError("zero leave-two-out radius")
                out.append(e/radius[:,None])
            return out

        u_s1,u_t1 = half_residuals(rows1,del1)
        u_s2,u_t2 = half_residuals(rows2,del2)
        forward = np.einsum("mi,mi->m",u_s1,u_t2)
        reverse = np.einsum("mi,mi->m",u_t1,u_s2)
        total += float(np.sum(h[s]**2*h[t]**2*(forward**2+reverse**2)))
    answer = N*N/(h2*(h2-1))*total
    if not np.isfinite(answer) or answer <= 0:
        raise FitError("nonpositive split trace")
    return float(answer)


def split_trace_literal(y, pr, D):
    T,N = y.shape
    total=0.
    for s in range(T-1):
        for t in range(s+1,T):
            rows=np.delete(np.arange(T),[s,t])
            # Use the same chronological split as the R package for odd T.
            midpoint = len(rows) // 2
            first,second=rows[:midpoint],rows[midpoint:]
            f1,f2=pr.factors[first],pr.factors[second]
            B1=np.linalg.solve(f1.T@f1,f1.T@y[first])
            B2=np.linalg.solve(f2.T@f2,f2.T@y[second])
            Us=[]
            for target,B in ((s,B1),(t,B2),(t,B1),(s,B2)):
                e=(y[target]-pr.factors[target]@B)/np.sqrt(D)
                Us.append(e/np.linalg.norm(e))
            total += pr.h[s]**2*pr.h[t]**2*((Us[0]@Us[1])**2+(Us[2]@Us[3])**2)
    return float(N*N/(pr.h2*(pr.h2-1))*total)


def sign_bootstrap_bias(residual, pr, B, rng, multipliers=None, tol=1e-7,
                        max_iter=500, batch_size=20):
    """Asset-time Rademacher bootstrap, exact batch version of primary Q."""
    T,N=residual.shape
    if multipliers is not None:
        multipliers=np.asarray(multipliers)
        if multipliers.shape != (B,T,N):
            raise ValueError("bootstrap multiplier shape mismatch")
    Q=[]
    for start in range(0,B,batch_size):
        end=min(start+batch_size,B)
        xi=(rng.integers(0,2,size=(end-start,T,N))*2-1
            if multipliers is None else multipliers[start:end])
        z=residual[None,:,:]*xi
        z -= pr.q @ (pr.q.T @ z)
        sq=z*z
        D=np.mean(sq,axis=1)
        D/=np.mean(D,axis=1)[:,None]
        for iteration in range(max_iter):
            radius2=np.sum(sq/D[:,None,:],axis=2)
            if np.any(radius2<=0):
                raise FitError("zero bootstrap radius")
            upd=N*np.mean(sq/radius2[:,:,None],axis=1)
            moments=upd/D
            upd/=np.mean(upd,axis=1)[:,None]
            eq=np.max(np.abs(moments-1))
            change=np.max(np.abs(np.log(upd/D)))
            D=upd
            if max(eq,change)<=tol:
                break
        else:
            raise FitError("bootstrap sign scales did not converge")
        e=z/np.sqrt(D[:,None,:])
        U=e/np.linalg.norm(e,axis=2)[:,:,None]
        sums=np.einsum("t,btn->bn",pr.h,U)
        off=np.sum(sums*sums,axis=1)-np.einsum("t,bt->b",pr.h**2,np.sum(U*U,axis=2))
        Q.extend(N/pr.h2*off)
    return float(np.mean(Q)),np.asarray(Q)


def spatial_median_fit(z, tol=1e-7, max_iter=500):
    T,N=z.shape
    theta=np.mean(z,axis=0)
    d=np.std(z,axis=0,ddof=1)
    if np.any(d<=0):
        raise FitError("nonpositive initial spatial median scale")
    d/=np.exp(np.mean(np.log(d)))
    for iteration in range(1,max_iter+1):
        e=(z-theta)/d
        r=np.linalg.norm(e,axis=1)
        if np.min(r)<=100*np.finfo(float).eps*np.median(r):
            raise FitError("collapsed spatial median radius")
        U=e/r[:,None]
        ir=1/r
        theta_new=np.sum(z*(ir/np.sum(ir))[:,None],axis=0)
        d_new=d*np.sqrt(N*np.mean(U*U,axis=0))
        d_new/=np.exp(np.mean(np.log(d_new)))
        change=max(np.max(np.abs((theta_new-theta)/d))/(np.median(r)/np.sqrt(N)),
                   np.max(np.abs(np.log(d_new/d))))
        theta,d=theta_new,d_new
        if change<tol:
            e=(z-theta)/d;r=np.linalg.norm(e,axis=1);U=e/r[:,None]
            eq1=np.sqrt(N)*np.max(np.abs(np.mean(U,axis=0)))
            eq2=np.max(np.abs(N*np.mean(U*U,axis=0)-1))
            if max(eq1,eq2)<tol:
                break
    else:
        raise FitError("spatial median did not converge")
    return {"theta":theta,"d":d,"r":r,"iterations":iteration,
            "equation":float(max(eq1,eq2))}


def pqr_scale(residual):
    T,N=residual.shape
    s,t=np.triu_indices(T,1)
    order=int(np.ceil(.5*len(s)))-1
    dist=np.abs(residual[t].T-residual[s].T)
    dist.partition(order,axis=1)
    q=dist[:,order].copy()
    if np.any(q<=0):
        raise FitError("nonpositive PQR scale")
    return q/np.sqrt(np.mean(q*q))


def in_fit(z,pr,d,tol=1e-8,max_iter=1000):
    T,N=z.shape
    h=pr.h
    alpha=h@z/pr.h2
    for iteration in range(1,max_iter+1):
        e=(z-h[:,None]*alpha)/d
        r=np.linalg.norm(e,axis=1)
        if np.min(r)<=100*np.finfo(float).eps*np.median(r):
            raise FitError("collapsed IN radius")
        ir2=1/(r*r)
        updated=np.sum(z*(h*ir2)[:,None],axis=0)/np.sum(h*h*ir2)
        change=np.max(np.abs((updated-alpha)/d))/(np.median(r)/np.sqrt(N))
        alpha=updated
        if change<tol:
            e=(z-h[:,None]*alpha)/d;r=np.linalg.norm(e,axis=1);ir2=1/r**2
            eq=np.sqrt(N)*np.median(r)*np.max(np.abs(np.mean((h*ir2)[:,None]*e,axis=0)))
            if eq<tol:
                break
    else:
        raise FitError("IN root did not converge")
    A=float(np.mean(h*h*ir2))
    pi=pr.q@(pr.q.T@(h*ir2))
    Bc=A+float(np.mean(r*r*pi*pi))
    zscore=np.sqrt(pr.df*N)*A*alpha/(d*np.sqrt(Bc))
    p,M=gumbel_p(np.max(zscore*zscore),N)
    return {"alpha":alpha,"d":d,"r":r,"A":A,"Bc":Bc,"z":zscore,
            "p":p,"M":M,"iterations":iteration,"equation":float(eq),
            "min_radius_ratio":float(np.min(r)/np.median(r))}


@lru_cache(maxsize=8)
def symmetric_sqrt(N,rho=.5):
    Sigma=rho**np.abs(np.arange(N)[:,None]-np.arange(N)[None,:])
    eigenval,eigenvec=np.linalg.eigh(Sigma)
    return (eigenvec*np.sqrt(eigenval))@eigenvec.T


def draw_factors(T,rng):
    c=np.array([.53,.19,.19]);a=np.array([.06,.19,.05])
    v0=np.array([.89,.62,.80]);b=np.array([.85,.74,.76]);g=np.array([.11,.19,.15])
    fprev=np.zeros(3);vprev=np.ones(3);zprev=np.zeros(3)
    out=[]
    for _ in range(T+50):
        v=v0+b*vprev+g*zprev*zprev
        z=rng.standard_normal(3)
        f=c+a*fprev+np.sqrt(v)*z
        out.append(f)
        fprev,vprev,zprev=f,v,z
    return np.asarray(out[-T:])


def draw_data(T,N,distribution,rng):
    f=draw_factors(T,rng)
    if distribution=="IC-t3":
        eps=rng.standard_t(3,size=(T,N))@symmetric_sqrt(N).T
    else:
        z=rng.standard_normal((T,N));innov=z*np.sqrt(.75);innov[:,0]=z[:,0]
        eps=lfilter([1.],[1.,-.5],innov,axis=1)
        if distribution=="t3":
            eps/=np.sqrt(rng.chisquare(3,size=T))[:,None]
        elif distribution=="mixture":
            eps*=np.where(rng.uniform(size=T)<.8,1.,3.)[:,None]/np.sqrt(2.6)
        elif distribution!="normal":
            raise ValueError("unknown error distribution")
    # Every tested statistic is invariant to adding an observed-factor fit.
    # The OBES source does not specify the loading values; these fixed values
    # are stated explicitly rather than attributed to the source paper.
    beta=np.ones((N,3))
    y=eps+f@beta.T
    return y,f


def signal_settings(T,N):
    settings=[("null",0,0.)]
    if T==60 and N in (100,200):
        settings += [("sparsity",s,.5) for s in range(1,47,5)]
    if T==60 and N==200:
        settings += [("strength",s,float(k/10)) for s in (2,10,25) for k in range(1,11)]
    return settings


