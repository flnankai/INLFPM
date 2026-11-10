"""Reproduce the manuscript panel streams and accumulate nine-method decisions in memory."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse, csv, zipfile, time
from pathlib import Path
from multiprocessing import Pool
from collections import Counter
import numpy as np
from scipy.special import ndtr, stdtr
from panel_math import ols_projected, diagonal_sign_fit, split_trace, sign_bootstrap_bias, sign_q
from pair_scales import pair_diagonals, signal_pvalues, StatFitError
from panel_math import prepare_projection,draw_data,signal_settings,cauchy_p,FitError
from calibration_math import Reference,sidak
DIST=['normal','t3','mixture','IC-t3']
METHOD=['LS-SUM','LS-MAX','LS-CC','SS-SUM','SS-MAX','SS-CC','IN-SUM','IN-MAX','IN-CC']

def baseline(y, pr, seed):
    """Rebuild the five unchanged p-values using the original bootstrap stream."""
    T, N = y.shape
    settings = signal_settings(T, N)
    z = pr.project(y)
    residual = ols_projected(z, pr)["residual"]
    try:
        bias, _ = sign_bootstrap_bias(residual, pr, 100, np.random.default_rng(seed))
    except FitError:
        bias = np.nan
    try:
        in_sum = signal_pvalues(z, pr, settings, pair_diagonals(residual))
    except StatFitError:
        in_sum = np.full(len(settings), np.nan)
    rows = np.full((len(settings), 9), np.nan)
    for j, (_, count, energy) in enumerate(settings):
        shift = np.zeros(N)
        if count:
            shift[:count] = np.sqrt(energy / count)
        shifted = y + shift
        projected = pr.project(shifted)
        ols = ols_projected(projected, pr)
        try:
            sign = diagonal_sign_fit(projected)
            trace = split_trace(shifted, pr, sign["D"])
            ss_sum = ndtr(-(sign_q(sign["U"], pr.h) - bias) / np.sqrt(2 * trace))
        except FitError:
            ss_sum = np.nan
        rows[j, :7] = [seed, *ols["p"], ss_sum, in_sum[j], pr.h2]
    return rows

def evaluate(task):
    T,N,dist,rep,old=task
    dc=DIST.index(dist)+1
    streams=np.random.SeedSequence([20261108,T,N,dc,rep]).spawn(2)
    rng=np.random.default_rng(streams[0])
    seed=int(np.random.default_rng(streams[1]).integers(1,2**32-1))
    archived = old is not None
    if archived and not np.all(old[:,0]==seed): raise AssertionError('archived seed mismatch')
    y,f=draw_data(T,N,dist,rng)
    pr=prepare_projection(f)
    if not archived:
        old = baseline(y, pr, seed)
    ref=Reference(y,pr)
    errors=Counter()
    try:
        infit=ref.inverse()
    except FitError as e:
        errors['IN '+str(e)]+=1
        infit=None
    center=2*np.log(N)-np.log(np.log(N))
    settings=signal_settings(T,N)
    counts=np.zeros((len(settings),9,2),int)
    for j,(_,s,delta) in enumerate(settings):
        row=old[j]
        shift=np.zeros(N)
        if s: shift[:s]=np.sqrt(delta/s)
        try:
            ps,previous=ref.spatial(ref.z+np.outer(pr.h,shift))
            if archived and not np.isclose(previous-center,row[7],atol=3e-7,rtol=3e-7):
                raise AssertionError(('archived SS statistic mismatch',T,N,dist,rep,j,previous-center,row[7]))
        except FitError as e:
            errors['SS '+str(e)]+=1
            ps=np.nan
        pi=np.nan
        if infit is not None:
            a,sc,adj,m=infit
            maximum=np.max(abs((a+shift)*sc))
            if archived and not np.isclose(maximum**2-center,row[8],atol=3e-7,rtol=3e-7):
                raise AssertionError(('archived IN statistic mismatch',T,N,dist,rep,j,maximum**2-center,row[8]))
            pi=sidak(2*stdtr(m,-maximum*adj),N)
        p=np.array([*row[1:5],ps,cauchy_p([row[4],ps],True)[0],row[5],pi,cauchy_p([row[5],pi])[0]])
        counts[j,:,0]=np.isfinite(p)
        counts[j,:,1]=p<.05
    return (T,N,dist),counts,errors

def prepare(archive):
    import pandas as pd
    keys=['T','N','distribution','replicate']
    cols=['bootstrap_seed',*METHOD[:4],'IN-SUM','h2','SS-MAX_stat','IN-MAX_stat']
    with zipfile.ZipFile(archive) as z:
        n=next(n for n in z.namelist() if n.endswith('simulation_raw.csv.gz'))
        with z.open(n) as f:
            df=pd.read_csv(f,compression='gzip',usecols=keys+['design','s','delta']+cols)
    df['design']=df['design'].fillna('null')
    df=df.sort_values(keys+['design','s','delta'])
    payload=df[cols].to_numpy()
    design=list(zip(df.design,df.s,df.delta))
    tasks=[]
    offset=0
    for (T,N,dist,rep), length in df.groupby(keys,sort=False).size().items():
        stop=offset+length
        if design[offset:stop]!=signal_settings(T,N): raise AssertionError('archive design keys')
        # The historical archive predates the unit-variance IC-t3 correction.
        archived_rows = None if dist == "IC-t3" else payload[offset:stop]
        tasks.append((T,N,dist,rep,archived_rows))
        offset=stop
    return tasks

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", required=True,
                        help="Run 1,000 replications per design cell")
    parser.add_argument("--workers", type=int, default=25)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--archive", type=Path,
                        help="Reuse unchanged p-values from the original research ZIP")
    args = parser.parse_args()
    workers = args.workers
    if workers < 1:
        parser.error("--workers must be positive")
    start=time.monotonic()
    if args.archive:
        tasks = prepare(args.archive)
    else:
        tasks = [(T, N, d, rep, None) for rep in range(1, 1001)
                 for T in (60, 120) for N in (100, 200, 400) for d in DIST]
    if len(tasks)!=24000: raise AssertionError('expected 24,000 base panels')
    # Balance panel costs across the worker pool.
    tasks.sort(key=lambda a:(a[3],a[0],a[1],a[2]))
    totals={};reps=Counter();fail=Counter()
    with Pool(workers) as pool:
        for n,(cell,counts,errors) in enumerate(pool.imap_unordered(evaluate,tasks,chunksize=4),1):
            if cell not in totals: totals[cell]=np.zeros_like(counts)
            totals[cell]+=counts
            reps[cell]+=1
            fail.update(errors)
            if n%500==0: print(n,'of',len(tasks),'base panels;',round(time.monotonic()-start),'s',flush=True)
    records=[]
    for (T,N,dist),arr in sorted(totals.items()):
        n=reps[T,N,dist]
        if n!=1000: raise AssertionError('incomplete cell')
        for (design,s,delta),cnt in zip(signal_settings(T,N),arr,strict=True):
            for name,(valid,reject) in zip(METHOD,cnt,strict=True):
                records.append([T,N,dist,design,s,delta,name,n,valid,reject,reject/n])
    if len(records)!=2016: raise AssertionError('incomplete result')
    with args.output.open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f)
        writer.writerow(['T','N','distribution','design','s','delta','method','attempted','valid','rejections','rate'])
        writer.writerows(records)
    print('Complete',len(records),'rows; fit exceptions:',dict(fail),flush=True)
if __name__=='__main__':
    main()


