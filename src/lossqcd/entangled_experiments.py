"""Reproduce the full Figure 2 CRE family; benchmark Appendix B.4 conditionals.

The initial conditional-detector comparison is calibrated by DIRECT Monte
Carlo at moderate ARL. It is not a reproduction of Figure 7, and it does not
reuse the iid scalar two-million-pulse thresholds for correlated increments.
"""
import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .models import coherent_kl,squeezed_kl
from .entangled import entangled_pair,entangled_kl_per_pulse,conditional_llrs
from .conditional_detection import calibrate_conditional,simulate_conditional
from .calibration import independent_ratio
from .simulation import summarize_arl,summarize_delay

PAPER='https://arxiv.org/html/2503.12276v2'


def figure2(root):
    nas=np.unique(np.r_[0,np.geomspace(1e-7,.1,100),np.linspace(.1,5,491)])
    sizes=[2,4,8,16,32,64,128,256]
    curves={'coherent':np.array([coherent_kl(100,a,.9,.85) for a in nas]),
            'squeezed_n1':np.array([squeezed_kl(100,a,.9,.85) for a in nas])}
    for n in sizes:
        curves[f'entangled_n{n}']=np.array([entangled_kl_per_pulse(100,a,.9,.85,n) for a in nas])
    np.savetxt(root/'figure2_all_curves.csv',np.column_stack([nas,*curves.values()]),delimiter=',',
        header=','.join(['Na',*curves.keys()]),comments='')
    fig,ax=plt.subplots(figsize=(10,5.5))
    ax.plot(nas,curves['coherent'],'r--',label='Coherent (Eq. 2)')
    ax.plot(nas,curves['squeezed_n1'],'r-',label='Independent squeezing (n=1)')
    for k,n in enumerate(sizes):
        ax.plot(nas,curves[f'entangled_n{n}'],ls='-.',color=plt.cm.Blues(.35+.6*k/7),label=f'Entangled n={n}')
    ax.set(xlabel='Augmentation energy Na (photons/pulse)',ylabel='Post/pre relative entropy (nats/pulse)',
        xlim=(0,5),ylim=(0,None),title='Figure 2 analytical reproduction: N=100, eta1=0.9, eta2=0.85')
    ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=9)
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(root/'figure2_entangled.png',dpi=180)
    fig.savefig(root/'figure2_entangled.svg');plt.close(fig)
    return {'block_sizes':sizes,'formula':'KL_block/n; seed squeezing asinh(sqrt(n*Na))',
        'Na0':[float(values[0]) for values in curves.values()],
        'Na5':{name:float(values[-1]) for name,values in curves.items()},
        'Na0_1':{'coherent':coherent_kl(100,.1,.9,.85),'squeezed':squeezed_kl(100,.1,.9,.85),
                 **{f'entangled_n{n}':entangled_kl_per_pulse(100,.1,.9,.85,n) for n in sizes}}}


def model_check(root):
    seed=2026102101
    pre,post=entangled_pair(n=8,Na=.1)
    x=post.sample(np.random.default_rng(seed),100000)
    covariance=np.cov(x,rowvar=False)
    subset=x[:200]
    def dense_logpdf(state):
        delta=subset-state.mean
        L=np.linalg.cholesky(state.covariance)
        q=np.sum(np.linalg.solve(L,delta.T)**2,axis=0)
        return -.5*(state.n*np.log(2*np.pi)+2*np.log(np.diag(L)).sum()+q)
    dense=dense_logpdf(post)-dense_logpdf(pre)
    llr=post.llr(subset,pre)
    stats={'n':8,'Na':.1,'seed':seed,'samples':len(x),
        'mean_max_abs_error':float(np.max(np.abs(x.mean(axis=0)-post.mean))),
        'covariance_max_abs_error':float(np.max(np.abs(covariance-post.covariance))),
        'dense_vs_reduced_llr_max_error':float(np.max(np.abs(dense-llr))),
        'conditional_chain_vs_joint_llr_max_error':float(np.max(np.abs(conditional_llrs(subset,pre,post).sum(axis=-1)-llr)))}
    if stats['dense_vs_reduced_llr_max_error']>1e-9 or stats['conditional_chain_vs_joint_llr_max_error']>1e-9:
        raise RuntimeError('independent dense/conditional likelihood check failed')
    np.savez_compressed(root/'covariance_check.npz',samples=x,expected=post.covariance,empirical=covariance)
    fig,axes=plt.subplots(1,2,figsize=(9,4),layout='constrained')
    for ax,values,title in zip(axes,[post.covariance,covariance],['Model covariance','Empirical covariance: 100,000 blocks']):
        from matplotlib.colors import TwoSlopeNorm
        im=ax.imshow(values,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.03,vcenter=0,vmax=.25))
        ax.set(title=title,xlabel='Pulse within block',ylabel='Pulse within block')
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.8,label='Covariance',ticks=[-.03,0,.05,.1,.15,.2,.25])
    fig.savefig(root/'entangled_covariance.png',dpi=180);plt.close(fig)
    return stats


def write_report(root,report):
    rows=[]
    baseline=report['models'][0]
    for model in report['models']:
        for nu in report['settings']['change_at']:
            d=model['delay'][str(nu)]
            b=baseline['delay'][str(nu)]
            ratio=({'ratio':1.,'se':0.,'approx_95_interval':[1.,1.]} if model is baseline else
                independent_ratio(b['conditional_delay'],b['conditional_delay_se'],
                                   d['conditional_delay'],d['conditional_delay_se']))
            rows.append({'label':model['label'],'n':model['n'],'change_at':nu,
                'phase_of_change':(nu-1)%model['n']+1,
                'threshold':model['calibration']['threshold'],'validation_arl':model['validation']['arl'],
                'validation_arl_se':model['validation']['arl_se'],
                'arl_matched_to_coherent':model['arl_match']['passed'],
                'conditional_delay':d['conditional_delay'],'conditional_delay_se':d['conditional_delay_se'],
                'false_alarm_probability':d['false_alarm_probability'],'survivors':d['survivors'],
                'coherent_delay_ratio':ratio['ratio'],'ratio_se':ratio['se'],
                'ratio_95_low':ratio['approx_95_interval'][0],'ratio_95_high':ratio['approx_95_interval'][1]})
    report['summary']=rows
    with (root/'conditional_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    labels=[m['label'] for m in report['models']]
    index=np.arange(len(labels))
    axes[0].errorbar(index,[m['validation']['arl'] for m in report['models']],
        yerr=[1.96*m['validation']['arl_se'] for m in report['models']],fmt='o',capsize=4)
    axes[0].axhline(report['settings']['target_arl'],color='black',ls='--')
    axes[0].set(xticks=index,xticklabels=labels,ylabel='No-change ARL (pulses)',title='Fresh direct Monte Carlo validation')
    for nu in report['settings']['change_at']:
        data=[row for row in rows if row['change_at']==nu]
        axes[1].errorbar(index,[x['conditional_delay'] for x in data],
            yerr=[1.96*x['conditional_delay_se'] for x in data],fmt='o-',capsize=3,label=f'Change at {nu}')
    axes[1].set(xticks=index,xticklabels=labels,ylabel='Conditional delay (pulses)',
        title='Appendix B.4 pulse-by-pulse conditional detector')
    axes[1].legend(fontsize=9)
    for ax in axes:
        ax.tick_params(axis='x',rotation=25);ax.grid(axis='y',alpha=.2)
    fig.suptitle(f"Initial conditional benchmark: N=100, Na={report['settings']['Na']:g}; bars: approximate 95% MC intervals")
    fig.tight_layout();fig.savefig(root/'conditional_detection.png',dpi=180)
    fig.savefig(root/'conditional_detection.svg');plt.close(fig)
    matched=sum(m['arl_match']['passed'] for m in report['models'])
    lines=['# Entangled probe model and conditional CUSUM implementation','',
        f'Source: [paper Eqs. (4)–(6), Appendix C.1–C.2 and Appendix B.4]({PAPER}#A2.SS4).','',
        '## Reproduced and implemented','',
        '- Full Figure 2 analytical curve family: ordinary coherent light, independent squeezing, and entangled blocks n=2,4,8,...,256, for N=100, eta1=0.9, eta2=0.85 and Na in [0,5]. The seed has n*Na photons. This is a formula-based reproduction, not a pixel/digitized-error comparison.',
        '- Correlated full-block Gaussian sampling, covariance, joint log density, log-likelihood ratio and per-pulse KL. Dense multivariate calculations check the reduction to the collective statistic Y=sum(X)/sqrt(n). n=1 reduces to independent squeezing; zero squeezing gives the ordinary coherent information per pulse.',
        '- The paper\'s Appendix B.4 pulse-by-pulse conditional-LLR detector. Conditioning history resets at block boundaries; Page evidence persists across blocks and reflects at zero after each pulse. A reflection does not erase the measured prefix. This variant is not claimed to satisfy iid-CUSUM optimality.',
        '- Within-block physical changes use the exact mixed-loss covariance, not independent pre/post samples within an entangled block. Full-block decisions are also provided as a separate reference implementation; the benchmark below uses pulse-by-pulse decisions.','',
        '## Initial conditional-detector benchmark','',
        f"Na={report['settings']['Na']:g}; target ARL={report['settings']['target_arl']:g} PULSES. Calibration uses fresh direct no-change simulations for each detector, followed by independent holdout ARL validation and independent changed-run ensembles. **{matched}/{len(report['models'])} detector settings passed** the predeclared {100*report['settings']['match_tolerance']:g}% approximate ARL equivalence criterion against coherent light.",'',
        'The entire pointwise approximate 95% ARL-ratio interval must fit within the equivalence band, and both point estimates must be close to target. Thresholds are never retuned on validation. Any unmatched detector remains visible and its comparison is exploratory only. No calibrated iid scalar two-million-pulse threshold is reused here: conditional increments are correlated, so those ARLs do not transfer.','',
        '| Detector | Threshold | Validated ARL ± SE (pulses) | Match passed |',
        '|---|---:|---:|---|']
    for m in report['models']:
        lines.append(f"| {m['label']} | {m['calibration']['threshold']:.6f} | {m['validation']['arl']:.2f} ± {m['validation']['arl_se']:.2f} | {m['arl_match']['passed']} |")
    for nu in report['settings']['change_at']:
        lines+=['',f'### Change on pulse {nu}','',
            '| Detector | Phase in block | Delay ± SE (pulses) | Early alarms | Surviving runs | Coherent / detector delay |',
            '|---|---:|---:|---:|---:|---:|']
        for row in [x for x in rows if x['change_at']==nu]:
            lines.append(f"| {row['label']} | {row['phase_of_change']} | {row['conditional_delay']:.3f} ± {row['conditional_delay_se']:.3f} | {100*row['false_alarm_probability']:.2f}% | {row['survivors']} | {row['coherent_delay_ratio']:.3f} |")
    lines+=['','## Derivation and time conventions','',
        'K_eta=I/4+(v_eta-1/4)uu^T, where u=ones(n)/sqrt(n) and v_eta=[1-eta+eta*exp(-2s)]/4. Orthogonal modes have the same variance 1/4 under both hypotheses, so they cancel from the homogeneous joint LLR. The per-pulse CRE is [(v2/v1)-1-log(v2/v1)]/(2n)+N*(sqrt(eta2)-sqrt(eta1))^2/(2v1).','',
        'For k previously measured pulses, c=(v_eta-1/4)/n gives conditional mean mu+c/(1/4+k*c)*(sum(previous X)-k*mu) and conditional variance (1/4)*(1/4+(k+1)c)/(1/4+k*c). These match the general Gaussian conditioning equations. Unreflected conditional LLRs sum to the full block LLR; the per-pulse Page reflections make the conditional detector distinct from a complete-block detector.','',
        'If transmissivities differ within a block, the covariance is I/4+q*w*w^T, where q=[exp(-2s)-1]/(4n) and w_i=sqrt(eta_i). Sequential generation uses Gaussian conditioning for that actual covariance. The detector score compares homogeneous pre/post conditional hypotheses; it is not a maximization over all within-block change hypotheses.','',
        'Pulses are one based: block 1 is pulses 1..n, change nu makes pulse nu the first changed pulse, and delay is T-nu. The table states the phase explicitly. Indexing and ARL settings are therefore explicit rather than inferred from the caption of the published Figure 7. Delay conditions on T>=nu; early alarms are reported, not counted as successful detections. nu=1 checks initialization without early alarms. A partially observed final block is allowed for the conditional detector.','',
        '## Evidence and limitations','',
        f"Calibration uses {report['settings']['calibration_trials']:,} runs/evaluation, holdout validation {report['settings']['validation_trials']:,}, and each delay ensemble {report['settings']['delay_trials']:,}. The horizon is {report['settings']['horizon']:,} pulses. No calibration/validation/surviving changed runs were censored; unavailable unrestricted estimates cause refusal rather than silent omission.",'',
        'The covariance_check.npz contains 100,000 independently sampled n=8 blocks, expected and empirical covariances. report.json records moment errors and dense/conditional likelihood checks. Separate model folders contain calibration histories, seed streams and raw validation/delay stopping times with event flags. CSV files retain plot coordinates and summaries. All uncertainty is pointwise approximate Monte Carlo uncertainty, conditional on chosen thresholds, not simultaneous coverage or optical-model uncertainty.','',
        'This completes the entangled Gaussian model, full Figure 2 analytical family, and an initial Appendix B.4 implementation. It is not a reproduction of Figure 7 at all of its settings, nor an independently calibrated two-million-pulse conditional-detector study. The earlier iid two-million-pulse results remain valid for their original detectors. Next: calibrate and independently validate the correlated conditional detector at the larger ARL, then sweep change phase and Na.','']
    (root/'RESULTS.md').write_text('\n'.join(lines))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('results_entangled'))
    p.add_argument('--seed',type=int,default=20261021)
    p.add_argument('--Na',type=float,default=.1)
    p.add_argument('--target-arl',type=float,default=1000.)
    p.add_argument('--block-sizes',type=int,nargs='+',default=[2,4,8,16])
    p.add_argument('--change-at',type=int,nargs='+',default=[1,101,1005])
    p.add_argument('--calibration-trials',type=int,default=20000)
    p.add_argument('--validation-trials',type=int,default=50000)
    p.add_argument('--delay-trials',type=int,default=30000)
    p.add_argument('--horizon',type=int,default=100000)
    p.add_argument('--match-tolerance',type=float,default=.05)
    args=p.parse_args(argv)
    if not np.isfinite(args.Na) or args.Na<0 or not np.isfinite(args.target_arl) or args.target_arl<=1:
        p.error('Na must be finite/nonnegative and target ARL finite/>1')
    if any(n<2 for n in args.block_sizes) or len(set(args.block_sizes))!=len(args.block_sizes):
        p.error('benchmark block sizes must be unique integers >=2')
    if min(args.calibration_trials,args.validation_trials,args.delay_trials)<2:
        p.error('ensemble sizes must be >=2')
    if args.horizon<1 or any(not 1<=nu<=args.horizon for nu in args.change_at) or len(set(args.change_at))!=len(args.change_at):
        p.error('change points must be unique pulses in the observation window')
    if not np.isfinite(args.match_tolerance) or not 0<args.match_tolerance<1:
        p.error('match tolerance must lie in (0,1)')
    args.output.mkdir(parents=True,exist_ok=False)
    report={'settings':vars(args)|{'output':str(args.output)},'paper':PAPER,
        'figure2':figure2(args.output),'model_check':model_check(args.output),'models':[],
        'python':sys.version,'platform':platform.platform(),
        'versions':{k:importlib.metadata.version(k) for k in ['numpy','matplotlib']},
        'conditional_two_million_arl_validated':False}
    choices=[('coherent',1,'Coherent'),('entangled',1,'Squeezed n=1')]+[
        ('entangled',n,f'Entangled n={n}') for n in args.block_sizes]
    streams=np.random.SeedSequence(args.seed).spawn(len(choices)*(2+len(args.change_at)))
    for i,(kind,n,label) in enumerate(choices):
        count=2+len(args.change_at)
        seeds=[int(s.generate_state(1,dtype=np.uint64)[0]) for s in streams[count*i:count*(i+1)]]
        print(f'{label}: calibrating conditional detector in pulses ...',flush=True)
        calibration=calibrate_conditional(kind,n,args.target_arl,args.calibration_trials,args.horizon,seeds[0],Na=args.Na)
        h=calibration['threshold']
        val=simulate_conditional(kind,n,Na=args.Na,threshold=h,trials=args.validation_trials,horizon=args.horizon,seed=seeds[1])
        stats=summarize_arl(val)
        if stats['arl'] is None:
            raise RuntimeError('censored holdout ARL; increase horizon')
        raw={'validation_stop':val.stop,'validation_event':val.event}
        delay={}
        for j,nu in enumerate(args.change_at):
            runs=simulate_conditional(kind,n,Na=args.Na,threshold=h,trials=args.delay_trials,
                horizon=args.horizon,change_at=nu,seed=seeds[2+j])
            delay[str(nu)]=summarize_delay(runs)
            if delay[str(nu)]['conditional_delay'] is None:
                raise RuntimeError('censored/no-survivor delay; unrestricted estimate unavailable')
            raw[f'nu{nu}_stop'],raw[f'nu{nu}_event']=runs.stop,runs.event
        if i==0:
            match={'passed':abs(stats['arl']/args.target_arl-1)<=args.match_tolerance,'ratio':1.}
        else:
            base=report['models'][0]['validation']
            ratio=independent_ratio(base['arl'],base['arl_se'],stats['arl'],stats['arl_se'])
            low,high=ratio['approx_95_interval'];tol=args.match_tolerance
            passed=(1-tol<=low and high<=1+tol and abs(stats['arl']/args.target_arl-1)<=tol
                and abs(base['arl']/args.target_arl-1)<=tol)
            match=ratio|{'passed':passed}
        model={'label':label,'kind':kind,'n':n,'calibration':calibration,'validation':stats,
               'delay':delay,'stream_seeds':seeds,'arl_match':match}
        report['models'].append(model)
        folder=args.output/f'model_{i:02d}';folder.mkdir()
        (folder/'report.json').write_text(json.dumps(model,indent=2)+'\n')
        np.savez_compressed(folder/'raw_runs.npz',**raw)
        print(f"  validated ARL={stats['arl']:.1f} +/- {stats['arl_se']:.1f}; matched={match['passed']}",flush=True)
    write_report(args.output,report)
    print('Full Figure 2 curves and initial conditional detector benchmark complete.',flush=True)


if __name__=='__main__':
    main()
