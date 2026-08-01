# -----------------------------------------------------------------------------
# make_fig4_panels.py — render the expression panels (marginal-vs-joint bar,
# forest, RNA dose-response, and WHITE-only volcano + dose-response).
# Author: Jonathan Bohlen
# (c) 2026 Jonathan Bohlen. Code accompanying Sengl et al. (2026),
# "Surveillance and selection of 5S ribosomal RNA genes in the human genome."
# Released under the MIT License; see LICENSE at the repository root.
# -----------------------------------------------------------------------------
"""Render the expression panels at 4x4in, 8pt, 1pt, 300 DPI.
Covers: marginal-vs-joint bar, forest, RNA dose-response, and WHITE-only volcano + dose-response."""
import os, pandas as pd, numpy as np, re, de_figstyle as S, matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
S.setup()
HD=f"{Path(os.environ.get('FIVES_DATA','data'))}/de_v2_handoff/meta"; FIG=str(Path(os.environ.get('FIVES_OUT','output'))/"Figure4_new"); E=f"{Path(os.environ.get('FIVES_DATA','data'))}/results/eqtl/extreme"
sym=pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/de_v2/ensg2symbol.tsv",sep="\t").set_index("ensg").symbol.to_dict()
mod=set(pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/de_v2/translation_dosage_module.tsv",sep="\t").symbol)
BLUE,RED="#4a90d9","#d62728"

def modz(scope):
    d=pd.read_csv(f"{HD}/meta_{scope}_cont.tsv",sep="\t"); d["symbol"]=d.ensg.map(sym)
    z=d[d.symbol.isin(mod)].z.clip(-10,10); return z.mean(), z.std()/np.sqrt(len(z))

# ---------- 1. marginal-vs-joint bar ----------
rows=[("DNA\nmarginal","cont_mut",BLUE),("RNA\nmarginal","cont_rna",RED),
      ("DNA\njoint","joint_mut",BLUE),("RNA\njoint","joint_rna",RED)]
vals=[(lab,*modz(sc),c) for lab,sc,c in rows]
fig,ax=plt.subplots(figsize=(S.PANEL,S.PANEL))
x=np.arange(4); alphas=[1,1,.55,.55]
for i,v in enumerate(vals):
    ax.bar(i,v[1],yerr=1.96*v[2],color=v[3],alpha=alphas[i],capsize=3,error_kw=dict(lw=1))
ax.axhline(0,color="k",lw=1); ax.axvline(1.5,color="grey",ls=":",lw=1)
ax.set_xticks(x); ax.set_xticklabels([v[0] for v in vals],fontsize=7.5)
ax.set_ylabel("translation module effect (mean z)",fontsize=8)
ax.set_title("Translation module effect: marginal vs joint dosage",fontsize=8)
for i,v in enumerate(vals): ax.text(i,v[1]-0.15,f"{v[1]:.2f}",ha="center",va="top",fontsize=7)
ax.tick_params(labelsize=8); fig.tight_layout(); S.save(fig,f"{FIG}/Fig4IJ_marginal_vs_joint_bar.pdf")

# ---------- 2. forest ----------
items=[("RNA marginal (cont_rna)","cont_rna",RED),("RNA joint | DNA (joint_rna)","joint_rna",RED),
       ("DNA marginal (cont_mut)","cont_mut",BLUE),("DNA joint | RNA (joint_mut)","joint_mut",BLUE)]
fig,ax=plt.subplots(figsize=(S.PANEL,S.PANEL))
for i,(lab,sc,c) in enumerate(items):
    b,se=modz(sc); y=len(items)-1-i
    ax.errorbar(b,y,xerr=1.96*se,fmt="o",color=c,ms=6,capsize=3,lw=1.2)
    ax.text(b,y+0.18,f"{b:.2f}",ha="center",fontsize=7,color=c)
ax.axvline(0,color="k",ls=":",lw=1)
ax.set_yticks(range(len(items))); ax.set_yticklabels([l for l,_,_ in items][::-1],fontsize=7)
ax.set_xlabel("translation module effect (mean z, 95% CI)",fontsize=8)
ax.set_title("Translation module: marginal vs joint",fontsize=8); ax.tick_params(labelsize=8)
ax.set_xlim(-4.2,0.6); fig.tight_layout(); S.save(fig,f"{FIG}/Fig4IJ_3panel_forest.pdf")

# ---------- 3+4. dose-response helper (audit_donor_resid clean score) ----------
DG=pd.read_pickle(f"{E}/audit_donor_resid.pkl")
syms=[s for s in mod if s in DG.columns]; score=DG[syms].mean(1)*100
dm=pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/de_v2/donor_metrics.tsv",sep="\t").set_index("donor")
sp=pd.read_csv(f"{Path(os.environ.get('FIVES_DATA','data'))}/eqtl_inputs/SubjectPhenotypesDS.txt",sep="\t").set_index("SUBJID")
white=set(sp.index[sp.RACE==3])
def resid_on(y,x): b=np.polyfit(x,y,1); return y-(b[0]*x+b[1])
def build(donors):
    D=pd.DataFrame({"score":score}).join((dm.rna_excess*100).rename("RNA")).join((dm.sum_wgs*100).rename("DNA")).dropna()
    if donors is not None: D=D[D.index.isin(donors)]
    D["score_ctrlDNA"]=resid_on(D.score.values,D.DNA.values); D["score_ctrlRNA"]=resid_on(D.score.values,D.RNA.values)
    return D
def bp(ax,D,xc,yc,color,title,nb=6):
    d=D.sort_values(xc).reset_index(); n=len(d); bins=[np.arange(int(i*n/nb),int((i+1)*n/nb)) for i in range(nb)]
    x=d[xc]; y=d[yc]
    ax.errorbar([x[b].mean() for b in bins],[y[b].mean() for b in bins],yerr=[y[b].sem() for b in bins],fmt="o-",color=color,ms=5,lw=1,capsize=2,elinewidth=1)
    r,p=stats.spearmanr(d[xc],d[yc]); ax.axhline(0,color="grey",ls=":",lw=1)
    ax.set_title(title,fontsize=8); ax.text(0.95,0.92,f"ρ={r:+.2f}\np={p:.0e}",transform=ax.transAxes,ha="right",va="top",fontsize=8,color=color); ax.tick_params(labelsize=8)

# 3. standalone RNA dose-response (all), 4x4
Dall=build(None)
fig,ax=plt.subplots(figsize=(S.PANEL,S.PANEL)); bp(ax,Dall,"RNA","score",RED,"RNA dose-response (all)")
ax.set_ylabel("translation program (%)",fontsize=8); ax.set_xlabel("aggregate gene-region RNA-VAF (%)",fontsize=8)
fig.tight_layout(); S.save(fig,f"{FIG}/Fig4J_RNA_doseresponse.pdf")

# 4. WHITE-only 2x2 dose-response, shared Y
Dw=build(white); print(f"white dose-response n={len(Dw)}")
fig,axes=plt.subplots(2,2,figsize=(2*S.PANEL,2*S.PANEL),sharey=True)
bp(axes[0,0],Dw,"DNA","score",BLUE,"DNA — single-variable"); bp(axes[0,1],Dw,"RNA","score",RED,"RNA — single-variable")
bp(axes[1,0],Dw,"DNA","score_ctrlRNA",BLUE,"DNA — joint (ctrl RNA)"); bp(axes[1,1],Dw,"RNA","score_ctrlDNA",RED,"RNA — joint (ctrl DNA)")
for a in axes[:,0]: a.set_ylabel("translation program (%)",fontsize=8)
axes[1,0].set_xlabel("gene-region DNA-VAF (%)",fontsize=8); axes[1,1].set_xlabel("gene-region RNA-VAF (%)",fontsize=8)
axes[0,0].set_xlabel("gene-region DNA-VAF (%)",fontsize=8); axes[0,1].set_xlabel("gene-region RNA-VAF (%)",fontsize=8)
fig.suptitle("Dose-response — WHITE only (RACE==3)",y=1.01,fontsize=8)
fig.tight_layout(); S.save(fig,f"{FIG}/FigS_Fig4IJ_doseresponse_WHITE.pdf")

# 5. gene-level trans-effect volcanoes (RNA, DNA, WHITE) — capped at 12, saturated dots shown at 12
cyto={e for e,s in sym.items() if re.match(r"^RP[LS]\d+[A-Z]?$",str(s)) or s in {"RPLP0","RPLP1","RPLP2","RPSA","FAU","UBA52"}}
mito={e for e,s in sym.items() if re.match(r"^MRP[LS]\d",str(s))}
CAP=12
def gene_volcano(scope,title,outfile):
    d=pd.read_csv(f"{HD}/meta_{scope}_cont.tsv",sep="\t"); d=d[d.n_tissue>=20].copy()
    raw=-np.log10(d.padj.clip(lower=1e-300)); d["nlp"]=np.minimum(raw,CAP); d["sat"]=raw>CAP; d["b"]=d.beta_meta*100
    fig,ax=plt.subplots(figsize=(S.PANEL,S.PANEL))
    ax.scatter(d.b,d.nlp,s=6,color="#ddd",alpha=.5,linewidths=0)
    mt=d[d.ensg.isin(mito)]; ax.scatter(mt.b,mt.nlp,s=28,color=BLUE,linewidths=0,label="mito-RP")
    rp=d[d.ensg.isin(cyto)]; ax.scatter(rp.b,rp.nlp,s=28,color=RED,linewidths=0,label="cyto-RP")
    # saturated dots: draw as upward triangles pinned at the cap so they read as "off the top"
    st=d[d.sat]; ax.scatter(st.b,[CAP]*len(st),s=20,color="#bbb",marker="^",linewidths=0,zorder=1)
    for e,c in [(mito,BLUE),(cyto,RED)]:
        s2=st[st.ensg.isin(e)]; ax.scatter(s2.b,[CAP]*len(s2),s=36,color=c,marker="^",linewidths=0,zorder=4)
    n=d[d.ensg=="ENSG00000181163"]
    if len(n): yy=min(n.nlp.iloc[0],CAP); ax.scatter(n.b,[yy],s=80,color="#ff8c00",edgecolor="k",lw=.6,zorder=5); ax.annotate("NPM1",(n.b.iloc[0],yy),fontsize=8,fontweight="bold",xytext=(6,-2),textcoords="offset points")
    ax.axvline(0,color="k",ls=":",lw=1); ax.axhline(-np.log10(0.05),color="grey",ls="--",lw=1)
    xl=np.nanpercentile(np.abs(d.b),99.7); ax.set_xlim(-xl,xl); ax.set_ylim(-0.5,CAP+0.6)
    ax.set_xlabel("Δ covariate-adj. expr. (%)",fontsize=8); ax.set_ylabel(f"−log10 FDR (capped {CAP}; ▲=saturated)",fontsize=8)
    ax.set_title(title,fontsize=8); ax.legend(frameon=False,fontsize=8,loc="upper left"); ax.tick_params(labelsize=8)
    fig.tight_layout(); S.save(fig,f"{FIG}/{outfile}")
gene_volcano("cont_rna","RNA-dosage trans-effect","Fig4H_volcano_RNAdose.pdf")
gene_volcano("cont_mut","DNA-dosage trans-effect","Fig4H_volcano_DNAdose.pdf")
gene_volcano("contW_rna","RNA-dosage trans-effect — WHITE only","FigS_Fig4H_volcano_RNAdose_WHITE.pdf")
print("done: bars, forest, RNA dose-response, WHITE dose-response + volcanoes (cap 12)")
