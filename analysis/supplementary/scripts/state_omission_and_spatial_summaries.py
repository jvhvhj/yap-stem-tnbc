# Purpose: Leave-one-cancer-state-out and paired-section summaries
# Inputs: module README and config/input_manifest.tsv.
# Outputs: preserved source filename conventions in a separate analysis workspace.
# Related manuscript: Supplementary Fig. S1–S4.
"""Requested evidence-gap audit only. No plotting and no writes to existing results.

Run with the project's scientific Python environment (scipy/h5py required).
The full-analysis reproduction gates must pass before state omissions are run.
"""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats
import h5py

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
DATA = OUT / "supporting_data"
DATA.mkdir(exist_ok=True)
SOURCES: dict[Path,str] = {}

def source(rel, role):
    p = Path(rel)
    if not p.is_absolute(): p = ROOT / p
    if not p.exists(): raise FileNotFoundError(p)
    SOURCES[p] = role
    return p

def tab(rel, role):
    return pd.read_csv(source(rel,role), sep="\t")

def write(df, name):
    df.to_csv(DATA/name,sep="\t",index=False,na_rep="NA",float_format="%.17g")

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(8*1024**2), b""): h.update(b)
    return h.hexdigest()

def load_functions(path,names,ns):
    """Load named function definitions only; NEVER execute a legacy pipeline body."""
    tree=ast.parse(path.read_text(encoding="utf-8-sig"))
    nodes=[x for x in tree.body if isinstance(x,(ast.FunctionDef,ast.AsyncFunctionDef)) and x.name in names]
    assert set(x.name for x in nodes)==set(names)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),"exec"),ns)
    return ns

cache_path=source("0716_yan2026_validation/inputs/yan_cancer_cell_scores.tsv.gz","ARTEMIS cell identifiers, author states, original scores, states and technical covariates")
h5path=source("scouter/af8c4fce-4c63-4671-b339-91a383cf36f6.h5ad","Authoritative ARTEMIS cell-level normalized expression and annotation")
a_script=source("Figure5_panelA_stage1/scripts/01_panelA_stage1_audit.py","Current frozen patient-level YAP-Stem raw and technical estimators")
p_script=source("0728_sensitivity_and_cnv_closure/scripts/03_yan_state_and_gate_sensitivity.py","Frozen Program146 expression extraction and median High-minus-Other effect definition")
prep=source("0716_yan2026_validation/scripts/01_prepare_yan_and_frozen_signatures.py","Author state metadata and original pooled YS_state assignment")
refA=tab("Figure5_panelA_stage1/PanelA_patient_level_complete.tsv","78-patient full-analysis coupling reference")
refP=tab("0728_sensitivity_and_cnv_closure/06_Yan_sensitivity/Yan_primary146_state_definition_sensitivity.tsv","77-testable-patient full-program reference")
support=tab("Figure5_N1_state_programme_testability/N1_patient_state_support.tsv","Authoritative 97 x 11 patient-state support grid")
meas=tab("Figure5_N1_state_programme_testability/N1_programme_measurability_audit.tsv","146 membership, 141 measurable and five unavailable")
coverage=tab("0728_sensitivity_and_cnv_closure/06_Yan_sensitivity/Yan_direction_gate_gene_coverage.tsv","Exact existing symbol mapping, no new alias search")
primaryfile=source("0728_sensitivity_and_cnv_closure/02_direction_gate_sensitivity/gate_6of8_primary_genes.txt","Existing exact 6/8 primary 146-gene list explicitly checked by the original Yan sensitivity script; recovers its stale PRIMARY path")
primary=primaryfile.read_text(encoding="utf-8-sig").split()
assert len(primary)==len(set(primary))==146
assert set(primary)==set(meas.loc[meas.record_type.eq("gene"),"gene"])
current_program=tab("Figure3_bottom_panels_HI_feasibility_and_plotting/Figure3_I_effect_consistency_source.tsv","Independent current strict_v2 146 membership cross-check")
assert set(primary)==set(current_program.loc[current_program.program_146.astype(str).str.lower().eq("true"),"gene"])

ns={"np":np,"pd":pd,"stats":stats}
ns=load_functions(a_script,["z_within","residualize","safe_spearman"],ns)
pns={"np":np,"pd":pd,"sparse":sparse,"h5py":h5py}
pns=load_functions(p_script,["contiguous_runs","extract_csr_rows_columns","residualize_within_patient"],pns)
cache=pd.read_csv(cache_path,sep="\t")
assert len(cache)==49275 and cache.cell_id.is_unique
states=list(support.sort_values("state_order").state.unique())
assert len(states)==11 and set(cache.cell_state)==set(states)|{"Unresolved"}
patients=sorted(cache.patient.unique(),key=lambda x:int(x[1:]))
counts=cache.groupby("patient").size()
eligible=counts[counts>=50].index
d=cache[cache.patient.isin(eligible)].copy().reset_index(drop=True)
assert len(eligible)==78 and len(d)==48885
assert set(eligible)==set(refA.patient_id)

print("A: checking original H5AD cell identity and all author states",flush=True)
def decode_col(h,key):
    obj=h[key]
    if isinstance(obj,h5py.Group):
        cats=obj["categories"].asstr()[:]
        codes=obj["codes"][:]
        return np.array([cats[c] if c>=0 else "NA" for c in codes],dtype=object)
    return obj.asstr()[:] if obj.dtype.kind in "OSU" else obj[:]

with h5py.File(h5path,"r") as h:
    ids=h["obs/_index"].asstr()[:]
    rr=pd.Index(ids).get_indexer(cache.cell_id)
    assert (rr>=0).all()
    for column,expected in [("donor_id",cache.patient),("cell_state",cache.cell_state)]:
        assert np.array_equal(decode_col(h,"obs/"+column)[rr],expected.values)
    assert np.all(decode_col(h,"obs/author_cell_type")[rr]=="Tumor")
    assert tuple(h["X"].attrs["shape"])==(427823,32354)
    gene_symbols=h["var/gene_symbols"].asstr()[:]
    rows=pd.Index(ids).get_indexer(d.cell_id)
    assert (np.diff(rows)>0).all()
    selected=coverage.set_index("requested_gene").loc[primary]
    selected=selected[selected.present.astype(str).str.lower().eq("true")]
    assert len(selected)==141
    columns=selected.h5ad_var_index.astype(int).to_numpy()
    assert np.array_equal(gene_symbols[columns],selected.resolved_symbol.to_numpy())
    present_meas=set(meas.loc[meas.record_type.eq("gene") & meas.measurable_in_yan.astype(str).str.lower().eq("true"),"gene"])
    assert set(selected.index)==present_meas

# Recover only the frozen arithmetic-mean score; no new genes or normalization.
# float32 expression extraction follows the exact existing sensitivity function.
print("A: recovering fixed 141-measurable-gene score from authoritative X",flush=True)
expr=pns["extract_csr_rows_columns"](h5path,rows,columns)
program_score=np.asarray(expr.mean(axis=1)).ravel()
assert len(program_score)==len(d) and np.isfinite(program_score).all()
d["program_score"]=program_score
del expr

support_rows=[]
for p in patients:
    dd=cache[cache.patient.eq(p)]
    for state in states:
        z=dd[dd.cell_state.eq(state)]
        support_rows.append(dict(patient_id=p,state=state,n_malignant_patient=len(dd),n_cells_state=len(z),
                                 n_high_state=int(z.YS_state.eq("High").sum()),
                                 n_other_state=int(z.YS_state.isin(["Intermediate","Low"]).sum()),
                                 primary_patient=p in eligible))
support_new=pd.DataFrame(support_rows)
ck=support_new.merge(support,on=["patient_id","state"],suffixes=("_new","_ref"),validate="one_to_one")
assert len(ck)==1067
for col in ["n_malignant_patient","n_cells_state","n_high_state","n_other_state"]:
    assert (ck[col+"_new"]==ck[col+"_ref"]).all()
write(support_new,"S3D_patient_state_counts.tsv")

def scenario_effects(omitted):
    z=d if omitted=="FULL" else d[~d.cell_state.eq(omitted)]
    z=z.copy().reset_index(drop=True)
    prog_adj=pns["residualize_within_patient"](z,z.program_score.to_numpy())
    result=[]
    for p in sorted(eligible,key=lambda x:int(x[1:])):
        pos=np.flatnonzero(z.patient.eq(p).to_numpy())
        dd=z.iloc[pos]
        hi=dd.YS_state.eq("High").to_numpy()
        oth=~hi
        qualified=len(dd)>=50
        raw=tech=effect=effectadj=np.nan
        if qualified:
            cov=np.column_stack([np.log1p(dd.nCount_RNA),dd.nFeature_RNA,dd.S_score,dd.G2M_score])
            raw=ns["safe_spearman"](dd.YAP_score,dd.Stemness_score)[0]
            tech=ns["safe_spearman"](ns["residualize"](dd.YAP_score,cov),ns["residualize"](dd.Stemness_score,cov))[0]
            if hi.any() and oth.any():
                effect=float(np.median(dd.program_score.to_numpy()[hi])-np.median(dd.program_score.to_numpy()[oth]))
                effectadj=float(np.median(prog_adj[pos][hi])-np.median(prog_adj[pos][oth]))
        result.append(dict(omitted_state=omitted,patient_id=p,n_cells_full=int(counts[p]),n_cells_remaining=len(dd),
                           n_cells_omitted=int(counts[p])-len(dd),n_high=int(hi.sum()),n_other=int(oth.sum()),
                           eligible_remaining_ge50=qualified,coupling_testable=bool(qualified and np.isfinite(raw) and np.isfinite(tech)),
                           rho_raw=raw,rho_technical=tech,program_testable=bool(np.isfinite(effect) and np.isfinite(effectadj)),
                           program_raw_high_minus_other=effect,program_technical_high_minus_other=effectadj,
                           program_low_group_support_lt5=bool(hi.sum()<5 or oth.sum()<5),
                           exclusion_reason="" if qualified else "Fewer than original minimum 50 remaining malignant cells",
                           high_definition="Original YS_state unchanged; no re-tertiling",
                           program_definition="Median of fixed 141-gene cell mean in High minus Other"))
    return pd.DataFrame(result)

full=scenario_effects("FULL")
qa=full.merge(refA,on="patient_id",suffixes=("_new","_ref"))
raw_diff=float(np.max(abs(qa.rho_raw-qa.raw_rho)))
tech_diff=float(np.max(abs(qa.rho_technical-qa.technical_rho)))
assert raw_diff<1e-12 and tech_diff<1e-12,(raw_diff,tech_diff)
ps=refP[refP.record_type.eq("patient") & refP.definition.eq("Primary pooled tertile")]
pdiffs={}
for model,col in [("raw","program_raw_high_minus_other"),("technical_adjusted","program_technical_high_minus_other")]:
    jj=full.merge(ps[ps.model.eq(model)],left_on="patient_id",right_on="patient",validate="one_to_one")
    assert np.array_equal(jj.program_testable,jj.testable)
    diff=float(np.nanmax(abs(jj[col]-jj.median_high_minus_other)))
    assert diff<=1e-7,(model,diff)
    pdiffs[model]=diff
assert full.program_testable.sum()==77
print("A: full baseline reproduction passed",raw_diff,tech_diff,pdiffs,flush=True)
effects=[full]
for state in states:
    effects.append(scenario_effects(state))
    print("A: omission complete",state,flush=True)
effects=pd.concat(effects,ignore_index=True)
write(effects,"S3D_leave_one_state_out_patient_effects.tsv")

summaries=[]
for omit,g in effects.groupby("omitted_state",sort=False):
    row=dict(omitted_state=omit,n_cells_remaining=int(g.n_cells_remaining.sum()),
             n_eligible_patients=int(g.eligible_remaining_ge50.sum()),
             patients_below50=";".join(g.loc[~g.eligible_remaining_ge50,"patient_id"]) or "NONE",
             program_low_group_support_lt5_n=int((g.program_testable & g.program_low_group_support_lt5).sum()))
    for col in ["rho_raw","rho_technical","program_raw_high_minus_other","program_technical_high_minus_other"]:
        vals=g[col].dropna().to_numpy()
        row.update({col+"_tested":len(vals),col+"_positive":int(sum(vals>0)),col+"_median":float(np.median(vals)),
                    col+"_q1":float(np.quantile(vals,.25)),col+"_q3":float(np.quantile(vals,.75)),
                    col+"_min":float(min(vals)),col+"_max":float(max(vals))})
    summaries.append(row)
summ=pd.DataFrame(summaries)
write(summ,"S3D_full_and_11_omissions.tsv")
write(pd.DataFrame([dict(raw_rho_max_error=raw_diff,technical_rho_max_error=tech_diff,
                         program_raw_max_error=pdiffs["raw"],program_adjusted_max_error=pdiffs["technical_adjusted"],
                         program_float32_reproduction_tolerance=1e-7,original_cell_scores_changed=False,
                         baseline_state_assignments_changed=False,full_reference_pass=True)]),"S3D_reference_reproduction_QC.tsv")

# B: paired coefficients are recovered, not re-estimated from spot data.
paired_path="Figure6_S1_GSE210616_patient_recurrent_spatial_coupling/Figure6_S1_paired_section_consistency.tsv"
paired=tab(paired_path,"21 original paired section coefficients; values are tumour+depth adjusted")
sec=tab("Figure6_S1_GSE210616_patient_recurrent_spatial_coupling/Figure6_S1_section_associations.tsv","Verifies exact adjusted coefficient identity by sample")
main=tab("Figure6A_final_reset/Figure6A_final_reset_source.tsv","Current Main Figure 6A paired-section coefficient cross-check")
source(".tmp/Figure6_S1/assemble_s1_outputs.R","Original paired-section assembly and rank-residual definition")
assert len(paired)==21 and paired.patient_id.is_unique and "P19" not in set(paired.patient_id)
for n in [1,2]:
    v=sec.set_index("sample_id").loc[paired[f"section_{n}_sample"],"rho_tumour_depth_adjusted"].to_numpy()
    assert np.max(abs(v-paired[f"section_{n}_rho"]))<1e-14
    m=main[main.record_type.eq("PAIRED_SECTION")].set_index("patient_id").loc[paired.patient_id,f"section_{n}_rho"].astype(float).to_numpy()
    assert np.max(abs(m-paired[f"section_{n}_rho"]))<1e-14
paired["delta_rho_section2_minus_section1"]=paired.section_2_rho-paired.section_1_rho
paired["absolute_delta_rho"]=abs(paired.delta_rho_section2_minus_section1)
write(paired,"S4B_paired_section_differences.tsv")
x=paired.section_1_rho.to_numpy(); y=paired.section_2_rho.to_numpy()
cov=np.mean((x-x.mean())*(y-y.mean()))
ccc=2*cov/(np.var(x)+np.var(y)+(x.mean()-y.mean())**2)
delta=y-x
pairsummary=dict(n_patients=21,concordant_positive_pairs=int(sum((x>0)&(y>0))),
                 lin_concordance_correlation=ccc,spearman_between_sections=float(stats.spearmanr(x,y).statistic),
                 pearson_between_sections=float(stats.pearsonr(x,y).statistic),
                 median_delta_rho=float(np.median(delta)),mean_delta_rho=float(np.mean(delta)),
                 median_absolute_delta_rho=float(np.median(abs(delta))),
                 absolute_delta_q1=float(np.quantile(abs(delta),.25)),absolute_delta_q3=float(np.quantile(abs(delta),.75)),
                 maximum_absolute_delta=float(max(abs(delta))),largest_difference_patient=paired.iloc[np.argmax(abs(delta))].patient_id,
                 new_p_values="NONE",coefficient_identity="TumorPurity+nCount adjusted partial rank correlation")
write(pd.DataFrame([pairsummary]),"S4B_concordance_summary.tsv")

# C: audit graph definitions and source portability only, do not choose new k/radius.
graph=tab("Figure6E_neighborhood_spatial_coorganization/Figure6E_spatial_graph_QC.tsv","Primary first-order hex-grid graph and 43 source pairs")
source(".tmp/Figure6_enrichment_gate/run_DE_gates.R","Exact primary adjacency offsets, row weights, Lee L formula, null, aggregation")
source(".tmp/Figure6_enrichment_gate/build_final_outputs.mjs","Copies formal D/E analysis outputs to publication analysis directories")
source("Figure6E_neighborhood_spatial_coorganization/Figure6E_readout.md","Primary first-order graph scientific readout")
source("Figure6E_neighborhood_spatial_coorganization/Figure6E_primary_test.tsv","Retained primary Lee L patient-level test")
graph["original_score_path_exists"]=graph.score_source.map(lambda p:Path(p).exists())
graph["original_coordinate_path_exists"]=graph.coordinate_source.map(lambda p:Path(p).exists())
graph["archived_coordinate_path"]=[str(ROOT/".tmp/Figure6_S0Q2/geo_objects"/r.sample_id/Path(r.coordinate_source).name) for r in graph.itertuples()]
graph["archived_coordinate_exists"]=graph.archived_coordinate_path.map(lambda p:Path(p).exists())
graph["alternative_graph_prespecified"]=False
graph["extension_status"]="NOT JUSTIFIED under current prespecification; primary score chunks not recovered"
write(graph,"S4C_graph_provenance_audit.tsv")

# D: copy all three effects and link each existing LOO value by omitted patient.
bs=tab("Figure6G1_BSW2_spatial_replication/Figure6G1_patient_associations.tsv","Nine BSW2 raw/depth/nFeature patient-level effects")
bsn=tab("Figure6G1_BSW2_spatial_replication/Figure6G1_nFeature_sensitivity.tsv","Independent nFeature table and cohort summary")
loo=tab("Figure6G1_BSW2_spatial_replication/Figure6G1_LOO.tsv","Existing nine leave-one-patient-out medians")
btest=tab("Figure6G1_BSW2_spatial_replication/Figure6G1_primary_test.tsv","Primary nine-patient median and exact Wilcoxon result")
source("Figure6G1_BSW2_spatial_replication/scripts/01_run_Figure6G1.R","Frozen BSW2 normalization, patient effect and LOO generation")
source("Figure6G1_BSW2_spatial_replication/Figure6G1_readout.md","sAA9 author-final exclusion and unresolved reason")
expected={"sAA1","sAA2","sAA6","sAA8","sEA1","sEA2","sEA5","sEA6","sEA7"}
assert len(bs)==9 and set(bs.patient_id)==expected and set(loo.omitted_patient)==expected
for row in loo.itertuples():
    assert abs(row.median_rho-np.median(bs.loc[bs.patient_id.ne(row.omitted_patient),"rho_depth_adjusted"]))<1e-14
assert np.max(abs(bs.set_index("patient_id").rho_nFeature_adjusted-
                  bsn[bsn.record_type.eq("PATIENT")].set_index("patient_id").rho_nFeature_adjusted))<1e-14
assert abs(np.median(bs.rho_depth_adjusted)-btest.median_rho.iloc[0])<1e-14
bsout=bs.merge(loo[["omitted_patient","median_rho","positive_n","n_patients"]],
              left_on="patient_id",right_on="omitted_patient",validate="one_to_one")
bsout=bsout.rename(columns={"median_rho":"loo_median_when_this_patient_omitted",
                            "positive_n":"loo_positive_n","n_patients":"loo_remaining_patients"})
write(bsout,"S4D_BSW2_existing_robustness.tsv")
bs_summary=[]
for col in ["rho_raw","rho_depth_adjusted","rho_nFeature_adjusted"]:
    v=bs[col]
    bs_summary.append(dict(model=col,n=9,positive_n=int(sum(v>0)),median=float(np.median(v)),
                           q1=float(np.quantile(v,.25)),q3=float(np.quantile(v,.75)),minimum=float(min(v)),maximum=float(max(v))))
write(pd.DataFrame(bs_summary),"S4D_BSW2_summary.tsv")

decisions=pd.DataFrame([
    dict(candidate="S3D",status="PASS",computation="Full reference reproduced, then 11 state omissions using original scores/states and existing estimators",
         claim_boundary="Descriptive omission robustness; non-independent overlapping omission datasets; no new P values; retain coverage losses and low High/Other counts"),
    dict(candidate="S4B",status="PASS",computation="21 existing adjusted pairs, delta and absolute delta; Lin CCC primary descriptive agreement, Spearman secondary association",
         claim_boundary="Same positive direction is not equivalent to high magnitude agreement; no new section correlation or patient test"),
    dict(candidate="S4C-extension",status="PARTIAL",computation="Primary first-order graph recovered; extension NOT JUSTIFIED; no alternative graph run",
         claim_boundary="No prespecified alternative k/radius/hop definition; cannot equate computability with scientific justification; original full score chunks absent at recorded locations"),
    dict(candidate="S4D",status="PASS",computation="All frozen nine-patient raw/depth/nFeature effects and nine LOO medians recovered and reconciled",
         claim_boundary="sAA9 absent from author-final cohort; exclusion reason unresolved, not documented QC failure")])
decisions.to_csv(OUT/"Supplementary_S3S4_gap_audit.tsv",sep="\t",index=False)

print("Hashing exact sources (including authoritative H5AD)",flush=True)
manifest=[]
for p,role in SOURCES.items():
    manifest.append(dict(file=str(p),role=role,bytes=p.stat().st_size,sha256=sha(p),modified=False))
manifest=pd.DataFrame(manifest)
assert manifest.loc[manifest.file.eq(str(cache_path)),"sha256"].iloc[0]=="eaf404d2d0863bc8888da526cda7c5f0450d4eb9c36fce2e305f21a3e552c5e9"
assert manifest.loc[manifest.file.eq(str(h5path)),"sha256"].iloc[0]=="601f90ce2329459681868cc72cfd8ee631625202aca4570169783cde84ae1063"
manifest.to_csv(OUT/"authoritative_sources.tsv",sep="\t",index=False)
details={"raw_full_max_difference":raw_diff,"technical_full_max_difference":tech_diff,
         "program_full_max_differences":pdiffs,"paired_summary":pairsummary,
         "states":states,"primary_unresolved_cells":int(d.cell_state.eq("Unresolved").sum()),
         "primary_program_genes":146,"measurable_genes":141,
         "original_graph_score_files_available":int(graph.original_score_path_exists.sum()),
         "archived_graph_coordinate_files_available":int(graph.archived_coordinate_exists.sum()),
         "bsw2_loo_min":float(loo.median_rho.min()),"bsw2_loo_max":float(loo.median_rho.max())}
(DATA/"audit_details.json").write_text(json.dumps(details,ensure_ascii=False,indent=2),encoding="utf-8")
print(decisions.to_string(index=False),flush=True)
print(summ[["omitted_state","n_eligible_patients","rho_raw_median","rho_raw_positive","rho_technical_median","rho_technical_positive",
            "program_raw_high_minus_other_tested","program_raw_high_minus_other_positive","program_raw_high_minus_other_median",
            "program_technical_high_minus_other_positive","program_technical_high_minus_other_median"]].to_string(index=False),flush=True)
print("PAIRED",pairsummary,flush=True)
print("AUDIT_DONE_NO_FIGURES",flush=True)
