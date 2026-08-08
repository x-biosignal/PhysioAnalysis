"""Generate the spm1d / rft1d numerical-parity reference for PhysioAnalysis (VAL-01).

Replaces the Worsley closed-form *surrogate* with fixtures captured from an
ACTUAL spm1d + rft1d run:
  * RFT critical thresholds  (rft1d.t.isf / rft1d.f.isf over an (alpha,df,Q,FWHM) grid)
  * SPM{t} two-sample + paired, and SPM{F} one-way ANOVA on fixed seeded data
    (t/F field, residual FWHM, critical height zstar, and cluster-level inference).

Writes /tmp/spm1d_ref/spm1d_reference.json; data-raw/spm1d_reference.R turns it
into tests/testthat/fixtures/spm1d-reference.rds. Data are stored so the R side
runs spmTTest/spmPairedTTest/spmAnova on the SAME input.

PhysioAnalysis stores fields as (Q nodes x N obs); spm1d uses (N obs x Q), so the
R fixture transposes. Two-tailed t uses alpha/2 (matching .rftCriticalT); F is
one-tailed.

Run:  micromamba run -n spm1d python data-raw/spm1d_reference.py
"""
import json, os
import numpy as np
import spm1d, rft1d, scipy
from scipy.ndimage import gaussian_filter1d

OUT = os.environ.get("OUT", "/tmp/spm1d_ref")
os.makedirs(OUT, exist_ok=True)
Q = 101

def smooth(m, fwhm):
    return gaussian_filter1d(m, sigma=fwhm / np.sqrt(8 * np.log(2)), axis=1, mode="nearest")

def clusters_of(inf):
    out = []
    for c in inf.clusters:
        e = [float(x) for x in c.endpoints]
        out.append({"endpoints": e, "P": float(c.P), "extent": float(c.extent)})
    return out

# ---- 1. RFT critical thresholds (rft1d) ------------------------------------
rft_t, rft_f = [], []
for df in (6, 10, 20, 40):
    for fwhm in (5.0, 10.0, 20.0):
        rft_t.append(dict(alpha=0.05, df=df, Q=Q, fwhm=fwhm,
                          resel=(Q - 1) / fwhm,
                          crit=float(rft1d.t.isf(0.05 / 2, df, Q, fwhm))))
for df1 in (2, 3):
    for df2 in (20, 40):
        for fwhm in (5.0, 10.0, 20.0):
            rft_f.append(dict(alpha=0.05, df1=df1, df2=df2, Q=Q, fwhm=fwhm,
                              resel=(Q - 1) / fwhm,
                              crit=float(rft1d.f.isf(0.05, (df1, df2), Q, fwhm))))

# ---- 1b. RFT cluster-extent probabilities (rft1d.p_cluster), one-tailed ------
# The parity target for .rftClusterPValue: rft1d's cluster-extent distribution
# at (extent-in-resels k, threshold u, df, resel count R).
clusterp_t, clusterp_f = [], []
for df in (10, 22, 40):
    for u in (2.5, 3.5, 4.5):
        for R in (5.0, 13.0, 20.0):
            fwhm = (Q - 1) / R
            for k in (0.1, 0.3, 0.6):
                clusterp_t.append(dict(df=df, u=u, R=R, fwhm=fwhm, k=k,
                                       p=float(rft1d.t.p_cluster(k, u, df, Q, fwhm))))
for (df1, df2) in ((2, 27), (3, 40)):
    for u in (6.0, 9.0):
        for R in (5.0, 13.0):
            fwhm = (Q - 1) / R
            for k in (0.1, 0.3):
                clusterp_f.append(dict(df1=df1, df2=df2, u=u, R=R, fwhm=fwhm, k=k,
                                       p=float(rft1d.f.p_cluster(k, u, (df1, df2), Q, fwhm))))

# ---- 2. SPM{t} two-sample (effect -> informative cluster p; and a null) -----
def ttest2_case(seed, nA, nB, bump, lo=45, hi=56):
    rng = np.random.default_rng(seed)
    yA = smooth(rng.standard_normal((nA, Q)), 8.0)
    yB = smooth(rng.standard_normal((nB, Q)), 8.0)
    if bump:
        yB[:, lo:hi] += bump
    t = spm1d.stats.ttest2(yA, yB)
    ti = t.inference(0.05, two_tailed=True, interp=False)
    df = int(nA + nB - 2)
    return dict(yA=yA.tolist(), yB=yB.tolist(), alpha=0.05,
                df=df, fwhm=float(t.fwhm), zstar=float(ti.zstar),
                rft_crit=float(rft1d.t.isf(0.05 / 2, df, Q, t.fwhm)),
                tfield=[float(x) for x in t.z], clusters=clusters_of(ti))

ttest2 = {"effect": ttest2_case(0, 12, 12, 0.7, lo=47, hi=54),
          "null":   ttest2_case(7, 12, 12, 0.0)}

# ---- 3. SPM{t} paired ------------------------------------------------------
def paired_case(seed, n, bump, lo=45, hi=56):
    rng = np.random.default_rng(seed)
    yA = smooth(rng.standard_normal((n, Q)), 8.0)
    yB = yA + smooth(rng.standard_normal((n, Q)), 8.0) * 0.6
    if bump:
        yB[:, lo:hi] += bump
    t = spm1d.stats.ttest_paired(yA, yB)
    ti = t.inference(0.05, two_tailed=True, interp=False)
    df = int(n - 1)
    return dict(yA=yA.tolist(), yB=yB.tolist(), alpha=0.05, df=df,
                fwhm=float(t.fwhm), zstar=float(ti.zstar),
                rft_crit=float(rft1d.t.isf(0.05 / 2, df, Q, t.fwhm)),
                tfield=[float(x) for x in t.z], clusters=clusters_of(ti))

paired = {"effect": paired_case(3, 14, 0.25, lo=49, hi=52)}

# ---- 4. SPM{F} one-way ANOVA ----------------------------------------------
def anova1_case(seed, per, bump, lo=45, hi=56):
    rng = np.random.default_rng(seed)
    groups = []
    for g in range(3):
        y = smooth(rng.standard_normal((per, Q)), 8.0)
        if g == 2 and bump:
            y[:, lo:hi] += bump
        groups.append(y)
    Y = np.vstack(groups)
    A = np.array([0] * per + [1] * per + [2] * per)
    F = spm1d.stats.anova1(Y, A)
    Fi = F.inference(0.05, interp=False)
    J, k = Y.shape[0], 3
    df1, df2 = k - 1, J - k                     # F-statistic df (matches the field)
    return dict(Y=Y.tolist(), labels=[int(a) for a in A], alpha=0.05,
                df1=df1, df2=df2,
                fwhm=float(F.fwhm), zstar=float(Fi.zstar),
                rft_crit=float(rft1d.f.isf(0.05, (df1, df2), Q, F.fwhm)),
                spm1d_df=[int(Fi.df[0]), int(Fi.df[1])],  # spm1d's own (differs: df2-1)
                Ffield=[float(x) for x in F.z], clusters=clusters_of(Fi))

anova1 = {"effect": anova1_case(9, 10, 0.4, lo=49, hi=52)}

ref = {"rft_t": rft_t, "rft_f": rft_f,
       "clusterp_t": clusterp_t, "clusterp_f": clusterp_f,
       "ttest2": ttest2, "paired": paired, "anova1": anova1, "Q": Q,
       "provenance": {"spm1d": spm1d.__version__, "rft1d": rft1d.__version__,
                      "numpy": np.__version__, "scipy": scipy.__version__,
                      "smoothing": "scipy.ndimage.gaussian_filter1d, mode=nearest, FWHM 8.0",
                      "note": "two-tailed t uses alpha/2; F is one-tailed; fields (Nobs x Q) here, transposed to (Q x N) in R"}}
json.dump(ref, open(os.path.join(OUT, "spm1d_reference.json"), "w"))
print("wrote", os.path.join(OUT, "spm1d_reference.json"))
print("versions:", ref["provenance"])
for k in ("effect", "null"):
    c = ttest2[k]["clusters"]
    print("ttest2", k, "zstar=%.4f nclust=%d" % (ttest2[k]["zstar"], len(c)),
          [round(x["P"], 4) for x in c])
print("paired effect: zstar=%.4f nclust=%d" % (paired["effect"]["zstar"], len(paired["effect"]["clusters"])),
      [round(x["P"], 4) for x in paired["effect"]["clusters"]])
print("anova1 effect: zstar=%.4f nclust=%d" % (anova1["effect"]["zstar"], len(anova1["effect"]["clusters"])),
      [round(x["P"], 4) for x in anova1["effect"]["clusters"]])
