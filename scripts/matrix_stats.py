"""Statistics on the interaction matrix shared by the figure scripts and the regenerated-matrix analyses: PERMANOVA of the phage host range,
NODF nestedness and BRIM modularity with null matrices of the same size and connectance."""
import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform


def design(factors, data):
    cols = [np.ones(len(data))]
    for f in factors:
        cols.append(pd.get_dummies(data[f], drop_first=True).values.astype(float))
    return np.column_stack(cols)


def permanova(Dm, X, Xr, perms=999, seed=0):
    """Pseudo-F test of full design X against reduced Xr on a distance matrix (Anderson 2001); returns R2, adjusted R2, p."""
    n = len(Dm); G = -0.5 * (np.eye(n) - 1 / n) @ (Dm ** 2) @ (np.eye(n) - 1 / n)
    hat = lambda Z: Z @ np.linalg.pinv(Z)
    Hf, Hr = hat(X), hat(Xr)
    def stat_(perm):
        Gp = G[np.ix_(perm, perm)]
        ssm = np.trace(Hf @ Gp @ Hf) - np.trace(Hr @ Gp @ Hr)
        sse = np.trace(Gp) - np.trace(Hf @ Gp @ Hf)
        return (ssm / (np.linalg.matrix_rank(X) - np.linalg.matrix_rank(Xr))) / (sse / (n - np.linalg.matrix_rank(X)))
    f0 = stat_(np.arange(n)); rng = np.random.default_rng(seed)
    p = (1 + sum(stat_(rng.permutation(n)) >= f0 for _ in range(perms))) / (perms + 1)
    tot = np.trace(G); r2 = (np.trace(Hf @ G @ Hf) - np.trace(Hr @ G @ Hr)) / tot
    rk = np.linalg.matrix_rank(X)
    return r2, 1 - (1 - (np.trace(Hf @ G @ Hf) / tot)) * (n - 1) / (n - rk), p


def phage_distance(M):
    """Jaccard distance between the binary lytic vectors of the phages (columns of M; missing = not lytic)."""
    return np.nan_to_num(squareform(pdist((M.fillna(0) > 0).values.T.astype(bool), "jaccard")))


def phage_models(Dm, meta, perms=999):
    """PERMANOVA of the phage host range on genus, isolation strain and both (meta has columns Genus and Host)."""
    one = np.ones((len(meta), 1)); Xg, Xh = design(["Genus"], meta), design(["Host"], meta)
    return {"genus": permanova(Dm, Xg, one, perms), "isolation strain": permanova(Dm, Xh, one, perms),
            "genus x isolation strain": permanova(Dm, np.column_stack([Xg, Xh[:, 1:]]), one, perms)}


def nodf(A):
    deg = A.sum(1); ov = A @ A.T
    def one(deg, ov):
        n = len(deg); i, j = np.triu_indices(n, 1)
        di, dj = deg[i], deg[j]; mn = np.minimum(di, dj)
        pn = np.where((di != dj) & (mn > 0), ov[i, j] / np.maximum(mn, 1), 0.0)
        return pn.sum(), n * (n - 1) / 2
    r = one(deg, ov); c = one(A.sum(0), A.T @ A)
    return 100 * (r[0] + c[0]) / (r[1] + c[1])


def brim(A, c, tries=3, rg=None):
    """Barber's bipartite modularity, optimised by BRIM (alternating assignment of rows and columns to c modules)."""
    m = A.sum(); k, dg = A.sum(1), A.sum(0); B = A - np.outer(k, dg) / m; best = -1
    for _ in range(tries):
        R = np.eye(c)[rg.integers(0, c, A.shape[1])]   # column memberships
        q_old = -1
        for _ in range(50):
            Sr = np.eye(c)[(B @ R).argmax(1)]
            R = np.eye(c)[(B.T @ Sr).argmax(1)]
            q = np.trace(Sr.T @ B @ R) / m
            if q - q_old < 1e-9: break
            q_old = q
        best = max(best, q)
    return best


def nested_modular(A, n_null=1000, seed=0):
    """NODF and BRIM modularity of the binary matrix A against n_null Bernoulli matrices of the same size and connectance."""
    rng = np.random.default_rng(seed)
    T0, dens = nodf(A), A.mean()
    nulls = np.array([nodf((rng.random(A.shape) < dens).astype(float)) for _ in range(n_null)])
    p_nodf = (1 + (nulls >= T0).sum()) / (n_null + 1)
    q_by_c = {c: brim(A, c, rg=rng) for c in range(2, 9)}
    cbest = max(q_by_c, key=q_by_c.get); Q0 = q_by_c[cbest]
    qn = np.array([brim((rng.random(A.shape) < dens).astype(float), cbest, tries=2, rg=rng) for _ in range(n_null)])
    return dict(T0=T0, nulls=nulls, p_nodf=p_nodf, cbest=cbest, Q0=Q0, qn=qn, p_q=(1 + (qn >= Q0).sum()) / (n_null + 1))
