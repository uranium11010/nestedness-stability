import numpy as np
import scipy.special as sc
from scipy.optimize import fsolve


NESTEDNESS_TYPES_MAP = {}
NESTEDNESS_TYPES_NAMES = {}

def register(label, name):
    def decorator(nestedness_func):
        NESTEDNESS_TYPES_MAP[label] = nestedness_func
        NESTEDNESS_TYPES_NAMES[label] = name
        return nestedness_func
    return decorator


def nestedness(L: np.ndarray, nestedness_type: str):
    return NESTEDNESS_TYPES_MAP[nestedness_type](L)


@register("nodf", "NODF")
def nestedness_nodf(L: np.ndarray):
    """
    L: adjacency matrix as 2D array; truthy/falsy values denote edges
    returns nestedness of `L` as float between 0 and 1, according to NODF
        (https://doi.org/10.1111/j.0030-1299.2008.16644.x)
    """
    L = make_nested(L.astype(bool))
    R, C = L.shape
    row_sums = L.sum(axis=1)
    col_sums = L.sum(axis=0)
    paired_overlap = 0.
    for i in range(R-1):
        for j in range(i+1, R):
            if row_sums[i] > row_sums[j]:
                paired_overlap += np.logical_and(L[i], L[j]).sum() / row_sums[j]
    for k in range(C-1):
        for l in range(k+1, C):
            if col_sums[k] > col_sums[l]:
                paired_overlap += np.logical_and(L[:,k], L[:,l]).sum() / col_sums[l]
    return paired_overlap / (R*(R-1)/2 + C*(C-1)/2)


@register("dipole_moment", r"$N_p$")
def nestedness_dipole_moment(L: np.ndarray):
    """
    L: adjacency matrix as 2D array; truthy/falsy values denote edges
    returns nestedness of `L` as float between 0 and 1, according to deviation
    from perfect nestedness following antidiagonal, with higher penalty
    farther from isocline
    """
    L = make_nested(L)
    R, C = L.shape
    score = 0
    for i in range(R):
        for j in range(C):
            dist = (i+0.5)/R + (j+0.5)/C - 1
            if L[i,j]:
                score -= dist
            else:
                score += dist
    return 3 * score / (R * C)


@register("antidiag_deviation", r"$N_{\text{ADD}}")
def nestedness_antidiag_deviation(L: np.ndarray):
    """
    L: adjacency matrix as 2D array; truthy/falsy values denote edges
    returns nestedness of `L` as float between 0 and 1, according to deviation
    from perfect nestedness following antidiagonal
    """
    L = make_nested(L)
    R, C = L.shape
    score = 0
    for i in range(R):
        for j in range(C):
            if ((i+0.5)/R + (j+0.5)/C <= 1) == L[i,j]:
                score += 1
    return score / (R * C)


@register("temperature", r"$N_T$")
def nestedness_temperature(L: np.ndarray):
    """
    L: adjacency matrix as 2D array; truthy/falsy values denote edges
    returns nestedness of `L` as float between 0 and 1, based on a simplified
    version of BINMATNEST (https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2699.2006.01444.x)
    """
    L = make_nested(L.astype(bool))
    first_not_full_i = int(np.sum(np.all(L, axis=1)).item())
    first_not_full_j = int(np.sum(np.all(L, axis=0)).item())
    L = L[max(1, first_not_full_i) - 1 :, max(1, first_not_full_j) - 1 :]
    R, C = L.shape

    def area_ipn_fill(p):
        superellipse_area = sc.gamma(1 + 1/p) ** 2 / sc.gamma(1 + 2/p)
        inner_rect_area = (1 - 1/R) * (1 - 1/C)
        return (1 - inner_rect_area) / 2 + inner_rect_area * superellipse_area

    fill = np.mean(L)
    if fill < (1 + (1 - 1/R) * (1 - 1/C)) / 2:
        p = fsolve(lambda p: area_ipn_fill(p) - fill, 1).item()
    else:
        p = np.inf

    def inside_ipn(i, j):
        if np.isinf(p):
            return True
        return (i / (R - 1)) ** p + (j / (C - 1)) ** p <= 1
    
    def normalized_dist_to_ipn(i, j):
        # y = i / (R - 1)
        # x = j / (C - 1)
        # t = fsolve(lambda t: (x + t) ** p + (y + t) ** p - 1, 0)
        # return abs(t) / (1 - abs(x - y) + 1/R + 1/C)
        y = (i + 0.5) / R
        x = (j + 0.5) / C
        if np.isinf(p):
            t = min(1 - 0.5/C - x, 1 - 0.5/R - y)
        else:
            t = fsolve(lambda t: ((x + t - 0.5/C) / (1 - 1/C)) ** p + ((y + t - 0.5/R) / (1 - 1/R)) ** p - 1, 0, factor=0.01).item()
        return abs(t) / (1 - abs(x - y))
    
    temp_score = 0
    for i in range(R):
        for j in range(C):
            if inside_ipn(i, j) != L[i,j]:
                temp_score += normalized_dist_to_ipn(i, j) ** 2
    
    return 1 - temp_score / (R * C * 0.04145)


@register("discrepancy", r"$N_d$")
def nestedness_discrepancy(L: np.ndarray):
    """
    L: adjacency matrix as 2D array; truthy/falsy values denote edges
    returns nestedness of `L` as float between 0 and 1 using
    Bruald-Sanderson discrepancy (https://link.springer.com/article/10.1007/s004420050784),
    normalized according to https://onlinelibrary.wiley.com/doi/full/10.1111/j.2006.0906-7590.04493.x
    """
    L = make_nested(L.astype(bool))
    R, C = L.shape
    row_sums = np.sum(L, axis=1)
    nested_L = np.array([[j < row_sums[i] for j in range(C)] for i in range(R)])
    return 1 - 4 * np.sum(~L & nested_L) / (R * C)


def make_nested(Fns):
    """
    order species in decreasing order of degree to make the topology look as nested as possible
    (used to assess degree of nestedness and detect complete topological nestedness)
    """
    R, C = Fns.shape
    rowNonzero = np.count_nonzero(Fns, axis=1)
    rNZ = []
    for i in range(R):
        rNZ.append((rowNonzero[i], i))
    colNonzero = np.count_nonzero(Fns, axis=0)
    cNZ = []
    for j in range(C):
        cNZ.append((colNonzero[j], j))
    rNZ.sort(key=lambda tup: tup[0], reverse=True)
    cNZ.sort(key=lambda tup: tup[0], reverse=True)

    F = np.empty(Fns.shape, dtype=Fns.dtype)
    for i in range(R):
        for j in range(C):
            F[i,j] = Fns[rNZ[i][1],cNZ[j][1]]

    return F


#determine whether the topology is completely nested
def is_nested(A): #A must be an output from makeNested(Fns)
    B = A.astype(bool)
    return bool(np.sum(np.append(np.full((1, B.shape[1]), True), B, 0) != np.append(B, np.full((1, B.shape[1]), False), 0)) == B.shape[1])


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', action='store_true')
    args = parser.parse_args()

    if args.benchmark:
        import time
        import matplotlib.pyplot as plt
        import seaborn as sns
        from foodweb import get_random_L_from_nested

        Ls = [get_random_L_from_nested(10, 10, s) for s in range(40) for _ in range(5)]
        nestedness_values = {nestedness_type: [] for nestedness_type in NESTEDNESS_TYPES_MAP}
        for nestedness_type in NESTEDNESS_TYPES_MAP:
            print(nestedness_type)
            start_time = time.time()
            for L in Ls:
                nestedness_values[nestedness_type].append(nestedness(L, nestedness_type))
            print('avg. time: {:.4f}ms'.format((time.time() - start_time) * 5))
        
        fig, axs = plt.subplots(len(nestedness_values), len(nestedness_values))
        fig.set_size_inches(12.8, 9.6)
        fig.tight_layout()
        nestedness_types = list(NESTEDNESS_TYPES_MAP)
        for i, ntype1 in enumerate(nestedness_types):
            for j, ntype2 in enumerate(nestedness_types):
                axs[i][j].scatter(nestedness_values[ntype1], nestedness_values[ntype2])
                axs[i][j].set_xlabel(ntype1)
                axs[i][j].set_ylabel(ntype2)
                sns.regplot(x=nestedness_values[ntype1], y=nestedness_values[ntype2],
                            ci=False, line_kws={'color': 'red'}, ax=axs[i][j])
        plt.show()
    else:
        Ls = [
            np.array([
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [1, 0, 1, 0, 1],
                ]),
            np.array([
                [0, 1, 1, 1, 1],
                [1, 0, 0, 1, 0],
                [1, 1, 0, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 0, 0, 1, 0],
                ]),
            np.array([
                [0, 1, 1, 1, 1],
                [1, 1, 1, 1, 1],
                [0, 0, 0, 1, 0],
                [0, 0, 1, 1, 0],
                [0, 1, 0, 1, 1],
                ]),
        ]
        for L in Ls:
            print(make_nested(L))
            for nestedness_type in NESTEDNESS_TYPES_MAP:
                print(f"{nestedness_type}: {nestedness(L, nestedness_type)}")