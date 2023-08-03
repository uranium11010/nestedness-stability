import numpy as np
import scipy.optimize as opt
import scipy.linalg
from nestedness import is_nested


def AB(S, sigma, tau):
    """
    get A(i) and B(j) from sigma and tau
    """
    A = np.empty((sigma[-1],), dtype=int)
    B = np.empty((tau[-1],), dtype=int)
    for k in range(S):
        A[sigma[k]:sigma[k+1]] = k+1
        B[tau[k]:tau[k+1]] = k+1

    return A, B


def get_E(S, sigma, tau):
    """
    get number of edges from sigma and tau
    """
    ans = 0
    for i in range(S):
        ans += sigma[S-i] * (tau[i+1] - tau[i])
    return ans


def get_random_L_fixed_E(R, C, E):
    """
    get random topology with given # edges
    """
    L = np.zeros((R,C), dtype=bool)
    while True:
        ind = np.random.choice(np.arange(R*C, dtype=int), (E,), False)
        for i in ind:
            L[i//C, i%C] = True
        if np.all(np.any(L, axis=1)) and np.all(np.any(L, axis=0)):
            break
    return L


def get_random_L_bernoulli(R, C, p):
    return np.random.uniform(size=(R, C)) < p


def get_random_L_nested(R, C, get_L=True):
    """
    get random topology given R,C
    """
    S = np.random.choice(np.arange(2, min(R,C)+1, dtype=int))
    sigma = np.concatenate((np.concatenate(([0], np.sort(np.random.choice(np.arange(1,R), size=(S-1,), replace=False)))), [R]))
    tau = np.concatenate((np.concatenate(([0], np.sort(np.random.choice(np.arange(1,C), size=(S-1,), replace=False)))), [C]))
    E = get_E(S, sigma, tau)
    if get_L:
        L = np.zeros(shape=(R, C), dtype=bool)
        for i in range(S):
            for j in range(S-i):
                L[sigma[i]:sigma[i+1], tau[j]:tau[j+1]] = True
        return S, sigma, tau, E, L
    else:
        return S, sigma, tau, E


def get_random_L_from_nested(R, C, swaps=0):
    while True:
        S, sigma, tau, E, L = get_random_L_nested(R, C, get_L=True)
        edge_locs = []
        non_edge_locs = []
        for i in range(R):
            for j in range(C):
                if L[i,j]:
                    edge_locs.append((i, j))
                else:
                    non_edge_locs.append((i, j))
        # print(0)
        # print(L)
        for t in range(swaps):
            edge1idx = np.random.randint(len(edge_locs))
            edge2idx = np.random.randint(len(non_edge_locs))
            L[edge_locs[edge1idx]], L[non_edge_locs[edge2idx]] = L[non_edge_locs[edge2idx]], L[edge_locs[edge1idx]]
            # print(edge_locs[edge1idx], non_edge_locs[edge2idx])
            edge_locs[edge1idx], non_edge_locs[edge2idx] = non_edge_locs[edge2idx], edge_locs[edge1idx]
            # print(t+1)
            # print(L)
        if np.all(np.sum(L, axis=0)) and np.all(np.sum(L, axis=1)):
            return L


def max_ent_model(SH, SP, d, e, S=None, sigma=None, tau=None, L=None, raise_bad=True):
    """
    get MaxEnt flow rates
    """
    if L is not None: #numerical solution
        def eqs(coef): #[a2, a3, ..., b1, b2, ...]
            D = coef[:(SH-1)] * np.dot(L, coef[(SH-1):])[1:] - d[1:]
            E = coef[(SH-1):] * np.dot(np.transpose(L), np.append(1, coef[:(SH-1)])) - e
            return np.append(D, E)

        sol, _, status, msg = opt.fsolve(eqs, np.append(d[1:], e), full_output=True)

        if raise_bad and status != 1:
            raise Exception(msg)

        A = np.append(1, sol[:(SH-1)])
        B = sol[(SH-1):]

        F = np.dot(np.transpose(np.array([A,])), np.array([B,])) * L #MODEL FLOW RATES
    
    else:
        #analytic solution
        A, B = AB(S, sigma, tau)
        F = np.zeros((SH, SP))
        for x in range(SH):
            for y in range(tau[S+1-A[x]]):
                F[x,y] = d[x]*e[y] \
                * np.prod([np.sum(d[:sigma[S-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+1-A[x])]) \
                / np.prod([np.sum(d[:sigma[S+1-t]]) + np.sum(e[:tau[t]]) - 1 for t in range(B[y], S+2-A[x])])

    if raise_bad and np.any(F < 0):
        raise Exception("Negative flow rates")

    return F


def get_random_F(R, C, S, sigma, tau, F_dist, e_d_dist, nested=True, max_ent=True, L=None, allow_neg=False):
    """
    generate random flow rate matrix
    """
    F = np.zeros((R, C))
    if max_ent:
        d = e_d_dist(size=(R,)) / R
        e = e_d_dist(size=(C,)) / C
        F = max_ent_model(R, C, d, e, S, sigma, tau, L)
        if not allow_neg and np.any(F < 0):
            return None
    else:
        if nested:
            A, _ = AB(S, sigma, tau)
            linkages = 0
            for i in range(R):
                linkages += tau[S+1-A[i]]
            for i in range(R):
                for j in range(tau[S+1-A[i]]):
                    F[i,j] = F_dist() * R*C / linkages
        else:
            F = F_dist(size=(R,C))
    return F


def get_random_F_config_model(L, dx, ey, SH, SP, N, fact=1000):
    """
    Use the configuration model to generate N random networks with nested topology given by L
    and aggregate flow rates given by dx and ey
    """
    assert is_nested(L)

    #multiply by fact to make dx, ey integers
    dx = np.rint(dx * fact).astype(int)
    ey = np.rint(ey * fact).astype(int)

    #generate stubs on resources and consumers (for later generation of multigraph)
    Hstubs = []
    for x in range(SH):
        for i in range(dx[x]):
            Hstubs.append(x)
    Pstubs = []
    for y in range(SP):
        for i in range(ey[y]):
            Pstubs.append(y)

    #create random networks
    Farrays = []
    for k in range(N):
        # print("random network #" + str(k))
        length = min(len(Hstubs), len(Pstubs))

        #create a multigraph by joining stubs
        success = False
        while not success:
            F = np.zeros((SH, SP))
            copy = [Pstubs[i] for i in range(length)]
            for s in range(length-1,-1,-1):
                good = False
                for t in range(s+1):
                    if L[Hstubs[s],copy[t]]:
                        good = True
                        break
                if not good:
                    break #no available stub on a predator to connect to prey; we're stuck so restart process of connecting stubs
                t = np.random.randint(s+1)
                
                #find a random available stub on a predator to connect to prey
                while not L[Hstubs[s],copy[t]]:
                    t = np.random.randint(s+1)
                F[Hstubs[s],copy[t]] += 1
                copy.pop(t)
            
            success = s == 0
        
        F /= fact #flow network from multigraph by dividing by fact
        Farrays.append(F)

    return Farrays


def get_random_F_MH_reparamed(L: np.ndarray, dx, ey, SH, SP, N, warmup=1000, sigma=None):
    """
    Use the the Metropolis-Hastings (MH) algorithm to generate N random networks
    with topology given by L and aggregate flow rates given by dx and ey
    Random steps are conducted on reparameterization of flow network such that new parameters
    are free to change without violating aggregate flow constraints
    """
    param_labels = []  # list of labels (i, j, k, l) of free parameters
    for i in range(SH-1):
        for j in range(i+1, SH):
            for k in range(SP-1):
                for l in range(k+1, SP):
                    if L[i,k] and L[i,l] and L[j,k] and L[j,l]:
                        param_labels.append((i, j, k, l))
    Fs = []
    F = max_ent_model(SH, SP, dx, ey, L=L)
    if sigma is None:
        sigma = 0.5 / np.sqrt(SH * SP)
    for n in range(N):
        for t in range(warmup):
            F_step = np.zeros(F.shape)
            params_step = np.random.normal(scale=sigma, size=len(param_labels))
            for param, param_step in enumerate(params_step):
                i, j, k, l = param_labels[param]
                F_step[i,k] += param_step
                F_step[i,l] -= param_step
                F_step[j,k] -= param_step
                F_step[j,l] += param_step
            new_F = F + F_step
            if np.all(new_F >= 0):
                F = new_F
        Fs.append(F)
    return Fs


def get_random_F_MH_projected(L: np.ndarray, dx, ey, SH, SP, N, warmup=100, delta=None):
    """
    Use the Metropolis-Hastings (MH) algorithm to generate N random networks
    with topology given by L and aggregate flow rates given by dx and ey
    Random steps are conducted on entries of flow matrix, and we project
    onto nullspace of constraint matrix to make sure constraints are satisfied
    """
    L = L.astype(bool)
    cum_num_edges = np.concatenate([[0], np.cumsum(L.sum(axis=1))])
    E = cum_num_edges[-1]
    constraints = np.zeros((SH+SP-1, E))
    for i in range(SH):
        constraints[i, cum_num_edges[i]:cum_num_edges[i+1]] = 1
    cum_edges_rows = np.concatenate([np.zeros((SH, 1)), np.cumsum(L, axis=1)], axis=1).astype(int)
    for j in range(SP-1):
        for i in range(SH):
            if L[i,j]:
                constraints[SH+j, cum_num_edges[i]+cum_edges_rows[i,j]] = 1
    projection = scipy.linalg.pinv(constraints) @ constraints - np.identity(E)
    Fs = []
    F = max_ent_model(SH, SP, dx, ey, L=L)
    # print("maxent:", F)
    F_flat = F[L]
    if delta is None:
        delta = 0.2
    sigma = F_flat * delta
    DeltaF = np.zeros(E)
    for n in range(N):
        for t in range(warmup):
            F_step = np.random.normal(scale=sigma)
            new_DeltaF = projection @ (DeltaF + F_step)
            if np.all(F_flat + new_DeltaF >= 0):
                DeltaF = new_DeltaF
        new_F_flat = F_flat + DeltaF
        new_F = np.zeros(F.shape)
        for i in range(SH):
            for j in range(SP):
                if L[i,j]:
                    new_F[i,j] = new_F_flat[cum_num_edges[i] + cum_edges_rows[i,j]]
        Fs.append(new_F)
    return Fs


def get_random_F_from_maxent(L: np.ndarray, dx, ey, SH, SP, N, sigma=0.7, project=False):
    """
    Use deviations from the MaxEnt solution to generate N random networks
    with topology given by L and aggregate flow rates given by dx and ey.
    Deviations are log-normal factors multiplying the MaxEnt flow matrix,
    and if `project` is True, then we project onto nullspace of constraint
    matrix to make sure constraints are satisfied; the process is repeated
    if there's a resultant negative flow.
    """
    L = L.astype(bool)
    Fs = []
    F = max_ent_model(SH, SP, dx, ey, L=L)
    # print("maxent:", F)
    if project:
        cum_num_edges = np.concatenate([[0], np.cumsum(L.sum(axis=1))])
        E = cum_num_edges[-1]
        constraints = np.zeros((SH+SP-1, E))
        for i in range(SH):
            constraints[i, cum_num_edges[i]:cum_num_edges[i+1]] = 1
        cum_edges_rows = np.concatenate([np.zeros((SH, 1)), np.cumsum(L, axis=1)], axis=1).astype(int)
        for j in range(SP-1):
            for i in range(SH):
                if L[i,j]:
                    constraints[SH+j, cum_num_edges[i]+cum_edges_rows[i,j]] = 1
        projection = scipy.linalg.pinv(constraints) @ constraints - np.identity(E)
        F_flat = F[L]
    for n in range(N):
        if project:
            while True:
                DeltaF = projection @ ((np.exp(np.random.normal(size=F_flat.shape, scale=sigma)) - 1) * F_flat)
                new_F_flat = F_flat + DeltaF
                if np.all(new_F_flat >= 0):
                    break
            new_F = np.zeros(F.shape)
            for i in range(SH):
                for j in range(SP):
                    if L[i,j]:
                        new_F[i,j] = new_F_flat[cum_num_edges[i] + cum_edges_rows[i,j]]
        else:
            new_F = np.exp(np.random.normal(size=F.shape, scale=sigma)) * F
        Fs.append(new_F)
    return Fs


if __name__ == "__main__":
    print(get_random_L_from_nested(4, 6, swaps=10))