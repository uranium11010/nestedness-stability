import argparse
from tqdm import tqdm
import math
import numpy as np
import pandas as pd
from scipy.stats import linregress
from foodweb import get_random_F, get_random_L_from_nested, get_random_F_config_model, get_random_F_MH_projected, get_random_F_from_maxent, max_ent_model
from nestedness import make_nested, nestedness, NESTEDNESS_TYPES_MAP

np.seterr(all='raise')

RANDOM_F_TYPES_MAP = {
    "config_model": get_random_F_config_model,
    "MH_projected": get_random_F_MH_projected,
    "from_maxent": get_random_F_from_maxent
}
RANDOM_ETA_TYPES_MAP = {
    "normal": lambda size: np.max(np.array(np.random.normal(size=size, scale=0.5) + 1, np.zeros(size)), axis=0),
    "log_normal": lambda size: np.exp(np.random.normal(size=size, scale=0.5))
}
STABILITY_TYPES_MAP = {
    "fraction_zero": lambda x: 1 - np.count_nonzero(x) / len(x),
    "fraction_nonzero": lambda x: np.count_nonzero(x) / len(x),
    "mean": np.mean
}
STABILITY_TYPES_NAMES = {
    "fraction_zero": "$p_{m = 0}$",
    "fraction_nonzero": "$p_{m > 0}$",
    "mean": r"$\overline m$",
}

def run_stability_simulation(args):
    topologies = []
    d = np.ones(args.R) / args.R
    e = np.ones(args.C) / args.C
    max_num_swaps = math.ceil(args.num_swaps_fraction * args.R * args.C)
    num_topologies = args.num_topologies
    topologies_per_swap = args.topologies_per_swap
    if num_topologies is None:
        num_topologies = max_num_swaps * topologies_per_swap
        print("No. topologies:", num_topologies)
    if topologies_per_swap is None:
        topologies_per_swap = num_topologies / max_num_swaps
    for i in range(num_topologies):
        while True:
            L = make_nested(get_random_L_from_nested(args.R, args.C, swaps=int(i // topologies_per_swap)))
            try:
                F = max_ent_model(args.R, args.C, d, e, L=L, raise_bad=True)
                topologies.append(L)
                if args.verbose >= 2:
                    print(i)
                    print('L', L)
                    print('F', F)
                break
            except:
                continue
    nestednesses = {nestedness_type: [] for nestedness_type in args.nestedness_types}
    stabilities = {stability_type: [] for stability_type in args.stability_types}
    edge_nums = []
    swap_nums = []
    for i, L in enumerate(tqdm(topologies)):
        L = np.array(L)
        if args.verbose >= 2:
            print('topology', L)
        eig_val_list = []
        max_imag_sqrt_list = []
        sample_Fs = []
        for F in RANDOM_F_TYPES_MAP[args.random_F_type](L, np.ones(shape=args.R) / args.R, np.ones(shape=args.C) / args.C, args.R, args.C, args.num_samples):
            sample_Fs.append(F)
            if args.verbose >= 3:
                print('random', F)
            F1 = F * RANDOM_ETA_TYPES_MAP[args.random_eta_type](F.shape)
            FTF = F1.T @ F if args.R >= args.C else F @ F1.T
            eig_val, eig_vec = np.linalg.eig(FTF)
            eig_val_list.extend(eig_val)
            max_imag_sqrt_list.append(np.max(np.imag(np.sqrt(eig_val.astype(complex)))))
        sample_Fs = np.array(sample_Fs)
        # print('mean random:', np.mean(sample_Fs, axis=0))
        # print('stdev random:', np.std(sample_Fs, axis=0))
        # PLOT EIGENVALUES IN COMPLEX PLANE
        # eig_val_list = np.array(eig_val_list)
        # x = np.real(eig_val_list)
        # y = np.imag(eig_val_list)
        # axs[i].set_xlim([0., 0.02])
        # y = np.ones(x.shape) * i
        # plt.scatter(x, y, label=str(i), s=10)
        # PLOT HISTOGRAM OF LARGEST REAL PART
        # axs[i].set_ylim([0., num_samples])
        # axs[i].hist(x, label="{:.4f}".format(nestedness), range=[0., 0.02], bins=100)
        # axs[i].legend()
        # REGRESS WITH NESTEDNESS
        for nestedness_type in args.nestedness_types:
            nestedness = NESTEDNESS_TYPES_MAP[nestedness_type](L)
            nestednesses[nestedness_type].append(nestedness)
            if args.verbose >= 1:
                print('{}: {:.4f}'.format(nestedness_type, nestedness))
        max_imag_sqrt_list = np.array(max_imag_sqrt_list)
        for stability_type in args.stability_types:
            stability = STABILITY_TYPES_MAP[stability_type](max_imag_sqrt_list)
            stabilities[stability_type].append(stability)
            if args.verbose >= 1:
                print('{}: {:.4f}'.format(stability_type, stability))
        edge_nums.append(np.sum(F.astype(bool)))
        swap_nums.append(int(i // topologies_per_swap))
    data_dfs = []
    regress_dfs = []
    idx = 0
    for nestedness_type in args.nestedness_types:
        for stability_type in args.stability_types:
            data_dfs.append(pd.DataFrame({
                "Nestedness": nestednesses[nestedness_type],
                "Instability": stabilities[stability_type],
                "Nestedness metric": nestedness_type,
                "Instability metric": stability_type,
                "No. swaps": swap_nums
            }))
            regress_result = linregress(nestednesses[nestedness_type], stabilities[stability_type])
            regress_dfs.append(pd.DataFrame({
                "rvalue": regress_result.rvalue,
                "pvalue": regress_result.pvalue,
                "Nestedness metric": nestedness_type,
                "Instability metric": stability_type,
                "slope": regress_result.slope,
                "intercept": regress_result.intercept,
                "stderr": regress_result.stderr,
                "intercept_stderr": regress_result.intercept_stderr
            }, [idx]))
            print(f"{nestedness_type} - {stability_type}")
            print(regress_result)
            idx += 1
    return pd.concat(data_dfs), pd.concat(regress_dfs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-R', type=int)
    parser.add_argument('-C', type=int)
    parser.add_argument('-Rs', nargs='+', type=int)
    parser.add_argument('-Cs', nargs='+', type=int)
    parser.add_argument('--num-topologies', type=int)
    parser.add_argument('--topologies-per-swap', type=int)
    parser.add_argument('--num-swaps-fraction', type=float)  # overriden if both --num-topologies and --topologies-per-swap are specified
    parser.add_argument('--num-samples', required=True, type=int)
    parser.add_argument('--random-F-type', required=True, choices=list(RANDOM_F_TYPES_MAP.keys()))
    parser.add_argument('--random-eta-type', required=True, choices=list(RANDOM_ETA_TYPES_MAP.keys()))
    parser.add_argument('--nestedness-types', nargs='+', choices=list(NESTEDNESS_TYPES_MAP.keys()))
    parser.add_argument('--stability-types', nargs='+', choices=list(STABILITY_TYPES_MAP.keys()))
    parser.add_argument('--seed', type=int)
    parser.add_argument('--output-path', required=True)
    parser.add_argument('--verbose', '-v', type=int, default=0)

    args = parser.parse_args()

    np.seterr(all='warn')

    if args.seed is not None:
        np.random.seed(args.seed)

    if args.R is not None:
        assert args.C is not None
        assert args.Rs is None and args.Cs is None

        data_df, regress_df = run_stability_simulation(args)
        data_df.to_csv(args.output_path + ".csv")
        regress_df.to_csv(args.output_path + "_regress.csv")
    
    else:
        assert args.C is None
        assert args.Rs is not None and args.Cs is not None

        regress_dfs = []
        for R in range(*args.Rs):
            for C in range(*args.Cs):
                print('R:', R, 'C:', C)
                args.R = R
                args.C = C
                _, regress_df = run_stability_simulation(args)
                regress_df['R'] = R
                regress_df['C'] = C
                regress_dfs.append(regress_df)
        regress_df = pd.concat(regress_dfs)
        regress_df.to_csv(args.output_path + "_batch_regress.csv")