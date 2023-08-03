import argparse
import numpy as np
import pandas as pd
import plotly.express as px

from stability import STABILITY_TYPES_MAP, STABILITY_TYPES_NAMES
from nestedness import NESTEDNESS_TYPES_MAP, NESTEDNESS_TYPES_NAMES

parser = argparse.ArgumentParser()
parser.add_argument('--input-path', required=True)
parser.add_argument('--save-fig-path', required=True)
parser.add_argument('-Rs', nargs='+', type=int)
parser.add_argument('-Cs', nargs='+', type=int)
parser.add_argument('--nestedness-type', choices=list(NESTEDNESS_TYPES_MAP.keys()))
parser.add_argument('--stability-type', choices=list(STABILITY_TYPES_MAP.keys()))
parser.add_argument('--nestedness-types', nargs='+', choices=list(NESTEDNESS_TYPES_MAP.keys()))
parser.add_argument('--stability-types', nargs='+', choices=list(STABILITY_TYPES_MAP.keys()))

args = parser.parse_args()

def get_rvalues_array(R_list, C_list, nestedness_type, stability_type):
    R2idx = {R: idx for idx, R in enumerate(R_list)}
    C2idx = {C: idx for idx, C in enumerate(C_list)}

    plot_df = pd.read_csv(args.input_path)
    plot_df = plot_df[(plot_df["Nestedness metric"] == nestedness_type) & (plot_df["Instability metric"] == stability_type)]

    rvalues = np.zeros((len(R_list), len(C_list))) * np.nan
    pvalues = np.zeros((len(R_list), len(C_list))) * np.nan
    for row_idx, row in plot_df.iterrows():
        i = R2idx.get(row['R'])
        j = C2idx.get(row['C'])
        if i is not None and j is not None:
            assert np.isnan(rvalues[i,j]) and np.isnan(pvalues[i,j])
            rvalues[i,j] = row['rvalue']
            pvalues[i,j] = row['pvalue']

    print(f"Largest p-value: {np.max(pvalues)}")

    return rvalues

if args.nestedness_type is not None:
    assert args.stability_type is not None and args.nestedness_types is None and args.stability_types is None

    R_list = list(range(*args.Rs))
    C_list = list(range(*args.Cs))
    rvalues = get_rvalues_array(R_list, C_list, args.nestedness_type, args.stability_type)

    fig = px.imshow(rvalues, y=R_list, x=C_list, range_color=[-1, 1], color_continuous_scale="rdbu")

    fig.update_xaxes(tick0=4, dtick=2, title_text="$S_C$", side="top") 
    fig.update_yaxes(tick0=4, dtick=2, title_text="$S_R$") 
    fig.update_coloraxes(colorbar_title_text="Correlation")
    # fig.update_yaxes(matches=None) 
    # fig.update_traces(marker={"size": 2})
    fig.write_image(args.save_fig_path)
    fig.write_html('plots/figure.html', auto_open=True)

else:
    assert args.stability_types is not None and args.nestedness_type is None and args.stability_type is None

    R_list = list(range(*args.Rs))
    C_list = list(range(*args.Cs))

    rvalues_arrays = []
    label_list = []

    for stability_type in args.stability_types:
        for nestedness_type in args.nestedness_types:
            print(f"{nestedness_type} - {stability_type}")
            mathify_str = lambda string: string[1:-1] if string[0] == string[-1] == '$' else rf"\text{{{string}}}"
            stability_name = mathify_str(STABILITY_TYPES_NAMES[stability_type])
            nestedness_name = mathify_str(NESTEDNESS_TYPES_NAMES[nestedness_type])
            label_list.append(rf"${stability_name}\text{{ vs. }}{nestedness_name}$")
            rvalues_arrays.append(get_rvalues_array(R_list, C_list, nestedness_type, stability_type))

    facet_row_spacing = 0.16
    fig = px.imshow(np.array(rvalues_arrays), y=R_list, x=C_list,
                    facet_col=0, facet_col_wrap=len(args.nestedness_types),
                    range_color=[-1, 1], color_continuous_scale="rdbu",
                    width=len(args.nestedness_types) * 300, height=len(args.stability_types) * 300,
                    facet_row_spacing=facet_row_spacing)

    fig.update_xaxes(showticklabels=True, tick0=4, dtick=2, title_text="$S_C$", title_standoff=0, side="top") 
    fig.update_yaxes(showticklabels=True, tick0=4, dtick=2, title_text="$S_R$", title_standoff=0) 
    fig.update_coloraxes(colorbar_title_text="Correlation")
    annotation_ys = [(1 + facet_row_spacing) / len(args.stability_types) * i for i in range(len(args.stability_types) - 1, -1, -1)]
    fig.for_each_annotation(lambda a: a.update(text=label_list[int(a.text.split("=")[1])], 
                                               y=annotation_ys[int(a.text.split("=")[1]) // len(args.nestedness_types)],
                                               yanchor="top")) 
    fig.write_image(args.save_fig_path)
    fig.write_html('plots/figure.html', auto_open=True)