import argparse
import pandas as pd
import plotly.express as px

from stability import STABILITY_TYPES_NAMES
from nestedness import NESTEDNESS_TYPES_NAMES

def format_axis_title(var_type, var_name):
    if var_name[0] == '$' and var_name[-1] == '$':
        stripped_var_name = var_name[1:-1]
        return rf"$\text{{{var_type}, }}{stripped_var_name}$"
    return f"{var_type}, {var_name}"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-path', required=True)
    parser.add_argument('--save-fig-path', required=True)
    parser.add_argument('--nestedness-type', choices=list(NESTEDNESS_TYPES_NAMES.keys()))
    parser.add_argument('--stability-type', choices=list(STABILITY_TYPES_NAMES.keys()))

    args = parser.parse_args()

    plot_df = pd.read_csv(args.input_path)

    if args.nestedness_type is None:
        assert args.stability_type is None

        fig = px.scatter(plot_df, x="Nestedness", y="Instability",
                        facet_row="Instability metric", facet_col="Nestedness metric",
                        # color="No. swaps",
                        trendline="ols")
        fig.update_xaxes(matches=None) 
        fig.update_yaxes(matches=None) 
        fig.update_traces(marker={"size": 2})
        fig.write_image(args.save_fig_path)
        fig.write_html('plots/figure.html', auto_open=True)

    else:
        assert args.stability_type is not None

        plot_df = plot_df[(plot_df["Nestedness metric"] == args.nestedness_type) & (plot_df["Instability metric"] == args.stability_type)]

        fig = px.scatter(plot_df, x="Nestedness", y="Instability",
                        labels = {"Nestedness": format_axis_title("Nestedness", NESTEDNESS_TYPES_NAMES[args.nestedness_type]),
                                "Instability": format_axis_title("Instability", STABILITY_TYPES_NAMES[args.stability_type])},
                        trendline="ols")
        fig.update_traces(marker={"size": 2})
        fig.write_image(args.save_fig_path)
        fig.write_html('plots/figure.html', auto_open=True)