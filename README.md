# Nestedness and stability in a MaxEnt energy flow model

This is one of two repositories associated with the paper
"Nestedness Promotes Stability in Maximum-Entropy Bipartite Food Webs" by Li and Harte, 2023.
The other repository is located [https://github.com/uranium11010/network-flow-model](here).

If you use our code in your work, please cite our paper:
```
@article{li2023nestedness,
  title={Nestedness Promotes Stability in Maximum-Entropy Bipartite Food Webs},
  author={Li, Zhening and Harte, John},
  year={2023}
}
```

## Dependencies

Code is implemented in Python 3. Dependencies are installed with the following command:
```shell
pip install numpy scipy pandas plotly tqdm statsmodels kaleido
```
Installation should take about a minute.

## Example usage

Code in this repository is used to analyze the relationship between the nestedness of a bipartite food web and stability of the ecological community
for various combinations of number of resources ($S_R$) and number of consumers ($S_C$).

Single $(S_R, S_C)$ pair: (takes about a minute to run)
```shell
python stability.py -R 10 -C 6 --num-swaps-fraction 0.4 --topologies-per-swap 25 --num-samples 1000 --random-F-type from_maxent --random-eta-type log_normal --nestedness-types dipole_moment nodf temperature discrepancy --stability-types mean fraction_nonzero --seed 0 --output-path output/10_6_040_25_1000
```
Scatter plot: (takes a few seconds to run)
```shell
python plot_scatter.py --input-path output/10_6_040_25_1000.csv --save-fig-path plots/10_6_040_25_1000_dipole_mean.pdf --nestedness-type dipole_moment --stability-type mean
```

Range of $(S_R, S_C)$ pairs: (takes about a day to run)
```shell
python stability.py -Rs 3 21 -Cs 3 21 --num-swaps-fraction 0.4 --topologies-per-swap 25 --num-samples 1000 --random-F-type from_maxent --random-eta-type log_normal --nestedness-types dipole_moment --stability-types mean --seed 0 --output-path output/3-21_3-21_040_25_1000
```
Plot correlations: (takes a few seconds to run)
```shell
python plot_correlation.py --input-path output/3-21_3-21_040_25_1000_batch_regress.csv --save-fig-path plots/3-21_3-21_040_25_1000_heatmap.pdf -Rs 3 21 -Cs 3 21 --nestedness-type dipole_moment --stability-type mean
```
