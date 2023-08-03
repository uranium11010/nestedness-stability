# Nestedness and stability in a MaxEnt energy flow model

Example usage:
Single $(S_R, S_C)$ pair:
```shell
python stability.py -R 10 -C 6 --num-swaps-fraction 0.4 --topologies-per-swap 25 --num-samples 1000 --random-F-type from_maxent --random-eta-type log_normal --nestedness-types dipole_moment nodf temperature discrepancy --stability-types mean fraction_nonzero --seed 0 --output-path output/10_6_040_25_1000
```
Scatter plot:
```shell
python plot_scatter.py --input-path output/10_6_040_25_1000.csv --save-fig-path plots/10_6_040_25_1000_dipole_mean.pdf --nestedness-type dipole_moment --stability-type mean
```

Range of $(S_R, S_C)$ pairs:
```shell
python stability.py -Rs 3 21 -Cs 3 21 --num-swaps-fraction 0.4 --topologies-per-swap 25 --num-samples 1000 --random-F-type from_maxent --random-eta-type log_normal --nestedness-types dipole_moment nodf --stability-types fraction_nonzero mean --seed 0 --output-path output/3-21_3-21_040_25_1000
```
Plot correlations:
```shell
python plot_correlation.py --input-path output/3-21_3-21_040_25_1000_batch_regress.csv --save-fig-path plots/3-21_3-21_040_25_1000_heatmap.pdf -Rs 3 21 -Cs 3 21 --nestedness-type dipole_moment --stability-type mean
```
```shell
python plot_correlation.py --input-path output/3-21_3-21_040_25_1000_batch_regress.csv --save-fig-path plots/3-21_3-21_040_25_1000_heatmap.pdf -Rs 3 21 -Cs 3 21 --nestedness-types dipole_moment nodf temperature discrepancy --stability-types mean fraction_nonzero
```