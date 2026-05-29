## Data
The datasets used in this repository are subsets extracted or calculated from the original dataset hosted on Dryad.

### Extracted Subsets
* `model_simulation/data/Q_data.csv`
* `model_simulation/data/C1_data.csv`
* `model_simulation/data/C2_data.csv`

### Calculated Summaries
The following files contain values calculated from `Q_data.csv`, `C1_data.csv`, and `C2_data.csv`, respectively:
* `model_simulation/data/Q_mean_CI_data.dat`
* `model_simulation/data/C1_mean_CI_data.dat`
* `model_simulation/data/C2_mean_CI_data.dat`

## Original Dataset
Boisseau, Romain P.; Vogel, David; Dussutour, Audrey (2016). Data from: Habituation in non-neural organisms: evidence from slime moulds [Dataset]. Dryad. https://doi.org/10.5061/dryad.51j89 (Licensed under CC0 1.0)

## Data Processing
* **For `Q_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Quinine` on days 1–6, `Bridge=Agar` on days 7–8, and `Bridge=Quinine` on day 9. Individuals with missing data on day 5 were excluded.
* **For `C1_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Agar` on days 1–5, `Bridge=Quinine` on day 6, `Bridge=Agar` on days 7–8, and `Bridge=Quinine` on day 9.
* **For `C2_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Agar` on days 1–9.