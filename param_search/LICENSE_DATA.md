## Data
The datasets used in this repository are subsets extracted from the original dataset hosted on Dryad.

### Extracted Subsets
* `data/Q_data.csv`
* `data/C1_data.csv`
* `data/C2_data.csv`

## Original Dataset
Boisseau, Romain P.; Vogel, David; Dussutour, Audrey (2016). Data from: Habituation in non-neural organisms: evidence from slime moulds [Dataset]. Dryad. https://doi.org/10.5061/dryad.51j89 (Licensed under CC0 1.0)

## Data Processing
* **For `Q_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Quinine` on days 1–6, `Bridge=Agar` on days 7–8, and `Bridge=Quinine` on day 9. Individuals with missing data on day 5 were excluded.
* **For `C1_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Agar` on days 1–5, `Bridge=Quinine` on day 6, `Bridge=Agar` on days 7–8, and `Bridge=Quinine` on day 9.
* **For `C2_data.csv`:** We extracted the daily `Crossing_time` for individuals exposed to `Bridge=Agar` on days 1–9.