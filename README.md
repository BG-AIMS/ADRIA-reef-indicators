# ADRIA-reef-indicators
Time series analyses for GBR health using ADRIA and CoralBlox model.
- Looking for reefs that decline early compared to their surrounding reefs. These could be considered as bellwether reefs that indicate oncoming decline.
- Looking for reefs that have a delay in decline compared to their surrounding reefs. These could be considered as more resilient to the surrounding stressors.

## Structure
``` code
ADRIA-reef-indicators/
├─ src/           # Analysis scripts
├─ outputs/
├─ figs/
├─ data/
├─ .gitignore
├─ config.toml    # configuration file for ADRIA
├─ Project.toml   # Julia project spec
└─ README.md      # this file
```
### Analysis scrips - `src`
- `1_initial_results.jl` : Loads RME domain and produces initial result set

- `2_initial_analysis_and_metrics.jl` : Produces coral cover metrics from `rs` result set

- `3_Subregion_correlation.jl` : Lagged correlation analysis at subregion scale

- `4_Bioregion_correlation.jl` : Lagged correlation analysis at bioregion scale

- `5_exploratory_analysis.jl` : Exploratory analysis of common reefs identified at subregion and bioregion scale.

## Methods
### Domain
Using GBR-wide domain from ReefModEngine data: "rme_ml_2024_01_08".

Currently using RCP 45 from domain.

#### Includes:
- Initial coral cover data from RME
- DHW/environmental data from RME
- Reef spatial data from RME (Canonical-Reefs geopackage)
- Connectivity data from (`connectivity team source`)

### Scales
Timeseries are grouped by reefs into their subregion (currently based on the closest port to each reef), and their bioregion (GBRMPA reefal bioregion data).
The timeseries are then compared for each reef within these groups with the overall group median timeseries.

### Lagged Correlation Analyses
Lagged correlation analysis is performed between each reef in a group and the group's median time series.
A strong-positive correlation at positive lags means a reef's changes are early compared to group-median. A strong-positive correlation at negative lags means a reef's changes are delayed compared to the group-median.

## Results
Results from lagged correlation analyses are used to explore characteristics of reefs that occur early/late compared to their surrounding group.
Reefs with timeseries that occur early/late at both subregion/bioregion scales are considered to be candidates for further analyses.

These candidate reefs are then explored for their common characteristics that differ from other surrounding reefs, such as:

- strength/number of outgoing and incoming larval connections
- initial coral cover values
- reef size
- reef depth
- proportion of functional groups occupying reefs (e.g. functional group that is dominant)
- environmental stressors compared to surrounding reefs (e.g. consistently higher DHW compared to surrounding reefs)