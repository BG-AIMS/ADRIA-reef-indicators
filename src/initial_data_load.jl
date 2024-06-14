using DataFrames
using DynamicCoralCoverModel
using ADRIA

# Local scale domain
dom = ADRIA.load_domain()

# GBR wide domain
dom = ADRIA.load_domain(RMEDomain, "","45")

# generate 128 sample scenarios
scens = ADRIA.sample(dom, 128)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(dom, scens, "45")

# ... or repeat scenario runs across multiple RCPs
rs = ADRIA.run_scenarios(dom, scens, ["45", "60", "85"])

# Analysis
