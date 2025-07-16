# lacertid-skull-NC
Script for the paper "Neural crest cell biology shapes lizard skull evolution across evolutionary time scales"

# Scripts (to be the run in order):

- Podarcis-IT.R : Analyses of shape variations in Podarcis muralis
- lacertids-phylogenetic-analyses.R : Preparing phylogenetic data on Lacertidae. Conduct morphological phylogenetic analyses in Lacertidae
- lacertids-phylogenetic-analyses-not-residuals.R : follow-up on modified-plot-functions-bayou.R and lacertids-phylogenetic-analyses.R to check the consistency of the results when using raw data instead of residuals from a regression of shape on size and origin. 
- Podarcis-distances-vs-Lacertid-rates-variances.R : regression of per-landmark Euclidean distances between extreme values of the nigriventris syndrome and landmark  variance or evolutuonary rates in Lacertidae. 
- BayesTraitEvolvability.R : phyologenetic modelling with Fabric model

# Data and metadata: 
- metadata-Podarcis : metadata file for P. muralis specimens
- metadata-lacertid : metadata file for Lacertidae specimens
- 2024-12-03_15_49_50 : SlicerMorph output folders containing the landmark configruations of the subset of data for P. muralis
- 2024-12-06_16_45_56 : SlicerMorph output folders containing the landmark configruations of the subset of data for Lacertidae
- 2024-12-19_10_33_32 : SlicerMorph output folders containing the landmark configruations of the entire data (P. muralis + Lacertidae)
- modules.txt : table of landmark names assiated with a neural crest- or a mesoderm origin
- 99814-mean-shape.ply reference-mesh-100429.ply : meshes of the specimens of P. muralis or Lacertidae to be used for generating the reference meshes
- squamates_Title_Science2024_ultrametric_constrained.tre : original tree for Title et al. (2024)
