# lightIonSummerSchool_2025

Get the repo:
```git clone https://github.com/KongTu/lightIonSummerSchool_2025.git```
```cd lightIonSummerSchool_2025```

## Setting up the environment

Fresh start with an eic_shell, if one hasn't done it yet:

```
wget --output-document install.sh http://get.epic-eic.org --no-check-certificate
	
bash install.sh
```
Run it:

```./eic-shell```

Run the miniTreeMaker:

first argument is inputfile
second argument is outputfile
third argument is 1: save MC, 0: no MC

```./run.sh root://dtn-eic.jlab.org//volatile/eic/EPIC/RECO/main/epic_craterlake/Test/Tutorials_June2025/EXCLUSIVE/DIFFRACTIVE_PHI_ABCONV/sartre1.39-1.0/eCa/coherent/bsat/18x137.5/q2_1to10/sartre1.39-1.0_coherent_phi_eCa_bsat_18x137.5_q2_1to10_ab.\*.eicrecon.edm4eic.root sartre_MiniTree_1M 1```
