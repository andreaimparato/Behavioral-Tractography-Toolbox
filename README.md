# Behavioral-Tractography-Toolbox

MATLAB toolbox for constructing and visualizing 3D Multi-Layer Networks (3D-MLN), and for computing **Behavioral Tractography (BT)** and **Behavioral Diffusion Analysis (BDA)** metrics from EMA-based time-series data.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Execution Guide](#execution-guide)
- [Scripts](#scripts)
- [Functions](#functions)
- [Inputs](#inputs)
- [MLNetwork](#mlnetwork)
- [BT2](#bt2)
- [Repository Authors](#repository-authors)

## Overview
This repository provides the full project structure required to run the analyses directly in MATLAB.

The repository is organized around three main components:
- **Scripts**: primary MATLAB scripts implementing the main analysis workflows.
- **Functions**: internal helper functions used by the scripts.
- **Inputs**: MATLAB input files required to run the analyses.

In addition, the repository contains a dedicated **BT2** folder for the clustering analyses described in the paper:

> *Shedding light on the dynamic interplay of positive and negative symptoms of psychosis with Behavioral Tractography*

This folder contains the commented script `BT2_Clustering.m` together with its required input files.

For questions regarding the repository structure or analysis pipeline, please contact the repository authors.

## Repository Structure

```text
Behavioral-Tractography-Toolbox/
├── Scripts/      # Main MATLAB analysis scripts
├── Functions/    # Internal helper functions
├── Inputs/       # Required input .mat and related files
├── BT2/          # Clustering analysis for the BT2 study
├── MLNetwork/    # Resources related to MLN visualization
└── README.md
```

## System Requirements

### Software
- **MATLAB version:** `R2024a`
- The code has been tested with MATLAB **R2024a**.
- Other MATLAB versions may work, but they have not been explicitly tested.

### Hardware
- A standard personal computer is sufficient.
- No specialized hardware is required.

### Operating System
The code is platform-independent and should run on any operating system supporting MATLAB R2024a, including:
- Windows
- macOS
- Linux

## Installation
No special installation procedure is required.

### 1. Install MATLAB
Make sure MATLAB **R2024a** is installed on your computer.

### 2. Clone or download this repository
You can either clone the repository with Git:

```bash
git clone https://github.com/andreaimparato/Behavioral-Tractography-Toolbox.git
```

or download it manually using the **Code > Download ZIP** option on GitHub.

### 3. Add repository folders to the MATLAB path
Open MATLAB, navigate to the repository directory, and add the relevant folders to your MATLAB path if needed.

## Execution Guide
1. Open MATLAB.
2. Navigate to the local repository folder.
3. Make sure all required input files are present in the `Inputs/` folder.
4. Run the desired script from the `Scripts/` folder.

### Expected execution time
- `MLN_Network_1.m`: approximately **2 minutes 30 seconds** on a standard desktop computer.
- `BT_BDA.m`: approximately **a few seconds** on a standard desktop computer.

## Scripts
This section describes the structure and purpose of the MATLAB scripts used in the analyses. For the theoretical rationale behind each step, please refer to the main text of the associated manuscript. Detailed comments are also included directly in the MATLAB code.

### `MLN_Network_1.m`
This script implements the analysis pipeline required for the construction and visualization of **3D Multi-Layer Networks (3D-MLN)** from time-series EMA data.

The input file loaded by the script contains EMA time-series data collected in the **Typically-Developing Replication** sample described in the main text.

The script performs the following main steps:
- Computes cross-sectional adjacency matrices using **mixed-effects linear models**.
- Performs **Network Dimensionality Reduction (NDR)** on the cross-sectional adjacency matrix.
- Tests NDR results through the relationship between **Euclidean distance** and **correlation strength**, as described in the main text.
- Constructs a **Multi-Layer Adjacency Matrix** that incorporates both cross-sectional and longitudinal correlations between EMA variables.
- Computes the metrics required for 3D-MLN visualization.
- Performs graph-theory analyses to compute **shortest paths** between EMA variables across temporal layers.
- Computes **Dynamic Betweenness Centrality** measures for each EMA variable.
- Compares Dynamic Betweenness estimates against randomly structured networks to identify EMA variables acting as **gateways** or **funnels** for psychological-contextual interactions.

For visualization, the script generates a **JSON** file containing the information required by the dedicated online 3D-MLN visualization platform:

- [MLNetwork visualization platform](https://dev.mlnetwork-diplab.ch/)

The platform is automatically launched at the end of the script.

### MLNetwork web platform for 3D Multi-Layer Network visualization
JSON files generated by `MLN_Network_1.m` can be uploaded to the dedicated 3D-MLN web platform for visualization.

The platform supports:
- 3D visualization and manipulation of the network (rotation, translation, zooming)
- Manual selection of network sub-portions
- Modification of visualization parameters, including node and edge color/size
- Export of modified JSON files
- Export of high-resolution figures and videos

A link to the 3D-MLN representation of the main networks is provided in the main text.

### `BT_BDA.m`
This script implements the analysis pipeline required to compute **Behavioral Tractography (BT)** and **Behavioral Diffusion Analysis (BDA)** metrics.

Please refer to the main text for the rationale and interpretation of these measures.

The input file loaded by the script contains:
- the **3D-MLN**, and
- the **cross-sectional coordinates** obtained from `MLN_Network_1.m`

The script performs the following main steps:
- Computes shortest paths connecting variables across temporal layers.
- Encodes each path as the 3D coordinates of the variables traversed.
- Decomposes each path into four XYZ coordinates corresponding to:
  - TL1 starting node
  - TL1 exit node
  - TL1 entry node
  - TL2 ending node
- Applies **k-means clustering** to these path coordinates to identify bundles of psychological-contextual pathways with similar 3D trajectories.
- Uses a consensus procedure across all samples to determine the optimal number of clusters/bundles.
- Estimates BDA at the symptom level by averaging the trajectories mediated by each symptom within each BT bundle.

The script also generates one **JSON** file per BT bundle for visualization on the same dedicated web platform:

- [MLNetwork visualization platform](https://dev.mlnetwork-diplab.ch/)

The platform is automatically launched at the end of the script.

## Functions
The scripts in this repository rely on both **internal** and **external** functions.

- **Internal functions** were developed by the project authors and are included in this repository.
- **External functions** are required for successful execution but must be downloaded separately.

These functions contribute to network analysis, data manipulation, and output formatting.

### Internal functions
1. **`BTW_SOURCE_TARGET_DIR`**  
   Calculates betweenness centrality for all nodes in the network, considering only **longitudinal shortest paths**.

2. **`Cross_corr_c`**  
   Uses a mixed-effects model framework to populate an adjacency matrix representing **cross-sectional** network connectivity between EMA items.

3. **`MIXED_CORR_MATRIX_c`**  
   Similar to `Cross_corr_c`, but computes both **cross-sectional** and **longitudinal** network connections.

4. **`jUpperTriMatToVec`**  
   Converts the upper triangular part of a matrix into a vector.

5. **`Contraire_jupperT`**  
   Reconstructs an upper triangular matrix from a vector.

6. **`Poseweight`**  
   Filters an adjacency matrix to retain only **positive** connections.

### External functions and toolboxes
7. **Brain Connectivity Toolbox (BCT)**  
   The scripts rely on functions from the Brain Connectivity Toolbox, a widely used toolbox for network analysis in neuroscience.  
   Download: [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home)

8. **`fdr_bh`**  
   Implements the Benjamini-Hochberg (1995) procedure to control the **false discovery rate (FDR)** in multiple hypothesis testing.  
   Download: [MATLAB Central - fdr_bh](https://ch.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh)

9. **`Prettyjson`**  
   Improves the readability and formatting of JSON structures.  
   Download: [MATLAB Central - prettyjson](https://ch.mathworks.com/matlabcentral/fileexchange/72667-prettyjson-m)

## Inputs
This section describes the files required for script execution. Please ensure that all necessary files are available in the `Inputs/` folder before running the scripts.

1. **`Matrix_1` (Multi-Layer Network)**  
   Adjacency matrix computed as an output of `MLN_Network_1.m`, representing both cross-sectional and longitudinal connectivity.

2. **`Coord_1` (Node 2D Coordinates)**  
   Two-dimensional spatial coordinates for each network node, obtained through Network Dimensionality Reduction applied to the cross-sectional adjacency matrix.

3. **`EMA_data` (Ecological Momentary Assessment Data)**  
   EMA data used throughout the analysis pipeline.

4. **`Names_3_v4` (Node names)**  
   Unique names or identifiers associated with each network node.

5. **`Cluster_HC`**  
   Clustering output generated by `BT_BDA.m`, where each cluster represents a bundle of psychological-contextual pathways sharing similar 3D trajectories.

### EMA items

| EMA Item Number | Question (ENG) | Question (FRA) | Name |
|---|---|---|---|
| 1 | Right now, I feel relaxed | En ce moment je me sens relaxé(e) | Lack_Relaxation |
| 2 | Right now, I feel lonely | En ce moment je me sens seul(e) | Loneliness |
| 3 | Right now, I feel anxious | En ce moment je me sens anxieux(se) | Anxiety |
| 4 | Right now, I feel happy, joyful | En ce moment je me sens content(e), joyeux(se) | Lack_Happiness |
| 5 | Right now, I feel irritated, angry | En ce moment je me sens irrité(e), en colère | Irritation |
| 6 | Right now, I feel troubled by sensory stimuli | En ce moment je me sens dérangé(e) par des stimulations sensorielles | Sensory_Issue |
| 7 | Right now, I feel excited | En ce moment je me sens excité(e) | Lack_Excitement |
| 8 | Right now, I feel sad | En ce moment je me sens triste | Sadness |
| 9 | Right now, I feel confident | En ce moment j’ai confiance en moi | Lack_Confidence |
| 10 | Right now, I feel like others don't like me | En ce moment j’ai l’impression que les autres ne m’aiment pas | Feeling_Rejected |
| 11 | Right now, I feel like I need to be cautious and that I'm not secure | En ce moment j’ai l’impression que je dois rester sur mes gardes, que je ne suis pas en sécurité | Feeling_Unsafe |
| 12 | Right now, I feel like my imagination is blending with reality | En ce moment j’ai l’impression que mon imagination se mélange avec la réalité | Confusion |
| 13 | Right now, I feel like I'm hearing or seeing things that others don't perceive | En ce moment j’ai l’impression d’entendre ou de voir des choses que les autres ne perçoivent pas | Hallucinations |
| 14 | Right now, I feel tired | En ce moment je me sens fatigué(e) | Feeling_Tired |
| 15 | Right now, I want to do a lot of things, I feel motivated | En ce moment j’ai envie de faire beaucoup de choses, je me sens motivé(e) | Lack_Motivation |
| 16 | Since the last beep, I've been physically active | Depuis le dernier bip, j’ai été physiquement actif(ve) | Lack_Physical_Activity |
| 17 | This activity is difficult | Cette activité est difficile | Finding_Activity_Difficult |
| 18 | I enjoy doing this activity | J’ai du plaisir à faire cette activité | Lack_Enjoing_Actvity |
| 19 | I'm focused on this activity | Je suis concentré sur cette activité | Lack_Concentration |
| 20 | Are you alone? | Es-tu seul? | Being_Alone |

## MLNetwork
The `MLNetwork/` folder contains the open-source code of the visualization tool used for 3D Multi-Layer Network representation.

This folder corresponds to the visualization platform referenced throughout the repository for interactive display and export of 3D-MLN networks.

## BT2
This folder is dedicated to the clustering component of the study:

> *Shedding light on the dynamic interplay of positive and negative symptoms of psychosis with Behavioral Tractography*

It contains the MATLAB script `BT2_clustering.m` together with the required input files.

The analysis pipeline is designed to provide a data-driven view of the correspondence between:
- multidimensional clinical patterns measured with the **SIPS interview**, and
- psychological phenomena measured with **EMA** in daily life.

The pipeline includes the following steps:
- Derivation of **SIPS-Dimension scores** using **Principal Component Analysis (PCA)**.
- Use of these SIPS dimensions in a **k-means clustering** procedure.
- Identification of participant subgroups that differentially express the multidimensional clinical patterns.

The outputs of the script include:
- PCA scores
- cluster assignments
- visualization of participants in PCA space with cluster labels

These outputs allow reproduction of the clustering results reported in the paper.

## Repository Authors
**Corrado Sandini**  
corrado.sandini@unige.ch

**Andrea Imparato**  
andrea.imparato@unige.ch
