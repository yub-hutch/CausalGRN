<table>
  <tr>
    <td><img src="man/figures/Logo.png" alt="CausalGRN Logo" width="300"/></td>
    <td><h1>CausalGRN: Deciphering Causal Gene Regulatory Networks from Single-Cell CRISPR Screens</h1></td>
  </tr>
</table>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**CausalGRN** is a scalable computational framework that infers causal gene regulatory networks (GRNs) and predicts cellular responses to unseen perturbations. It is designed to translate the complex outputs from large-scale single-cell CRISPR screens with scRNA-seq readouts into reliable causal insights.

## Overview

Large-scale single-cell CRISPR screens provide critical data to map causal GRNs. However, analyzing this data to extract reliable causal relationships is a major challenge. CausalGRN addresses this by:

1.  **Mitigating Spurious Partial Correlations**: It employs a novel adaptive thresholding correction to reduce the impact of pervasive spurious partial correlations found in sparse scRNA-seq data, enabling a more robust inference of the network's undirected skeleton.
2.  **Orienting the Network**: It orients the graph using observed perturbation outcomes from CRISPR screens.
3.  **Predicting Perturbation Effects**: The resulting directed GRN can be used to predict the downstream effects of novel, unseen perturbations via network propagation.

Across both simulations and diverse experimental datasets, CausalGRN substantially outperforms existing approaches in network reconstruction accuracy and in predicting the effects of unseen perturbations.

## Installation

You can install the development version of CausalGRN from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("yub-hutch/CausalGRN")
```

## Example Usage

Here is a toy example demonstrating a complete CausalGRN workflow. We will:
1.  Simulate a dataset with wild-type cells and three different gene knockdowns.
2.  Use a subset of the data (wild-type and one knockdown) to infer a causal GRN and compare the result to other common methods.
3.  Use the inferred GRN to build a predictive model.
4.  Use the model to predict the effects of the other two knockdowns that were not used for training.

### 0. Setup and Data Simulation

First, let's load the necessary libraries and define our simulation parameters. We will create a dataset with wild-type (WT) cells and cells with knockdowns for genes 'A', 'B', and 'C' based on the ground truth network `A -> B -> C`.

```r
# --- 0. SETUP & DATA SIMULATION ---
library(dplyr)
library(igraph)
library(CausalGRN)

# Define all simulation parameters upfront
a <- b <- 1
sd <- 2
s <- -4
nwt <- 1e5
npt <- 1e4
knockdown_efficacy <- 0.9 # 90% knockdown efficiency

# Define the ground truth graph for this simulation
ground_truth <- igraph::make_graph(~ A -+ B, B -+ C)
E(ground_truth)$weight <- c(a, b)

# --- Simulate Data ---
set.seed(123)

# Generate Wild-Type (WT) Data
x_latent_wt <- rnorm(nwt, 0, sd)
y_latent_wt <- rnorm(nwt, a * x_latent_wt + s, sd)
z_latent_wt <- rnorm(nwt, b * (y_latent_wt - s), sd)
wt_counts <- cbind(
  A = rpois(nwt, exp(x_latent_wt)),
  B = rpois(nwt, exp(y_latent_wt)),
  C = rpois(nwt, exp(z_latent_wt))
)

# Generate Perturb-A Data
x_latent_kdA <- rnorm(npt, 0, sd) + log(1 - knockdown_efficacy)
y_latent_kdA <- rnorm(npt, a * x_latent_kdA + s, sd)
z_latent_kdA <- rnorm(npt, b * (y_latent_kdA - s), sd)
kdA_counts <- cbind(
  A = rpois(npt, exp(x_latent_kdA)),
  B = rpois(npt, exp(y_latent_kdA)),
  C = rpois(npt, exp(z_latent_kdA))
)

# Generate Perturb-B Data
x_latent_kdB <- rnorm(npt, 0, sd)
y_latent_kdB <- rnorm(npt, a * x_latent_kdB + s, sd) + log(1 - knockdown_efficacy)
z_latent_kdB <- rnorm(npt, b * (y_latent_kdB - s), sd)
kdB_counts <- cbind(
  A = rpois(npt, exp(x_latent_kdB)),
  B = rpois(npt, exp(y_latent_kdB)),
  C = rpois(npt, exp(z_latent_kdB))
)

# Generate Perturb-C Data
x_latent_kdC <- rnorm(npt, 0, sd)
y_latent_kdC <- rnorm(npt, a * x_latent_kdC + s, sd)
z_latent_kdC <- rnorm(npt, b * (y_latent_kdC - s), sd) + log(1 - knockdown_efficacy)
kdC_counts <- cbind(
  A = rpois(npt, exp(x_latent_kdC)),
  B = rpois(npt, exp(y_latent_kdC)),
  C = rpois(npt, exp(z_latent_kdC))
)

# --- Combine and Prepare Inputs ---
count <- rbind(wt_counts, kdA_counts, kdB_counts, kdC_counts)
group <- c(rep('WT', nwt), rep('A', npt), rep('B', npt), rep('C', npt))
Y <- scale(log1p(count), center = TRUE, scale = TRUE)
```

Now that we have our simulated data, we can proceed with the CausalGRN workflow.

### 1. Inferring a Causal GRN

For network inference, we will pretend we only have access to the wild-type data and the data from the perturbation of gene 'A'. We will use this subset to infer the network structure with `CausalGRN`.

```r
# --- 1. Infer GRN from a subset of data (WT and Perturb-A) ---
cat("Running CausalGRN to infer the network...\n")
train_idx <- which(group %in% c('WT', 'A'))

# Infer the graph using the training data
skel <- infer_skeleton(count[train_idx, ], Y[train_idx, ], alpha = 0.05, min_abspcor = 0, ncores = 1)
stat <- calc_perturbation_effect(Y[train_idx, ], group[train_idx], ncores = 1)
inferred_graph <- infer_causalgrn(skel$graph, stat, alpha = 0.05, max_order = 2)

cat("CausalGRN inferred graph:\n")
print(inferred_graph)
```

The plot below shows how the output of CausalGRN compares to other GRN inference methods. For the simple chain `A -> B -> C`, CausalGRN correctly identifies the causal structure, while other methods may infer incorrect edges or directions.

<p align="center">
  <img src="man/figures/illustration.png" alt="CausalGRN Illustration" width="500"/>
</p>

<details>
<summary>Click to see the code for running other common GRN inference methods for comparison</summary>

Note that the following methods are not part of the core CausalGRN algorithm but are included in the package as wrappers for convenient benchmarking. Here is how you could run them on the same dataset.

```r
# --- Compare with other common GRN inference methods ---
cat("Running other GRN methods for comparison...\n")

# Prepare data subsets for baseline methods
# Observational methods are typically run on wild-type data
wt <- Y[group == 'WT', ]
# Interventional methods can use a list of perturbed data
pts <- list(
  A = Y[group == 'A', ]
)

# Observational method: PC Algorithm (on WT data)
graph_pc <- run_pc(wt, alpha = 0.05)

# Observational method: GES (on WT data)
graph_ges <- run_ges(wt)

# Interventional method: GIES (uses perturbation data)
graph_gies <- run_gies(wt, pts)

# Observational method: Lasso (on WT data)
graph_lasso_wt <- run_lasso(wt, ncores = 1)

# Interventional method: Lasso (on all data)
graph_lasso_all <- run_lasso(Y, ncores = 1)

# Observational method: GENIE3 (on WT data)
graph_genie3 <- run_genie3(wt, ncores = 1)

# GRNBoost2 (requires Python environment)
# For this example, we assume pre-computed results
# graph_grnboost2 <- ...
```
</details>

### 2. Predicting Effects of Unseen Perturbations

Now, using the `inferred_graph` from the previous step, we will fit an expression model. We will train the model on the same data we used for inference (WT and Perturb-A) and then predict the effects for the held-out, unseen perturbations of 'B' and 'C'.

```r
# --- 2. Fit model and predict effects for unseen perturbations ---
cat("Fitting expression model...\n")
B_fit <- fit_expression_model(
  Y[train_idx, ],
  group[train_idx],
  graph = inferred_graph,
  ncores = 1,
  method = 'lm'
)

# --- Predict effects for B and C knockdown ---
cat("Predicting effects for held-out perturbations...\n")
# We need the mean expression of the perturbed gene in the knockdown cells,
# and the mean expression in the WT cells.
wt_expressions <- colMeans(Y[group == 'WT', ])
knockdown_expressions <- c(
  'B' = mean(Y[group == 'B', 'B']),
  'C' = mean(Y[group == 'C', 'C'])
)

# Predict the delta (change from WT) for all genes
pred_effects <- predict_standard_effect(B_fit, knockdown_expressions, wt_expressions)

cat("Predicted effects matrix:\n")
print(pred_effects)


# --- 3. Visualize the predicted effects ---
# Convert matrix to data frame for plotting and create the heatmap
plot_df <- as.data.frame(as.table(pred_effects))
names(plot_df) <- c("Perturbation", "Gene", "Effect")

ggplot2::ggplot(plot_df, ggplot2::aes(x = Gene, y = Perturbation, fill = Effect)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Predicted Effects of Gene Knockdowns", fill = "Effect")
```

The plot below is the output of this prediction example. It shows the predicted effects of knocking down genes 'B' and 'C'. As expected from the ground truth graph `A -> B -> C`, knocking down 'B' affects both 'B' and 'C', while knocking down 'C' only affects 'C'.

<p align="center">
  <img src="man/figures/Predicted_effect.png" alt="Predicted Effect Plot" width="400"/>
</p>

## Citation

If you use CausalGRN in your research, please cite our paper:

> **CausalGRN: deciphering causal gene regulatory networks from single-cell CRISPR screens**
>
> (Further citation details will be added upon publication)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
