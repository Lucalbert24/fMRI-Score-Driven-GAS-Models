# Score-Driven Models for fMRI BOLD Time Series

An empirical study of **Generalized Autoregressive Score (GAS)** models applied to resting-state **fMRI BOLD time series**, with particular emphasis on robustness to heavy-tailed and non-Gaussian observations.

This project was developed as part of the **Advanced Time Series** coursework at the **University of Bologna**.

## Overview

Functional Magnetic Resonance Imaging (fMRI) produces complex time-series data that may exhibit serial dependence, heteroskedasticity, outliers, and departures from Gaussianity.

This project investigates how **score-driven models** can adapt to these characteristics by comparing Gaussian and Student-t specifications on two BOLD time series with very different distributional properties.

The analysis combines:

- exploratory time-series analysis;
- normality diagnostics;
- ACF and PACF analysis;
- a GAS mean model implemented from scratch;
- Gaussian and Student-t score updates;
- maximum-likelihood estimation;
- dynamic location and scale modelling with the `GAS` R package;
- residual diagnostics;
- AIC and BIC model comparison.

## Dataset

The analysis uses resting-state fMRI data from the **NKI-Rockland Sample**.

The dataset is organised as a four-dimensional array:

```text
70 ROIs × 404 time points × 24 subjects × 2 sessions
