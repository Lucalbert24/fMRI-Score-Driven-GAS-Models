# Score-Driven Models for fMRI BOLD Time Series

An empirical study of **Generalized Autoregressive Score (GAS)** models applied to resting-state **fMRI BOLD time series**, with particular emphasis on robustness to heavy-tailed and non-Gaussian observations.

This project was developed as part of the **Advanced Time Series** coursework at the **University of Bologna**.

## Overview

Functional Magnetic Resonance Imaging (fMRI) produces complex time-series data that may exhibit serial dependence, heteroskedasticity, outliers, and departures from Gaussianity.

This project investigates how **score-driven models** can adapt to these characteristics by comparing Gaussian and Student-t specifications on two BOLD time series with markedly different distributional properties.

The analysis includes:

- exploratory time-series analysis;
- descriptive statistics;
- normality diagnostics;
- ACF and PACF analysis;
- a score-driven mean model implemented from scratch;
- Gaussian and Student-t score updates;
- maximum-likelihood estimation;
- dynamic location and scale modelling using the `GAS` package in R;
- residual diagnostics;
- model comparison using AIC and BIC.

---

## Dataset

The analysis uses resting-state fMRI data from the **NKI-Rockland Sample**.

The dataset is organised as a four-dimensional array:

`70 ROIs × 404 time points × 24 subjects × 2 sessions`

where:

- 70 represents the number of brain Regions of Interest (ROIs);
- 404 represents the number of equally spaced time observations;
- 24 represents the number of subjects;
- 2 represents the scan-rescan sessions.

Only the first scan session is used in the analysis.

A systematic Shapiro-Wilk normality screening across all subjects and ROIs was used to identify two series with strongly different distributional behaviour.

The analysis focuses on **ROI 69**:

- **Subject 18 — ROI 69:** approximately Gaussian;
- **Subject 19 — ROI 69:** strongly non-Gaussian, heavy-tailed and negatively skewed.

The selected series provide a useful comparison for studying how the assumed conditional distribution affects score-driven parameter updates.

---

## Methodology

### 1. Exploratory Analysis

The selected BOLD time series are first analysed using:

- descriptive statistics;
- Shapiro-Wilk normality tests;
- Jarque-Bera tests;
- Q-Q plots;
- time-series plots;
- autocorrelation functions;
- partial autocorrelation functions.

The diagnostics show a strong contrast between the two selected series.

Subject 18 behaves approximately as a Gaussian series, while Subject 19 exhibits strong negative skewness, extremely high kurtosis and large outlying observations.

---

## 2. Manual Score-Driven Mean Model

A simplified score-driven model is implemented directly in R in order to illustrate the core mechanism of the GAS framework.

The conditional mean evolves according to:

\[
\mu_{t+1}
=
\alpha
+
\beta(\mu_t-\alpha)
+
\gamma u_t
\]

where \(u_t\) is the score of the conditional log-likelihood.

The model treats the scale parameter as constant and dynamically updates only the conditional mean.

Two score specifications are considered.

### Gaussian Score

Under the Gaussian specification, the score with respect to the conditional mean is proportional to the prediction error:

\[
u_t = y_t - \mu_t
\]

The Gaussian score is unbounded, meaning that extreme observations can have a large influence on the next parameter update.

### Student-t Score

Under the Student-t specification, the score is:

\[
u_t =
\frac{(\nu+1)(y_t-\mu_t)}
{\nu\sigma^2+(y_t-\mu_t)^2}
\]

where \(\nu\) represents the degrees of freedom.

Unlike the Gaussian score, the Student-t score limits the influence of very large prediction errors.

This provides the main robustness mechanism studied in the project.

The model parameters are estimated by numerical maximum likelihood using R's `nlminb()` optimisation function.

---

## 3. Full GAS Models

The analysis is then extended using the R package `GAS`.

The general GAS recursion is:

\[
\theta_{t+1}
=
\omega
+
A s_t
+
B\theta_t
\]

where \(s_t\) is a scaled score of the conditional log-likelihood.

Identity score scaling is used in this project.

Unlike the simplified implementation, the full GAS specification allows both:

- conditional location \(\mu_t\);
- conditional scale \(\sigma_t\);

to evolve dynamically over time.

Two conditional distributions are compared:

- **GAS-Normal**
- **GAS-Student-t**

The Student-t degrees-of-freedom parameter is estimated as a static parameter rather than being dynamically updated.

---

## Main Results

### Gaussian Series — Subject 18

Subject 18 is extremely close to Gaussian.

The exploratory analysis gives approximately:

- Shapiro-Wilk p-value: **0.9933**
- Jarque-Bera p-value: **0.9789**
- skewness: **-0.023**
- excess kurtosis: **0.004**

Both normality tests therefore fail to reject Gaussianity.

The manual Gaussian and Student-t score-driven models produce almost identical parameter dynamics.

In the Student-t specification, the estimated degrees of freedom are approximately:

\[
\hat{\nu} \approx 185
\]

A Student-t distribution with such a large number of degrees of freedom is effectively indistinguishable from a Gaussian distribution.

The information criteria for the full GAS models are:

| Model | Log-Likelihood | AIC | BIC |
|---|---:|---:|---:|
| GAS-Normal | -521.789 | **1055.58** | **1079.59** |
| GAS-Student-t | -521.845 | 1057.69 | 1085.70 |

The additional Student-t parameter therefore provides no meaningful improvement.

For this series, the **Gaussian GAS specification is sufficient**.

---

### Non-Gaussian Series — Subject 19

Subject 19 shows substantially different behaviour.

The descriptive statistics indicate:

- strong negative skewness;
- excess kurtosis of approximately **52.5**;
- Shapiro-Wilk p-value below **3.86 × 10⁻²⁴**;
- large extreme observations.

The Gaussian assumption is therefore strongly rejected.

### Manual Model

For the manual score-driven mean model:

| Model | Log-Likelihood | AIC |
|---|---:|---:|
| Gaussian score | -897.287 | 1802.574 |
| Student-t score | -731.915 | **1473.830** |

The Student-t model therefore provides a very large improvement in likelihood and information criteria.

The estimated degrees of freedom are:

\[
\hat{\nu} \approx 3.52
\]

which confirms strongly heavy-tailed behaviour.

### Full GAS Model

The same conclusion emerges from the package-based GAS models.

| Model | Log-Likelihood | AIC | BIC |
|---|---:|---:|---:|
| GAS-Normal | -710.154 | 1432.31 | 1456.32 |
| GAS-Student-t | -704.636 | **1423.27** | **1451.28** |

The estimated degrees of freedom are approximately:

\[
\hat{\nu} \approx 2.83
\]

confirming extremely heavy-tailed conditional behaviour.

A particularly clear difference appears in the estimated dynamic scale.

Under the GAS-Normal model, an extreme negative observation produces a very large spike in the filtered scale, reaching values above 1500.

The GAS-Student-t specification reacts much less aggressively because extreme observations receive lower weight through the Student-t score.

The resulting scale dynamics are therefore smoother and more robust.

---

## Residual Diagnostics

Residual adequacy is evaluated using Ljung-Box tests on both:

- standardised residuals;
- squared standardised residuals.

For the non-Gaussian series, both full GAS specifications remove the main serial dependence.

For the GAS-Student-t model, the Ljung-Box test on squared residuals produces a p-value of approximately:

\[
p \approx 0.97
\]

indicating no evidence of remaining autocorrelation in the conditional variance.

Residual diagnostics must nevertheless be considered together with the likelihood and information criteria.

Although both specifications can remove serial dependence, the Student-t model provides a substantially better distributional description of the heavy-tailed series.

---

## Key Findings

The analysis illustrates an important feature of score-driven models:

> The behaviour of the parameter update depends directly on the assumed conditional distribution.

For approximately Gaussian observations, the Gaussian score provides an adequate and parsimonious specification.

For strongly heavy-tailed observations, the Gaussian score can react excessively to extreme values.

The Student-t score instead automatically downweights large innovations, resulting in:

- more robust parameter updates;
- smoother filtered scale dynamics;
- higher likelihood;
- lower AIC and BIC;
- better handling of extreme observations.

The project therefore demonstrates the practical importance of choosing an appropriate conditional distribution when applying dynamic score-driven models.

---

## Repository Structure

    .
    ├── code/
    │   └── code-fMRI-ROI-time-series.Rmd
    │
    ├── data/
    │   └── data-fMRI-ROI-time-series.RData
    │
    ├── report/
    │   └── report-fMRI-ROI-time-series.pdf
    │
    └── README.md

---

## Files

### `code/code-fMRI-ROI-time-series.Rmd`

Complete R Markdown source containing:

- data loading and preparation;
- exploratory analysis;
- normality testing;
- ACF and PACF analysis;
- manual GAS model implementation;
- Gaussian and Student-t score functions;
- numerical maximum-likelihood estimation;
- GAS package model estimation;
- residual diagnostics;
- information criteria;
- figures and tables.

### `data/data-fMRI-ROI-time-series.RData`

Resting-state fMRI ROI time-series dataset used by the analysis.

### `report/report-fMRI-ROI-time-series.pdf`

Complete academic report containing the theoretical framework, methodology, empirical results, diagnostic analysis and conclusions.

---

## Software

The project was developed in **R**.

The main packages used include:

    library(tseries)
    library(stats)
    library(numDeriv)
    library(e1071)
    library(GAS)
    library(FinTS)

The full GAS models are primarily estimated using:

    UniGASSpec()
    UniGASFit()
    getFilteredParameters()
    residuals()

The simplified score-driven model is implemented manually and estimated using numerical maximum likelihood.

---

## Reproducing the Analysis

Clone the repository:

    git clone https://github.com/YOUR_USERNAME/fMRI-Score-Driven-GAS-Models.git

Move into the repository:

    cd fMRI-Score-Driven-GAS-Models

Open:

`code/code-fMRI-ROI-time-series.Rmd`

in RStudio.

Install any required packages that are not already available.

The dataset should be loaded from the `data` directory using the relative path:

    load("../data/data-fMRI-ROI-time-series.RData")

Using a relative path rather than an absolute machine-specific path allows the analysis to be reproduced after cloning the repository.

The R Markdown document can then be executed or knitted to reproduce the statistical analysis and report outputs.

---

## Authors

- Youssef Abdelsamie
- Luca Alberti
- Matteo Canton
- Gabriele Toro

**University of Bologna**  
Advanced Time Series  
2026

---

## References

The project builds primarily on the Generalized Autoregressive Score framework and its implementation in R.

Main references include:

- Ardia, D., Boudt, K., & Catania, L. (2019). *Generalized Autoregressive Score Models in R: The GAS Package*. Journal of Statistical Software, 88(6), 1–28.

- Creal, D., Koopman, S. J., & Lucas, A. (2013). *Generalized Autoregressive Score Models with Applications*. Journal of Applied Econometrics, 28(5), 777–795.

- Harvey, A. C. (2013). *Dynamic Models for Volatility and Heavy Tails*. Cambridge University Press.

- Ljung, G. M., & Box, G. E. P. (1978). *On a Measure of Lack of Fit in Time Series Models*. Biometrika, 65(2), 297–303.

- Engle, R. F. (1982). *Autoregressive Conditional Heteroscedasticity with Estimates of the Variance of United Kingdom Inflation*. Econometrica, 50(4), 987–1007.

- Shapiro, S. S., & Wilk, M. B. (1965). *An Analysis of Variance Test for Normality*. Biometrika, 52(3–4), 591–611.

- Jarque, C. M., & Bera, A. K. (1987). *A Test for Normality of Observations and Regression Residuals*. International Statistical Review, 55(2), 163–172.

The complete bibliography is available in the project report.

---

## Academic Context

This project was originally developed as coursework for the **Advanced Time Series** course at the **University of Bologna**.

The objective was to select Gaussian and non-Gaussian fMRI time series, estimate appropriate score-driven models, compare them with Gaussian specifications and evaluate the resulting model behaviour and diagnostics.

---

## Disclaimer

This repository contains an academic statistical modelling project.

The fMRI data are used to study statistical time-series modelling techniques and the robustness properties of Generalized Autoregressive Score models.

The results should not be interpreted as medical, neurological, psychiatric or clinical conclusions.

The dataset is organised as a four-dimensional array:

```text
70 ROIs × 404 time points × 24 subjects × 2 sessions
