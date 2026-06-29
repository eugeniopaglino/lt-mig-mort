# Relational Mortality Models

This repository provides tools for fitting Generalized Relational Additive Models (GARM) for mortality models using a flexible Generalized Additive Model (GAM) approach.

The primary utility of these scripts is to convert abridged, sparse, or noisy empirical mortality data into a smooth, continuous, single-year mortality schedule. It achieves this by anchoring the empirical data to a known standard demographic schedule and using a penalized spline to smooth the localized deviations.

---

## The Shared Framework: Likelihood and Identifiability

We model discrete death counts ($D_x$) occurring over a volume of exposure ($E_x$), expressed in person-years.

By modeling counts directly, we naturally handle age intervals with zero deaths, and we correctly weight observations based on their population size.

$${}_nD_x \sim \text{NB}({}_n\lambda_x, \theta)$$ 
$$\log({}_n\lambda_x) = \log(\mu_a) +\log({}_nE_x)$$
$$\log(\mu_a) = \alpha + \beta \log(m_a^{\text{std}}) + s(a)$$

Models of this form face a fundamental identifiability challenge: the standard mortality schedule and the smoothing spline are highly collinear over age $x$.

We solve this using a single, robust procedure: by applying an additional shrinkage penalty to the null space of the spline (`select = TRUE` in `mgcv::gam`), the model performs data-driven variable selection. The optimizer relies on the standard curve for the global trend. If the flexible spline $s(a)$ is not necessary to explain the data, the penalty seamlessly shrinks its terms to zero, preventing it from "stealing" variance from the demographic standard. 

---

## Core Functions

### Mortality Modeling (`mortality_models_counts.R`)
The core of this approach is `fit_abridged_sim_lt()`.

  * **`singularize_single_step()`**: The internal engine that fits the GAM and generates posterior simulations using multivariate normal draws from the estimated coefficients.
  * **`fit_abridged_sim_lt()`**: The primary user-facing function. It takes empirical data and standard schedules, fits the model, and performs simulation-based inference. This function returns a list containing two distinct summaries:
        
    1.  `abridged_summary`: A life table summary following your specified `max_abridged_age` (e.g., 85+).
    2.  `single_year_summary`: A continuous, granular life table summary extending to age 110+.

#### Handling Data Artifacts
The model implements a `skip_ages` argument within the modeling functions. Users can pass a vector `c(min, max)` or a list of vectors `list(c(0, 5), c(20, 25))` to flag localized reporting artifacts or age-heaping. The models treat these bands as missing data and interpolate over them using the standard curve and the spline.

---

## Testing and Examples

We provide testing scripts to validate model performance and provide usage templates:

**`test_count_models.R`**: Demonstrates the comprehensive modeling and evaluation workflow.

  * **Integrated Visualization:** The script evaluates model fit by plotting the **continuous modeled trend** (mean line with 95% credible interval ribbons) directly against the **empirical abridged data** (points). This allows for immediate visual comparison between the modeled single-year schedule and the input abridged data.
  * **Stress-Testing:** Includes tests for sparse data (injecting structural zeros) and data requiring the `skip_ages` argument to interpolate over artifacts.
  * **Application:** Includes a demonstration using Russia 2000 data.

**`simple_usage_script.R`**: A streamlined template for basic usage.

  * **Simplified Workflow:** Provides a minimal, end-to-end example of how to format empirical data midpoints, apply the standard schedule, fit the GAM model, and directly extract and export the finalized abridged life table summary. This is ideal for quickly processing data without the overhead of diagnostic plotting.

---

## Utility Functions: Classical Life Table Construction

**File:** `lt_functions.R`

This file provides two functions for computing standard life table columns ($q_x, l_x, d_x, L_x, T_x, e_x$):

  * **`lt_abridged(nMx)`**: Constructs a life table for classical abridged age groups (0, 1-4, 5-9, 10-14, ..., 90+).
  * **`lt_single(Mx)`**: Constructs a life table for complete, single-year age groups (0, 1, 2, ..., 100+).

*Note: These utilities use standard demographic $a_x$ approximations (e.g., 0.5 for most ages, specific adjustments for infant mortality). For highly specific regional datasets, manual adjustment of $a_x$ vectors may be required.*