# S3 Pipe-Friendly Design Implementation Summary

## Overview
The brsm package now supports pipe-friendly function composition through S3 method dispatching. This enables users to work with `brsm_fit` objects directly without needing to explicitly manage `factor_names` and other metadata.

## Architecture

### Core Design Pattern
Each main analysis function now follows this pattern:

```r
# 1. Generic function with S3 dispatch
function_name <- function(object, ...) {
  UseMethod("function_name")
}

# 2. S3 method for brsm_fit objects
function_name.brsm_fit <- function(object, factor_names = NULL, ...) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names  # Extract from metadata
  }
  function_name.default(
    as_brsm_draws(object), 
    factor_names = factor_names, 
    ...
  )
}

# 3. S3 method for data frames (backward compatibility)
function_name.default <- function(object, factor_names = NULL, ...) {
  # Original implementation
}
```

### Benefits
- **Metadata Extraction**: `factor_names` and other parameters are automatically extracted from `brsm_fit` objects
- **Pipe-Friendly**: Enables natural composition: `fit |> congruence_parameters() |> rope_congruence()`
- **Backward Compatible**: Existing code using data frames + explicit `factor_names` continues to work
- **Type Safe**: Functions validate inputs and provide clear error messages

## Implemented Functions

### 1. **congruence_parameters()**
Computes RSA congruence parameters (a1-a5) from posterior draws.

**S3 Methods:**
- `.brsm_fit()`: Extracts factor_names, calls default method
- `.brmsfit()`: Direct brmsfit support for legacy compatibility
- `.default()`: Data frame input with explicit factor_names

**Example:**
```r
fit <- fit_brsm(data, response = "y", factor_names = c("x1", "x2"))
params <- congruence_parameters(fit)  # No need for factor_names!
```

### 2. **stationary_point()**
Identifies and characterizes stationary points of the response surface.

**S3 Methods:**
- `.brsm_fit()`: Extract metadata and dispatch to default method
- `.default()`: Core implementation with validation

**Example:**
```r
sp <- fit |> stationary_point()
```

### 3. **classify_stationarity_point()**
Classifies stationary points (maximum, minimum, saddle point, etc.).

**S3 Methods:**
- `.brsm_fit()`: Metadata extraction
- `.default()`: Classification logic

**Example:**
```r
classified <- fit |> classify_stationarity_point()
```

### 4. **posterior_ridge_analysis()**
Performs Bayesian ridge analysis at specified radii.

**S3 Methods:**
- `.brsm_fit()`: Metadata extraction
- `.default()`: Ridge analysis computation

**Example:**
```r
ridge <- fit |> posterior_ridge_analysis(radii = seq(0, 3, 0.5))
```

### 5. **credible_optimum_region()**
Computes credible regions around the optimum.

**S3 Methods:**
- `.brsm_fit()`: Metadata extraction
- `.default()`: Credible region computation

**Example:**
```r
region <- fit |> credible_optimum_region(probs = c(0.025, 0.975))
```

### 6. **steepest_ascent()**
Computes steepest ascent paths from posterior draws.

**S3 Methods:**
- `.brsm_fit()`: Metadata extraction with default start = c(x1=0, x2=0)
- `.default()`: Path computation

**Example:**
```r
path <- fit |> steepest_ascent(start = c(x1=0, x2=0), n_steps = 20)
```

### 7. **as_brsm_draws()**
Helper function to convert various input types to standardized draws format.

**S3 Methods:**
- `.brsm_fit()`: Extracts draws from the internal brmsfit object
- `.brmsfit()`: Converts brmsfit to data frame with brsm column naming
- `.default()`: Direct data frame input with column validation

**Example:**
```r
draws <- as_brsm_draws(fit)  # Automatic format conversion
```

### 8. **Auxiliary Functions**
These functions leverage the S3 dispatch of core functions:

- **rope_congruence()**: Uses `congruence_parameters()` internally
- **summarize_congruence()**: Uses `congruence_parameters()` internally

```r
# These work automatically through S3 dispatch of congruence_parameters()
rope <- fit |> rope_congruence(rope = c(-0.1, 0.1))
summary <- fit |> summarize_congruence()
```

## Usage Examples

### Basic Workflow
```r
library(brsm)

# Fit a model
fit <- fit_brsm(
  data = my_data,
  response = "y",
  factor_names = c("x1", "x2")
)

# Use pipe-friendly composition
analysis <- fit |>
  congruence_parameters() |>   # S3 dispatch extracts factor_names
  rope_congruence(rope = c(-0.1, 0.1))

stationarity <- fit |> 
  stationary_point() |>
  classify_stationarity_point()

ridge <- fit |> posterior_ridge_analysis()
```

### Backward Compatibility
```r
# Old style still works
draws <- as.data.frame(fit$fit)  # Get raw draws
params <- congruence_parameters(
  draws, 
  factor_names = c("x1", "x2")  # Explicit specification
)
```

### With Data Preparation
```r
# Prepare data with coding
prepared <- prepare_brsm_data(
  data,
  factor_names = c("x1", "x2"),
  method = "zscore"
)

# Fit and analyze
fit <- fit_brsm(prepared, response = "y", factor_names = c("x1", "x2"))
params <- fit |> congruence_parameters()  # Still works!
```

## Implementation Details

### S3 Method Registration
All S3 methods are properly exported via `@export` roxygen tags:
```r
#' @rdname function_name
#' @export
function_name.brsm_fit <- function(...) { ... }
```

### Error Handling
- **Missing factor_names**: When data frame is passed without `factor_names`, error message is clear
- **Invalid brsm_fit**: Validation ensures `$factor_names` is properly structured
- **Type mismatches**: S3 dispatch provides sensible defaults

### Performance Considerations
- **No overhead**: S3 dispatch adds negligible performance cost
- **Lazy evaluation**: Code only evaluates what's needed
- **Memory efficient**: Maintains existing data structures

## Testing

Comprehensive test suite included in `tests/testthat/`:

- **test-congruence_parameters.R**: Tests for all S3 methods and variations
- **test-integration_workflow.R**: Full pipeline tests
- **test-fit_brsm_congruence.R**: Constrained model fits
- **test_s3_dispatch.R**: Explicit S3 dispatch verification

**Running tests:**
```bash
# Local environment
devtools::test()

# With brms fitting tests enabled
BRSM_RUN_BRMS_TESTS=true devtools::test()
```

## Design Rationale

### Why S3?
- **Simplicity**: Minimal learning curve for R developers
- **Flexibility**: Easy to extend with new object types
- **Package compatibility**: Works seamlessly with brms and tidyverse
- **Documentation**: Roxygen handles both generic and methods documentation

### Optimization Strategy
1. **Metadata Extraction**: `brsm_fit` objects carry all necessary context
2. **Lazy Defaults**: `factor_names = NULL` triggers extraction, not override
3. **Named Defaults**: Steepest ascent uses `start = NULL` to detect auto-defaults
4. **Backward Compatibility**: Old code continues to work unmodified

## Migration Guide

### For Existing Code
```r
# Old style (still works)
params <- congruence_parameters(draws, factor_names = c("x1", "x2"))

# New style (recommended)
fit <- fit_brsm(...)
params <- fit |> congruence_parameters()
```

### Best Practices
1. Use `fit_brsm()` to create a `brsm_fit` object
2. Pass the fit object directly to analysis functions
3. Let S3 dispatch handle metadata management
4. Use pipes for readable composition

## Future Extensions

The S3 architecture enables:
- Adding custom `brsm_fit` subclasses for specialized workflows
- Supporting additional input types (e.g., `cmdstanr` models)
- Integrating with other Bayesian analysis packages
- Creating domain-specific wrappers

## See Also

- [fit_brsm()] for creating brsm_fit objects
- [prepare_brsm_data()] for data preprocessing
- [as_brsm_draws()] for draws conversion
- [pipe] for pipe operator documentation
