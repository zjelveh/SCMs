# SCMs R Package - Project Context

## 🚨 CRITICAL DEVELOPMENT PRINCIPLE - HIGHEST PRIORITY 🚨

**FAIL HARD, NOT SOFT**
- **NO FALLBACKS**: Never implement fallback mechanisms or graceful degradation without explicit user approval
- **NO APPROXIMATIONS**: Never substitute approximations for exact functionality (e.g., correlation-based "interactions" when user requests true SHAP interactions)
- **HARD FAILURES ONLY**: When functionality cannot be implemented as requested, fail with clear error messages
- **SOLVE THE ROOT CAUSE**: Fix underlying issues rather than working around them
- **ASK PERMISSION**: If exact functionality is impossible, explicitly ask user for permission to use alternatives
- **PRIORITY**: Maximizing long-term robustness over short-term functionality
- **NO BACKWARD COMPATIBILITY**: Do not write backward compatible code. This repo has not been released so we don't need to worry about breaking changes. Remove any existing backward compatibility code.

## Project Overview
This is an R package for Synthetic Control Methods (SCM) with advanced analytics including specification curve analysis, CatBoost+SHAP integration, and comprehensive inference methods.
The highlight of the package is the specification curve where the feature for each specification are colored by their shapley value. 
The shapley value represents how much that feature is relevant to a predicted treatment effect.

## Code Review Focus Areas
- **Code Organization**: Ensure functions are in logical files, remove redundancy
- **Efficiency**: Identify verbose or inefficient code patterns
- **Refactoring**: Explore whether it makes sense to break supporting_functions.R into two or more separate files?
- **Documentation**: Verify docstrings match function behavior, update man pages
- **Package Integrity**: Ensure DESCRIPTION dependencies are accurate and minimal
- **File Management**: Remove debugging/temporary files and unused functions
- **Backward Compatibility Removal**: Identify and remove all backward compatibility code, fallback mechanisms, and legacy support since this package has not been released

## Key Package Structure
- **Core SC Methods**: `scest.R`, `estimate_sc.R`, `inference_*.R`
- **Data Processing**: `create_scm_dataset.R`, `covariate_processing.R`, `data_processing.R`
- **Specification Curve**: `spec_curve.R`, `spec_curve_analysis.R`
- **Machine Learning**: `catboost_shap_analysis.R` (note: moved from xgboost)
- **Visualization**: `plot_spec_curve.R`, `visualization_advanced.R`
- **Supporting**: `supporting_functions.R`, `outcome_models.R`

## Documentation Standards
- All exported functions must have complete roxygen2 documentation
- Examples should be runnable and demonstrate key functionality
- Parameter descriptions should be precise and include expected data types

## Dependencies to Verify
- Remove unused imports from DESCRIPTION
- Ensure all used packages are listed
- Check version requirements are appropriate

## SCPI Compatibility Requirements
- All synthetic control implementations must maintain compatibility with SCPI methodology
- Use C matrix approach for constant terms: `min ||A - Z*w - C*γ||_V^2` 
- Never implement row-based matrix modifications for constants
- Optimization framework must handle C matrix alongside A and Z matrices

## Specification Curve Architecture
- SHAP values must align perfectly with specifications via full_spec_id
- P-values and statistical inference must reflect only displayed/filtered specifications
- No internal approximations - require external SHAP computation or explicit user permission
- Filtering must occur BEFORE p-value calculations to ensure statistical validity

## Error Handling Standards
- All validation must fail immediately with specific error messages
- Error messages must point to root causes and provide resolution paths
- Never mask underlying issues with fallback mechanisms
- Validation errors must include required vs. actual data structure information
- Missing required functionality must stop execution with clear user instructions

## JSS Submission Roadmap - Efficient Ordering Strategy

### Phase 1: Foundation (Items 1-3)
1. **CRAN Checks First** - Reveals fundamental structural issues before investment
2. **Input Validation** - Establishes error handling foundation affecting API design  
3. **S3 Architecture** - Design user-facing API before documenting it

### Phase 2: Documentation (Items 4-5)
4. **Package-Level Documentation** - Forces design thinking, reveals API issues
5. **Function Documentation** - Document stable API without interface change risk

### Phase 3: Quality Assurance (Items 6-7)  
6. **Unit Tests** - Test stable interfaces to avoid rework
7. **Regression Tests** - Ensure consistent behavior on settled codebase

### Phase 4: Performance (Items 8-9)
8. **Profile & Optimize** - Optimize stable code, not changing architecture
9. **Memory Management** - Build on optimized, stable computational core

### Phase 5: User Experience (Items 10-11)
10. **Vignettes** - Teach workflow on mature, stable package
11. **Tidy Methods** - Polish ecosystem integration on solid foundation

### Phase 6: Presentation (Items 12-14)
12. **pkgdown Website** - Showcase mature package capabilities
13. **README** - Present stable, polished package to users
14. **Reproducible Examples** - Demonstrate full capabilities of finished product

**Key Principle**: Minimize rework by establishing stable foundations before building higher-level features.

## Build Commands
- Only run these when necessary.
- `devtools::document()` - Update documentation
- `devtools::check()` - Package checks
- `devtools::install()` - Install package

IMPORTANT: 
- Prioritize 100% success versus graceful failure. Avoid graceful failure if it conflicts the goal of 100% success. 
- Do not circumvent issues just to make the code run. Solve the issue.
- Do not make local fixes if a global fix would lead to a more robust codebase
- **CRITICAL**: Do not implement approximations or workarounds without explicit user approval. When the user requests specific functionality (e.g., "true SHAP interactions"), only implement that exact functionality or explicitly ask for permission to use alternatives. Do not silently substitute correlation-based approximations for true interaction effects.


