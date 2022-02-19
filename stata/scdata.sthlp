{smcl}
{* *!version 0.1 2022-01-25}{...}
{viewerjumpto "Syntax" "scdata##syntax"}{...}
{viewerjumpto "Description" "scdata##description"}{...}
{viewerjumpto "Options" "scdata##options"}{...}
{viewerjumpto "Examples" "scdata##examples"}{...}
{viewerjumpto "Stored results" "scest##stored_results"}{...}
{viewerjumpto "References" "scdata##references"}{...}
{viewerjumpto "Authors" "scdata##authors"}{...}

{title:Title}

{p 4 8}{cmd:scdata} {hline 2} Data Preparation for Synthetic Control Methods.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:scdata } {it:features} {ifin} 
{cmd:,} 
{cmd:id(}{it:idvar}{cmd:)} 
{cmd:time(}{it:timevar}{cmd:)}
{cmd:outcome(}{it:outcomevar}{cmd:)}
{cmd:treatment(}{it:treatmentvar}{cmd:)}  
{cmd:dfname(}{it:string}{cmd:)}
[{cmd:covadj(}{it:string}{cmd:)}
{cmd:anticipation(}{it:#}{cmd:)}
{cmd:cointegrated}
{cmd:constant}]{p_end}

{synoptset 28 tabbed}{...}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:scdata} prepares the data to be used by {help scest:scest} or {help scpi:scpi} to implement estimation and inference procedures for Synthetic Control (SC) methods. 
It allows the user to specify the outcome variable, the features of the treated unit to be matched, and covariate-adjustment feature by feature. 
The command follows the terminology proposed in {browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Cattaneo, Feng, and Titiunik (2021)}. 
The command is a wrapper of the companion Python package. As such, the user needs to have a running version of Python with the package installed. A tutorial on how to install Python and link it to Stata
can be found {browse "https://nppackages.github.io/scpi/":here}.{p_end}

{p 8 8} Companion {browse "www.r-project.org":R} and {browse "https://www.python.org/":Python} packages are described in 
{browse "https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2022_scpi.pdf":Cattaneo, Feng, Palomba and Titiunik (2022)}{p_end}

{p 8 8} Companion commands are: {help scest:scest} for point estimation, {help scpi:scpi} for inference procedures, and {help scplot:scplot} for SC plots.{p_end}

{p 4 8}Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:{p_end}

{p 8 8}{browse "https://nppackages.github.io/scpi/":https://nppackages.github.io/scpi/}{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Variables}

{p 4 8}{cmd:id(}{it:idvar}{cmd:)} specifies the variable containing the identifier for each unit.{p_end}

{p 4 8}{cmd:time(}{it:timevar}{cmd:)} specifies the variable containing the time period of each observation.{p_end}

{p 4 8}{cmd:outcome(}{it:outcomevar}{cmd:)} specifies the outcome variable of interest. Note that {it:outcomevar} may not be among the {it:features} specified.{p_end}

{p 4 8}{cmd:treatment(}{it:treatmentvar}{cmd:)} specifies the treatment indicator.{p_end}

{dlgtab:Estimand}

{p 4 8}{cmd:covadj(}{it:#}{cmd:)} specifies the variable to be used for adjustment for each feature. If the user wants
     to specify the same set of covariates for all features, a string should be provided according to the following format: {opt covadj("['cov1', 'cov2']")}. 
     If instead a different set of covariates per feature has to be specified, then the following format should be used {opt covadj("[['cov1', 'cov2'],['cov1','cov3']]")}. 
     Note that in this latter case the number of sub-lists must be equal to the number of {it:features}. Moreover, the order of the 
     sub-lists matters, in the sense that the first sub-list is interpreted as the set of covariates for the first 
     feature, and so on. Finally, the user can specify 'constant' and 'trend' as covariates even if they are not 
     present in the loaded dataset.

{p 4 8}{cmd:anticipation(}{it:#}{cmd:)} specifies the number of periods of potential anticipation effects. Default is {cmd:anticipation(0)}.{p_end}

{p 4 8}{cmd:cointegrated} if specified indicates that there is a belief the features form a cointegrated system. {p_end}

{p 4 8}{cmd:constant} if specified includes a constant term across features.{p_end}

{dlgtab:Others}

{p 4 8}{cmd:dfname(}{it:string}{cmd:)} specifies the name of the Python object that is saved and that will be passed to {help scest:scest} or {help scpi:scpi}.{p_end}

    {hline}


{marker examples}{...}
{title:Example: Cattaneo, Feng and Titiunik (2021) Germany Data}

{p 4 8}Setup{p_end}
{p 8 8}{cmd:. use scpi_germany.dta}{p_end}

{p 4 8}Prepare data{p_end}
{p 8 8}{cmd:. scdata gdp, dfname("python_scdata") id(country) outcome(gdp) time(year) treatment(status) cointegrated}{p_end}


{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:scdata} stores the following in {cmd:e()}:

{synoptset 25 tabbed}{...}

{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(J)}}number of donors{p_end}
{synopt:{cmd:e(KM)}}total number of covariates used for adjustment{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(features)}}name of features{p_end}
{synopt:{cmd:e(outcomevar)}}name of outcome variable{p_end}
{synopt:{cmd:e(constant)}}logical indicating the presence of a common constant across features{p_end}
{synopt:{cmd:e(cointegrated_data)}}logical indicating cointegration{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(A)}}pre-treatment features of the treated unit{p_end}
{synopt:{cmd:e(B)}}pre-treatment features of the control units{p_end}
{synopt:{cmd:e(C)}}covariates used for adjustment{p_end}
{synopt:{cmd:e(P)}}predictor matrix{p_end}

{marker references}{...}
{title:References}

{p 4 8}Cattaneo, M. D., Feng, Y., and Titiunik, R. 2021. 
{browse "https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf":Prediction intervals for synthetic control methods}. 
{it:Journal of the American Statistical Association}, 116(536), 1865-1880.{p_end}

{p 4 8}Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. 2022. 
{browse "https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2022_scpi.pdf":scpi - Uncertainty Quantification for Synthetic Control Estimators}.{p_end}

{marker authors}{...}
{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University, Princeton, NJ.
{browse "mailto:cattaneo@princeton.edu":cattaneo@princeton.edu}.{p_end}

{p 4 8}Yingjie Feng, Tsinghua University, Beijing, China.
{browse "mailto:fengyj@sem.tsinghua.edu.cn":fengyj@sem.tsinghua.edu.cn}.

{p 4 8}Filippo Palomba, Princeton University, Princeton, NJ.
{browse "mailto:fpalomba@princeton.edu":fpalomba@princeton.edu}.

{p 4 8}Rocio Titiunik, Princeton University, Princeton, NJ.
{browse "mailto:titiunik@princeton.edu":titiunik@princeton.edu}.{p_end}