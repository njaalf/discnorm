discnorm
==========

This package contains an implementation of a the bootstrap test for underlying non-normality proposed by Foldnes and Gronneberg (Structural Equation Modeling, 2019). Also contains an adjusted polychoric estimator proposed by Gronneberg and Foldnes (Psychological Methods, 2022).


How to install
--------------


You can install:

-   the stable release on CRAN:

    ``` r
    install.packages("discnorm")
    ```

-   the latest development version:

    ``` r
    devtools::install_github("njaalf/discnorm")
    ```

------------------------------------------------------------------------

Package overview
----------------


The package offers function bootTest() which tests an ordinal data frame for underlying normality.
A function catLSadj() is provided that computes the adjusted polychoric correlations based on user-provided non-normal marginals.

References
----------
Njål Foldnes & Steffen Grønneberg (2019) Pernicious Polychorics: The Impact and Detection of Underlying Non-normality, Structural Equation Modeling: A Multidisciplinary Journal, DOI: 10.1080/10705511.2019.1673168

Steffen Grønneberg & Njål Foldnes (2022) Factor Analyzing Ordinal Items Requires Substantive Knowledge of Response Marginals, Psychological Methods, DOI: 10.1037/met0000495
