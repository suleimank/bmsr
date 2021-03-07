# Bayesian multi-source regression
The bmsr package implements joint regression from multiple data sources in a Bayesian framework. The package provides implementation for both single-task and multi-task regression. The model is implemented using STAN and interface is provided using R programming language. Options for training the model using both NUTS sampler and variational inference are provided. The package is structured for ease of use and the included demo shows the model execution on real-life as well as simulated datasets.

## Installation
To install the package, run the following command in R:

``` r
devtools::install_github("suleimank/bmsr")
```

## Reference Manual
The repository has the R-package [Reference Manual](bmsr_1.0.0.pdf).

## Demo
To run demo of Bayesian multi-source regression on simulated data, execute the following in R:

``` r
demo('demo_bmsr')
```
![demo bmsr features](/images/demo_bmsr.png)

To run demo of Bayesian multi-task regression on simulated data, execute the following in R:

``` r
demo('demo_bmsmtr')
```


To run demo of Bayesian multi-source regression to predict biomarkers of lapatinib response in breast cancer using GDSC and CCLE datasets accessed and harmonized via the [PharmacoGx](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) package, execute the following in R:

``` r
demo('demo_bmsr_lapatinib_gdsc_ccle')
```

To run demo of Bayesian multi-source regression to predict biomarkers of venetoclax response in AML using the "Tavor" and "Beat AML" OHSU datasets  accessed from [ORCESTRA](https://www.orcestra.ca/) and harmonized via the [PharmacoGx](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) package, execute the following in R:

``` r
demo('demo_bmsr_venetoclax_tavor_ohsu')
```
![demo bmsr features](/images/demo_bmsr_venetoclax_tavor_ohsu.png)

## Contents
The repository contains the following:

``` bash
.
├── R
│   └── ...
├── demo
│   └── ...
└── inst
    └── stan
        └── ...
```

* BMSR: Bayesian multi-source regression
	* STAN code for Bayesian multi-source regression (./inst/stan/bmsr.stan)
	* R script for training and predicting from the BMSR method (./R/bmsr.R)
 
* BMSMTR: Bayesian multi-source multi-task regression
	* STAN code for Bayesian multi-source multi-task regression (./inst/stan/bmsmtr.stan)
	* R script for training and predicting from the BMSMTR method (./R/bmsr.R)

* Demo:
	* Demo code to run the multi-source method on simulated data (./demo/demo\_bmsr.R)
	* Demo code to run the multi-task method on simulated data (./demo/demo\_bmsmtr.R)
	* Demo code to run the multi-source method to predict biomarkers of lapatinib response in breast cancer using CCLE and GDSC datasets via PharmacoGx (./demo/demo\_bmsr\_lapatinib\_gdsc\_ccle.R)
	* Demo code to run the multi-source method to predict biomarkers of venetoclax response in AML using Tavor and OHSU datasets via PharmacoGx (./demo/demo\_bmsr\_venetoclax\_tavor\_ohsu.R)	

## Source Code
To access a stand-alone repository with raw STAN and R source code, browse the [sourcecode](https://github.com/suleimank/bmsr-code).


## Citation
Cite as: 

Brian S. White, Suleiman A. Khan, Mike J Mason, Muhammad Ammad-ud-din, Swapnil Potdar, Disha Malani, Heikki Kuusanmäki, Brian J. Druker, Caroline A Heckman, Olli Kallioniemi, Stephen E Kurtz, Kimmo Porkka, Cristina E. Tognon, Jeffrey W. Tyner, Tero Aittokallio, Krister Wennerberg, Justin Guinney,
__Bayesian multi-source regression and monocyte-associated gene expression predict BCL-2 inhibitor resistance in acute myeloid leukemia__,
To Appear, (2021)

## References

Anthony Mammoliti, Petr Smirnov, Minoru Nakano, Zhaleh Safikhani, Chantal Ho, Gangesh Beri, Benjamin Haibe-Kains, __ORCESTRA: a platform for orchestrating and sharing high-throughput pharmacogenomic analyses__, BioRxiv, https://doi.org/10.1101/2020.09.18.303842

Petr Smirnov, Zhaleh Safikhani, Nehme El-Hachem, Dong Wang, Adrian She, Catharina Olsen, Mark Freeman, Heather Selby, Deena M A Gendoo, Patrick Grossmann, Andrew H Beck, Hugo J W L Aerts, Mathieu Lupien, Anna Goldenberg, Benjamin Haibe-Kains, 2016. __PharmacoGx: An R Package for Analysis of Large Pharmacogenomic Datasets__, Bioinformatics  32 (8): 1244–46.

Sigal Tavor, Tali Shalit, Noa Chapal Ilani, Yoni Moskovitz, Nir Livnat, Yoram Groner, Haim Barr, Mark D Minden, Alexander Plotnikov, Michael W Deininger, Nathali Kaushansky, Liran I Shlush. 2020. __Dasatinib Response in Acute Myeloid Leukemia Is Correlated with FLT3/ITD, PTPN11 Mutations and a Unique Gene Expression Signature__, Haematologica 105 (12): 2795–2804.

Jeffrey W. Tyner, Cristina E. Tognon, Daniel Bottomly, Beth Wilmot, Stephen E. Kurtz, Samantha L. Savage, Nicola Long, et al. 2018. __Functional Genomic Landscape of Acute Myeloid Leukaemia__, Nature 562 (7728): 526–31.



## FAQ
1. Do I need to compile the stan code after installation?

	The stan-code is compiled automatically during package installation. You do not need to compile it.

2. How long does it takes to run?

	The runtime depends on the data size and number of sampling iterations you run. The demo_bmsr should run in 1-3 minutes depending on the machine.


