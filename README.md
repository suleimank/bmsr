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


To run demo of Bayesian multi-source regression on real data from GDSC, and CCLE datasets obtained from [PharmacoGx](https://bioconductor.org/packages/release/bioc/html/PharmacoGx.html) package, execute the following in R:

``` r
demo('demo_bmsr_lapatinib_gdsc_ccle')
```

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
	* Demo code to run the multi-source method on real data (./demo/demo\_bmsr\_lapatinib\_gdsc\_ccle.R)
	
## Source Code
To access a stand-alone repository with raw STAN and R source code, browse the [sourcecode](https://github.com/suleimank/bmsr-code).


## Citation
Cite as: 

Brian S. White, Suleiman A. Khan, Mike J Mason, Muhammad Ammad-ud-din, Swapnil Potdar, Disha Malani, Heikki Kuusanmäki, Brian J. Druker, Caroline A Heckman, Olli Kallioniemi, Stephen E Kurtz, Kimmo Porkka, Cristina E. Tognon, Jeffrey W. Tyner, Tero Aittokallio, Krister Wennerberg, Justin Guinney,
_Bayesian multi-source regression and monocyte-associated gene expression predict BCL-2 inhibitor resistance in acute myeloid leukemia_,
To Appear, (2021)


## FAQ
1. Do I need to compile the stan code after installation?

	No, the stan-code is compiled during package installation automatically. You do not need to compile it.

2. How long does it takes to run?

	The run-time depends on the data-size and number of sampling iterations you run. The demo_bmsr should run in 2-3 minutes depending on the machine.


