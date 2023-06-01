<!-- badges: start -->

<!-- For more info: https://usethis.r-lib.org/reference/badges.html -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

# Emperor Goose Harvest Strategy

## Overview

This code is provided to support [emperor goose (*Anser canagicus*)](https://en.wikipedia.org/wiki/Emperor_goose) harvest management efforts of the [Alaska Migratory Bird Co-management Council](https://www.alaskamigratorybirds.com/) (AMBCC) and [Pacific Flyway Council](https://www.pacificflyway.gov/). In 2016 the first AMBCC management plan for emperor goose was approved and hunting seasons were opened in the spring and fall of 2017. These were considered experimental and management strategies were to be revised after three years of observing harvest and population response data. This project present code used to model population response to harvest, predict consequences to various management actions, and estimate uncertainty associated with predictions.  

The main component of this project is a Bayesian integrated population model (theta.logistic.R) that combines harvest data with goose population survey data to estimate demographic parameters of a mathematical population model. This model and the full parametric uncertainty is then used to predict future responses of the population, expected harvest, and other management objectives under various harvest policies (simulate.pop.R). Finally, the population simulation is used to final the optimal management policy given a single objective or multiple objective function (threshold.opt.R).  

This is a work in progress and all products are under active development. The code is being provided for those that wish to contribute and to allow for maximum transparency in the modelling process.  

## Installation

You must have a current version of R ($>$v4.2) and [JAGS (
$\geq$ v4.3.1)](https://sourceforge.net/projects/mcmc-jags/) to run the code presented here.  Follow the instructions under [Usage](#usage) below to use the code and provide feedback.

## Usage

Contact the [project maintainer](emailto:erik_osnas@fws.gov) for data needed to run the model and produce posterior distributions. Population survey data are available in the R package [AKaerial](https://github.com/USFWS/AKaerial/) and harvest data is available with permission from the Alaska Department of Fish and Game. 

## Getting help

Contact the project maintainer for help with this code.  

## Contribute

Contact the project maintainer for information about contributing to this project. Submit a [GitHub Issue](https://github.com/USFWS/Emperor-Goose-Harvest-Strategy/issues) to report a bug or request a feature or enhancement.
