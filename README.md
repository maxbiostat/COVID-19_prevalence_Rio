## Estimating COVID-19 prevalence in Rio de Janeiro from the EPICOVID serological survey (and case data).

This repository contains code to reproduce the prevalence estimation analyses in Bastos, Carvalho & Gomes (2021). 

The contents of `data/` are an amalgamation of data from the [EPÌCOVID](http://www.rs.epicovid19brasil.org/) and [Brasil.io](https://brasil.io/home/) projects, for serological survey and case/death data, respectively. 
The `stan` folder contains Stan code for fitting (i) simple prevalence correction model where the true number of cases projected from the prevalence is a generated quantity and (ii) a slightly less simple model that takes case data and prior knowledge about underreporting to estimate prevalence and detection probability jointly from survey and case data.
