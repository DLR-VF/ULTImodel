# ULTIMO

[![License: BSD](https://img.shields.io/badge/License-BSD-green.svg)](link-to-license-github)
[![Documentation Status](link-to-readthedocs)](link-to.readthedocs)

ULTIMO &mdash; A universal transport distribution model.

## Introduction

ULTIMO is a distribution model that helps to spatially distribute road-based transport for countries, including 
border-crossing travel. It is set up using open data like [OSM](https://openstreetmap.org).
The software includes modules for network generation, trip generation and trip distribution based on two main inputs:

* Georeferenced traffic analysis zones (TAZ) for the respective region
* Target value for national transport volume (i.e. person-kilometres or tonne-kilometres)

*insert schema*

## Background

ULTIMO was initially developed in order to spatially distribute road transport emissions. The goal was to introduce a 
bottom-up approach using traffic flows from a transport model and emission factors to determine spatially distributed 
emissions (see e.g. [Mathias et al., 2020](https://doi.org/10.1016/j.trd.2020.102536)). For a large-scale implementation, a simplified approach was needed that:

* produces reliable results for multiple countries
* uses data that is widely available

With ULTIMO, it is possible to model travel movements in central Europe based on a NUTS-3 cell structure with total 
travel volumes given per country. The distribution was calibrated and validated for Germany, it is however possible to 
change the parameters for applications to other world regions.

For more information, see [Thomsen & Seum, 2020](https://aetransport.org/past-etc-papers/conference-papers-2021?abstractId=7202&state=b).