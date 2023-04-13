# ULTImodel

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/DLR-VF/ULTImodel/blob/master/LICENSE)
[![PyPI version](https://badge.fury.io/py/ultimodel.svg)](https://badge.fury.io/py/ultimodel)
[![Documentation Status](https://readthedocs.org/projects/ultimodel/badge/?version=latest)](https://ultimodel.readthedocs.io/en/latest/?badge=latest)
[![Cite-us](https://img.shields.io/badge/doi-10.5281%2Fzenodo.7817425-blue)](https://doi.org/10.5281/zenodo.7826486)

ULTImodel &mdash; A universal transport distribution model written in Python.

## Introduction

ULTImodel is a distribution model that helps to spatially distribute road-based transport for countries, including border-crossing travel. It is set up using open data like [OSM](https://openstreetmap.org).
The software includes modules for network generation, trip generation and trip distribution based on two main inputs:

* Georeferenced traffic analysis zones (TAZ) for the respective region
* Target value for national transport volume (i.e. person-kilometres or tonne-kilometres)

![Prim_Sec](ultimodel-mkdocs/docs/images/readme_visual_fr.png "Results of distribution and secondary model")

## Background

ULTImodel was initially developed in order to spatially distribute road transport emissions. The goal was to introduce a 
bottom-up approach using traffic flows from a transport model and emission factors to determine spatially distributed 
emissions (see e.g. [Mathias et al., 2020](https://doi.org/10.1016/j.trd.2020.102536)). For a large-scale implementation, a simplified approach was needed that:

* produces reliable results for multiple countries
* uses data that is widely available

With ULTImodel, it is possible to model travel movements in central Europe based on a NUTS-3 cell structure with total 
travel volumes given per country. The distribution was calibrated and validated for Germany, it is however possible to 
change the parameters for applications to other world regions.

For more information, see:
 
Thomsen, Nina und Seum, Stefan (2021) [Using Open Data for Spatial Transport Emission Modelling](https://aetransport.org/past-etc-papers/conference-papers-2021?abstractId=7202&state=b). 
European Transport Conference ETC 2021, 13.-15. Sep. 2021, online.

## Citation

Please cite as *German Aerospace Center (DLR): ULTImodel. https://doi.org/10.5281/zenodo.7826486*
