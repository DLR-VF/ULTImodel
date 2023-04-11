# ULTImodel

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://github.com/DLR-VF/ULTImodel/blob/master/LICENSE)
[![PyPI version](https://badge.fury.io/py/ultimodel.svg)](https://pypi.python.org/pypi/ultimodel)
[![Documentation Status](https://readthedocs.org/projects/ultimodel/badge/?version=latest)](https://ultimodel.readthedocs.io/en/latest/?badge=latest)
[![Cite-us](https://img.shields.io/badge/doi-10.5281%2Fzenodo.7817425-blue)](https://doi.org/10.5281/zenodo.7817425)
 
**ULTImodel** &mdash; A universal transport distribution model written in Python.

## Description
**ULTImodel** is a distribution model that helps to spatially distribute road-based transport for countries, including border-crossing travel. It is set up using open data like [OSM](https://openstreetmap.org).
The software includes modules for network generation, trip generation and trip distribution based on two main inputs:

* Georeferenced traffic analysis zones (TAZ) for the respective region
* Target value for national transport volume (i.e. person-kilometres or tonne-kilometres)

![Prim_Sec](ultimodel-mkdocs/docs/images/readme_visual_fr.png "Results of distribution and secondary model")

## Installation

The __current version__ is [ultimodel-1.0.0](https://github.com/DLR-VF/ULTImodel/releases/tag/1.0.0).

You may __install ULTImodel__ by executing the following

__pip__
```console
python -m pip install ultimodel
```
__conda__
```console
conda install -c conda-forge ultimodel
```

You may __download a copy or fork the code__ at [ULTImodel&apos;s github page](link-to-github).

Besides, you may __download the current release__ here:

* [ultimodel-1.0.0.zip](https://github.com/DLR-VF/ULTImodel/archive/refs/tags/1.0.0.zip)
* [ultimodel-1.0.0.tar.gz](https://github.com/DLR-VF/ULTImodel/archive/refs/tags/1.0.0.tar.gz)


## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Authors and acknowledgment
**ULTImodel** was developed by Nina Thomsen.

We want to thank the following persons for the help during **ULTImodel's** development: Lars Hedemann, Simon Metzler, Gernot Liedtke, Christian Winkler, Tudor Mocanu, and Daniel Krajzewicz.

## License
**ULTImodel** is licensed under the [MIT license](https://github.com/DLR-VF/ULTImodel/blob/master/LICENSE).

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
