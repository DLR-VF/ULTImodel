# =========================================================
# setup.py
# @author Nina Thomsen
# @date 20.03.2023
# @copyright Institut fuer Verkehrsforschung,
#            Deutsches Zentrum fuer Luft- und Raumfahrt
# @brief setup module for ULTIMO
# =========================================================

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='ultimo',
    version='1.0',
    author='German Aerospace Center - DLR (Nina Thomsen)',
    author_email='nina.thomsen@dlr.de',
    description='Universal transport distribution model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://gitlab.dlr.de/ultimo/ultimo-release',
    project_urls = {
        "Documentation": 'readthedocs',
        "Source": 'github',
        "Bug Tracker": "https://github.com/"
    },
    license='BSD',
    packages=['ultimo'],
    install_requires=['pandas', 'numpy', 'geopandas', 'shapely', 'osmnx', 'networkx', 'tqdm', 'requests', 'sklearn', 'scipy']
)
