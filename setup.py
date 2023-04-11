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

'''
import required packages from requirements file
with open("requirements.txt") as f:
    INSTALL_REQUIRES = [line.strip() for line in f.readlines()]'''

setuptools.setup(
    name='ultimodel',
    version='1.0.0',
    author='German Aerospace Center - DLR (Nina Thomsen)',
    author_email='nina.thomsen@dlr.de',
    description='Universal transport distribution model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/DLR-VF/ULTImodel',
    project_urls = {
        "Documentation": 'https://ultimodel.readthedocs.io/',
        "Source": 'https://github.com/DLR-VF/ULTImodel',
        "Bug Tracker": "https://github.com/DLR-VF/ULTImodel/issues "
    },
    license='MIT',
    packages=['ultimodel'],
    install_requires=['pandas', 'numpy', 'geopandas', 'shapely', 'osmnx', 'networkx', 'tqdm', 'requests', 'scikit-learn', 'scipy']
)