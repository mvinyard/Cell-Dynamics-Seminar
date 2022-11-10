
# -- import packages: ----------------------------------------------------------------------
import setuptools, re, os, sys


# -- run setuptools: ------------------------------------------------------------------------
setuptools.setup(
    name="cd_seminar",
    version="0.0.0",
    python_requires=">3.7.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="https://github.com/mvinyard/Cell-Dynamics-Seminar",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="Materials for the UC Berkeley CCB Cell Dynamics Seminar",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy<1.22",
        "pandas>=1.1.2",
        "vinplots>=0.0.73",
        "licorice_font>=0.0.3",
        "webfiles>=0.0.1",
        "nb_black",
        "torch>=1.12",
        "cellrank>=1.5.1",
        "pickle5",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.7",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
