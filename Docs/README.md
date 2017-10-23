# Overview

This explains how to generate the documentation for Warpx, and contribute to it.

## Generating the documentation

### Installing the requirements

Install the Python requirements for compiling the documentation:
```
pip install sphinx sphinx_rtd_theme
```
In addition, you need to install `pywarpx` (i..e the Python package for WarpX),
so that Sphinx can automatically build the documentation for this package.

### Compiling the documentation

`cd` into this directory and type
```
make html
```

### Cleaning the documentation

In order to remove all of the generated files, use:
```
make clean
```

## Contributing to the documentation

The documentation is built with several tools:

- The overall framework is [Sphinx](http://www.sphinx-doc.org/en/stable/intro.html). Sphinx reads the files that are located in the folder `source`, and which are written in [reStructuredText](http://docutils.sourceforge.net/rst.html) format. These files describe the structure of the documentation, and where to insert automatic code documentation. Please read the [Sphinx documentation](http://www.sphinx-doc.org/en/stable/intro.html) for more information.

- For Python code (e.g. the WarpX Python interface), Sphinx can generate automatic code documentation by itself (it imports the Python package directly, and extracts the docstrings).

- For Latex files, you need to convert them into [reStructuredText](http://docutils.sourceforge.net/rst.html) using [pandoc](http://pandoc.org/), and to commit them in the folder `source`. This is automatically done by the `Makefile`, for the theory section of the WarpX documentation.
