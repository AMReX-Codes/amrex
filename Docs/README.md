# Overview

This explains how to generate the documentation for Warpx, and contribute to it.

## Generating the documentation

### Installing the requirements

Install the Python requirements for compiling the documentation:
```
pip install sphinx sphinx_rtd_theme breathe
```
and install Doxygen. Doxygen can be installed under Linux Debian by using
```
apt-get install doxygen
```
and under MacOSX by using Homebrew
```
brew install doxygen
```

### Compiling the documentation

`cd` into this directory and type
```
make html
```

When executing this command, Doxygen parses the C++ code to find relevant annotations. Then Sphinx uses the Doxygen results.

You can then browse the documentation by opening the file `build/html/index.html` with a web browser.

### Cleaning the documentation

In order to remove all of the generated files, use:
```
make clean
```

## Contributing to the documentation

The documentation is built with several tools:

- The overall framework is [Sphinx](http://www.sphinx-doc.org/en/stable/intro.html). Sphinx reads the files that are located in the folder `source`, and which are written in [reStructuredText](http://docutils.sourceforge.net/rst.html) format. These files describe the structure of the documentation, where to insert automatic code documentation. Please read the [Sphinx documentation](http://www.sphinx-doc.org/en/stable/intro.html) for more information.

- For Python code (e.g. the WarpX Python interface), Sphinx can generate automatic code documentation by itself (it imports the Python package directly, and extracts the docstrings).

- For C++ code (e.g. the WarpX C++ interface), Sphinx couples with [Doxygen](http://www.stack.nl/~dimitri/doxygen/), using the [Breathe](https://breathe.readthedocs.io/en/latest/) package. More precisely, Doxygen extracts the structure of the C++ code and the corresponding inline documentation and creates a corresponding XML file. Breathe parses this XML file for Sphinx. Please read the [Doxygen documentation](http://www.stack.nl/~dimitri/doxygen/manual/index.html) and [Breathe documentation](https://breathe.readthedocs.io/en/latest/) for more information.

- For Latex files, you need to convert them into [reStructuredText](http://docutils.sourceforge.net/rst.html) using [pandoc](http://pandoc.org/), and to commit them in the folder `source`.
