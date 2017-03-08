# Overview

This explains how to generate the documentation for Warpx.

## Installing the requirements

Install the Python requirements for compiling the documentation:
```
pip install sphinx breathe
```
and install Doxygen. Doxygen can be installed under Linux Debian by using
```
apt-get install doxygen
```
and under MacOSX by using Homebrew
```
brew install doxygen
```

## Compiling the documentation

`cd` into this directory and type
```
make html
```
