---
layout: default
title: Installation
permalink: /installation/
---
### Dependencies
- cmake version 3.9 or higher.
- boost version 1.66 or higher.

### Compiling/Building
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make pepsirf
```

Alternatively, one can use the included build script.
```
chmod +x build.sh
./build.sh
```
Both of these options will create a ```pepsirf``` executable in the build directory.

### Running Tests
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make pepsirf_test

./pepsirf_test
```
#### Note: Building for HPC Clusters
Some additional arguments may need to be specified, depending upon how your system's
module system is managed. For example on NAU's Monsoon cluster:
```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DBOOST_ROOT=/packages/boost/1.66.0-gcc-6.2.0 ..
```
### Usage
PepSIRF is a module-based program, and each module has its own set of arguments.
```
USAGE: pepsirf [ --help | module_name <module_args*> ]
--help, -h displays this message, while 'pepsirf module_name --help' will display the help for the module module_name.
```

### Generating Code Documentation
To create documentation for the code run the following:
```
cd doc
doxygen doxygen.conf
```
From here you can either use the html or LaTeX versions.
If you choose to use the LaTeX version, run the following:
```
cd latex
make
```
This creates a file named ```refman.pdf```.
