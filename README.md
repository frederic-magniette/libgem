Libgem version 1.0
=================

Libgem is a library for tracks fitting with EMs algorithms. The whole project is freely available in LGPL.

The library itself contains a lot of numerical features : weighted gaussian distributions, line and circle regressions, vector computations, EM algorithm and a graphical interface based on gnuplot or sdl depending on the dimension.
The library is written in C for performances.

Four high level algorithms have been developped and are given as executables :

LFit : fit a fixed number of lines in arbitrary dimensions
MFit : fit lines or circles in a sequential way
UFit : fit combined mix of lines and circles in a sequential way
SPFit : fit a single spiral

Developped by Frédéric Magniette

Installation
------------

Install dependencies : SDL1.2 development package and gnuplot with wxt support

Compile the libgem and executables by simply typing : make

LFit usage
----------

lfit.exe datafile dim nb_lines seed output

- datafile is the path of the data file. This datafile has to be in gnuplot format (numbers in columns separated by tab or space)
- dim is an integer and represents the dimension of the space
- nb_lines is an integer and represents the number of lines to fit
- seed is an integer and represents the seed of the random generator. The random is used to initialize the initial distributions.
- convcrit is a float and represents the convergence criterion of the EM algorithm. 0.00001 is a typical value.
- output is a string and represents the output mode. It can take following values : 

- none : no graphical output.
- summary : only the result of the fit is graphically represented.
- detail : all the steps of the algorithm are graphically represented.

examples of good fit in 1D, 2D and 3D : 
- ./lfit.exe data/demo_1d.txt 1 2 0 0.00001 detail
- ./lfit.exe data/cross_noised_1.txt 2 2 0 0.00001 detail
- ./lfit.exe data/demo_3d.txt 3 2 2 0.00001 detail

example of bad fit because of bad initialization : 
- ./lfit.exe data/demo_3d.txt 3 2 0 0.00001 detail

MFit usage
----------

mfit.exe datafile dim type[line|circle] convcrit scalecrit output

- datafile is the path of the data file. This datafile has to be in gnuplot format (numbers in columns separated by tab or space)
- dim is an integer and represents the dimension of the space
- type is a string and reprensents the type of shape to fit. it can take following values : 
- line
- circle
- convcrit is a float and represents the convergence criterion of the EM algorithm. 0.001 is a typical value.
- scalecrit is a float and represents the scale criterion of the mfit algorithm. The precision of the algorithm but also its performance depends on this parameter. 40 is a reasonable value, giving good precision and not too expensive in term of performance.
- output is a string and represents the output mode. It can take following values : 

- none : no graphical output.
- summary : only the result of the fit is graphically represented.
- detail : all the steps of the algorithm are graphically represented.

examples of linear fits in 1D, 2D and 3D : 
- ./mfit.exe data/demo_1d.txt 1 line 0.00001 20 detail
- ./mfit.exe data/demo_2d.txt 2 line 0.00001 20 detail
- ./mfit.exe data/demo_3d.txt 3 line 0.00001 20 detail

examples of circular fits in 2D : 
- ./mfit.exe data/demo_mc1.txt 2 circle 0.001 40 detail
- ./mfit.exe data/demo_mc2.txt 2 circle 0.01 40 detail

UFit usage
----------

ufit.exe datafile dim convcrit scalecrit output

- datafile is the path of the data file. This datafile has to be in gnuplot format (numbers in columns separated by tab or space)
- dim is an integer and represents the dimension of the space
- convcrit is a float and represents the convergence criterion of the EM algorithm. 0.001 is a typical value.
- scalecrit is a float and represents the scale criterion of the mfit algorithm. The precision of the algorithm but also its performance depends on this parameter. 40 is a reasonable value, giving good precision and not too expensive in term of performance.
- output is a string and represents the output mode. It can take following values : 
- none : no graphical output.
- summary : only the result of the fit is graphically represented.
- detail : all the steps of the algorithm are graphically represented.

example in 2D : 
- ./ufit.exe data/demo_ufit.txt 2 0.001 40 summary

SPFit usage
-----------

spfit.exe data_filename

- datafile is the path of the data file. This datafile has to be in gnuplot format (numbers in columns separated by tab or space)

example : 
- ./spfit.exe data/demo_spiral.txt

Data
----

The data folder contains demo data files. 


You can also compiling tools for generating data by typing : make.

random_line.exe generates a random line given dimension, number of points you want, output filename, random seed and coordinates of limit points.

random_circle.exe generates a random circle given dimension, number of points you want, output filename, random seed and coordinates of limit points.

angle_line.exe generates a line with fixed angle with x axis (works only in 2D)

noise.exe can noise a datafile




