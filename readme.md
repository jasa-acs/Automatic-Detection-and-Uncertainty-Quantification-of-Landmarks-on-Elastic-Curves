# Automatic Detection and Uncertainty Quantification of Landmarks on Elastic Curves

# Author Contributions Checklist Form

## Data

### Abstract

Three main datasets are provided, covering almost all of the included examples within the main
manuscript. First is the Matlab file toycurve.mat, which provides the simulated curves used
throughout Section 5 of the manuscript. The other two datasets cover the application-based
examples of Section 6. In particular, MPEG7closed.mat is the MPEG-7 dataset (also available
online via the URL in the manuscript), consisting of numerous common shapes extracted from
computer vision scenes. The file mice.mat contains outlines of the second thoracic vertebrae of
mice split into three groups; this data was found originally in the ‘shapes’ package in R, with
detailed documentation provided there.

### Availability 

There are no restrictions on any of the provided data. MPEG-7 and the mice vertebrae data are
both available in the public domain.

### Description 

The simulated curve data was created by the authors, and will be posted to a data repository if
the manuscript is accepted, along with the webpages of at least one of the authors. The file is a
.mat file. The file ToyCurveLabels.txt describes each of the simulated curves included in the
corresponding .mat file.

The MPEG-7 Core Experiment CE-Shape-1 Test set is originally an image database which is
commonly used in computer vision to demonstrate improvements in shape matching algorithms
due to the wide variety of objects depicted in the images. The original data is publicly available
at http://www.dabi.temple.edu/~shape/MPEG7/dataset.html (the Shape Similarity Project from
members of the Department of Computer and Information Sciences at Temple University) and
consists of a zipped folder of .jpeg images. The MPEG-7 data provided with this manuscript is a
.mat file consisting only of the extracted shapes from the images of the original file. The
corresponding file MPEG7closedLabels.txt describes the contents of the .mat file, including
identification of the various object classes.

The mice vertebrae data contains the second thoracic vertebrae outlines of mice split into three
groups by body weight: Large, Small, and Control. Originally provided by Paul O’Higgins (Hull-
York Medical School) and David Johnson (Leeds), this dataset is now publicly available through
the ‘shapes’ package in R. Documentation for the ‘shapes’ package is found at https://cran.r-
project.org/web/packages/shapes/shapes.pdf. This package was developed by Ian Dryden to
complement his co-authored book with Kanti Mardia titled “Statistical Shape Analysis.” The
original application of this data used landmark-based shape representations for analysis; thus,
the file mice.mat includes the vertebrae outlines in several forms. The form used in this
manuscript consists of N = 61 discretization points for each shape, and can be identified using
the MiceLabels.txt file. Group labels for the vertebrae are also available within this file.

## Code

### Abstract 

We are providing the relevant code needed to perform model inference for a fixed and an
unknown number of landmarks using the computing methods described in the manuscript. In
particular, there are separate files for MCMC inference of landmark locations (given a fixed
number of landmarks), as well as RJMCMC inference of landmark locations and their quantity.
Both programs have been written to accept open and closed curve samples. Necessary sub-
programs are also included, along with detailed comments describing various lines of code.

### Description

The included code is written in Matlab, and will be placed in both a code repository and on the
authors’ websites after the manuscript is accepted. All written code was developed by the
authors, and consist of 11 .m files.

### Optional Information

In order to run the provided code, any reasonably recent version of Matlab is suitable. The file
Findbestk.m is the only one which requires an additional comment: parallel computing is used in
this program (to run many MCMC chains at once), through the use of a parfor loop. In the event
that this is not feasible, simply replace the parfor loop with a for loop.

## Instructions for use

### Reproducibility 

All examples, figures, and tables from the main manuscript should be reproducible provided the
given code, except for the brain substructure application of Section 6.3. The included file
ExampleCode.m provides step-by-step commands (with comments) for the examples provided
in Section 5. Some commands for examples in Section 6 are also provided. In order to use and
reproduce this work, one simply needs to load data in as listed in the ExampleCode.m file, and
then run various lines to produce results found throughout Sections 5 and 6.
