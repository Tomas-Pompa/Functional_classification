## Classification of Functional Data with vertical differences using SVM

This document serves as a supporting material for the article:

**Classification of Functional Data with vertical differences using SVM**,

which is the extension of authors master thesis with the official project assignment as follows.

*Assignment*: In many applications, measured data represent values of a function. Therefore, if the situation allows, it is advantageous to work with them in a functional form, i.e., as elements of an infinite-dimensional space. The thesis will build upon the student’s previous bachelor’s work. The aim is to generalize machine learning methods for functional data and to describe the properties of such approaches. The obtained results will be demonstrated on simulated or real data.

### Support Vector Machines for Functional Data

The purpose of this document is to apply knowledge of the Support Vector Machine (SVM) method for multivariate data to functional-type data, i.e., infinite-dimensional objects. To achieve this, we will primarily use the transformation (reduction) of objects from infinite to finite dimension, followed by the application of known procedures from finite dimensions. Several possible approaches will be demonstrated.

Another objective will be to compare various methods for the classification of functional data on real and simulated datasets. We will focus primarily on simulated data and, in addition to comparing the methods against each other through a simulation study, we will also examine the dependency of classification success on the parameters used in data generation (we will be interested in the variance around the generating curves and the variance of vertical shifts). We will also look into the impact of interval discretization on the effectiveness of classification methods, which is one approach to applying finite-dimensional methods to functional data.

The classification methods considered include:

  - $K$ nearest neighbors (KNN),

  - logistic regression (both standard (LR) and its functional modification (LR_fda)),

  - linear (LDA) and quadratic (QDA) discriminant analysis,

  - decision trees (DT),

  - random forests (RF), and

  - Support Vector Machines: here we will consider many variants, all of which are based on the principle of filtering (dimension reduction) or direct aplying the SVM on functional data.

We will go through each method, starting with simulated data, and then proceed to construct the Support Vector Machine method for functional data.

The primary package in `R` for working with functional objects is `fda`. Other useful packages include `MASS`, `e1071`, `fda.usc`, `refund`, and more.

In the application section of this document, we will examine dataset `tecator`. 
