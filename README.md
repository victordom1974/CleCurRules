# CleCurExpRule - Clenshaw-Curtis Quadrature with Exponential Weight

This MATLAB project implements Clenshaw-Curtis quadrature rules for computing integrals of the form:  

$$I = \int_0^b f(s) \exp(z s) {\rm d}s$$

## Installation
Clone or download the repository and add the files to your MATLAB path:

```matlab
addpath(genpath('CleCurExpRule'))
```

## Functions Overview

### Core Functions

* CleCurExpRule – Computes the integral using the product Clenshaw-Curtis rule.

### Core Functions - vpa version

* CleCurExpRule_vpa – VPA version of CleCurExpRule. Only for testing purposes. 

### Auxiliar functions

* computeWeights – Computes the quadrature weights for the Clenshaw-Curtis method.
* computeRhoMax – Computes rho_m(z) for sufficiently large m 
* idctI – Implements the Type-I inverse discrete cosine transform.
* idctII – Implements the Type-II inverse discrete cosine transform.
* thomas_algorithm - solver for tridiagonal linear systems

### Auxiliar functions - vpa versions

* computeWeights – Computes the quadrature weights for the Clenshaw-Curtis method.
* computeRhoMax – Computes rho_m(z) for sufficiently large m 
* idctI – Implements the Type-I inverse discrete cosine transform.
* idctII – Implements the Type-II inverse discrete cosine transform.
* mifft_vpa – Performs the inverse FFT using variable-precision arithmetic.

### Training and test area folder

Some script files for testing purposes.


## Usage 
Several examples: the integral

$$ I = \int_0^2 \frac{\cos(5\pi x)}{4+\sin(4 \pi x)}\exp((-20+15i) s) {\rm d}s $$

can be (numerically) computed with 

```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
CleCurExpRule(f,-20+15i,'NumberOfNodes',64)
```
or
```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
t = linspace(0,pi,65);
y = 1+cos(t(:));  
CleCurExpRule(f(t),-20+15i)
```

Other examples: integral in [0,5]

```matlab
CleCurExpRule(f,-20+15i,'NumberOfNodes',64,'EndPoint',5)
```

## More information:
For more information:  
    https://www.arxiv.org/abs/2503.08169
    
Author: Victor Dominguez  
Contact: victor.dominguez@unavarra.es  
Date: 11 March 2025

## Copyright

Copyright (C) 2025 Victor Dominguez

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.
