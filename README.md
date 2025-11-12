# CleCurExpRule – Clenshaw–Curtis Quadrature for Exponentially Weighted Integrals

This MATLAB project implements Clenshaw-Curtis quadrature rules for computing integrals of the form:

$$ I = \int_0^b f(s) \exp(z s) ds $$

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

### Auxiliary functions

* computeWeights – Computes the quadrature weights for the Clenshaw-Curtis method.
* computeRhoMax – Computes rho_m(z) for sufficiently large m 
* idctI – Implements the Type-I inverse discrete cosine transform.
* idctII – Implements the Type-II inverse discrete cosine transform.
* thomas_algorithm - solver for tridiagonal linear systems

### Auxiliary functions - vpa versions

* computeWeights_vpa – Computes the quadrature weights for the Clenshaw-Curtis method.
* computeRhoMax_vpa – Computes rho_m(z) for sufficiently large m 
* idctI_vpa – Implements the Type-I inverse discrete cosine transform.
* idctII_vpa – Implements the Type-II inverse discrete cosine transform.
* ifft_vpa – Performs the inverse FFT using variable-precision arithmetic.

### Training and test area folder

Some script files for testing purposes.


## Usage 
Several examples: the integral
$$ I = \int_0^2 \frac{cos(5\pi x)}{4+\sin(4 \pi x)}\exp((-20+15i) s) ds $$

can be (numerically) computed with 

```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
I = CleCurExpRule(f,-20+15i,'NumberOfNodes',64)

% I =
%   0.0076 + 0.0024i
```
or
```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
t = linspace(0,pi,65);
y = 1+cos(t(:));  
I = CleCurExpRule(f(t),-20+15i)
% I = 
% -0.0007 + 0.0014i
```
**Notice** that the values of $f$ at Chebyshev nodes are provided in reverse 
order in the second case 

Other examples: integral in [0,5]

```matlab
I = CleCurExpRule(f,-20+15i,'NumberOfNodes',64,'EndPoint',5)
%
% I =
%   0.0076 + 0.0024i
```
The implementation supports vector functions. When the input is (n+1) x m matrix assume that each column 
corresponds to pointwise values of m-component vector function at 
the n+1 nodes of Chebyshev nodes (in reverse order!): 

```matlab
f = @(x) [cos(5*pi*x)./(4+sin(4*pi*x))...
          sin(5*pi*x)./(4+cos(4*pi*x))];
I1 = CleCurExpRule(f,-20+15i,'NumberOfNodes',64,'EndPoint',3)
% Equivalently
nNodes = 64; b = 3;  
t = linspace(0,pi,nNodes+1); t = t(:);
t = (1+cos(t))*b/2;
I2 = CleCurExpRule(f(t),-20+15i,'NumberOfNodes',64,'EndPoint',3)
%
% I1 =
%    0.0076 + 0.0024i   0.0025 + 0.0039i
%
% I2 =
%    0.0076 + 0.0024i   0.0025 + 0.0039i
```


For **nNodes even** an error estimate is provided just comparing the result
of the quadrature rule with **nNodes /2**:
```matlab
f = @(x) [cos(5*pi*x)./(4+sin(4*pi*x))...
          sin(5*pi*x)./(4+cos(4*pi*x))];
[I2,errorEstimate] = ...
 CleCurExpRule(f,-20+15i,'NumberOfNodes',96,'EndPoint',3)
%
% I2 =
%   0.0076 + 0.0024i   0.0025 + 0.0039i
%
% errorEstimate =
%   1.0e-07 *
%
%   0.1982 - 0.0185i  -0.1112 + 0.0003i
%
```

## More information:
For more information:  
    https://www.arxiv.org/abs/2503.08169
    
Author:  Victor Dominguez  
Contact: victor.dominguez@unavarra.es  
Date:    12 November 2025

---
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

