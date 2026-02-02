# CleCurExpRule — Clenshaw–Curtis Quadrature for Exponentially Weighted Integrals

This MATLAB project implements Clenshaw-Curtis quadrature rules for computing integrals of the form:

$$
I = \int_0^b f(s) \exp(z s) \, ds
$$

where $f$ is a smooth function and $z \in \mathbb{C}$ may lead to highly oscillatory or exponentially decaying integrands.

---

## Installation

Clone or download the repository and add the files to your MATLAB path:

```matlab
addpath(genpath(pwd))
```
Note: this project relies on a structured directory layout; therefore, all
subfolders must be added to the MATLAB path before use.


---

## Functions Overview

### Core Functions

- **`CleCurExpRule`** — Computes the integral using the product Clenshaw-Curtis rule.

### Core Functions - VPA Version

- **`CleCurExpRule_vpa`** — Variable-precision arithmetic (VPA) version of `CleCurExpRule`. For testing purposes only.

### Auxiliary Functions

- **`computeWeights`** — Computes the quadrature weights for the Clenshaw-Curtis method.
- **`computeRhoMax`** — Computes $\rho_m(z)$ for sufficiently large $m$.
- **`idctI`** — Implements the Type-I inverse discrete cosine transform.
- **`idctII`** — Implements the Type-II inverse discrete cosine transform.
- **`thomas_algorithm`** — Solver for tridiagonal linear systems. Supports double, complex, and variable-precision (VPA) arithmetic.

### Auxiliary Functions - VPA Versions

- **`computeWeights_vpa`** — VPA version of `computeWeights`.
- **`computeRhoMax_vpa`** — VPA version of `computeRhoMax`.
- **`idctI_vpa`** — VPA version of `idctI`.
- **`idctII_vpa`** — VPA version of `idctII`.
- **`ifft_vpa`** — Performs the inverse FFT using variable-precision arithmetic.

### Training and Test Files

The repository includes various test scripts for validation and experimentation.

---

## Usage Examples

### Example 1: Basic Integration

Compute the integral:

$$
I = \int_0^2 \frac{\cos(5\pi x)}{4+\sin(4\pi x)} \exp((-20+15i) s) \, ds
$$

**Method 1: Using function handle**

```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
I = CleCurExpRule(f, -20+15i, 'NumberOfNodes', 64)

% Output:
% I = 0.0076 + 0.0024i
```

**Method 2: Using pre-computed values at Chebyshev nodes**

```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
t = linspace(0, pi, 65);
y = 1 + cos(t(:));  
I = CleCurExpRule(f(t), -20+15i)

% Output:
% I = -0.0007 + 0.0014i
```

> **Note:** When providing pre-computed values at Chebyshev nodes, they must be given in **reverse order** (as in the second example).

---

### Example 2: Custom Integration Interval

To integrate over the interval $[0, 5]$:

```matlab
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));
I = CleCurExpRule(f, -20+15i, 'NumberOfNodes', 64, 'EndPoint', 5)

% Output:
% I = 0.0076 + 0.0024i
```

---

### Example 3: Vector-Valued Functions

The implementation supports vector-valued functions. When the input is an $(n+1) \times m$ matrix, each column corresponds to pointwise values of an $m$-component vector function at the $n+1$ Chebyshev nodes (in reverse order).

```matlab
f = @(x) [cos(5*pi*x)./(4+sin(4*pi*x)), ...
          sin(5*pi*x)./(4+cos(4*pi*x))];

I1 = CleCurExpRule(f, -20+15i, 'NumberOfNodes', 64, 'EndPoint', 3)

% Equivalently, using pre-computed values:
nNodes = 64; 
b = 3;  
t = linspace(0, pi, nNodes+1); 
t = t(:);
t = (1 + cos(t)) * b / 2;

I2 = CleCurExpRule(f(t), -20+15i, 'NumberOfNodes', 64, 'EndPoint', 3)

% Output:
% I1 = 0.0076 + 0.0024i   0.0025 + 0.0039i
% I2 = 0.0076 + 0.0024i   0.0025 + 0.0039i
```

---

### Example 4: Error Estimation

For **even values of `nNodes`**, an error estimate is provided by comparing the result with the quadrature rule using `nNodes/2`:

```matlab
f = @(x) [cos(5*pi*x)./(4+sin(4*pi*x)), ...
          sin(5*pi*x)./(4+cos(4*pi*x))];

[I2, errorEstimate] = CleCurExpRule(f, -20+15i, 'NumberOfNodes', 96, 'EndPoint', 3)

% Output:
% I2 = 0.0076 + 0.0024i   0.0025 + 0.0039i
%
% errorEstimate =
%   1.0e-07 *
%   0.1982 - 0.0185i  -0.1112 + 0.0003i
```

---

## Code Structure

- **`CleCurExpRule.m`**  
  Main entry point. Computes the Clenshaw–Curtis exponential quadrature rule.

- **`computeWeights.m`, `computeRhoMax.m`**  
  Core routines defining the quadrature scheme.

- **`private/`**  
  Internal numerical utilities (DCTs, tridiagonal solvers). Not intended to be called directly by users.

- **`vpa/`**  
  Variable-precision backend (symbolic arithmetic), mainly for validation and reference computations.

---

## References

For more information, see:  
**arXiv preprint:** [https://www.arxiv.org/abs/2503.08169](https://www.arxiv.org/abs/2503.08169)

---

## Contact

**Author:** Victor Domínguez  
**Email:** victor.dominguez@unavarra.es  
**Date:** 2nd February 2026

---

## Copyright

Copyright © 2026 Victor Domínguez

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
