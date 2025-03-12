# CleCurExpRule - Clenshaw-Curtis Quadrature with Exponential Weight

This MATLAB project implements Clenshaw-Curtis quadrature rules for computing integrals of the form:

int(f(s) exp(z s), s, 0, b) 

## Installation

Clone the repository and add the files to your MATLAB path:
## Functions Overview

### Core Functions
- `CleCurExpRule` – Computes the integral using the standard Clenshaw-Curtis rule.

### Variable-Precision (VPA) Versions
- `CleCurExpRule_vpa` – VPA version of `CleCurExpRule`. Only for testing purposes
- `computeWeights_vpa` – VPA version of `computeWeights`.
- `computeRhoMax_vpa` – VPA version of `computeRhoMax`.

### Auxiliary Functions

- `computeWeights` – Computes the quadrature weights.
- `computeRhoMax` – Computes stability-related parameters.
- `idctI`, `idctII` – Implements inverse discrete cosine transforms.

### Auxiliary Functions (VPA) Versions

- `computeWeights_vpa` – VPA version of `computeWeights`.
- `computeRhoMax_vpa` – VPA version of `computeRhoMax`.
- `mifft_vpa` – Performs inverse FFT with variable-precision.
- `idctI_vpa`, `idctII_vp` – Implements inverse discrete cosine transforms.
- `thomas_algorithm` solver for tridiagonal linear system

## Usage

To compute an integral in [0,2] using the Clenshaw-Curtis rule:

```matlab
f = @(x) 4*cos(5*pi*s)./(4+sin(4*pi*s)
CleCurRule(f,-50+20i,'NumberOfNodes', 48)
```

Or alternatively
```matlab
m = 48
t =linspace(0,pi,m+1);
y = 1+cos(t(:))
CleCurRule(f(t),-50+20i)
```
Other examples: integral in [0,4]
```matlab  
CleCurRule(f,-210+122i,'NumberOfNodes', 65,'EndPoint',4)
```

## More Information:

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## More Information 

See the paper: 

https://www.arxiv.org/abs/2503.08169

## Author

**Víctor Domínguez**  
Email: victor.dominguez@unavarra.es  
ORCID: 0000-0002-6095-619X (https://orcid.org/0000-0002-6095-619X)
