## Title
Improving DomainMatrix class and Matrices module

## Table of Contents
- [About Me](#about-me)
  - [Personal Information](#personal-information)
  - [Personal Background](#personal-background)
  - [Programming Skills](#programming-skills)

- [Contributions to SymPy](#contributions-to-sympy)
  - [Pull Requests](#pull-requests)
  - [Issues](#issues)

- [The Project](#the-project)
  - [Problem & Motivation](#problem--motivation)
  - [Theory ](#theory)
    - [Bareiss Algorithm](#Bareiss-Algorithm)
    - [Sparse Fraction-Free Algorithm](#Sparse-Fraction-Free-Algorithm)
    - [Dense v/s Sparse](#Dense-v/s-Sparse-Representation)
    - [Cholesky Decomposition](#Cholesky-Decomposition)
    - [QR Decomposition](#QR-Decomposition)

- [Timeline](#timeline)

- [References](#references)
<p>&nbsp;</p>

## About Me

### Personal Information
Heading | Details
--------|--------
**Name**                  |   Kartik Sethi
**University**            |   [Indian Institute of Technology, Guwahati](http://www.iitg.ac.in/)
**Email**                 |   Kartiks31416@gmail.com
**GitHub/Gitter**     |   [Kartik Sethi](https://github.com/ks147)
**LinkedIn**              |   [Kartik Sethi](https://www.linkedin.com/in/kartik-sethi-iit-guwahati/)
**Timezone**              |   IST (UTC+5:30)




### Personal Background
I am Kartik Sethi, currently enrolled in a four-year B.Tech program at Indian Institute of Technology, Guwahati. I am majoring in Mathematics and Computing. I have prior experience as a SDE intern at IBM and also as a quantative researcher in a startup. I am currently in my 8th semester of college, I have an extensive math background as a result of my degree and have taken the following courses in recent times.

- [Matrix Computations](https://github.com/ks147/MA423_Matrix_Computations)
- Abstract Algebra (Covering Domains, Rings, Fields etc.)
- Real and Complex Analysis
- [Numerical Analysis of ODEs and PDEs](https://github.com/ks147/MA322_Scientific-Computing)


### Programming Skills
* I work on _Ubuntu 20.04.2.0 LTS (Focal Fossa)_. I prefer _PyCharm_ as my
text editor for working on Python projects. 
* I have been programming for about 5 years. I am proficient in _C/C++_ and _Python_. My [GitHub
Page](https://github.com/ks147) has details about all the projects I've taken up.

<p>&nbsp;</p>

## Contributions to Sympy

### Pull Requests
- Merged
  - [#21015](https://github.com/sympy/sympy/pull/21015) Added function zeros to DomainMatrix class
  - [#20922](https://github.com/sympy/sympy/pull/20922) Implemented function returning Wilkinson Matrix
  - [#20911](https://github.com/sympy/sympy/pull/20911) Implemented function returning upper and lower triangular part of a matrix 
  - [#20829](https://github.com/sympy/sympy/pull/20829) Implemented Upper Hessenberg Decomposition
  - [#20761](https://github.com/sympy/sympy/pull/20761) Implemented Singular Value Decomposition
  - [#20651](https://github.com/sympy/sympy/pull/20651) Implemented Matrix Student's t-distribution
  - [#20607](https://github.com/sympy/sympy/pull/20607) Implemented Logit Normal distribution
- Open
  - [#20992](https://github.com/sympy/sympy/pull/20992) Added Documentation to DomainMatrix Class
  - [#21060](https://github.com/sympy/sympy/pull/21060) Improving `Matrix.LUdecomposition()` using `DomainMatrix.lu()` 


### Issues
- [#20961](https://github.com/sympy/sympy/issues/20961) Complete Singular Value Decomposition
- [#20990](https://github.com/sympy/sympy/issues/20990) Domain Matrix convert_to function unable to convert from QQ to RR 
## The Project

### Problem and Motivation
The aim is to improve the DomainMatrix class and use the functionality of DomainMatrix class to improve algorithms implemented in [Matrix](https://docs.sympy.org/latest/modules/matrices/index.html) module of Sympy. I was motivated by issue [#20987](https://github.com/sympy/sympy/issues/20987). I feel I have the experience and expertise to improve this module.

List of things I intend to implement/improve
- Add `DomainMatrix.__getitem__` for slicing and element access
  using this function add other utilities like:
  - `DomainMatrix.diag()`, which returns diagonal elements of a matrix
  - `DomainMatrix.upper_triangular()` and `DomainMatrix.lower_triangular()` which returns the upper and lower triangular part of a Matrix
  - `DomainMatrix.permute()`, which returns a permutation of rows/columns of the given input matrix
  
- Implement Bareiss algorithm for `DomainMatrix.det`
- Implement Sparse Fraction Free algorithm
- Add a mechanism for deciding between dense and sparse representation.
- Implement Cholesky and QR decomposition
- Figure out methods that are faster for DomainMatrix and transport those methods to Matrix module to improve speed of computation.
- Improve Documentation of `DomainMatrix`, `DDM` class and `Domains`

#### Bareiss Algorithm[1]
Given a matrix of integers, we wish to compute the determinant using a method
that does not introduce fractions. This fraction-free algorithms make perfect sense for a symbolic math library where we usually want to avoid working with floats. Fraction-free algorithms like Bareiss ensure that any divisions that occour during computation are exact.

The algorithm is given by

![Bareiss Algo](/images/Bareiss_algo.png)

#### Sparse Fraction Free Algorithm[2][3]
The algorithm is a fraction-free Gaussian elimination for Sparse Matrices. It is essentially modified Bareiss' algorithm for sparse matrices.

![SFF_1](/images/SFF_1.png)
![SFF_2](/images/SFF_2.png) 

#### Dense v/s Sparse Representation
Some of the matrix algorithms scale very poorly for large matrices, with time complexity ~ O(n^3). There is a way to speed up the calculation if we know the matrix is sparse i.e. most of the entries are 0. There is no fixed criteria for a matrix to be classified as sparse but usually when the number of non-zero entries is ~ O(m + n), where m is number of columns and n = number of rows.

#### Cholesky Decomposition[4]
The Cholesky decomposition of a Hermitian positive-definite matrix A, is a decomposition of the form

A = L L*

where L is a lower triangular matrix with real and positive diagonal entries, and L* denotes the conjugate transpose of L.

There are various methods for calculating the Cholesky decomposition. The computational complexity of commonly used algorithms is O(n^3) in general.The algorithms described below all involve about (1/3)n^3 FLOPs (n^3/6 multiplications and the same number of additions) for real matrices and (4/3)n^3 FLOPs for complex matrices.

The Cholesky–Banachiewicz algorithm starts from the upper left corner of the matrix L and proceeds to calculate the matrix row by row.

```python

for (i = 0; i < dimensionSize; i++) {
    for (j = 0; j < (i + 1); j++) {
    float sum = 0;
    for (k = 0; k < j; k++)
        sum += L[i][k] * L[j][k];

    if (i == j)
        L[i][j] = sqrt(A[i][i] - sum);
    else
        L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
    }
}
```

The Cholesky–Crout algorithm starts from the upper left corner of the matrix L and proceeds to calculate the matrix column by column.

```python
for (j = 0; j < dimensionSize; j++) {
    float sum = 0;
    for (k = 0; k < j; k++) {
        sum += L[j][k] * L[j][k];
    }
    L[j][j] = sqrt(A[j][j] - sum);

    for (i = j + 1; i < dimensionSize; i++) {
        sum = 0;
        for (k = 0; k < j; k++) {
            sum += L[i][k] * L[j][k];
        }
        L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
    }
}
```


#### QR Decomposition[5][7]
QR decomposition for a matrix A is given by,

For square matrix A
A = QR, where Q is an orthogonal matrix and R is upper triangular

For rectangular matrix A, 
A = QR, where Q is left-orthogonal

There are two broad algorithms for finding the QR decomposition, Gram-Schmidt method and Householder reflector method, both of these take around the same number of operations/flop-count. The time complexity for both is ~ O(n^3). Currently Sympy uses Gram-Schmidt method which works better for symbolic matrices. It is possible that for `DomainMatrix` Householder method is faster. There is an optimization for Hessenberg matrices that can find out the QR decomposition in O(n^2). I implemented and tested this optimization for `Matrix` module but the time reduction was not significant. I believe for `DomainMatrix` this optimization could be useful.

![CGS](/images/CGS_1.png)
![CGS](/images/CGS_2.png)
![CGS](/images/MGS_1.png)

A matrix is said to be _upper hessenberg_ if it has nonzeros only in the upper triangle and first subdiagonal. The algorithm for finding QR decomposition for an upper Hessenberg matrix is given below.

```matlab
[H,Q] = hessred(A)

Compute the Hessenberg decomposition H = Q’*A*Q using
Householder transformations.

function [H,Q] = hessred(A)

n = length(A);
Q = eye(n);  Orthogonal transform so far
H = A;
Transformed matrix so far

for j = 1:n-2

  % -- Find W = I-2vv’ to put zeros below H(j+1,j)
  u = H(j+1:end,j);
  u(1) = u(1) + sign(u(1))* norm(u);
  v = u/norm(u);

  % -- H := WHW’, Q := QW
  H(j+1:end,:) = H(j+1:end,:)-2*v*(v’*H(j+1:end,:));
  H(:,j+1:end) = H(:,j+1:end)-(H(:,j+1:end)*(2*v))*v’;
  Q(:,j+1:end) = Q(:,j+1:end)-(Q(:,j+1:end)*(2*v))*v’;

end

end
``` 

## Timeline
I can easily give 7-8 hours to work on Sympy everyday throughout the summer. I will be mostly free throughout May - August. I'll also continue to contribute throughout the month of April as I have several ideas for new PRs and some pending PRs I want to improve and helpfully get merged.
### Community Bonding
I will use this time to learn more about Sympy's DomainMatrix and Domain classes and also try and figure out which all functions in `Matrix` module could benefit from using `DomainMatrix`.

### Week 1,2,3 (June 7 - June 27)
I will implement the following
- `DomainMatrix.__getitem__` and some other important functions like `DomainMatrix.diag()`, `DomainMatrix.upper_triangular()`, `DomainMatrix.lower_triangular()`.
- Figure out and implement a criteria for assigning Sparse or Dense representation to a Matrix
- Bareiss fraction free algorithm for getting determinant
### Week 4,5 (June 28 - July 11)
- Sparse Fraction Free algorithm

### Week 6,7 (July 11 - July 25)
- Implement QR decomposition and optimization for Hessenberg Matrices
- Implement Cholesky Decomposition

### Week 8, 9 (July 26 - August 8)
- Improve the Documentation for `DomainMatrix` class,`DDM` class and `Domains` 

### Week 10 (August 8 - August 15)
- Buffer week

## References 
- [FRACTION-FREE METHODS FOR DETERMINANTS, Deanna Richelle Leggett][1]
- [Sparse Fraction Free Algorithm pull request][2]
- [Fraction free Guassian elimination for Sparse Matrices][3]
- [https://en.wikipedia.org/wiki/Cholesky_decomposition][4]
- [Fundamentals of Matrix Computations (Second Edition) Chapter 3, DAVID S. WATKINS][5]
- [https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-10-21.pdf][7]
- [compare-gram-schmidt-and-householder-orthogonalization-algorithms][6]

[1]: https://aquila.usm.edu/cgi/viewcontent.cgi?article=1001&context=masters_theses
[2]: https://github.com/sympy/sympy/pull/9133
[3]: https://www.eecis.udel.edu/~saunders/papers/sffge/it5.ps
[4]: https://en.wikipedia.org/wiki/Cholesky_decomposition
[5]: https://epubs.siam.org/doi/10.1137/1035117
[6]: https://blogs.mathworks.com/cleve/2016/07/25/compare-gram-schmidt-and-householder-orthogonalization-algorithms/?s_tid=answers_rc2-1_p4_BOTH
[7]: https://www.cs.cornell.edu/~bindel/class/cs6210-f16/lec/2016-10-21.pdf
