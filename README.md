# diy-regression
Armed with nothing other than ANSI C (c99), here's some code to perform OLS regression and report coefficients and their t-stats.

# example
Look in `test_matrix.c` for an example regression set up. Creating the target variable vector and a covariate matrix is easiest (for small test cases) with the `make_matrix` variadic function:

```c
matrix *y = make_matrix(10, 1,
                         -0.05911,   
                          8.11098065,
                          2.24198061,
                         -8.0514162,
                         -8.07322677,
                         10.44922017,
                          8.82230702,
                         -3.18227219,
                          4.66138728, 
                        -21.56386711);

matrix *x = make_matrix(10, 2,
                        1.0, -0.18380953,
                        1.0,  0.54366253,
                        1.0,  0.14948514,
                        1.0, -0.99259025,
                        1.0, -0.82037301,
                        1.0,  0.79391838,
                        1.0,  0.91688518,
                        1.0, -0.51775223,
                        1.0,  0.45836676,
                        1.0, -2.18044347);

```
Note: this data was generated randomly for a Python example, then copied over to C -- this way I could compare the regression results with a baseline implementation in `pandas`.

Allocate a data structure that will hold the output results of regressing `y` on `x`:

```c
regression_fit *rfit = new_regression_fit();
```

Finally, use the function `lstsq` to actually populate `rfit` with the results.

```c
lstsq(y, x, rfit);
```

When done, use `print_regression_fit` to display a summary. Compile the code and run:

```bash
$ gcc -std=c99 test_matrix.c matrix.c -o test_matrix -g -L. -lm 
$ ./test_matrix 
-------------------------Summary of Regression Analysis-------------------------

Number of Observations:         10
Number of Degrees of Freedom:   2

R-squared:         0.9876
Adj R-squared:     0.9861

-----------------------Summary of Estimated Coefficients------------------------
     Variable       Coef      t-stat     
--------------------------------------------------------------------------------
     Column_0     1.2025      3.1898     
     Column_1     10.1868      25.2242     
---------------------------------End of Summary---------------------------------

```

# motivation
The [formula for OLS regression](https://en.wikipedia.org/wiki/Ordinary_least_squares#Estimation) is deceptively simple. The bare bones matrix functionality needed is: matrix creation, matrix transposition, matrix multiplication, and matrix inversion. 

It's a good exercise to prove you can write all of these things yourself. I chose simple algorithms, like Gaussian elimination with pivoting to perform matrix inversion. I did not do fancy adjustments for precision optimization. I did not put extra effort into making the printing facilities polished. It's just a bare-bones implementation.

But it's very rewarding to say that every single step, from the definition of a matrix data structure to the algorithm for matrix inversion to the computation of regression coefficient t-stats, was programmed directly by me, using nothing outside of the C standard library.

Further, to carry this all the way to obtaining t-stats for the regression coefficients will force you to confront the OLS formulation at a low level. I remember needing to do this many times in graduate school and it's good to review from time to time. For example, since you will need to manually compute the standard error for each coefficient, you have to be comfortable understanding the covariance matrix of the estimated coefficient vector. This makes it easier to understand extensions, like Seemingly Unrelated Regression in which possible correlation between error terms means that computing linear combinations of estimated coefficients (such as an average) requires the off-diagonal covariance terms from the coefficient covariance matrix.

# lessons
I have not had the chance to use C a great deal in professional experiences, and even when I have, it's mostly been through Cython. Writing bigger C programs always gives a good learning opportunity and I present a few lessons below:

1. Use `malloc` when allocating a new struct pointer. You may think you can simply declare a struct pointer, but the problem is that if it's not initialized, it's a "junk" pointer, and so if you attempt to manually access struct fields later to initialize, you might be accessing invalid memory. This means you also need to manually free your struct pointers.

2. There's no great solution for using return error codes when your return types have wide domains. If something should return a non-negative integer, then you could return a signed integer and let the sign represent error state. But what if your function returns a pointer, or a double? You could return NULL or a sentinel value, but this is sort of ugly and it only gives you a single sentinel value, so no way of indicating multiple error types. You can use this in conjunction with setting an error state with the error.h API, but you may not want to rope in all of that machinery. Another option is to create custom data structures for each different result type, extra fields for error states or messages. Of course, then you're committing to creating and maintaining structs for every function call, which could be tedious or make for an unreasonably complicated API.

3. stdarg-based variadic arguments do not do any type promotion or implicit casting. If you give a literal of "1" and ask for it to be a "double", it fails. And `va_arg` fails silently and moves on to the next item in the argument list, so you can get unexpected behavior. You have to provide exactly the right types. If you end up debugging a complicated variadic function, one simple but helpful tip is to just print the arguments out after `va_arg` obtains them. If things disappear from what you expect, it's likely a type error. 

