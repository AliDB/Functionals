
#Functionals (Functions of Functions)


In this tutorial, we are going to show how funcionals (functions of functions) can be used in our programming. For that reason we have implemented diff function. This tutorial was part of my project for Statistical Computing course with Prof. Harner.

Newton's method for finding a root of a differentiable function *f* takes a guess *y* and computes hopfully an improved guess as:
$$ y - \\frac{f(y)}{Df(y)}$$
 where *D**f* denotes the derivative of *f*.

#### 1.  Create a function called `newton_search` with four arguments: `f`, `df`, `guess`, `conv` (the convergence criterion).

``` r
newton_search <- function(f, df, guess, conv=0.001) {
# Put your function body here.
  # f= sin;df= cos; guess=3.2; conv=.00001
  
  makeNewGuess <- function(guess){
      newGuess <- guess - f(guess)/df(guess)
      newGuess
  }
  
  iter <- 0;
  while( abs(makeNewGuess(guess) - guess) >= conv && iter < 100 ){
    guess <- makeNewGuess(guess)
    iter <- iter +1;
    if (iter==100) print(" we couln't reach to the convergence criterion in 100 attempts")
  }
  guess
}

## There are different types of stopping criterion. 
#1) when $x_{n+1}$ - $x_n$ is sufficiently small. 
#2) when $\abs(F(x_n))$ sufficiently small in some sense for some  $x_{n}$
# I chose the first one here. I assumed that our function has simple format.
```

Hint: Define a local functions (or helper function) to compute the improvement and then test for convengence.

#### 1.  Use this function to find the root of *s**i**n*(*x*) near 3 using the actual symbolic derivative. The exact answer is *π*.

``` r
newton_search(sin,cos,3,.001)
```

    ## [1] 3.142547

#### 1.  In general you may not be able to compute the derivative exactly. Use the symmetric difference quotient to approxiate the derivative of *f* at *x* numerically by the defintion:
    $$ Df \\approx \\frac{f(x + h) - f(x - h)}{2h} $$
     for small *h*.

#### efine a function `make_derivative` with arguments `f` and `h`. The result returned should be a function closure that remembers both `f` and `h`.

``` r
make_derivative <- function(f, h) {
    # Put your function body here.
  make_derivative_closure <<- function(x){
    df <- (  (f(x+h)-f(x-h)) / (2*h) )
    df  
  }
  make_derivative_closure
}
```

#### 1.  Find the root of *s**i**n*(*x*) near 3 using numerical derivatives.

``` r
# Put your R code here.
rootSine <- make_derivative(sin,h=1e-6)

root <- 3 - sin(3) / rootSine(3) ; root
```

    ## [1] 3.142547

#### 1.  The log-likelihood of the gamma distribution with scale parameter 1 can be written as:
    (*α* − 1)*s* − *n*log*Γ*(*α*)
     where *α* is the shape parameter and *s* = ∑log*X*<sub>*i*</sub> is the sufficient statistic.

Randomly draw a sample of *n* = 30 with a shape parameter of *α* = 4.5. Using `newton_search` and `make_derivative`, find the maximun likelihood estimate of *α*. Use the moment esitmator of *α*, i.e., $\\bar X$ as the intial guess. The log-likelihood function in R is:

``` r
x <- rgamma(n=30, shape=4.5)

gllik <- function() {
  s <- sum(log(x))
  n <- length(x)
  s/n
  function(a) {
    (a - 1) * s - n * lgamma(a)
  }
}

guess <- mean(x); guess; 
```

    ## [1] 4.638775

``` r
LL <- gllik();
#LL2 <- LL(4.5)

df1 <- make_derivative(LL,h=1e-6)
df2 <- make_derivative(df1 ,h=1e-6)

newton_search(df1,df2,guess,conv=.001)
```

    ## [1] 4.569934

``` r
# Put more R code here.
```

Hint: You must apply `newton_search` to the first and second derivatives (derived numerically using `make_derivative`) of the log-likelihood. Your answer should be near 4.5.
