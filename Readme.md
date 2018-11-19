# Modified Bayesian Change Detection

Enclosed is modified Bayesian change detection (BCD) code for 1 dimensional time series data.  This code was used for label correction and event detection in electricity disaggregation (see paper citation below).  

## Code modifications

There are two main modifications from the original BCD code provided by Adams/MacKay.

1.  Reduction of quadratic space/time complexity to linear.

2.  Change from online changepoint detection to offline to increase accuracy.

## Running the code

### Prerequisites

Matlab (tested on 2016b.  Should also work on other versions.  Should work on Octave as well, but hasn't been tested...)

### Example

There is an included example script named 'Example.m.'  This can be run with the included time series 'BreadMaker.mat.'  To run it, simply execute the script in the Matlab interpreter.

```
Example.m
```

### Main function

The main function is 'Get_ChangePoints_Windowed.m,' which requires a number of parameters.  The domain-agnostic settings we used were provided by the original authors and not tuned for any specific data stream.  They include:

Starting statistical values (see Adams/MacKay 2007 paper below):

```
params.mu0    = 0;
params.kappa0 = 1;
params.alpha0 = 1;
params.beta0  = 1;
```

Function characteristics (here we use a constant hazard function with base probability 0.005):

```
lambda        = 200;
hazard_func   = @(r) constant_hazard(r, lambda);
```

And finally, the window segmenting size:

```
params.WINDOW_LENGTH = 1000;
```

The window size is necessary to reduce the time and space complexity from quadratic to linear time for sufficiently large time series (see Section 3.3 of Valovage-AAMAS2017 paper below).  

To run the function with a different time series, say a 1-dimensional array called 'my_time_series', pass the variables into the function along with the time series and it will find the changepoints:

```
ChangePoint_Indexes = Get_ChangePoints_Windowed(params.mu0, ...
                                                params.kappa0, ...
                                                params.alpha0, ...
                                                params.beta0, ...
                                                lambda, ...
                                                hazard_func, ...
                                                my_time_series, ...
                                                size(my_time_series,1), ...
                                                params.WINDOW_LENGTH);
```

## Built With

* [Bayesian Change Detection](http://hips.seas.harvard.edu/content/bayesian-online-changepoint-detection) - Original Bayesian change detection code.

## Authors

* **Mark Valovage** - *Initial work*

* **Akshay Shekhawat** - *Assistance with testing/github posting*

## Paper citations

If you use this code please cite the following papers (.bib files included):

* M. Valovage and M. Gini. Label correction and event detection for electricity disaggregation. In Proceedings of the 2017 International Conference on Autonomous Agents & Multiagent Systems, 2017.
* R. P. Adams and D. J. MacKay. Bayesian online changepoint detection. arXiv preprint arXiv:0710.3742, 2007.