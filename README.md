# malan: MAle Lineage ANalysis

An R package (<https://www.r-project.org/>) to perform **MA**le **L**ineage **AN**alysis 
by simulating genealogies backwards and 
imposing short tandem repeats (STR) mutations forwards. 
Intended for forensic Y chromosomal STR (Y-STR) haplotype analyses. 
Numerous analyses are possible, e.g. number of matches and meiotic distance to matches.

## Installation

You first need `R` (<https://www.r-project.org/>). 
Then you can install `malan` from GitHub by using the `remotes` package (<https://CRAN.R-project.org/package=remotes>):

``` r
# install.packages("remotes")
remotes::install_github("mikldk/malan")
```

### For Mac OS users

Some Mac OSX configurations may have installation problems (for example missing `gfortran` libraries). In such cases, it may help to use Conda package management system. See also the discussion in [issue #14](https://github.com/mikldk/malan/issues/14).

## Getting started

See documentation included in package (vignettes and manual) at <https://mikldk.github.io/malan/>. The introduction vignette is available at <https://mikldk.github.io/malan/articles/introduction.html>.

You can also get an overview of the included vignettes by the following `R` command:

```r
vignette(package = "malan")
```

To read a vignette, type:

```r
vignette("introduction", package = "malan")
```

### Running tests

Note that to also install the tests, you need to install the package as follows:

``` r
# install.packages("remotes")
remotes::install_github("mikldk/malan", INSTALL_opts="--install-tests")
```

You can now run the tests:

``` r
library(malan)
library(testthat)
test_package('malan')
```

## Contribute, issues, and support

Please use the issue tracker at <https://github.com/mikldk/malan/issues> 
if you want to notify us of an issue or need support.
If you want to contribute, please either create an issue or make a pull request.

## Dependencies

This package depends on R (<https://www.r-project.org/>) and the following R packages: 
`Rcpp`, `RcppProgress`, `RcppArmadillo`, `igraph`, `tibble`, `magrittr`, `dplyr`, and `tidygraph`.

## References

Andersen MM, Balding DJ (2017). *How convincing is a matching Y-chromosome profile?*. 
PLoS Genet 13(11): e1007028. <https://doi.org/10.1371/journal.pgen.1007028>.

## Disclaimer

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## License

License: GPL-2.

## Badges

The Journal of Open Source Software:

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00684/status.svg)](https://doi.org/10.21105/joss.00684)

Zenodo: 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1241769.svg)](https://doi.org/10.5281/zenodo.1241769)

Travis CI:

[![Travis-CI Build Status](https://travis-ci.org/mikldk/malan.svg?branch=master)](https://travis-ci.org/mikldk/malan)


