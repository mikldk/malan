# malan: MAle Lineage ANalysis

An R package (<https://www.r-project.org/>) to perform **MA**le **L**ineage **AN**alysis 
by simulating genealogies backwards and 
imposing short tandem repeats (STR) mutations forwards. 
Intended for forensic Y chromosomal STR (Y-STR) haplotype analyses. 
Numerous analyses are possible, e.g. number of matches and meiotic distance to matches.

See documentation included in package (vignettes and manual) at <https://mikldk.github.io/malan/>.

## Installation

You first need `R` (<https://www.r-project.org/>). 
Then you can install `malan` from GitHub by using the `remotes` package (<https://CRAN.R-project.org/package=remotes>):

``` r
# install.packages("remotes")
remotes::install_github("mikldk/malan")
```

## References

Andersen MM, Balding DJ (2017). *How convincing is a matching Y-chromosome profile?*. 
PLoS Genet 13(11): e1007028. <https://doi.org/10.1371/journal.pgen.1007028>.

## Disclaimer

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## License

License: GPL (>= 2).

## Build status

Travis CI:

[![Travis-CI Build Status](https://travis-ci.org/mikldk/malan.svg?branch=master)](https://travis-ci.org/mikldk/malan)

