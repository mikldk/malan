#' MAle Lineage ANalysis
#'
#' Simulating genealogies backwards and imposing STR mutations forwards.
#' 
#' See vignettes and manual for documentation.
#' 
#' Disclaimer:
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
#' EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
#' MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. 
#' IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING 
#' THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, 
#' TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
#' CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#'
#' @references
#'  Andersen MM, Balding DJ (2017) How convincing is a matching Y-chromosome profile? 
#'  PLoS Genet 13(11): e1007028. \doi{10.1371/journal.pgen.1007028}.
"_PACKAGE"


# To get RcppArmadillo fastLm imported into NAMESPACE such 
# that arma:: etc. can be used in C++.
#' @importFrom RcppArmadillo fastLm
RcppArmadilloFastLm <- function(...) {
  a <- RcppArmadillo::fastLm(...)
  return(a)
}

