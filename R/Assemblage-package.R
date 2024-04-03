#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

usethis::use_description(fields = list(Title='Assemblage Regression',Version='1.0.1',License ='GPLv3 License',Authors='Philippe Goulet Coulombe, Karin Klieber, Christophe Barrette, and Maximilian Göbel', Maintainer= 'PGC <goulet_coulombe.philippe@uqam.ca> & CB <christophebarrette@gmail.com>' ,Description='The concept of assemblage originates from the research paper titled "Maximally Forward-Looking Core Inflation" (2024) authored by Philippe Goulet Coulombe, Karin Klieber, Christophe Barrette, and Maximilian Göbel. Assemblage Regression represents a specialized form of generalized nonnegative ridge regression, designed to optimize the weights of subcomponents to maximize the predictive capability of the aggregate.',"Authors@R" =NULL) )

usethis::use_package("glmnet")
usethis::use_package("pracma")
usethis::use_package("CVXR")
usethis::use_package("foreach")
usethis::use_package("doParallel")
usethis::use_package("stats")
usethis::use_package("methods")
usethis::use_package("Matrix")
usethis::use_package("iterators")
usethis::use_package("datasets")
usethis::use_package("base")
usethis::use_package("utils")
#devtools::build_readme()
#usethis::use_news_md()
#usethis::use_readme_rmd()
