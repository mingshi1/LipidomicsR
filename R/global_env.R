
#' Global environment settings of LipidomicsR
#' @description Loading environment for LipidomicsR. Please do not call it directly.
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom graphics legend title
#' @importFrom stats dist median na.omit p.adjust prcomp sd shapiro.test t.test
#' @importFrom utils combn fix menu read.csv read.table write.csv
#' @importFrom tidyverse tidyverse_packages
#' @return No return value, called for loading environment.
lipidomicsR_env <- function() {
  grDevices::colorRampPalette(colors = 'black')
  grDevices::dev.off()
  grDevices::pdf()
  graphics::legend()
  graphics::title()
  stats::dist()
  stats::median()
  stats::na.omit()
  stats::p.adjust()
  stats::prcomp()
  stats::sd()
  stats::shapiro.test()
  stats::t.test()
  utils::combn()
  utils::fix()
  utils::menu()
  utils::read.csv()
  utils::read.table()
  utils::write.csv()
  tidyverse::tidyverse_packages()
}
