#' produce a static view of target boxplot
#' @importFrom png readPNG
#' @importFrom grid grid.raster
#' @export
plotDemo = function() {
 im = png::readPNG(system.file("png/ormdl3Demo.png", package="BiocRnaHap"))
 grid::grid.raster(im)
}
