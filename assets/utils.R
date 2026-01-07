# function to get threshold based on number of median absolute deviations
mad_threshold <- function(x, nmads, allow_negative_values = FALSE) {
  x <- x[!is.na(x)]
  threshold <- c(
    lower = (median(x) - (nmads * mad(x))),
    upper = (median(x) + (nmads * mad(x)))
  )
  
  if (!allow_negative_values & threshold["lower"] < 0) {threshold["lower"] <- 0}
  
  return(threshold)
}

# function to get point density for plotting
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}