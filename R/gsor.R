##' Open a GUESS smart output netcdf file
##'
##' Open a GUESS smart output netcdf file
##' @title Open a GUESS smart output netcdf file
##' @param filename 
##' @param ... additional parameters passed to \code{nc_open}
##' @return an object of class gso
##' @import ncdf4
##' @import dplyr
##' @import lubridate
##' @export
gso_open <- function(filename, ...) {
  nc <- nc_open(filename, ...)
  
  if (nc$format == "NC_FORMAT_64BIT")
  {
    format <- "nc3"
  }
  else if (nc$format == "NC_FORMAT_NETCDF4")
  {
     format <-"nc4"
  }
  else
  {
    stop("Unsuppored netcdf filetype")
  }
  
  delim = ""
  
  if (format == "nc4") 
  {
    lon <- ncvar_get(nc, "Base/Longitude")
    lat <- ncvar_get(nc, "Base/Latitude")
    pfts <- ncvar_get(file_new, "Base/Pfts")
    delim = "/"
  }
  
  else  
  {
    lon <- ncvar_get(nc, "Base_Longitude")
    lat <- ncvar_get(nc, "Base_Latitude")
    pfts <- unname(unlist(ncatt_get(nc, "Base_Pfts")))
    delim = "_"
  }
  
  cellid <- 1:length(lon)

  
  if (grepl("Monthly", filename, fixed = T))
  {
    resolution = "monthly"
  }
  
  else if(grepl("Annually", filename, fixed = T))
  {
    resolution = "annually"
    timevarstr = paste("Base", "Time", sep = delim)
    starttimestr = ncatt_get(nc, timevarstr)$units
    starttime =unlist(strsplit(starttimestr, split="years since ", fixed=TRUE))[2]
    timevals = ncvar_get(nc, timevarstr)
    time <- (dmy(starttime) + years(timevals)) %>% year
  }
  
  else  
  {
    stop("currently only annual or monthly resolution implemented!")
  }
  
  varnames = names(nc$var)

  o <- list(
    netcdf = nc,
    resolution = resolution,
    format = format,
    time = time,
    lon = lon,
    lat = lat,
    pfts = pfts,
    cellid = tibble(cellid, lon, lat)
  )
  class(o) <- "gso"
  o
}

##' Break a long line for pretty printing
##' @keywords internal
break_print_line <- function(line) {
  n <- nchar(line)
  crit <- 50
  s <- ""
  while (n > crit) {
    whitespaces <- lapply(strsplit(line, ''), function(x) which(x ==
                                                                  " "))[[1]]
    breakat <- which(whitespaces > crit)[1]
    breakpos <- whitespaces[breakat]
    s <- paste0(s, substr(line, 1, breakpos), "\n\t\t")
    line <- substr(line, (breakpos + 1), n)
    n <- nchar(line)
    if (is.na(n)) n <- 0
  }
  s <- paste0(s, line)
  s
}
  
##' @param x an object of class gso
##' @param ... ignored
##' @export
print.gso <- function(x, ...) {
  cat("GUESS smart output object\n")
  cat("File:\t\t", x$netcdf$filename, "\n")
  cat("Resolution:\t", x$resolution, "\n")
  cat("Pfts:\t\t", x$pfts, "\n")
  cat("Variables:\t", break_print_line(paste(names(x$netcdf$var),
                                             collapse = " ")), "\n")
}

##' Extract variables from GUESS smart output into a tibble
##'
##' .. content for \details{} ..
##' @title Extract variables from GUESS smart output into a tibble
##' @param gso a gso object as return from gso_open
##' @param vars character vector with names of variable to extract
##' @param newnames optional character vector for renaming the variables 
##' @param silent supress status updates
##' @return a tibble
##' @importFrom reshape2 melt
##' @import ncdf4
##' @import dplyr
##' @export
gso_getvar <- function(gso, vars, newnames = NULL, silent = FALSE) {
  if (is.null(newnames)) newnames <- vars
  if (length(newnames) != length(vars)) {
    stop("Length of new names does not match number of variables.")
  }
  if (gso$resolution != "annual") {
    stop("currently only annual resolution implemented!")
  }
  for (i in seq_along(vars)) {
    newname <- newnames[i]
    if (!silent) {
      cat("Extracting", vars[i], "as", newname, "...\n")
    }
    variable <- ncvar_get(gso$netcdf, vars[i])
    dimnames(variable) <- list(pft = gso$pfts,
                               year = gso$time,
                               cellid = gso$cellid$cellid)
    variable_m <- as_tibble(melt(variable, value.name = newname)) %>%
      left_join(gso$cellid, by = "cellid") %>%
      select(lon, lat, year, pft, !!newname, -cellid)
    if (i == 1) {
      M <- variable_m
    } else {
     if (!silent) {
       cat("Merging...\n")
     } 
      M <- left_join(M, variable_m, by = c("lon", "lat", "year", "pft"))
    }
  }
  M
}