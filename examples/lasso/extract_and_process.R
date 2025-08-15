## From Victor VderM ##
## Rec'd by email with the output/clim_extract.rds file on 12 Aug 2025 ##
## Edits by Lizzie started on 15 August 2025 ##

# Extract climate data for Lizzie
rm(list = ls())
# wd <- '/home/victor/projects/misc/lizzie'
wd <- '/Users/lizzie/Documents/git/projects/misc/miscmisc/bayesianflows/examples/lasso'

# Raw daily data come from climexp.knmi.nl
# ERA5 Europe (0.25 degree)

process_raw <- FALSE

if(process_raw){
  library(terra)
  climdir <- file.path(wd, 'data', 'era5_europe')
  
  # Read sites (from a file I found in the repo)
  sites <- data.frame(longitude = 8.71667, latitude = 50)
  sites_vect <-  vect(sites, geom=c("longitude", "latitude"))
  
  years <- 1950:2024
  vars <- list('tmean' = 't2m', 'tmin' = 'tmin', 'tmax' = 'tmax', 'tdew' = 'tdew', 
               'prec' = 'tp', 'evap' =  'evap', 'wind' = 'sfcWind')
  
  clim_df <- data.frame()
  for(v in vars){
    era5_r <- rast(file.path(climdir, paste0('era5_',v,'_daily_eu.nc')))
    era5_dat <- extract(era5_r, sites_vect, ID = FALSE)
    
    # note: we would need to tweak the below if we had several sites
    clim_df <- rbind(
      clim_df,
      data.frame(lat = sites$latitude, lon = sites$longitude, date = gsub(" .*","", time(era5_r)), var = names(vars[vars==v]), value = t(era5_dat))
    )
  }
  
  saveRDS(clim_df, file.path(wd, 'output', 'clim_extract.rds'))
}else{
  clim_df <- readRDS(file.path(wd, 'output', 'clim_extract.rds'))
}

# convert to degree celsius
clim_df[clim_df$var %in% c( "tmean","tmin","tmax","tdew"), 'value'] <- clim_df[clim_df$var %in% c( "tmean","tmin","tmax","tdew"), 'value'] - 273.15 

# define the variables you want to compute here (dates: from-to / range: use only for GDD and very simple chill units, but could be use for something else?)
# be careful, maybe it does not work for very specific cases
# but it should work most of the time! 
vars_def <-
  rbind(c(var = c('tmin'), fun = c('min'), from = '11-01', to = '03-31', range = NA, name = 'min_tmin_NDJFM'),
        c(var = c('tmean'), fun = c('mean'), from = '03-01', to = '05-31', range = NA, name = 'mean_MAM'),
        c(var = c('tmean'), fun = c('sum'), from = '03-01', to = '05-31', range = '5-31', name = 'gdd_MAM'),
        c(var = c('prec'), fun = c('sum'), from = '03-01', to = '05-31', range = NA, name = 'prec_MAM'),
        c(var = c('prec'), fun = c('sum'), from = '01-01', to = '12-31', range = NA, name = 'totalprec'),
        c(var = c('tmean'), fun = c('sum'), from = '10-01', to = '02-28', range = '0-7', name = 'chill_ONDJF'))

newclim_df <- data.frame()
for(i in 1:nrow(vars_def)){
  
  clim_here <- clim_df[clim_df$var == vars_def[i,'var'],]
  clim_here$year <- lubridate::year(clim_here$date)
  
  # should we look in previous year?
  if(as.Date(paste0('1950-', vars_def[i,'from'])) > as.Date(paste0('1950-', vars_def[i,'to']))){
    seqdates <- lapply((min(clim_here$year)+1):max(clim_here$year), 
      function(y){as.character(seq.Date(as.Date(paste0(y-1, '-', vars_def[i,'from'])), 
        as.Date(paste0(y, '-', vars_def[i,'to'])), by = 1))})
    clim_here <- clim_here[clim_here$date %in% unlist(seqdates), ]
    seqdates_prev <- lapply(min(clim_here$year):(max(clim_here$year)-1), 
      function(y){as.character(seq.Date(as.Date(paste0(y, '-', vars_def[i,'from'])), 
        as.Date(paste0(y, '-', '12-31')), by = 1))})
    clim_here[clim_here$date %in% unlist(seqdates_prev), 'year'] <- clim_here[clim_here$date %in% unlist(seqdates_prev), 'year'] + 1  # trick for aggregate later
  }else{
    seqdates <- lapply(min(clim_here$year):max(clim_here$year), 
      function(y){as.character(seq.Date(as.Date(paste0(y, '-', vars_def[i,'from'])), 
        as.Date(paste0(y, '-', vars_def[i,'to'])), by = 1))})
    clim_here <- clim_here[clim_here$date %in% unlist(seqdates), ]
  }
  
  # GDD case
  if(!is.na(vars_def[i,'range']) & grepl(pattern = 'gdd', vars_def[i,'name'], ignore.case = TRUE)){
    # apply lower and upper bounds
    bounds <- as.numeric(stringr::str_split(vars_def[i,'range'], '-')[[1]])
    if(bounds[1] >= bounds[2]){stop('Error with value range?')}
    clim_here$value <- ifelse(clim_here$value < bounds[1], bounds[1], ifelse(clim_here$value > bounds[2], bounds[2], clim_here$value))-bounds[1]
  }
  
  # chilling case
  if(!is.na(vars_def[i,'range']) & grepl(pattern = 'chill', vars_def[i,'name'], ignore.case = TRUE)){
    # apply lower and upper bounds
    bounds <- as.numeric(stringr::str_split(vars_def[i,'range'], '-')[[1]])
    if(bounds[1] >= bounds[2]){stop('Error with value range?')}
    clim_here$value <- ifelse(clim_here$value < bounds[1], 0, ifelse(clim_here$value > bounds[2], 0, 1)) # one chilling unit for every day inside the range
  }
  
  newvar <- aggregate(value ~ lat + lon + year, data = clim_here, FUN = match.fun(vars_def[i,'fun']))
  newvar$var <- vars_def[i,'name']
  newclim_df <- rbind(newclim_df, newvar)
}

## Add leafout date and daylength

# add daylength
library(geosphere) 

moreclim_df <- subset(clim_df, var=="tmean")
names(moreclim_df)[names(moreclim_df)=="value"] <- "tmean"

moreclim_df$daylength <- NA
moreclim_df$doy <- as.numeric(format(as.Date(moreclim_df$date, format="%Y-%m-%d") , "%j"))

for(i in 1:nrow(moreclim_df)){
  moreclim_df$daylength[i] <- daylength(50, moreclim_df$doy[i])
}

moreclim_df$gddtemp <- ifelse(moreclim_df$tmean>0, moreclim_df$tmean, 0)
moreclim_df$year <- as.numeric(format(as.Date(moreclim_df$date, format="%Y-%m-%d") , "%Y"))
# simulate leafout 
fstar <- 150
lodf <- data.frame(year=unique(moreclim_df$year), 
  loday=rep(NA, length(unique(moreclim_df$year))),
  gdd=rep(NA, length(unique(moreclim_df$year))), 
  daylength=rep(NA, length(unique(moreclim_df$year))))

for(yearhere in unique(moreclim_df$year)) {
  thisyear <- moreclim_df[which(moreclim_df$year==yearhere),]
  leafoutdate <- min(which(cumsum(thisyear[["gddtemp"]]) > fstar))
  lodf$loday[which(lodf$year==yearhere)] <- leafoutdate
  lodf$gdd[which(lodf$year==yearhere)] <- cumsum(thisyear[["gddtemp"]])[leafoutdate]
  lodf$daylength[which(lodf$year==yearhere)] <- thisyear$daylength[leafoutdate]
}

## Put everything together to write out ...
newclim_wide <- reshape(newclim_df, idvar=c("lat", "lon", "year"), timevar="var", direction = "wide")
names(newclim_wide) <- c( "lat", "lon", "year", 
  "tminwinter",
  "tmeanspring",
  "gddspring",
  "precspring",
  "totalprec",
  "chillwinter")
simdat <- merge(newclim_wide, lodf, by="year")


write.csv(simdat, file.path(wd, "output/simuldateddat.csv"), row.names=FALSE)
