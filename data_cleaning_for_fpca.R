library(ncdf4)


# Setup parameters
data_dir = "~/OneDrive/STA237A/"
data_file = "lai_Lmon_CCSM4_rcp45_r2i1p1_200601-210012.nc"
data_path = paste(data_dir, data_file, sep="")
attr_name = "lai" # leaf area index

# Load the.nc
ncin = nc_open(data_path)

# Extract lat, lon, time
lon = ncvar_get(ncin, "lon")
lat = ncvar_get(ncin, "lat")
loc.af = c(
  min(which(lat > 3)),
  min(which(lon > 25))
)
t = ncvar_get(ncin, "time")
tunits = ncatt_get(ncin, "time", "units")

# Extract attribute and data
attr = ncatt_get(ncin, attr_name)
attr_data = ncvar_get(ncin, attr_name)
d = attr_data # shorthand for the array of space X time
# NOTE: d is indexed d[lon, lat, time], but the natural location indexing is
# loc = [lat, lon] so to index d by a location you must do
# d[loc[2], loc[1], time] rather than d[loc[1], loc[2], time]

## Plot the global LAI at time 0
# extract variables to plot
x = lon # 288 longitudes
y = lat # 192 latitudes
z = d[,,1] # all LAI values across the planet at time 1

# remove zero vectors and na vectors in time-series
# list[[lon idx]][[lat idx]] calles leaf index over the time ()
listall = lapply(1:length(lon), function(i) lapply(1:length(lat), function(j) d[i,j,]))

# first want to remove the NA values
tmp1 = sapply(1:length(lon), function(i) sapply(1:length(lat), function(j) sum(is.na(listall[[i]][[j]]))))
# 0 indicating nonNA -> choose the location having 0 which indicates non NA value vector
# tmp1 is 192 by 288 matrix
non0idx = which(tmp1 == 0) #it returns vector, need to convert matrix coordinates (1 -> [1,1], 2)
# Generally (1) number %/% total row number = column, (2) number %% total row number = row 
# if (2) is zero, row is the ending row, column is (1), otherwise, row is (2) and column is (1) + 1
vec2mat = function(vec, numR)
  #Vector: index pointing a certain matrix (numR by numC)
  #numR : number of rows
{
  tmprow = vec %% numR
  tmpcol = vec %/% numR
  tmppair = cbind(tmprow, tmpcol)
  tmp2col = ifelse(tmppair[,1] == 0, tmppair[,2], tmppair[,2] + 1)
  tmp2row = ifelse(tmppair[,1] == 0, numR, tmppair[,1])
  pairs = cbind('Row' = tmp2row, 'Col' = tmp2col)
  return(pairs)
}

idxmat1 = vec2mat(non0idx, 192)
chosenlist = sapply(1:dim(idxmat1)[1], function(i) listall[[idxmat1[i,2]]][[idxmat1[i,1]]])
plot(chosenlist[,1596])
data1 = cbind('lon' = lon[idxmat1[,2]], 'lat' = lat[idxmat1[,1]], t(chosenlist))

#now take off only zero observations
zeroidx = sapply(1:dim(data1)[1], function(i) sum(data1[i,-(1:2)] != 0))
data2 = data1[which(zeroidx != 0),]
names(data2) = c('lon', 'lat', t)
summary(data2[,1:2])
write.csv(data2, file = "leafidx.csv")
save(data2, t, file = "leafidx.RData")
