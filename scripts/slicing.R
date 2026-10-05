#### Native climatologies to monsoon slices ####

slice_climatologies = function(data_dir,products,locations,
                               n_lat=361,n_long=720,range_km=750,
                               lambda_smooth=1250,days_wrap=60){
  wendland_kernel = function(d,r){
    (d<r)*(1/3)*((1-d/r)^6)*(35*(d/r)^2 + 18*(d/r) + 3)
  }

  distChordal = function(x,y,r=6371){
    r*acos(round(sin(x[,1])*sin(y[,1])+cos(x[,1])*cos(y[,1])*cos(x[,2]-y[,2]),10))
  }

  smooth_ts = function(y,lambda=1250,days.wrap=60){
    #wrap the year to connect December and January when smoothing
    len = length(y)
    ind.wrap = c((len-days.wrap+1):len,1:len,1:days.wrap)
    ind.unwrap = days.wrap + 1:len
    glmgen::trendfilter(y[ind.wrap],k=2,lambda=lambda,thinning=FALSE)$beta[ind.unwrap]
  }

  #### Nearest-neighbour regridding ####

  xyz = function(lon,lat){
    cbind(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat))
  }

  nearest_grid = function(lon,lat,coords,rectangular=TRUE){
    lon = as.numeric(lon)*(pi/180)
    lat = as.numeric(lat)*(pi/180)
    plon = coords$dlong*(pi/180)
    plat = coords$dlat*(pi/180)
    if(!rectangular){
      neighbours = FNN::get.knnx(xyz(lon,lat),xyz(plon,plat),k=8)
      ids = neighbours$nn.index
      ids[neighbours$nn.dist>neighbours$nn.dist[,1]+1e-14] = .Machine$integer.max
      return(as.integer(apply(ids,1,min)))
    }
    nx = length(lon)
    ny = length(lat)
    stopifnot(all(diff(lon)>0),all(diff(lat)>0) || all(diff(lat)<0))
    plon[plon<lon[1]] = plon[plon<lon[1]]+2*pi
    plon[plon>lon[1]+2*pi] = plon[plon>lon[1]+2*pi]-2*pi
    ii = pmax(1L,findInterval(plon,c(lon,lon[1]+2*pi),left.open=TRUE))
    ii[ii==nx] = 0L
    direction = sign(lat[ny]-lat[1])
    jj = pmax(1L,findInterval(plat*direction,lat*direction,left.open=TRUE))
    inside = plat>=min(lat) & plat<=max(lat)
    src = rep(.Machine$integer.max,length(plat))
    distance = rep(Inf,length(plat))
    query = xyz(plon,plat)
    clon = cos(lon %% (2*pi))
    slon = sin(lon %% (2*pi))
    clat = cos(lat)
    slat = sin(lat)

    #match CDO's regular-grid search
    for(j in -1:1){
      for(i in -1:1){
        y = jj+j+1L
        x = (ii+i) %% nx+1L
        valid = which(inside & y>=1L & y<=ny)
        a = (y[valid]-1L)*nx+x[valid]
        d = (query[valid,1]-clat[y[valid]]*clon[x[valid]])^2 +
            (query[valid,2]-clat[y[valid]]*slon[x[valid]])^2 +
            (query[valid,3]-slat[y[valid]])^2
        d = sqrt(readBin(writeBin(d,raw(),size=4),double(),n=length(d),size=4))
        take = d+1e-12<distance[valid] |
               (a<src[valid] & abs(d-distance[valid])<1e-12)
        src[valid[take]] = a[take]
        distance[valid[take]] = d[take]
      }
    }

    #handling for locations outside the native latitude extent
    for(k in which(!inside)){
      edge = if(abs(plat[k]-lat[1])<abs(plat[k]-lat[ny])) 1:2 else (ny-1):ny
      ids = as.vector(outer(1:nx,(edge-1)*nx,'+'))
      d = unlist(lapply(edge,function(y){
        suppressWarnings(acos(cos(plat[k])*cos(lat[y])*
          (cos(plon[k])*cos(lon)+sin(plon[k])*sin(lon)) + sin(plat[k])*sin(lat[y])))
      }))
      src[k] = ids[which.min(d)]
    }
    stopifnot(all(src>=1L & src<=nx*ny))
    as.integer(src)
  }

  stopifnot(n_lat>=3,n_long>=4,n_lat==as.integer(n_lat),n_long==as.integer(n_long),
            days_wrap>=1,days_wrap<365,range_km>0)
  coords = expand.grid(dlat=seq(90,-90,length.out=n_lat),
                        dlong=seq(0,360-360/n_long,length.out=n_long))
  grid = round(coords*(pi/180),10)
  lat_w = cos(coords$dlat*pi/180)
  slice_lats = unique(locations$dlat)
  slice_longs = unique(locations$dlong)
  stopifnot(identical(locations,expand.grid(dlat=slice_lats,dlong=slice_longs)))
  n = nrow(coords)
  n_slice = nrow(locations)

  #recycle computation by rotating latitude kernels along the integration grid
  phases = round(slice_longs %% (360/n_long),10)
  kernels = lapply(unique(phases),function(phase){
    lapply(slice_lats,function(lat){
      centre = matrix(round(c(lat,phase)*(pi/180),10),nrow=1)
      k = lat_w*wendland_kernel(distChordal(centre,grid),range_km)
      stopifnot(sum(k)>0)
      k = k/sum(k)
      keep = which(k!=0)
      list(index=keep,value=k[keep])
    })
  })
  rows = values = vector('list',n_slice)
  for(i in seq_len(n_slice)){
    lat = match(locations$dlat[i],slice_lats)
    lon = match(locations$dlong[i],slice_longs)
    k = kernels[[match(phases[lon],unique(phases))]][[lat]]
    shift = round((slice_longs[lon]-phases[lon])/(360/n_long))*n_lat
    rows[[i]] = (k$index-1+shift) %% n+1
    values[[i]] = k$value
  }
  convo_matrix = Matrix::sparseMatrix(i=unlist(rows),
                   j=rep(seq_len(n_slice),lengths(rows)),x=unlist(values),
                   dims=c(n,n_slice))
  fine_index = which(Matrix::rowSums(convo_matrix)!=0)
  convo_matrix = convo_matrix[fine_index,,drop=FALSE]
  fine_coords = coords[fine_index,]
  kernel_bytes = as.numeric(object.size(convo_matrix))
  grid_bytes = sum(vapply(list(coords,grid,fine_coords),function(x) as.numeric(object.size(x)),numeric(1)))
  rm(kernels,rows,values,k,coords,grid,lat_w)

  slices = list()
  memory = data.frame()
  for(model in products){
    message('Slicing: ',model)
    nc = ncdf4::nc_open(file.path(data_dir,paste0(model,'.nc')))
    stopifnot(identical(as.integer(ncdf4::ncvar_get(nc,'day_of_year')),1:365),
              ncdf4::ncatt_get(nc,'pr','units')$value=='mm/day')
    lat = as.numeric(ncdf4::ncvar_get(nc,'latitude'))
    lon = as.numeric(ncdf4::ncvar_get(nc,'longitude'))
    rectangular = !'cell' %in% names(nc$dim)
    src = nearest_grid(lon,lat,fine_coords,rectangular)
    native_index = sort(unique(src))
    native = matrix(0,365,length(native_index))

    #read one native daily field at a time and keep only the required cells
    for(day in 1:365){
      if(rectangular){
        field = ncdf4::ncvar_get(nc,'pr',start=c(1,1,day),count=c(length(lon),length(lat),1))
      }else{
        field = ncdf4::ncvar_get(nc,'pr',start=c(1,day),count=c(length(lat),1))
      }
      native[day,] = field[native_index]
    }
    ncdf4::nc_close(nc)
    stopifnot(all(is.finite(native)))
    smooth = apply(native,2,smooth_ts,lambda=lambda_smooth,days.wrap=days_wrap)
    up = smooth[,match(src,native_index),drop=FALSE]
    sliced = as.matrix(up %*% convo_matrix)
    dim(sliced) = c(365,length(slice_lats),length(slice_longs))
    slices[[model]] = sliced
    memory = rbind(memory,data.frame(dataset=model,native_cells=length(native_index),
               native_MB=as.numeric(object.size(native))/1e6,
               smooth_MB=as.numeric(object.size(smooth))/1e6,
               rescaled_MB=as.numeric(object.size(up))/1e6,
               slices_MB=as.numeric(object.size(sliced))/1e6))
    rm(native,smooth,up,sliced,field)
    gc()
  }
  message('Integration-grid coordinates: ',round(grid_bytes/1e6,2),' MB; kernel matrix: ',
          round(kernel_bytes/1e6,2),' MB; supported grid cells: ',length(fine_index))
  print(memory,row.names=FALSE)
  list(slices=slices,coords=locations,memory=memory)
}
