#### Setup ####

library(tidyverse)
library(Matrix)
library(ncdf4)
library(glmgen)
library(FNN)
library(fdasrvf)
library(viridis)
library(ggrepel)
library(patchwork)
library(sf)
library(maps)
library(rworldmap)
library(rworldxtra)

mods = c('NorESM2-MM','AWI-CM-1-1-MR','CESM2-WACCM',
         'ICON-ESM-LR','IPSL-CM5A2-INCA','NESM3')
data_dir = 'data/'
figure_dir = 'output/'
n_lat = 361
n_long = 720
range_km = 750
lambda_smooth = 1250
days_wrap = 60
set.seed(1)
dir.create(figure_dir,showWarnings=FALSE)

#### Slicing ####

source('scripts/slicing.R')
locations = expand.grid(dlat=seq(29.5,15.5,-1),dlong=68:88)
climatologies = slice_climatologies(data_dir,c('GPCP',mods),locations,
                                    n_lat,n_long,range_km,lambda_smooth,days_wrap)

#### Monsoon phase analysis and figures ####

source('scripts/monsoon.R')
source('scripts/plotting.R')
message('Retained slices: ',round(as.numeric(object.size(climatologies$slices))/1e6,2),' MB')
message('Timing-bias arrays: ',round(as.numeric(object.size(mods_tb))/1e6,3),' MB')
