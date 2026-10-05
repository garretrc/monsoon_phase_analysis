#### Plotting functions ####

cdate = function(t){
  t_date = as.Date(t, origin = '2024-12-31')
  t_date
}

phase_report = function(days){
  c(paste(round(abs(days[1]),1),'days',ifelse(days[1]<0,'early','late')),
    paste(round(abs(days[2]),1),'days',ifelse(days[2]<0,'early','late')))
}

efda = function(c1=gpcp_ts,n1='GPCP',
                c2=gpcp_ts,n2='GPCP',
                c3=gpcp_ts,n3='GPCP',
                mt_a=0.5,
                mt_r=0.5,
                d_lims = c(as.Date('2025-05-15'),as.Date('2025-11-01')),
                offset=0){
  t_arrive = min(which(c1>=mt_a*max(c1)))
  t_peak = which.max(c1)
  t_retreat = max(which(c1>=mt_r*max(c1)))
  tc = 'black'
  time=1:365
  res = efda_warp(c1,c2,time, 
                  lambda = 0, pen = "roughness", method = 'DP',nbhd_dim=7,
                  t_arrive,t_peak,t_retreat)
  j1 = 364*res[[1]][5]+1
  a1 = 364*res[[1]][6]+1
  o1 = 364*res[[1]][7]+1
  labels = phase_report(c(j1-t_arrive,o1-t_retreat))
  res2 = efda_warp(c1,c3,time, 
                  lambda = 0, pen = "roughness", method = 'DP',nbhd_dim=7,
                  t_arrive,t_peak,t_retreat)
  j2 = 364*res2[[1]][5]+1
  a2 = 364*res2[[1]][6]+1
  o2 = 364*res2[[1]][7]+1
  labels2 = phase_report(c(j2-t_arrive,o2-t_retreat))
  p1 = ggplot()+
          geom_abline(aes(intercept=0,slope=1),color=tc,linetype='dashed',linewidth=lw*0.5)+
          geom_line(aes(x=cdate(time),y=cdate(res[[3]]*364+1)),linewidth=1.25,color="#9E7AD6")+
          geom_vline(aes(xintercept = cdate(t_arrive)),color=tc,linetype='dotted',linewidth=lw*0.5)+
          geom_vline(aes(xintercept = cdate(t_retreat)),color=tc,linetype='dotted',linewidth=lw*0.5)+
          geom_segment(aes(yend = cdate(t_arrive), y = cdate(res[[3]]*364+1)[t_arrive],  
                           xend = cdate(t_arrive), x = cdate(t_arrive)),
                       arrow=arrow(type='closed',length=unit(0.15,"cm")),
                       color='#FD760F',linewidth=1)+
          geom_segment(aes(yend = cdate(t_retreat), y = cdate(res[[3]]*364+1)[t_retreat],  
                           xend = cdate(t_retreat), x = cdate(t_retreat)),
                       arrow=arrow(type='closed',length=unit(0.15,"cm")),
                       color='#FD760F',linewidth=1)+
          theme_bw()+
          theme(aspect.ratio = 1)+
          labs(x=paste('Date in',n1),y=paste('Date in',n2))+
          geom_label_repel(aes(x=cdate(t_arrive),y=cdate(t_arrive+(j1-t_arrive)/2),label=labels[1]),
                          nudge_x = 55,nudge_y=-7)+
          geom_label_repel(aes(x=cdate(t_retreat),y=cdate(t_retreat+(o1-t_retreat)/2),label=labels[2]),
                          nudge_x = -45, nudge_y=7)+
          coord_cartesian(xlim=d_lims,ylim=d_lims)+
          scale_x_date(date_breaks = "1 month", date_labels = "%b%e")+
          scale_y_date(date_breaks = "1 month", date_labels = "%b%e")
  c2w = res[[2]]
  c3w = res2[[2]]
  p2 = ggplot(mapping=aes(x = cdate(time)))+
          geom_vline(aes(xintercept=cdate(t_arrive)), linetype='dotted', color = tc,linewidth=lw*0.5)+
          geom_vline(aes(xintercept=cdate(t_retreat)), linetype='dotted', color = tc,linewidth=lw*0.5)+
          geom_line(aes(y=c1,color=n1),linewidth=lw)+
          geom_line(aes(y=c2,color=n2),linewidth=lw)+
          geom_line(aes(y=c2w,color=n2),linewidth=lw*0.75,linetype='dashed')+
          geom_segment(aes(y=c2w[t_arrive],yend=c2w[t_arrive],  x=cdate(j1),xend=cdate(t_arrive)),
                       arrow=arrow(type='closed',length=unit(0.15,"cm")),color='#FD760F',linewidth=1)+
          geom_segment(aes(y=c2w[t_retreat],yend=c2w[t_retreat],x=cdate(o1),xend=cdate(t_retreat)),
                       arrow=arrow(type='closed',length=unit(0.15,"cm")),color='#FD760F',linewidth=1)+
          scale_color_manual(values=c("black","#9E7AD6","#009E73") |> `names<-`(c(n1,n2,n3)))+
          labs(y='precipitation (mm)',x='Date',color=NULL)+
          theme_bw()+
          theme(aspect.ratio=1/1.75,
                legend.position = c(0.84,0.9),
                legend.box.background = element_rect(color="black", size=1),
                legend.margin=margin(0,1,0,0))+
          coord_cartesian(xlim=d_lims)+
          scale_x_date(date_breaks = "1 month", 
                       date_labels = "%b%e")
  
  
  list(p1,p2)
}

#### Monsoon maps ####

plot_ismr_date = function(field, coords, value='Date',zlims=FALSE,flip=FALSE,guide='colorbar'){
  cols = c("x", "y")
  if(flip==TRUE){
    cols = c('y','x')
  }
  m_df = data.frame(
    x = coords$dlong,
    y = coords$dlat
  )
  m_df[,value] = as.numeric(field)
  if(length(zlims)!=2){
    min_val = min(field)
    max_val = max(field)
  }else{
    min_val = zlims[1]
    max_val = zlims[2]
  }
  if(value=='Temperature'){
    units = '°C'
  }else{
    units=value
  }
  ggplot(data=m_df %>% filter(value!=0))+
    geom_tile(aes(x=x,y=y,fill=.data[[value]]),color=NA)+
    geom_tile(aes(x=x,y=y,color=.data[[value]]),fill=NA,linewidth=0.01)+
    geom_path(data=mp_ismr,mapping=aes(x=long,y=lat,group=group))+
    coord_quickmap(xlim = range(m_df$x),ylim = range(m_df$y))+
    scale_fill_viridis(labels = function(x){format(as.Date(x,origin='2024-12-31'), "%b %d")},
                                  breaks = round(seq(min_val,max_val,length.out=5)),
                                  option=date_palette, limits = c(min_val,max_val),
                                  guide=guide, na.value=NA)+
    scale_color_viridis(labels = function(x){format(as.Date(x,origin='2024-12-31'), "%b %d")},
                                   breaks = round(seq(min_val,max_val,length.out=5)),
                                   option=date_palette, limits = c(min_val,max_val),
                                   guide='none', na.value=NA)+
    theme_minimal()+
    labs(fill=units,x=NULL,y=NULL)
}

plot_ismr_tb = function(field, mask, coords, value='Date',zlims=FALSE,flip=FALSE,guide='colorbar'){
  cols = c("x", "y")
  if(flip==TRUE){
    cols = c('y','x')
  }
  m_df = data.frame(
    x = coords$dlong,
    y = coords$dlat
  )
  m_df[,value] = as.numeric(field)
  m_df[,"mask"] = as.numeric(mask)
  if(length(zlims)!=2){
    min_val = min(m_df[,value][m_df$mask!=0])
    max_val = max(m_df[,value][m_df$mask!=0])
  }else{
    min_val = zlims[1]
    max_val = zlims[2]
  }
  m_df[,value][m_df[,value]<min_val] = min_val
  m_df[,value][m_df[,value]>max_val] = max_val
  if(value=='Temperature'){
    units = '°C'
  }else{
    units=value
  }
  ggplot(data=m_df %>% filter(mask!=0))+
    geom_tile(aes(x=x,y=y,fill=.data[[value]]),color=NA)+
    geom_tile(aes(x=x,y=y,color=.data[[value]]),fill=NA,linewidth=0.01)+
    geom_path(data=mp_ismr,mapping=aes(x=long,y=lat,group=group))+
    coord_quickmap(xlim = range(m_df$x),ylim = range(m_df$y))+
    scale_fill_gradientn(colours=bias_colours, values=bias_colour_values, 
                                                breaks = seq(-42,42,14),
                                                labels = paste0(abs(seq(-6,6,2)),c('+','','','','','','+'),' weeks',c(rep(' early',3),'',rep(' late',3))),
                                                limits = c(min_val,max_val),
                                                guide=guide, na.value=NA)+
    scale_color_gradientn(colours=bias_colours, values=bias_colour_values,
                                                 breaks = seq(-42,42,14),
                                                 labels = paste0(abs(seq(-6,6,2)),c('+','','','','','','+'),' weeks',c(rep(' early',3),'',rep(' late',3))),
                                                 limits = c(min_val,max_val),
                                                 guide='none', na.value=NA)+
    theme_minimal()+
    labs(fill='Timing bias',color='Timing bias',x=NULL,y=NULL)
}

#### Monsoon region ####

ismr_coords = climatologies$coords
india_border = map_data('world') %>% filter(region=='India')
india_poly = india_border %>% filter(group==839) %>%
  st_as_sf(coords=c('long','lat'),crs=4326) %>%
  summarise(geometry=st_combine(geometry)) %>%
  st_cast('POLYGON')
ismr_coords_sf = st_as_sf(ismr_coords,coords=c('dlong','dlat'),crs=4326)
india_mask = 1:nrow(ismr_coords) %in% st_intersects(india_poly,ismr_coords_sf)[[1]]
mask_ismr = matrix(india_mask,length(unique(ismr_coords$dlat)),length(unique(ismr_coords$dlong)))

data('countriesHigh')
mp1 = fortify(countriesHigh)
mp1$group = as.numeric(mp1$group)
mp2 = mp1
mp2$long = mp2$long+360
mp2$group = mp2$group+max(mp2$group)+1
mp_ismr = rbind(mp1,mp2)

#cividis colors with hyperbolic scale
date_palette = 'cividis'
bias_colours = viridis(257,option='cividis')
palette_positions = seq(-1,1,length.out=257)
bias_colour_values = (1+sinh(palette_positions*asinh(3))/3)/2

#### Figure 5: GPCP onset and retreat ####

range_ismr = function(x){range(as.numeric(x*mask_ismr)[as.numeric(x*mask_ismr)>0])}
po = plot_ismr_date(onset*mask_ismr,ismr_coords,zlims=range_ismr(onset))+
  labs(title='Monsoon onset',y='Latitude',x='Longitude')+
  scale_y_continuous(breaks=c(15,20,25,30),labels=paste0(c(15,20,25,30),'°N'))+
  scale_x_continuous(breaks=c(70,75,80,85),labels=paste0(c(70,75,80,85),'°E'))
pr = plot_ismr_date(retreat*mask_ismr,ismr_coords,zlims=range_ismr(retreat))+
  labs(title='Monsoon retreat',x='Longitude')+
  scale_y_continuous(breaks=c(15,20,25,30),labels=paste0(c(15,20,25,30),'°N'))+
  scale_x_continuous(breaks=c(70,75,80,85),labels=paste0(c(70,75,80,85),'°E'))
figure5 = po | pr
ggsave(paste0(figure_dir,'onset_retreat.pdf'),figure5,width=8.5,height=3,device=cairo_pdf)
ggsave(paste0(figure_dir,'onset_retreat.png'),figure5,width=8.5,height=3,dpi=160,bg='white')

#### Figure 6: Alignment at 22.5 N, 78 E ####

high_mod = 'IPSL-CM5A2-INCA'
low_mod = 'ICON-ESM-LR'
high_ismr = climatologies$slices[[high_mod]]
low_ismr = climatologies$slices[[low_mod]]
slice_y = which(unique(ismr_coords$dlat)==22.5)
slice_x = which(unique(ismr_coords$dlong)==78)
gpcp_ts = gpcp_ismr[,slice_y,slice_x]
high_ts = high_ismr[,slice_y,slice_x]
low_ts = low_ismr[,slice_y,slice_x]
lw = 1.5
high_plots = efda(gpcp_ts,'GPCP',high_ts,high_mod)
low_plots = efda(gpcp_ts,'GPCP',low_ts,low_mod)
figure6 = ((high_plots[[2]]+labs(x=NULL)|high_plots[[1]]+labs(x=NULL))+
             plot_layout(widths=c(1.75,1),heights=c(1,1)))/
          ((low_plots[[2]]|low_plots[[1]])+plot_layout(widths=c(1.75,1),heights=c(1,1)))
ggsave(paste0(figure_dir,'model_slices.pdf'),figure6,width=626/72,height=468/72,device=cairo_pdf)
ggsave(paste0(figure_dir,'model_slices.png'),figure6,width=626/72,height=468/72,dpi=160,bg='white')

#### Figure 7: Onset and retreat timing bias ####

bias_grid = function(component){
  panels = lapply(1:length(mods),function(i){
    ret = plot_ismr_tb(mods_tb[[i]][component,,],mask_ismr,ismr_coords,zlims=c(-42,42))+
      labs(subtitle=mods[i])+
      scale_y_continuous(breaks=c(15,20,25,30),labels=paste0(c(15,20,25,30),'°N'))+
      scale_x_continuous(breaks=c(70,75,80,85),labels=paste0(c(70,75,80,85),'°E'))
    if(i %in% c(2,3,5,6)){ret = ret+theme(axis.text.y=element_blank())}
    if(i %in% c(1,2,3)){ret = ret+theme(axis.text.x=element_blank())}
    ret
  })
  wrap_plots(panels,ncol=3,nrow=2)+plot_layout(guides='collect') &
    theme(legend.position='right') &
    guides(fill=guide_colorbar(barheight=unit(0.8,'npc')))
}
figure7_onset = bias_grid(1)
figure7_retreat = bias_grid(3)
ggsave(paste0(figure_dir,'onset_bias.pdf'),figure7_onset,width=8,height=4,device=cairo_pdf)
ggsave(paste0(figure_dir,'onset_bias.png'),figure7_onset,width=8,height=4,dpi=160,bg='white')
ggsave(paste0(figure_dir,'retreat_bias.pdf'),figure7_retreat,width=579/72,height=4,device=cairo_pdf)
ggsave(paste0(figure_dir,'retreat_bias.png'),figure7_retreat,width=579/72,height=4,dpi=160,bg='white')

