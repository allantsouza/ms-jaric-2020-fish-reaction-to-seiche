source(file = 'data_loading.R', local = T)

#Plots
ggplot(thermocline.sub , aes(x = timestamp_utc, y = th_depth, col = pos_name))+
  geom_line()+
  geom_smooth()



ggplot(detections.m.sub, aes(x = therm_minus_balance_cat, y = fish_minus_balance, group = interaction(therm_minus_balance_cat, day_night), col = day_night))+
  geom_boxplot()+
  facet_wrap(~up_tag_sn, scales = "free")



#After Souza's modifications####
setwd('~/Documents/Czechia/Ivan/plots/')

#Plot of fish depth and thermocline through time
ggplot(detections_std, aes(x=timestamp, y=fishdepth))+
  geom_point(aes(col=day_night))+
  scale_y_reverse()+
  geom_line(aes(x=timestamp, y=balanced_therm_depth))+
  facet_wrap(~fishid)+
  labs(x='Time', y='Fish depth (m)', title='Tench (Chabarovice reservoir)')+
  ggsave(filename = 'fish_depth_time.tiff', device = 'tiff', width = 12, height = 8, dpi = 600, units = 'cm', compression = 'lzw')

ggplot(detections_std, aes(y=fish_minus_therm_depth, x=seiche_state, col=day_night))+
  geom_boxplot()


ggplot(detections_std, aes(x=thermocline_depth_column, y=fishdepth, col=day_night))+
  geom_smooth()#+
  #geom_point(alpha=0.2)

ggplot(detections.m.sub, aes(x=dd_timestamp_utc, y=dd_depth, col=thermpart))+
  geom_line()+
  facet_wrap(~thermpart)
