## Fish vertical velocities by time of day
We use the data from real detections to calculate the difference in depth between two successive detections (every minute). The absolute values of the velocity are therefore given.

Calculate fish speeds
``` r
speeds<-detections%>%group_by(fishid)%>%mutate(diff=det_depth-lag(det_depth,default=first(det_depth)))
speeds <- mutate(speeds, speed = diff/60)
```

Drop roach level
``` r
speeds<- droplevels(speeds[!speeds$species == 'roach',])
speeds$Species<-as.factor(speeds$species)
```
Group speed values by fishid, species and diel period (otherwise there are too many points)
``` r
speeds_mean <- speeds %>%
                      group_by(fishid,species,diel_period) %>%
                      summarise_at(vars(speed), mean)
```
Convert speeds to absolute values (negative speeds are reversed changes in depth - going up)
``` r
speeds_mean_abs <- speeds_mean %>%
                group_by(fishid,species,diel_period) %>%
                summarise_at(vars(speed), abs)
```
Violin plot of vertical speeds
``` r
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
               "#E956B4","#B4E956","#22485D","#000000","#DDF0FA","#485D22","#742B5A","#E9E356","#EF223E",
               "#22EF58","#169D39","#9D1579","#9D1579","#149D38","#6E149D","#9D1479","#B985D5","#439D14",
               "#1B149D","#9D5114","#B65440","#40A1B6","#77CFC7","#B69B40","#979C13","#15DB1E","#DB155B",
               "#9C1343","#D6EB4C","#F1E284","#84D4F1","#F1A184","#604CEB","#AFA5F5","#E2E6FB","#EEDF8F")
```
``` r
scientific <- function(x){
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}
ggplot(speeds_mean_abs, aes(x=species, y=speed, fill=species)) +
  geom_violin(trim=FALSE, fill = "white")+
  geom_boxplot(width=0.1) + scale_y_continuous(label=scientific) +
    scale_fill_manual(values = c(cbPalette[4], cbPalette[9], cbPalette[15], cbPalette[27]),
                       name="Species", labels=c("Northern pike", "Rudd", "Tench","Wels catfish")) +
    facet_wrap(~diel_period) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          plot.title = element_text( size = 18, family="Helvetica"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text( size = 17, family="Helvetica"),
          axis.text.y=element_text(size = 16, family="Helvetica"),
          axis.title.y = element_text( size = 18,margin = unit(c(0, 9, 0, 0), "mm"), family="Helvetica"),
          legend.title = element_blank(),
          legend.text = element_text(color = "black", size = 17, family="Helvetica"),
          legend.key.width=unit(1, "cm"),
          legend.key.height=unit(1, "cm"),
          strip.text = element_text(size = 18, face = "bold", family="Helvetica")) +
          theme(strip.background =element_rect(fill="white")) +
          theme(strip.text = element_text(colour = 'black')) + ylab("Vertical speed (m/s)")
```

## Seiche speeds

Calculate seiche speeds
``` r
speeds_seiche = detections %>%
                             filter(is_valid_seiche == TRUE) %>%
                             mutate(startindex = ifelse(
                               test = row_number(dets_ts) == 1, yes = T, no = F)) %>%
                             mutate(fishid = as_factor(fishid))
```
``` r
speeds_seiche<-speeds_seiche%>%group_by(fishid)%>%mutate(diff_seiche=amplitude-lag(amplitude,default=first(amplitude)))
speeds_seiche <- mutate(speeds_seiche, speed = diff_seiche/60)
```
Drop roach level
``` r
speeds_seiche<- droplevels(speeds_seiche[!speeds_seiche$species == 'roach',])
speeds_seiche$species<-as.factor(speeds_seiche$species)
```
Separate up and down seiche
``` r
speeds_seiche$up_down_seiche <- sapply(speeds_seiche$amplitude, function(amplitude){
                               if(amplitude < 0){"down"
                               } else if(amplitude > 0){"up"
                               } else{"up"}})
```
``` r
speeds_seiche<-as.data.frame(speeds_seiche)
speeds_seiche$up_down_seiche<-as.factor(speeds_seiche$up_down_seiche)
speeds_seiche$speed <- abs(speeds_seiche$speed)
```

Upward and downward seiche velocities over time (gam)
``` r
ggplot(speeds_seiche, aes(x=date4, y=speed, color=up_down_seiche, group=up_down_seiche)) +
  stat_smooth(size=1.5, method = "gam", level = 0.95, fullrange = TRUE, se = TRUE) +
    scale_x_date(breaks = "1 months", minor_breaks = "1 month", date_labels = "%b %y", expand = c(0,0), labels = scales::number_format(accuracy = 0.2)) +
    scale_y_continuous(label=scientific_10) +
    scale_color_manual(values = c("dodgerblue4","firebrick")) +
    theme_bw() +  ggtitle("") + xlab("Date") + ylab("Seiche speed (m/s)") + labs(color = "seiche") +
    guides(fill=guide_legend(override.aes = list(fill="white",size=1.2))) +
              theme(panel.grid.major = element_blank(),
              plot.title = element_text( size = 18, family="Helvetica"),
              axis.title.x = element_text( size = 18,margin = unit(c(5, 0, 0, 0), "mm"), family="Helvetica"),
              axis.text.x = element_text(angle=0, hjust=1, size = 17 ),
              axis.text = element_text( size = 17, family="Helvetica"),
              axis.text.y=element_text(size = 17, family="Helvetica"),
              axis.title.y = element_text( size = 18,margin = unit(c(0, 9, 0, 0), "mm"), family="Helvetica"),
              legend.position ='right',
              legend.justification='center',
              legend.direction='vertical',
              legend.title = element_text(color = "black", size = 18, face = "plain", family="Helvetica"),
              legend.text = element_text(color = "black", size = 17, family="Helvetica"),
              legend.key.width=unit(1, "cm"),
              legend.key.height=unit(1, "cm"),
              strip.text = element_text(size = 18, face = "bold", family="Helvetica")) +
              theme(strip.background =element_rect(fill="white")) +
              theme(strip.text = element_text(colour = 'black'))
```
Upward and downward seiche velocities as a function of the amplitude of the thermocline
``` r
ggplot(data = speeds_seiche) +
  geom_point(mapping = aes(x = amplitude, y = speed, colour = up_down_seiche))  +
   theme_bw() + ggtitle("") +
   xlab("Amplitude") +
   ylab("Seiche speed (m/s)")+
   guides(fill=guide_legend(override.aes = list(fill="white",size=1.2))) +
   theme(panel.grid.major = element_blank(),
     plot.title = element_text( size = 18, family="Helvetica"),
     axis.title.x = element_text( size = 18,margin = unit(c(5, 0, 0, 0), "mm"), family="Helvetica"),
     axis.text.x = element_text(angle=0, hjust=1, size = 17 ),
     axis.text = element_text( size = 17, family="Helvetica"),
     axis.text.y=element_text(size = 17, family="Helvetica"),
     axis.title.y = element_text( size = 18,margin = unit(c(0, 9, 0, 0), "mm"), family="Helvetica"),
     legend.position ='none',
     legend.justification='center',
     legend.direction='vertical',
     legend.title = element_text(color = "black", size = 18, face = "plain", family="Helvetica"),
     legend.text = element_text(color = "black", size = 17, family="Helvetica"),
     legend.key.width=unit(1, "cm"),
     legend.key.height=unit(1, "cm"),
     strip.text = element_text(size = 18, face = "bold", family="Helvetica")) +
     theme(strip.background =element_rect(fill="white")) +
     theme(strip.text = element_text(colour = 'black'))
```
