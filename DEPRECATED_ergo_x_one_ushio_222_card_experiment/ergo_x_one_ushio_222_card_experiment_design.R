# ergo x one experiment design script
# CRUDELY measure light field with Ushio 222 cards https://www.ushio.com/products/infection-prevention/dose222-far-uv-c-indicator-cards/##

# library(tidyverse)
# needs function uv_can_model from ergo_x_one_ushio_222_card_experiment, so run that script first!

# data collected on my kitchen island with janky setup. 
# Illustrative but could easily be off by a factor of 2!
d<-read.csv('ergo_x_one_ushio_222_card_experiment_data.R.csv') %>% as_tibble()

# discretize angles for plotting later
d$theta_factor = factor(d$theta,levels=sort(unique(d$theta)))

# NOTE! My tabletop geometry was clearly off.  When I put the device flush on
# the card, I get roughly 1500 mu-W/cm^2 as expected. But the intercept for finite r
# empirically is around 1000 mu-W/cm^2.
# two likely sources of error
#  1) my cards were off-vertical-axis of the device. The amoung required is ~30 degrees
#     up close, which feels like more than it could've possibly been, but this surely contributes.
#     1 cm vertically made similar relative differences for Joey Fox https://twitter.com/joeyfox85/status/1577071654311989254
#  2) My color timing was off. 10 mJ has maximum disiminability for my eye, but it still 
#     wasn't easy to stop at the right time, and I had variability. 
#     If I was biased slow, that could also contribute. The ambient light for the r=0
#     measurement was brighter and more direct, and so this can be part of it.
# Either way, I'll fit the shape to the observed data excluding r=0, but will simulate scaled up.

# fit full model
m <- nls(power ~ p*(1-(1/sqrt(1+(R/r)^2))) * cos(theta/(T/2)*pi/2),
         data=d %>% filter(r!=0),
         start=list(R=6,p=1000,T=135))
summary(m)  # note that I fit the power and the cone angle low relative to the manufacturer spec!
            # I'm pretty sure this has more to do with errors in my data than with the device.
                # for example, the smaller cone angle and smaller peak power are consistent with
                # my targets were off-vertical-axis, probably low. 
            # The smaller effective R estimate is not surprising to me, but could also be related in error.
            # Reason it isn't surprising is the mirror isn't as deep, and the cone is wider,
            # So I expect that the characteristic length scale to enter 1/r^2 is shorter for this than the Lily.
            # which is to say, the x one "looks more point-source-y."

d$fit = predict(m,newdata=d)

ggplot(d) +
  geom_point(aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=fit,group=theta_factor,color=theta_factor),linetype='solid') +
  geom_point(aes(x=r,y=fit,group=theta_factor,color=theta_factor),shape=0) +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.intensity_model.png',height=3.5,width=6)

# residual is worse for this than UV Can Lily but not surprising since the data are worse too!
idx = d$r!=0
sqrt(sum((d$power[idx]-d$fit[idx])^2)/sum(idx))
ggplot(d %>% filter(idx)) +
  geom_point(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('residual [mu-W / cm^2')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.intensity_model.residuals.png',height=3.5,width=6)


coef(m)
# power_xone=1500 # mu-W/cm2  https://www.ergo-healthtech.com/xone
# cone_angle_xone = 135
power_xone=1054 # observed
cone_angle_xone = 100 #observed
eff_R_xone=4.1 # observed

pred_df = expand.grid(r=seq(0,50,by=1),theta=sort(unique(d$theta)))
pred_df$power = uv_can_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = TRUE,
                             cone_angle=cone_angle_xone,power=power_xone,eff_R=eff_R_xone)
pred_df$power2 = uv_can_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = FALSE,
                              cone_angle=cone_angle_xone,power=power_xone,eff_R=eff_R_xone)
pred_df$theta_factor=factor(pred_df$theta, levels=sort(unique(pred_df$theta)))

# plot fit and interpolation
ggplot() +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  geom_line(data=pred_df,aes(x=r,y=power2,group=theta_factor,color=theta_factor),linetype='dashed') +
  geom_point(data=d,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')

# final model with interpolation
ggplot() +
  geom_point(data=d,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.intensity_model.png',height=3.5,width=6)


##########################################
### start looking at the light cone
#########################################


# make intensity heat plot by gridding cone

cone_field=expand.grid(x=seq(0,50,by=1),z=seq(0,50,by=1))
cone_field$r = sqrt((cone_field$x-25)^2 + cone_field$z^2)
cone_field$theta=360/(2*pi)*atan((cone_field$x-25)/(cone_field$z))
cone_field$theta[cone_field$r==0]=0 # deal with limit at r=0
cone_field$power = uv_can_model(r=cone_field$r,theta=cone_field$theta,
                                cone_angle=cone_angle_xone,power=power_xone,eff_R=eff_R_xone) # estimated parameters
                                # cone_angle=135,power=1500,eff_R=eff_R_xone) # manufacturer-specified parameters

# linear intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=power)) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100,300),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_continuous(name='mu-W / cm^2') 
ggsave('ergo_x_one_ushio_222_card_experiment.observed.intensity_cone.png',height=3.5,width=5)

# log intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmax(power,1))) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100,300),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='mu-W / cm^2',trans='log10',na.value='black')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.intensity_cone_log10.png',height=3.5,width=5)


######################################
#### let's estimate some killing times
######################################

# energy required to reduce virus to exp(-1) of original concentration.
# equivalent to 1 additional Air Change (AC) in well-mixed ventilation
ac_equiv_energy = 1/12.4*1e3 # mu-J-cm2 https://twitter.com/famulare_mike/status/1573514127448100865

# 10x reduction is 2.3 AC equivalent
# 100x reduction is 4.6 AC equivalent

##### time required to get integrated power for 1 AC equivalent neutralization
cone_field$t_ac_equiv = ac_equiv_energy/cone_field$power

# log killing time for 1 eAC
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to 1 equivalent air change (68% reduction)')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.killing_time_1_equiv_air_change.png',height=3.5,width=5)

# log killing time for 50%
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*0.69,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*0.69, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to kill 50%')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.killing_time_50_percent.png',height=3.5,width=5)


# log killing time for 90% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 90%')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.killing_time_90_percent.png',height=3.5,width=5)

# log killing time for 99% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3*2,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3*2, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 99%')
ggsave('ergo_x_one_ushio_222_card_experiment.observed.killing_time_99_percent.png',height=3.5,width=5)




##########################################
### Rerun with manufacturer expected specs, as they are possibly more correct
### than my biased measurements
#########################################


# make intensity heat plot by gridding cone

cone_field=expand.grid(x=seq(0,50,by=1),z=seq(0,50,by=1))
cone_field$r = sqrt((cone_field$x-25)^2 + cone_field$z^2)
cone_field$theta=360/(2*pi)*atan((cone_field$x-25)/(cone_field$z))
cone_field$theta[cone_field$r==0]=0 # deal with limit at r=0
cone_field$power = uv_can_model(r=cone_field$r,theta=cone_field$theta,
                                cone_angle=135,power=1500,eff_R=eff_R_xone) # manufacturer-specified parameters

# linear intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=power)) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100,300),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_continuous(name='mu-W / cm^2') 
ggsave('ergo_x_one_ushio_222_card_experiment.expected.intensity_cone.png',height=3.5,width=5)

# log intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmax(power,1))) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100,300),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='mu-W / cm^2',trans='log10',na.value='black')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.intensity_cone_log10.png',height=3.5,width=5)


######################################
#### let's estimate some killing times
######################################

##### time required to get integrated power for 1 AC equivalent neutralization
cone_field$t_ac_equiv = ac_equiv_energy/cone_field$power

# log killing time for 1 eAC
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to 1 equivalent air change (68% reduction)')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.killing_time_1_equiv_air_change.png',height=3.5,width=5)

# log killing time for 50%
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*0.69,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*0.69, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to kill 50%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.killing_time_50_percent.png',height=3.5,width=5)


# log killing time for 90% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 90%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.killing_time_90_percent.png',height=3.5,width=5)

# log killing time for 99% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3*2,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3*2, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 99%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.killing_time_99_percent.png',height=3.5,width=5)


#########################
######### expected field on bigger domain, reflecting large dining table

# make intensity heat plot by gridding cone

cone_field=expand.grid(x=seq(0,100,by=1),z=seq(0,75,by=1))
cone_field$r = sqrt((cone_field$x-50)^2 + cone_field$z^2)
cone_field$theta=360/(2*pi)*atan((cone_field$x-50)/(cone_field$z))
cone_field$theta[cone_field$r==0]=0 # deal with limit at r=0
cone_field$power = uv_can_model(r=cone_field$r,theta=cone_field$theta,
                                cone_angle=135,power=1500,eff_R=eff_R_xone) # manufacturer-specified parameters

# log intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmax(power,1))) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(1,3,10,30,100,300),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='mu-W / cm^2',trans='log10',na.value='black')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.table.intensity_cone_log10.png',height=3.5,width=5)


######################################
#### let's estimate some killing times
######################################

##### time required to get integrated power for 1 AC equivalent neutralization
cone_field$t_ac_equiv = ac_equiv_energy/cone_field$power

# log killing time for 1 eAC
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to 1 equivalent air change (68% reduction)')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.table.killing_time_1_equiv_air_change.png',height=3.5,width=5)

# log killing time for 50%
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*0.69,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*0.69, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7")+
  ggtitle('time to kill 50%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.table.killing_time_50_percent.png',height=3.5,width=5)


# log killing time for 90% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 90%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.table.killing_time_90_percent.png',height=3.5,width=5)

# log killing time for 99% reduction
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(t_ac_equiv*2.3*2,100))) +
  stat_contour2(aes(x=x,y=z,z = t_ac_equiv*2.3*2, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='seconds',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('time to kill 99%')
ggsave('ergo_x_one_ushio_222_card_experiment.expected.table.killing_time_99_percent.png',height=3.5,width=5)

