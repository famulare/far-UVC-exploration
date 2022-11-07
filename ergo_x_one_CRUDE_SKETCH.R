# estimating simple model of Ergo X One irradiance field
# Mike Famulare
# 
# THIS IS A CRUDE SKETCH PULLED FROM A PICTURE ON TWITTER NOT REALLY INTENDED FOR THIS PURPOSE!
# THE MODELING IS VERY CRUDE. I HAVE NO PROFESSIONAL EXPERTISE IN FAR UV-C STERILIZATION.
# THIS IS NOT PROFESSIONAL ENGINEERING!!! JUST AN INTERESTED PERSON LEARNING WITH YOU.
# READER BEWARE!!!!!!

library(tidyverse)
library(metR) # for better contour plots

# Ergo X One irradiance data from Joey Fox
# https://twitter.com/joeyfox85/status/1589708278409220096
# NO BLAME ATTACHED TO HIM! I PULLED THIS FROM HIS SKETCH WITHOUT ASKING.

# there is an assymetry, but it's not big, and only noticeable close to the device, so I'm using averages
ergo = data.frame(r=c(0,5,10,15,20,30,40,50,60,80,
                      10,15,20,30,40,50,60,
                      50,30,15,
                      10,15,20,30,40,50,60),
                    theta=c(0,0,0,0,0,0,0,0,0,0,
                            36,32,29,23,23,18,21,
                            50,39,33,
                            -42,-45,-51,-32,-43,-35,-36),# angles from normal, assuming symmetry
                    power=c(1611,391,146,76,47,21,13,8,5,3,
                            53,46,36,18,12,8,5,
                            5,5,3,
                            87,22,11,18,4,5,3))# mu-W/cm2

# discretize angles for plotting later
ergo$theta_factor = factor(5*round(abs(ergo$theta)/5),levels=sort(unique(5*round(abs(ergo$theta)/5))))


# axial fit first to get a feel for the problem. 

# geometry of axial e-field from disk https://www.physics.udel.edu/~watson/phys208/exercises/kevan/efield1.html
# goes like 1/r^2 at long radius
m <- nls(power ~ p*(1-(1/sqrt(1+(R/(r-r0))^2))),
         data=ergo %>% filter(theta==0),
         start=list(R=6,p=500,r0=0))
summary(m)
coef(m) # marketed intensity at opening is 490+ https://uv-can.com/products/lily-handheld-personal-far-uv-disinfection-light
        # small and positive r0 looks like overfitting and makes little difference, so I will leave it out later on.
        # empirical "disk radius" of ~5 cm is much larger than the device face.
        # this is a function of the mirror geometry and placement of the light source in it,
        # but we don't know much about those and I'm not about to guess and ray trace.
ergo$fit = predict(m,newdata=ergo)

ggplot(ergo %>% filter(theta==0)) +
  geom_point(aes(x=r,y=power)) +
  geom_line(aes(x=r,y=fit),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
  

# off-axis fit by assuming axial r-dependence and simple angle dependence
  # cosine angle dependence guessed because it's the lowest order term in any taylor series with the right symmetry
  # also not crazy to think a large lamp positioned inside a small reflector is not unlike
  # a lambertian source within the cone https://en.wikipedia.org/wiki/Lambert%27s_cosine_law, 
  # (which is also not totally inconsistent with the "disk-like" source field observed empirically)
# better would probably spherical harmonics but this is just for sketching! 
# and I've had zero luck hunting down "flashlight" beam physics formulas
  # also note this works since all the off-angle data are far enough into the 
  # far field to be in the 1/r^2 part (up to what looks like measurement floor censoring),
  # and so we don't need to model the changing r-dependence to fit the data 
  # (although we come back to that later!). 
  # In anology with a charged rod of finite aspect,
  # http://online.cctt.org/physicslab/content/phyapc/lessonnotes/Efields/EchargedRods.asp
  # off-axis will enter 1/r^2 faster. 
m <- nls(power ~ p*(1-(1/sqrt(1+(R/(r-r0))^2))) * cos(theta/(T/2)*pi/2),
         data=ergo,
         start=list(R=6,p=500,T=60,r0=0))
summary(m)

coef(m)
# note that cone angle fits at 54 degrees. I suspect my angle guesses from the picture are a bit off,
# and will use 60 degrees in modeling later, after looking at fit.

ergo$fit = predict(m,newdata=ergo)

ggplot(ergo) +
  geom_point(aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=fit,group=theta_factor,color=theta_factor),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
ggsave('ergo_x_one_CRUDE_SKETCH.intensity_model.png',height=3.5,width=6)


# residuals are well-behaved except 3 points at large angles and short-ish distances
sqrt(sum((ergo$power-ergo$fit)^2)/nrow(ergo))
ggplot(ergo) +
  geom_point(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('residual [mu-W / cm^2')
ggsave('ergo_x_one_CRUDE_SKETCH.intensity_model.residuals.png',height=3.5,width=6)


###########################################################
## flesh out prediction set and code up model as a function
###########################################################
coef(m)
ergo_x_one_model = function(r,theta,
                        eff_R=4.3,cone_angle=115,power=1940,r0=-0.75,
                        near_field_hack=TRUE){
  
  axial_power = power * (1-(1/sqrt(1+(eff_R/(r-r0))^2)))
  
  if(near_field_hack){
    # hack for finite area of opening, because near the light, everything is axial.
    # when within eff_R, interpolate to theta=0 power.
    
    # interpolate with exponential smoothing, chosen so that interpolation is "turned off" by eff_R
    # Note I have no justification for this functional form other than smoothness.
    # see more discussion in comments around line ~50
    p_r_theta = axial_power * (exp(-5*r/eff_R) + cos(theta*(2/cone_angle)*pi/2)*(1-exp(-5*r/eff_R)))
    
  }else{
    # using fit to only data from r > eff_R. This is surely more wrong close to the light than the hack.
    p_r_theta= axial_power * cos(theta*(2/cone_angle)*pi/2)
  }
  
  # remove stuff outside the cone
  outside_idx = abs(theta)>(cone_angle/2)
  p_r_theta[outside_idx]=0
  
  return(p_r_theta)
}
  
pred_df = expand.grid(r=seq(0,80,by=1),theta=as.numeric(as.character(levels(ergo$theta_factor))))
pred_df$power = ergo_x_one_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = TRUE)
pred_df$power2 = ergo_x_one_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = FALSE)
pred_df$theta_factor=factor(pred_df$theta, levels=sort(unique(pred_df$theta)))

# plot fit and interpolation
ggplot() +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  geom_line(data=pred_df,aes(x=r,y=power2,group=theta_factor,color=theta_factor),linetype='dashed') +
  geom_point(data=ergo,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')

# final model with interpolation
ggplot() +
  geom_point(data=ergo,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
ggsave('ergo_x_one_CRUDE_SKETCH.intensity_model.png',height=3.5,width=6)


##########################################
### start looking at the light cone
#########################################


# make intensity heat plot by gridding cone

cone_field=expand.grid(x=seq(0,100,by=1),z=seq(0,99,by=1))
cone_field$r = sqrt((cone_field$x-50)^2 + cone_field$z^2)
cone_field$theta=360/(2*pi)*atan((cone_field$x-50)/(cone_field$z))
cone_field$theta[cone_field$r==0]=0 # deal with limit at r=0
cone_field$power = ergo_x_one_model(r=cone_field$r,theta=cone_field$theta)

# linear intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=power)) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_continuous(name='mu-W / cm^2') 
ggsave('ergo_x_one_CRUDE_SKETCH.intensity_cone.png',height=3.5,width=5)

# log intensity
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmax(power,1))) +
  stat_contour2(aes(x=x,y=z,z = power, label = stat(level)),
                color='white',breaks=c(3,10,30,100),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='mu-W / cm^2',trans='log10',na.value='black')
ggsave('ergo_x_one_CRUDE_SKETCH.intensity_cone_log10.png',height=3.5,width=5)


######################################
#### let's estimate some killing times
######################################

# energy required to reduce virus to exp(-1) of original concentration.
# equivalent to 1 additional Air Change (AC) in well-mixed ventilation
ac_equiv_energy = 1/12.4*1e3 # mu-J/cm2 https://twitter.com/famulare_mike/status/1573514127448100865

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
ggsave('ergo_x_one_CRUDE_SKETCH.killing_time_1_equiv_air_change.png',height=3.5,width=5)

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
ggsave('ergo_x_one_CRUDE_SKETCH.killing_time_50_percent.png',height=3.5,width=5)


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
ggsave('ergo_x_one_CRUDE_SKETCH.killing_time_90_percent.png',height=3.5,width=5)

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
ggsave('ergo_x_one_CRUDE_SKETCH.killing_time_99_percent.png',height=3.5,width=5)

## SAFETY
# time until maximum exposure per https://twitter.com/joeyfox85/status/1589708284021600256
eye_max_exposure= 160*1e3 #mu-J/cm².
skin_max_exposure= 479*1e3 #mu-J/cm².

##### time required to get integrated power for safety thresholds
cone_field$h_eye_max = eye_max_exposure/cone_field$power/3600
cone_field$h_skin_max = skin_max_exposure/cone_field$power/3600


# log time to max daily eye exposure
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(h_eye_max,24))) +
  stat_contour2(aes(x=x,y=z,z = h_eye_max, label = stat(level)),
                color='white',breaks=c(1,3,8,12),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='hours',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('hours until maximumum daily eye exposure')
ggsave('ergo_x_one_CRUDE_SKETCH.max_eye_exposure.png',height=3.5,width=5)

# log time to max daily skin exposure
ggplot(data=cone_field) +
  geom_tile(aes(x=x,y=z,fill=pmin(h_skin_max,24))) +
  stat_contour2(aes(x=x,y=z,z = h_skin_max, label = stat(level)),
                color='white',breaks=c(1,3,8,12,16,24),skip = 0) +
  theme_bw() + 
  coord_equal() +
  scale_x_continuous(name='position [cm]',expand=c(0,0))+
  scale_y_continuous(name='height [cm]',expand=c(0,0))+
  scale_fill_gradient(name='hours',trans='log10',na.value='black',high = "#132B43", low = "#56B1F7") + 
  ggtitle('hours until maximumum daily skin exposure')
ggsave('ergo_x_one_CRUDE_SKETCH.max_skin_exposure.png',height=3.5,width=5)

