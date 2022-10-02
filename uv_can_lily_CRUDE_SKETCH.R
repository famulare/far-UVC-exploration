# estimating simple model of UV-Can Lily intensity field
# Mike Famulare
# 
# THIS IS A CRUDE SKETCH PULLED FROM A PICTURE ON TWITTER NOT REALLY INTENDED FOR THIS PURPOSE!
# THE MODELING IS VERY CRUDE. I HAVE NO PROFESSIONAL EXPERTISE IN FAR UV-C STERILIZATION.
# THIS IS NOT PROFESSIONAL ENGINEERING!!! JUST AN INTERESTED PERSON LEARNING WITH YOU.
# READER BEWARE!!!!!!

library(tidyverse)
library(metR) # for better contour plots

# UV can Lily irradiance data from Joey Fox
# https://twitter.com/joeyfox85/status/1575882708525682688/photo/1
# NO BLAME ATTACHED TO HIM! I PULLED THIS FROM HIS SKETCH WITHOUT ASKING.

uv_can = data.frame(r=c(0,2,5,10,15,20,30,40,50,
                        10,15,20,30,40,50,
                        10,10,15,15,15,
                        20,20,30,30,
                        40,50),
                    theta=c(0,0,0,0,0,0,0,0,0,
                            20,15,15,22,15,12,
                            25,15,22,15,8, # these are - angles, but I'm assuming symmetry
                            20,10,15,10,# these are - angles, but I'm assuming symmetry
                            20,10),# these are - angles, but I'm assuming symmetry
                    power=c(497,370,185,70,35,21,10,7,4, # mu-W/cm2
                            39,28,13,4.5,5,3.2,
                            3,41,6,16,30,
                            5,15,3.5,6,
                            3.5,3.3))

# discretize angles for plotting later
uv_can$theta_factor = factor(uv_can$theta,levels=sort(unique(uv_can$theta)))


# axial fit first to get a feel for the problem. 
# (Note this code doesn't reflect how I got here. I wasn't committing while playing.)

# geometry of axial e-field from disk https://www.physics.udel.edu/~watson/phys208/exercises/kevan/efield1.html
# goes like 1/r^2 at long radius
m <- nls(power ~ p*(1-(1/sqrt(1+(R/(r-r0))^2))),
         data=uv_can %>% filter(theta==0),
         start=list(R=6,p=500,r0=0))
summary(m)
coef(m) # marketed intensity at opening is 490+ https://uv-can.com/products/lily-handheld-personal-far-uv-disinfection-light
        # small and positive r0 looks like overfitting and makes little difference, so I will leave it out later on.
        # empirical "disk radius" of ~6 cm is much larger than the device face.
        # this is a function of the mirror geometry and placement of the light source in it,
        # but we don't know much about those and I'm not about to guess and ray trace.
uv_can$fit = predict(m,newdata=uv_can)

ggplot(uv_can %>% filter(theta==0)) +
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
m <- nls(power ~ p*(1-(1/sqrt(1+(R/(r))^2))) * cos(theta/(T/2)*pi/2),
         data=uv_can,
         start=list(R=6,p=500,T=60))
summary(m)

coef(m)
# note that cone angle fits at 52 degrees. I suspect my angle guesses from the picture are a bit off,
# and will use 60 degrees in modeling later, after looking at fit.

uv_can$fit = predict(m,newdata=uv_can)

ggplot(uv_can) +
  geom_point(aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=fit,group=theta_factor,color=theta_factor),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')

# residuals are well-behaved: sd = 6 mu-W/cm2 and none above 20 mu-W (where also disk assumption isn't good enough) 
sqrt(sum((uv_can$power-uv_can$fit)^2)/nrow(uv_can))
ggplot(uv_can) +
  geom_point(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  geom_line(aes(x=r,y=power-fit,group=theta_factor,color=theta_factor)) +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('residual [mu-W / cm^2')
ggsave('uv_can_lily_CRUDE_SKETCH.intensity_model.residuals.png',height=3.5,width=6)


###########################################################
## flesh out prediction set and code up model as a function
###########################################################
coef(m)
uv_can_model = function(r,theta,
                        eff_R=6.2,cone_angle=52,power=500,
                        near_field_hack=TRUE){
  disk_r = 1-(r/sqrt(eff_R^2+r^2))
  
  if(near_field_hack){
    # hack for finite area of opening, because near the light, everything is axial.
    # when within eff_R, interpolate to theta=0 power.
    
    # interpolate with exponential smoothing, chosen so that interpolation is "turned off" by eff_R
    # Note I have no justification for this functional form other than smoothness.
    # see more discussion in comments around line ~50
    p_r_theta = (power * (1-(r/sqrt(eff_R^2+r^2)))) * (exp(-5*r/eff_R) + cos(theta*(2/cone_angle)*pi/2)*(1-exp(-5*r/eff_R)))
    
  }else{
    # using fit to only data from r > eff_R. This is surely more wrong close to the light than the hack.
    p_r_theta= power * (1-(r/sqrt(eff_R^2+r^2))) * cos(theta*(2/cone_angle)*pi/2)
  }
  
  # remove stuff outside the cone
  outside_idx = abs(theta)>(cone_angle/2)
  p_r_theta[outside_idx]=0
  
  return(p_r_theta)
}
  
pred_df = expand.grid(r=seq(0,50,by=1),theta=sort(unique(uv_can$theta)))
pred_df$power = uv_can_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = TRUE)
pred_df$power2 = uv_can_model(r=pred_df$r,theta=pred_df$theta,near_field_hack = FALSE)
pred_df$theta_factor=factor(pred_df$theta, levels=sort(unique(pred_df$theta)))

# plot fit and interpolation
ggplot() +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  geom_line(data=pred_df,aes(x=r,y=power2,group=theta_factor,color=theta_factor),linetype='dashed') +
  geom_point(data=uv_can,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')

# final model with interpolation
ggplot() +
  geom_point(data=uv_can,aes(x=r,y=power,group=theta_factor,color=theta_factor)) +
  geom_line(data=pred_df,aes(x=r,y=power,group=theta_factor,color=theta_factor),linetype='solid') +
  scale_y_continuous(trans='log10') +
  # scale_x_continuous(trans='log10') +
  scale_color_discrete(name='degrees off axis') + 
  theme_bw() + xlab('distance [cm]') + ylab('power [mu-W / cm^2]')
ggsave('uv_can_lily_CRUDE_SKETCH.intensity_model.png',height=3.5,width=6)


##########################################
### start looking at the light cone
#########################################


# make intensity heat plot by gridding cone

cone_field=expand.grid(x=seq(0,50,by=1),z=seq(0,50,by=1))
cone_field$r = sqrt((cone_field$x-25)^2 + cone_field$z^2)
cone_field$theta=360/(2*pi)*atan((cone_field$x-25)/(cone_field$z))
cone_field$theta[cone_field$r==0]=0 # deal with limit at r=0
cone_field$power = uv_can_model(r=cone_field$r,theta=cone_field$theta,cone_angle=60,power=490) # manufacturer-specified parameters

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
ggsave('uv_can_lily_CRUDE_SKETCH.intensity_cone.png',height=3.5,width=5)

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
ggsave('uv_can_lily_CRUDE_SKETCH.intensity_cone_log10.png',height=3.5,width=5)


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
ggsave('uv_can_lily_CRUDE_SKETCH.killing_time_1_equiv_air_change.png',height=3.5,width=5)

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
ggsave('uv_can_lily_CRUDE_SKETCH.killing_time_90_percent.png',height=3.5,width=5)

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
ggsave('uv_can_lily_CRUDE_SKETCH.killing_time_99_percent.png',height=3.5,width=5)

