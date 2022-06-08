#%% 
G           = 6.67e-11
Mb          = 4.0e30                    # black hole
Ms          = 2.0e30                    # sun
Me          = 5.972e24                  # earth        
Mm          = 6.39e23                   # mars
Mc          = 6.39e20                   # unknown comet
AU          = 1.5e11
daysec      = 24.0*60*60

e_ap_v      = 29290                     # earth velocity at aphelion
m_ap_v      = 21970                     # mars velocity at aphelion
commet_v    = 7000

gravconst_e = G*Me*Ms
gravconst_m = G*Mm*Ms
gravconst_c = G*Mc*Ms

# setup the starting conditions
# earth
xe,ye,ze    = 1.0167*AU,0,0
xve,yve,zve = 0,e_ap_v,0

# mars
xm,ym,zm    = 1.666*AU,0,0
xvm,yvm,zvm = 0,m_ap_v,0

# comet
xc,yc,zc    = 2*AU,0,0
xvc,yvc,zvc = 0,commet_v,3000

# sun
xs,ys,zs    = 0,0,0
xvs,yvs,zvs = 0,0,0

t           = 0.0
dt          = 1*daysec # every frame move this time

xelist,yelist,zelist = [],[],[]
xslist,yslist,zslist = [],[],[]
xmlist,ymlist,zmlist = [],[],[]
xclist,yclist,zclist = [],[],[]

# start simulation
while t<5*365*daysec:
    ################ earth #############
    # compute G force on earth
    #rx,ry,rz = xs - xe, ys - ye, zs - ze
    rx,ry,rz = xe - xs, ye - ys, ze - zs
    modr3_e = (rx**2+ry**2+rz**2)**1.5
    fx_e = -gravconst_e*rx/modr3_e
    fy_e = -gravconst_e*ry/modr3_e
    fz_e = -gravconst_e*rz/modr3_e
    
    # update quantities how is this calculated?  F = ma -> a = F/m
    xve += fx_e*dt/Me
    yve += fy_e*dt/Me
    zve += fz_e*dt/Me
    
    # update position
    xe += xve*dt
    ye += yve*dt 
    ze += zve*dt
    
    # save the position in list
    xelist.append(xe)
    yelist.append(ye)
    zelist.append(ze)
    
    ################ Mars ##############
    # compute G force on mars
    rx_m,ry_m,rz_m = xm - xs, ym - ys, zm - zs
    modr3_m = (rx_m**2+ry_m**2+rz_m**2)**1.5
    fx_m = -gravconst_m*rx_m/modr3_m
    fy_m = -gravconst_m*ry_m/modr3_m
    fz_m = -gravconst_m*rz_m/modr3_m
    
    xvm += fx_m*dt/Mm
    yvm += fy_m*dt/Mm
    zvm += fz_m*dt/Mm
    
    # update position
    xm += xvm*dt
    ym += yvm*dt
    zm += zvm*dt
    
    # add to list
    xmlist.append(xm)
    ymlist.append(ym)
    zmlist.append(zm)
    
    ################ comet ##############
    # compute G force on comet
    rx_c,ry_c,rz_c = xc - xs, yc - ys, zc - zs
    modr3_c = (rx_c**2+ry_c**2+rz_c**2)**1.5
    fx_c = -gravconst_c*rx_c/modr3_c
    fy_c = -gravconst_c*ry_c/modr3_c
    fz_c = -gravconst_c*rz_c/modr3_c
    
    xvc += fx_c*dt/Mc
    yvc += fy_c*dt/Mc
    zvc += fz_c*dt/Mc
    
    # update position
    xc += xvc*dt
    yc += yvc*dt 
    zc += zvc*dt
    
    # add to list
    xclist.append(xc)
    yclist.append(yc)
    zclist.append(zc)
    
    ################ the sun ###########
    # update quantities how is this calculated?  F = ma -> a = F/m
    xvs += -(fx_e+fx_m)*dt/Ms
    yvs += -(fy_e+fy_m)*dt/Ms
    zvs += -(fz_e+fz_m)*dt/Ms
    
    # update position
    xs += xvs*dt
    ys += yvs*dt 
    zs += zvs*dt
    xslist.append(xs)
    yslist.append(ys)
    zslist.append(zs)
    
    # update dt
    t +=dt

print('data ready')
#print(xalist,yalist)

#ploy it in ipy
from matplotlib.pyplot import xlabel
import numpy as np
import ipyvolume as ipv

# configure style
fig = ipv.figure(width=800,height=600)
ipv.pylab.style.axes_off()
ipv.pylab.xyzlim(vmin=-2*AU,vmax=2*AU)
ipv.pylab.style.use('dark')

# sun
s_p_x,s_p_y,s_p_z = np.array(xslist).reshape(-1,1),np.array(yslist).reshape(-1,1),np.array(zslist).reshape(-1,1)
s_p = ipv.pylab.scatter(s_p_x,s_p_y,s_p_z,size=5,marker='sphere',color='yellow')

# earth 
e_p_x,e_p_y,e_p_z = np.array(xelist).reshape(-1,1),np.array(yelist).reshape(-1,1),np.array(zelist).reshape(-1,1)
e_p = ipv.pylab.scatter(e_p_x,e_p_y,e_p_z,size=3,marker='sphere',color='blue')

e_l_x,e_l_y,e_l_z = [],[],[]
for i in range(len(xelist)):
    elxt = np.empty(len(xelist)); elxt[:]=np.NaN
    elxt[:i+1] = xelist[:i+1]
    e_l_x.append(elxt)
    elyt = np.zeros(len(yelist))
    elyt[:i+1] = yelist[:i+1]
    e_l_y.append(elyt)
    elzt = np.zeros(len(zelist))
    elzt[:i+1] = zelist[:i+1]
    e_l_z.append(elzt)
e_l_x = np.array(e_l_x)
e_l_y = np.array(e_l_y)
e_l_z = np.array(e_l_z)
e_l = ipv.pylab.plot(e_l_x,e_l_y,e_l_z,size=1,color='white')

# mars
m_p_x,m_p_y,m_p_z = np.array(xmlist).reshape(-1,1),np.array(ymlist).reshape(-1,1),np.array(zmlist).reshape(-1,1)
m_p = ipv.pylab.scatter(m_p_x,m_p_y,m_p_z,size=2,marker='sphere',color='red')

m_l_x,m_l_y,m_l_z = [],[],[]
for i in range(len(xmlist)):
    mlxt = np.empty(len(xmlist)); mlxt[:]=np.NaN
    mlxt[:i+1] = xmlist[:i+1]
    m_l_x.append(mlxt)
    mlyt = np.zeros(len(ymlist))
    mlyt[:i+1] = ymlist[:i+1]
    m_l_y.append(mlyt)
    mlzt = np.zeros(len(zmlist))
    mlzt[:i+1] = zmlist[:i+1]
    m_l_z.append(mlzt)
m_l_x = np.array(m_l_x)
m_l_y = np.array(m_l_y)
m_l_z = np.array(m_l_z)
m_l = ipv.pylab.plot(m_l_x,m_l_y,m_l_z,size=1,color='white')

# comet
c_p_x,c_p_y,c_p_z = np.array(xclist).reshape(-1,1),np.array(yclist).reshape(-1,1),np.array(zclist).reshape(-1,1)
c_p = ipv.pylab.scatter(c_p_x,c_p_y,c_p_z,size=2,marker='sphere',color='gray')

c_l_x,c_l_y,c_l_z = [],[],[]
for i in range(len(xclist)):
    clxt = np.empty(len(xclist)); clxt[:]=np.NaN
    clxt[:i+1] = xclist[:i+1]
    c_l_x.append(clxt)
    clyt = np.zeros(len(yclist))
    clyt[:i+1] = yclist[:i+1]
    c_l_y.append(clyt)
    clzt = np.zeros(len(zclist))
    clzt[:i+1] = zclist[:i+1]
    c_l_z.append(clzt)
c_l_x = np.array(c_l_x)
c_l_y = np.array(c_l_y)
c_l_z = np.array(c_l_z)
c_l = ipv.pylab.plot(c_l_x,c_l_y,c_l_z,size=1,color='white')

ipv.animation_control([s_p,e_p,e_l,m_p,m_l,c_p,c_l],interval=20)
ipv.show()