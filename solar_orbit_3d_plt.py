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


#%% 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import animation

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.axis('auto')

axis_size = 2.5
ax.set_xlim(-axis_size*AU,axis_size*AU)
ax.set_ylim(-axis_size*AU,axis_size*AU)
ax.set_zlim(-axis_size*AU,axis_size*AU)

# ax.set_aspect('auto')
# ax.grid()
datadict = {}
dataset_s = [xslist,yslist,zslist]
dataset_e = [xelist,yelist,zelist]
dataset_m = [xmlist,ymlist,zmlist]
dataset_c = [xclist,yclist,zclist]
datadict['s'] = dataset_s
datadict['e'] = dataset_e
datadict['m'] = dataset_m
datadict['c'] = dataset_c

vis_dict = {}
# sun
line_s,     = ax.plot([0],[0],[0],'-g',lw=1)
point_s,    = ax.plot([AU],[0],[0], marker="o", markersize=7, markeredgecolor="yellow", markerfacecolor="yellow")
text_s      = ax.text(AU,0,0,'Sun')
vis_dict['s'] = [line_s,point_s,text_s]

# earth
line_e,     = ax.plot([0],[0],[0],'-g',lw=1)
point_e,    = ax.plot([AU],[0],[0], marker="o", markersize=4, markeredgecolor="blue", markerfacecolor="blue")
text_e      = ax.text(AU,0,0,'Earth')
vis_dict['e'] = [line_e,point_e,text_e]

# mars 
line_m,     = ax.plot([0],[0],[0],'-g',lw=1)
point_m,    = ax.plot([AU],[0],[0], marker="o", markersize=4, markeredgecolor="red", markerfacecolor="red")
text_m      = ax.text(AU,0,0,'Mars')
vis_dict['m'] = [line_m,point_m,text_m]

# comet
line_c,     = ax.plot([0],[0],[0],'-g',lw=1)
point_c,    = ax.plot([AU],[0],[0], marker="o", markersize=4, markeredgecolor="black", markerfacecolor="black")
text_c      = ax.text(AU,0,0,'comet')
vis_dict['c'] = [line_c,point_c,text_c]

def update(num,data_dict,vis_dict):
    # sun 
    dataset_s               = data_dict['s']
    line_s,point_s,text_s   = vis_dict['s'][0],vis_dict['s'][1],vis_dict['s'][2]
    line_s.set_data_3d(dataset_s[0][:num],dataset_s[1][:num],dataset_s[2][:num])
    point_s.set_data_3d(dataset_s[0][num],dataset_s[1][num],dataset_s[2][num])
    text_s.set_position((dataset_s[0][num],dataset_s[1][num],dataset_s[2][num]))
    
    # earth 
    dataset_e               = data_dict['e']
    line_e,point_e,text_e   = vis_dict['e'][0],vis_dict['e'][1],vis_dict['e'][2]
    line_e.set_data_3d(dataset_e[0][:num],dataset_e[1][:num],dataset_e[2][:num])
    point_e.set_data_3d(dataset_e[0][num],dataset_e[1][num],dataset_e[2][num])
    text_e.set_position((dataset_e[0][num],dataset_e[1][num],dataset_e[2][num]))
    
    # mars
    dataset_m               = data_dict['m']
    line_m,point_m,text_m   = vis_dict['m'][0],vis_dict['m'][1],vis_dict['m'][2]
    line_m.set_data_3d(dataset_m[0][:num],dataset_m[1][:num],dataset_m[2][:num])
    point_m.set_data_3d(dataset_m[0][num],dataset_m[1][num],dataset_m[2][num])
    text_m.set_position((dataset_m[0][num],dataset_m[1][num],dataset_m[2][num]))
    
    # mars
    dataset_c               = data_dict['c']
    line_c,point_c,text_c   = vis_dict['c'][0],vis_dict['c'][1],vis_dict['c'][2]
    line_c.set_data_3d(dataset_c[0][:num],dataset_c[1][:num],dataset_c[2][:num])
    point_c.set_data_3d(dataset_c[0][num],dataset_c[1][num],dataset_c[2][num])
    text_c.set_position((dataset_c[0][num],dataset_c[1][num],dataset_c[2][num]))

ani = animation.FuncAnimation(
    fig
    ,update
    ,len(xelist)
    ,fargs=(datadict, vis_dict)
    ,interval=1
)

plt.show()