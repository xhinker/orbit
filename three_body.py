#%% 
G           = 6.67e-11
Ma          = 3.0e30                    # black hole
Mb          = 2.0e30                    # sun
Mc          = 3.0e30                    # star c
AU          = 1.5e11
daysec      = 24.0*60*60

gravconst_ab = G*Ma*Mb
gravconst_bc = G*Mb*Mc
gravconst_ca = G*Mc*Ma

# setup the starting conditions
# a
xa,ya,za    = -2*AU,0,0
xva,yva,zva = 10000,-10000,0

# b
xb,yb,zb    = 0,2*AU,0
xvb,yvb,zvb = -10000,-10000,0

# c
xc,yc,zc    = 3*AU,0,0
xvc,yvc,zvc = -10000,10000,0

t           = 0.0
dt          = 0.2*daysec # every frame move this time

xalist,yalist,zalist = [],[],[]
xblist,yblist,zblist = [],[],[]
xclist,yclist,zclist = [],[],[]

# start simulation
while t<3*365*daysec:
    # Common numbers 
    rx_ab,ry_ab,rz_ab = xb - xa, yb - ya, zb - za
    modr3_ab = (rx_ab**2+ry_ab**2+rz_ab**2)**1.5

    rx_bc,ry_bc,rz_bc = xc-xb,yc-yb,zc-zb
    modr3_bc = (rx_bc**2+ry_bc**2+rz_bc**2)**1.5

    rx_ca,ry_ca,rz_ca = xa-xc,ya-yc,za-zc
    modr3_ca = (rx_ca**2+ry_ca**2+rz_ca**2)**1.5
    ################ a #############
    # compute G force on a
    fx_a_b = gravconst_ab*rx_ab/modr3_ab
    fy_a_b = gravconst_ab*ry_ab/modr3_ab
    fz_a_b = gravconst_ab*rz_ab/modr3_ab
    
    fx_a_c = -gravconst_ca*rx_ca/modr3_ca
    fy_a_c = -gravconst_ca*ry_ca/modr3_ca
    fz_a_c = -gravconst_ca*rz_ca/modr3_ca

    # update quantities how is this calculated?  F = ma -> a = F/m
    xva += (fx_a_b+fx_a_c)*dt/Ma
    yva += (fy_a_b+fy_a_c)*dt/Ma
    zva += (fz_a_b+fz_a_c)*dt/Ma
    
    # update position
    xa += xva*dt
    ya += yva*dt 
    za += zva*dt
    
    # save the position in list
    xalist.append(xa)
    yalist.append(ya)
    zalist.append(za)

    ################ b #############
    # compute G force on a
    fx_b_a = -gravconst_ab*rx_ab/modr3_ab
    fy_b_a = -gravconst_ab*ry_ab/modr3_ab
    fz_b_a = -gravconst_ab*rz_ab/modr3_ab
    
    fx_b_c = gravconst_bc*rx_bc/modr3_bc
    fy_b_c = gravconst_bc*ry_bc/modr3_bc
    fz_b_c = gravconst_bc*rz_bc/modr3_bc

    # update quantities how is this calculated?  F = ma -> a = F/m
    xvb += (fx_b_a+fx_b_c)*dt/Mb
    yvb += (fy_b_a+fy_b_c)*dt/Mb
    zvb += (fz_b_a+fz_b_c)*dt/Mb
    
    # update position
    xb += xvb*dt
    yb += yvb*dt 
    zb += zvb*dt
    
    # save the position in list
    xblist.append(xb)
    yblist.append(yb)
    zblist.append(zb)
    
    ################ c #############
    # compute G force on a
    fx_c_a = gravconst_ca*rx_ca/modr3_ca
    fy_c_a = gravconst_ca*ry_ca/modr3_ca
    fz_c_a = gravconst_ca*rz_ca/modr3_ca
    
    fx_c_b = -gravconst_bc*rx_bc/modr3_bc
    fy_c_b = -gravconst_bc*ry_bc/modr3_bc
    fz_c_b = -gravconst_bc*rz_bc/modr3_bc

    # update quantities how is this calculated?  F = ma -> a = F/m
    xvc += (fx_c_a+fx_c_b)*dt/Mc
    yvc += (fy_c_a+fy_c_b)*dt/Mc
    zvc += (fz_c_a+fz_c_b)*dt/Mc
    
    # update position
    xc += xvc*dt
    yc += yvc*dt 
    zc += zvc*dt
    
    # save the position in list
    xclist.append(xc)
    yclist.append(yc)
    zclist.append(zc)
    
    # update dt
    t +=dt

print('data ready')

#%% plot the final 
# import matplotlib.pyplot as plt
# plt.plot(xalist,yalist,'-g',lw=2,c='red')
# plt.plot(xblist,yblist,'-g',lw=2,c='blue')
# plt.plot(xclist,yclist,'-g',lw=2,c='black')
# plt.axis('equal')
# plt.show()

#%%
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib
matplotlib.rcParams['animation.embed_limit'] = 2**128
matplotlib.use("TkAgg") # for mac
from IPython.display import HTML

fig, ax = plt.subplots(figsize=(10,10))
ax.set_aspect('equal')
ax.grid()

line_a,     = ax.plot([],[],'-g',lw=1)
point_a,    = ax.plot([AU], [0], marker="o", markersize=4, markeredgecolor="blue", markerfacecolor="blue")
text_a      = ax.text(AU,0,'A')

line_b,     = ax.plot([],[],'-g',lw=1)
point_b,    = ax.plot([1.5*AU], [0], marker="o", markersize=3, markeredgecolor="red", markerfacecolor="red")
text_b      = ax.text(1.666*AU,0,'B')

line_c,     = ax.plot([],[],'-g',lw=1)
point_c,    = ax.plot([2*AU], [0], marker="o", markersize=2, markeredgecolor="black", markerfacecolor="black")
text_c      = ax.text(2*AU,0,'C')


axdata,aydata = [],[]                   # earth track
bxdata,bydata = [],[]                   # sun track
cxdata,cydata = [],[]                   # mars track

print(len(xalist))

def update(i):
    axdata.append(xalist[i])
    aydata.append(yalist[i])
    
    bxdata.append(xblist[i])
    bydata.append(yblist[i])
    
    cxdata.append(xclist[i])
    cydata.append(yclist[i])
    
    line_a.set_data(axdata,aydata)
    point_a.set_data(xalist[i],yalist[i])
    text_a.set_position((xalist[i],yalist[i]))
    
    line_b.set_data(bxdata,bydata)
    point_b.set_data(xblist[i],yblist[i])
    text_b.set_position((xblist[i],yblist[i]))
    
    line_c.set_data(cxdata,cydata)
    point_c.set_data(xclist[i],yclist[i])
    text_c.set_position((xclist[i],yclist[i]))

    ax.axis('equal')
    ax.set_xlim(-4*AU,4*AU)
    ax.set_ylim(-4*AU,4*AU)
    #print(i)
    return line_a,point_a,text_a,line_b,point_b,text_b,line_c,point_c,text_c

anim = animation.FuncAnimation(fig,func=update,frames=len(xalist),interval=1,blit=True)
plt.show()
