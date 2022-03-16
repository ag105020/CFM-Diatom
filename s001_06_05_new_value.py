'''
Cleaned from s001 06 01
'''

from pylab import *
from sf import *
from fs import *

rcParams['image.cmap'] = 'RdYlBu_r'

########################################################

Vpower=2        #(dimensionless) power factor for V
V=10**(Vpower)*10**(-18)       #(m3) volume of the cell
R=(3*V/(4*pi))**(1/3)   #(m) radius of the cell
dr=R/100
da=dr
r0=arange(dr,R+dr/10,dr)
r1=r0.reshape(size(r0),1)       #converting a row into a column
r2=r1[::-1]              #reverse the column
padsize=1
pad0=zeros((padsize,1))
pad=pad0+nan
r=vstack((pad,r2))
                        #add a pad to look a figure nice
r_number=arange(0,size(r),1)
a=zeros((size(r),size(r)))      

for i in r_number:
    ri=r[i,0]
    ai=arange(dr,sqrt(R**2-ri**2)+dr/10,dr)
    a[i,:len(ai)]=a[i,:len(ai)]+ai

a[a<dr]=nan
dv=(pi*a**2-pi*(a-da)**2)*dr

#=======
Qr_other=0.216*(V*10**18)**0.939*10**(-12)/12      #(molC cell-1) carbon per cell for other types of phytoplankton (unit conversion 194-31)
#=======

#============================================
# Obtaining Chl_v  refering to a702_06_04
#============================================
Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll
V_syn=4/3*pi*1.63e-6**3
Q=10**(-12)/12/V_syn*V  #(molC/cell) biomass C per cell (164-18)(from Healey 1985)

rho=Qr_other/V     #(molC m-3) biomass C per volume in the cell 
Chlmax=0.048*1000    #(ug chlorophyll a mg C-1) (see 157-36 for conversion)
ChlC=Chlmax/3    #(ug chlorophyll a mg C-1) approximate chlorohpyll (see 157-36 for conversion)
ChlC0=ChlC/(1/12/1000*Mchl/55*10**6)      #(mol C chl mol C -1) Chlorophyll to carbon ratio (see 157-36 for conversion)
Chlv=ChlC0*rho          #(mol C chl m-3) Chlorophyl concentration per volume cell (value check done with a702_06_04)

#=================================
# Photosynthesis parameter
#=================================
Pmax0=7*1.9*1.27                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990)
Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion)
O=0.00863
T=1
I0=200    #(uEmol m-2 s-1)
kchlv=2450*2   #( (chlCmol/m3)-1 m-1) extinguishing coefficient by chlorophyll
kphyto0=kchlv   #(chlmol/m3)-1 m-1) chlorophyll normalized light absorption coefficient by phytoplankton (look at a520_21) 
kphyto=kphyto0*55/1000/Mchl   #(m2 (mg chl a)-1) chlorophyll normalized light absorption coefficient by phytoplankton (unite same as Hickman 2010) (from a520_21)
kchl=kchlv*Chlv         #(m-1) extinguishing coefficient by chlorophyll
k=kchl          #(m-1) total extinguishing coefficient by chlorophyll
Zup=sqrt(R**2-a**2)-r       #(m) path length of the upper hemisphere    (192-7)
Zbottom=sqrt(R**2-a**2)+r   #(m) path length off the bottom hemisphere    (192-9)
Iup=I0*exp((-1)*k*Zup)      #(uEmol m-2 s-1) Light intensity distribution of the upper hemisphere    (192-6)
Ibottom=I0*exp((-1)*k*Zbottom)      #(uEmol m-2 s-1) light intensity distribution of the bottom hemisphere (192-6)
Pup=dv*Chlv*Pmax*(1-exp(-O*T*Iup))      #(mol C s-1) carbon fixation rate off the upper hemisphere (192-6)
Pbottom=dv*Chlv*Pmax*(1-exp(-O*T*Ibottom))      #(mol C s-1) carbon fixtion rate of the bottom hemisphere (192-6) 

#================================================================
#For plotting P-I curve
#================================================================
Iplot=arange(0,500)             #(uEmol m-2 s-1) light intensity for P-I curve plot
Pplot=Pmax0*(1-exp(-O*T*Iplot))      #(g C g Chl-1 h-1) carbon fixation rate per chlorophyll 

#========================================================================================
#This part creates circular views of the distribution of the local carbon fixation rate
#========================================================================================
Pup_right=copy(Pup)
Pup_left=copy(Pup[:,::-1])
Pbottom_right=copy(Pbottom)
Pbottom_left=copy(Pbottom[:,::-1])
P_circle=vstack((hstack((Pup_left,Pup_right)),hstack((Pbottom_left,Pbottom_right))))

#======================================================================================
#This part creates circular views of the distribution of the local light intensity
#======================================================================================
Iup_right=copy(Iup)
Iup_left=copy(Iup[:,::-1])
Ibottom_right=copy(Ibottom[::-1,:])
Ibottom_left=copy(Ibottom[::-1,::-1])
I_circle=vstack((hstack((Iup_left,Iup_right)),hstack((Ibottom_left,Ibottom_right))))
I_circle_masked=ma.masked_where(isnan(I_circle),I_circle)

#====================================================================================
#This part creates the x and y location vectors
#====================================================================================
padsize_number=arange(0,padsize,1)
padsize_number_reverse=padsize_number[::-1]

xaxis_right=a[size(r)-1,:]

for i in padsize_number:
    xaxis_right[-padsize+i]=xaxis_right[-padsize-1+i]+da    #convert nan to appropriate values

xaxis_left=-1*xaxis_right[::-1]
xaxis=hstack((xaxis_left,xaxis_right))
xaxis_um=xaxis*10**6

yaxis_up=copy(r)
for i in padsize_number_reverse:
    yaxis_up[i]=yaxis_up[i+1]+dr            #convert nan to appropriate values

yaxis_bottom=-1*copy(yaxis_up[::-1])
yaxis0=vstack((yaxis_up,yaxis_bottom))
yaxis=yaxis0[:,0]
yaxis_um=yaxis*10**6
#====================================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#5.Plot  (from 605_06_01_12)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

Dpi = 450
#===================================================================
figure(0)
pcolormesh(xaxis_um,yaxis_um,I_circle_masked,vmin=0,vmax=I0)
title('Cell V = '+ str(int(V*10**18))+' ($\mu$m$^{3}$)',y=1.015)
ylabel('($\mu$m)')    #copied from 73
xlabel('($\mu$m)')
cbar = colorbar()
cbar.ax.set_title('Light',y=1.03,fontsize='23')
legend(loc='upper left', borderaxespad=0.6, fontsize=25)
Edge = 3.5
xlim(-Edge,Edge)
ylim(-Edge,Edge)
sf('Light attenuation',int(V*10**18),300)
#===================================================================

show()
