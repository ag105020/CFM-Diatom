'''
Created on Nov 23, 2015
Note 
about s003

In this part, we will develop a model for growth rate predction 
based on the carbon and silica content in the cell.

Introducing 3 dimension

@author: Keisuke
'''

from pylab import *
from Steph_data_01 import *
from sf001_01_Energy_calculation import *
from sf2 import *

rcParams.update({'mathtext.default': 'regular' })

#=============================================
#for creating 3rd dimension cell volume array
#=============================================

E3=evalue()
E=E3.E
print(E)

V_power_max=8
V_power_min=1
dv_power=0.1
V_power0=arange(V_power_min,V_power_max+dv_power/10,dv_power)
V_power=V_power0.reshape(size(V_power0),1,1)
V=10**V_power*(10**(-6))**3      #(m3) Volume of the cell
#=======
Qr_other=0.216*(V*10**18)**0.939*10**(-12)/12      #(molC cell-1) carbon per cell for other types of phytoplankton (unit conversion 194-31)
#=======

#=======
R=(3*V/(4*pi))**(1/3)           #(m) Cell radius
rsize=100          #size of dr vector
dr=R/rsize        #(m) r step
asize=rsize   #size of da vector
da=copy(dr)     #(m) a step

R_number=arange(0,size(R),1)
r=zeros((size(R),rsize,1))


for i in R_number:
    
    r0=arange(dr[i],R[i]+dr[i]/10,dr[i])
    r1=r0.reshape(size(r0),1)       #converting a row into a column
    r2=r1[::-1]              #reverse the column
    r[i]=r2                 #inserting r values

r_number=arange(0,rsize,1)
#a=zeros((size(r),size(r)))
a=zeros((size(R),rsize,asize))

for i in R_number:
    for j in r_number:
        rj=r[i,j,0]
        aj=arange(dr[i],sqrt(R[i]**2-rj**2)+dr[i]/10,dr[i])
        a[i,j,:len(aj)]=a[i,j,:len(aj)]+aj

a[a<dr]=nan
dv=(pi*a**2-pi*(a-da)**2)*dr

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
#check this calculation kphyto#####
kphyto=kphyto0*55/1000/Mchl   #(m2 (mg chl a)-1) chlorophyll normalized light absorption coefficient by phytoplankton (unite same as Hickman 2010) (from a520_21)
print('kphyto',kphyto)
kchl=kchlv*Chlv         #(m-1) extinguishing coefficient by chlorophyll
k=kchl          #(m-1) total extinguishing coefficient by chlorophyll

Zup=sqrt(R**2-a**2)-r       #(m) path length for the upper half (192-7)
Zbottom=sqrt(R**2-a**2)+r   #(m) path length for the upper bottom half of the sphere (192-7)
Iup=I0*exp((-1)*k*Zup)      #(uEmol m-2 s-1) light intensity for the upper half
Ibottom=I0*exp((-1)*k*Zbottom)  #(uEmol m-2 s-1) light intensity for the bottom half

Pup=dv*Chlv*Pmax*(1-exp(-O*T*Iup))  
Pbottom=dv*Chlv*Pmax*(1-exp(-O*T*Ibottom))

Ptot=zeros((size(R),1,1))       #(molC cell-1 s-1) Photosynthesis rate per cell
for i in R_number:
    Ptot[i]=nansum(Pup[i])+nansum(Pbottom[i])

#E=0.7

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
#Carbon size and silica consideration
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
ATP=50      #(kJ molATP-1) 

#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

Ysi_c=0.163         #(molSi molC-1) Si to C ratio in the cell (194-33) (
Yres_si_c=1/6       #(molC molSi) Carbonconsumption to Si uptake ratio (194-33) (Werner 1977)

Qr_diatom=0.288*(V*10**18)**0.811*10**(-12)/12     #(molC cell-1) carbon per cell for diatom (unit conversion 194-31)
Qr_other=0.216*(V*10**18)**0.939*10**(-12)/12      #(molC cell-1) carbon per cell for other types of phytoplankton (unit conversion 194-31)

ChlC0_diatom=ChlC0*Qr_other/Qr_diatom   #(mol C chl mol C -1) Chlorophyll to carbon ratio for diatom
ChlC0_plot=ChlC0/V*V                #(mol C chl mol C -1) Chlorophyll to carbon ratio for other phytoplankton for plot

Mu=Ptot/(Q*(1+E))     #(s-1) growth rate
Mu_diatom=Ptot/(Qr_diatom*(1+E+Ysi_c*Yres_si_c))    #(s-1) growth rate for diatom (194-33)
Mu_other=Ptot/(Qr_other*(1+E))      #(s-1) growth rate of other phytoplankton (194-33)
print('Si',Ysi_c*Yres_si_c)
Muday=Mu*86400          #(d-1) growth rate per day
Muday_diatom=Mu_diatom*86400     #(d-1) growth rate per day for diatom (194-33)
Muday_other=Mu_other*86400      #(d-1) growth rate per day for other phytoplankton (194-33)

#=====================================================
#For Cs (carbohydrate consumption rate) stack plot
#=====================================================
Csbio_diatom=Mu_diatom*Qr_diatom/V                #(molC s-1 m-3) Carbohydrate consumption for biomass for diatom (194-34)
Csresp_bio_diatom=Mu_diatom*Qr_diatom*E/V              #(molC s-1 m-3) Carbohydrate consumption for respiraiton supporting biomass production for diatom (194-34)
Csresp_Si_diatom=Mu_diatom*Qr_diatom*Ysi_c*Yres_si_c/V     #(molC s-1 m-3) Carbohydrate coonsumption for respration supporting silica uptake (194-34) * uptake through the cell membrane and uptake into the silica deposition vessicle

Csbio_other=Mu_other*Qr_other/V                #(molC s-1 m-3) Carbohydrate consumption for biomass for other phytoplankton (194-34)
Csresp_bio_other=Mu_other*Qr_other*E/V              #(molC s-1 m-3) Carbohydrate consumption for respiraiton supporting biomass production for other phytoplankton (194-34)

#============================
#for calcifier
#============================
Calcifier_CO2_limiting_factor=0.55
Muday_calcifier=Muday*Calcifier_CO2_limiting_factor

#==========================================
#for dino : ***Overwriting some values***
#==========================================
Dino_chl_factor=0.33   #Chlorophyll concentration relative to diatom
Chlv_dino=Chlv*Dino_chl_factor   #(mol C chl m-3) Chlorophyl concentration per volume cell for dino (value check done with a702_06_04)
kchl=kchlv*Chlv_dino         #(m-1) extinguishing coefficient by chlorophyll
k=kchl          #(m-1) total extinguishing coefficient by chlorophyll
Zup=sqrt(R**2-a**2)-r       #(m) path length for the upper half (192-7)
Zbottom=sqrt(R**2-a**2)+r   #(m) path length for the upper bottom half of the sphere (192-7)
Iup=I0*exp((-1)*k*Zup)      #(uEmol m-2 s-1) light intensity for the upper half
Ibottom=I0*exp((-1)*k*Zbottom)  #(uEmol m-2 s-1) light intensity for the bottom half
Pup=dv*Chlv_dino*Pmax*(1-exp(-O*T*Iup))  
Pbottom=dv*Chlv_dino*Pmax*(1-exp(-O*T*Ibottom))
Ptot=zeros((size(R),1,1))
for i in R_number:
    Ptot[i]=nansum(Pup[i])+nansum(Pbottom[i])
Mu=Ptot/(Q*(1+E))     #(-s) growth rate
Muday_dino=Mu*86400          #(-d) growth rate per day

#================================================================
#For plotting P-I curve
#================================================================
Iplot=arange(0,500)             #(uEmol m-2 s-1) light intensity for P-I curve plot
Pplot=Pmax0*(1-exp(-O*T*Iplot))      #(g C g Chl-1 h-1) carbon fixation rate per chlorophyll 
#================================================================

#==================================
#test with non self shading case
#==================================
Iex=50
PIchl=Pmax*(1-exp(-O*Iex*T)) #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
Chlcell=Chlv*V          #(mol C chl cell-1) Chlorophyl concentration 
PIex=PIchl*Chlcell      #(mol C s-1 cell-1 s-1) photosynthesis rate
Muex=PIex/(Q*(1+E))     #(s-1) growth rate
Muexday=Muex*86400      #(d-1) growth rate


#========================================================
#Getting experimental data
#========================================================

data_diatom = genfromtxt('Data_diatom.csv', delimiter=',').T
data_others = genfromtxt('Data_others.csv', delimiter=',').T

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#5.Plot  (from 605_06_01_12)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rcParams.update({'font.size': 20,
                 'lines.markersize':10,
                 'lines.markeredgewidth':1.5,
                 'lines.linewidth':2.5})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'font.serif': 'Times New Roman'})
rcParams.update({'figure.autolayout': True})
rcParams.update({'patch.edgecolor':'none'})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
lineW = 1.5
rcParams.update({'axes.linewidth':lineW}) 
rcParams.update({'xtick.major.width':lineW})  #
rcParams.update({'xtick.minor.width':lineW})  #
rcParams.update({'ytick.major.width':lineW})  #
rcParams.update({'ytick.minor.width':lineW})  #
rcParams.update({'xtick.major.size':6})  #
rcParams.update({'xtick.minor.size':3})  #
rcParams.update({'ytick.major.size':6})  #
rcParams.update({'ytick.minor.size':3})  #
rcParams['xtick.major.pad']='8'
rcParams.update({'lines.markeredgewidth': 1}) #default = 1

FolderPath = "C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\01\\Size-growth rate\\"
#===================================================================
figure('P-I curve')
plot(Iplot,Pplot)
title('P-I curve',y=1.015)
ylabel('($\mu$g C $\mu$g Chl$^{-1}$ h$^{-1}$)')    #copied from 73
xlabel('($\mu$Emol m$^{-2}$ s$^{-1}$)')

Dpi = 450
#===================================================================
#===================================================================
figure(3)
plot(10**data_diatom[0],data_diatom[1],'o',markersize=10,color='none',markeredgecolor='red',markeredgewidth=2.5)
plot(10**data_others[0],data_others[1],'o',markersize=10,color='none',markeredgecolor='cyan',markeredgewidth=2.5)

plot(V[:,0,0]*10**18,Muday_diatom[:,0,0], linewidth=3.5, color='red', label='Diatom')
plot(V[:,0,0]*10**18,Muday_other[:,0,0], linewidth=3.5, color='cyan', label='Others')

xlabel('Cell V ($\mu$m$^{3}$)')    #copied from 73
ylabel('$\mu$ (d$^{-1}$)')         #copied from 73
xscale('log')
#xlim(left=10)
legend(loc='upper right', borderaxespad=0.5, fontsize=25)
sf('','Fig3',Dpi) 

figure(4)
plot(V[:,0,0]*10**18,Qr_diatom[:,0,0]*1e12*12,color='green',label='Diatom')
plot(V[:,0,0]*10**18,Qr_other[:,0,0]*1e12*12,color='orange',label='Other')
xscale('log')
yscale('log')
xlabel('Cell V ($\mu$m$^{3}$)')
ylabel('Cell C (pg cell$^{-1}$)')
legend(loc=4,fontsize=25)
sf('','C-V',Dpi)

xlimmin=0
xlimmax=1.2000001
rcParams['xtick.major.pad']='10'

UC = 86400*1e-18 #Unit conversion from s-1 m-3 to d-1 um-3)
figure(7)
Names = ["Biosynthesis","Silica deposition"]
Colors = ['#007BFF','red']
for i in arange(size(Names)):
    plot([],[],linewidth=10,color=Colors[-i-1],label=Names[-i-1])
stackplot(Muday_diatom[:,0,0],(Csbio_diatom[:,0,0]+Csresp_bio_diatom[:,0,0])*UC,Csresp_Si_diatom[:,0,0]*UC,colors=Colors)
print(Csresp_Si_diatom[:,0,0]/(Csbio_diatom[:,0,0]+Csresp_bio_diatom[:,0,0]))

# xlim(left=0.4)
# xlim(xmax=xlimmax)
#ylim(ymax=0.36*UC)
#ylim(top=0.3500001)
xlim(0,2)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$C_s$ (molC d$^{-1}$ $\mu$m$^{-3}$)')
title('Diatom',y=1.015)
legend(loc=2,fontsize=25)
sf('','Fig4a',Dpi)

figure(71)
Names = ["Biosynthesis","Si deposition"]
Colors = ['#007BFF','red']
for i in arange(size(Names)):
    plot([],[],linewidth=10,color=Colors[-i-1],label=Names[-i-1])
stackplot(Muday_diatom[:,0,0],Csresp_Si_diatom[:,0,0]*UC,(Csbio_diatom[:,0,0]+Csresp_bio_diatom[:,0,0])*UC,colors=Colors)
print(Csresp_Si_diatom[:,0,0]/(Csbio_diatom[:,0,0]+Csresp_bio_diatom[:,0,0]))

# xlim(left=0.4)
# xlim(xmax=xlimmax)
# ylim(ymax=0.36*UC)
xlim(0,2)
#ylim(top=0.3500001)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$C_s$ (molC d$^{-1}$ $\mu$m$^{-3}$)')
title('Diatom',y=1.015)
legend(loc=2,fontsize=25)
sf('','Fig4a_flip',Dpi)




figure(8)
stackplot(Muday_other[:,0,0],(Csbio_other[:,0,0]+Csresp_bio_other[:,0,0])*UC,colors=Colors)
# xlim(xmin=xlimmin)
# xlim(xmax=xlimmax)
# ylim(top=0.3600001*UC)
xlim(0,2)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$C_s$ (molC d$^{-1}$ $\mu$m$^{-3}$)')
title('Other',y=1.015)
sf('','Fig4b',Dpi)

figure(9)
plot(V[:,0,0]*10**18,Ptot[:,0,0]*86400)
xlabel('Cell volume ($\mu$m$^{3}$)')
ylabel('F$_{Pho}$ (molC d$^{-1}$ cell$^{-1}$)')
xscale('log')
sf('','Ph_log',Dpi)

figure(10)
plot(V[:,0,0]*10**18,Ptot[:,0,0]*86400)
xlabel('Cell volume ($\mu$m$^{3}$)')
ylabel('F$_{Pho}$ (molC d$^{-1}$ cell$^{-1}$)')
sf('','Ph',Dpi)

figure(11)
plot(V[:,0,0]*10**18,Ptot[:,0,0]*86400)
xlabel('Cell V ($\mu$m$^{3}$)')
ylabel('F$_{Pho}$ (molC d$^{-1}$ cell$^{-1}$)')
xscale('log')
yscale('log')
sf('','Ph_log_log',Dpi)

V0 = log10(V[:,0,0]*10**18)
Pho = log10(Ptot[:,0,0]*86400)

print((Pho[-1]-Pho[-2])/(V0[-1]-V0[-2]))

figure(12)
plot(V[:,0,0]*10**18,Qr_diatom[:,0,0])
xlabel('Cell volume ($\mu$m$^{3}$)')
ylabel('Q$_{C}$ (molC cell$^{-1}$)')
xscale('log')
sf('','Qc_log',Dpi)

figure(13)
plot(V[:,0,0]*10**18,Qr_diatom[:,0,0])
xlabel('Cell volume ($\mu$m$^{3}$)')
ylabel('Q$_{C}$ (molC cell$^{-1}$)')
xscale('log')
yscale('log')
sf('','Qc_log_log',Dpi)

figure(14)
plot(V[:,0,0]*10**18,Qr_diatom[:,0,0])
xlabel('Cell volume ($\mu$m$^{3}$)')
ylabel('Q$_{C}$ (molC cell$^{-1}$)')
sf('','Qc',Dpi)

figure(15)
plot(V[:,0,0]*10**18,Ptot[:,0,0]*86400/(V[:,0,0]*10**18))
xlabel('Cell V ($\mu$m$^{3}$)')
ylabel('F$_{Pho}$/V (molC d$^{-1}$ $\mu$m$^{-3}$)')
xscale('log')
yscale('log')
sf('','Ph_V_log_log',Dpi)

show()

#======================================================
rcParams.update({'font.size': 35,
                 'lines.markersize':10,
                 'lines.markeredgewidth':1.5,
                 'lines.linewidth':2.5})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'xtick.major.pad': 15})
rcParams.update({'font.serif': 'Times New Roman'})
rcParams.update({'figure.autolayout': True})
rcParams.update({'patch.edgecolor':'none'})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'None'})
lineW = 3
rcParams.update({'axes.linewidth':lineW}) 
rcParams.update({'xtick.major.width':lineW})  #
rcParams.update({'xtick.minor.width':lineW})  #
rcParams.update({'ytick.major.width':lineW})  #
rcParams.update({'ytick.minor.width':lineW})  #
rcParams.update({'xtick.major.size':6})  #
rcParams.update({'xtick.minor.size':3})  #
rcParams.update({'ytick.major.size':6})  #
rcParams.update({'ytick.minor.size':3})  #
rcParams['xtick.major.pad']='8'
rcParams.update({'lines.markeredgewidth': 1}) #default = 1

figure(77)
Names = ["Biosynthesis","Si deposition"]
Colors = ['#007BFF','red']
for i in arange(size(Names)):
    plot([],[],linewidth=10,color=Colors[-i-1],label=Names[-i-1])
stackplot(Muday_diatom[:,0,0],Muday_diatom[:,0,0]*0,Csresp_Si_diatom[:,0,0]*UC,colors=Colors)
xlabel('$\mu$ (d$^{-1}$)')
ylabel('$C_s$ (molC d$^{-1}$ $\mu$m$^{-3}$)')
xlim(1e-10,2)
ylim(0,10e-16)
sf('','Inset',Dpi)
#print(Csresp_Si_diatom)

show()

