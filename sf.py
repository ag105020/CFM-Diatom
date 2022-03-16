'''
Created on May 18, 2014
This one reads dpi
@author: Keisuke
'''
from pylab import * 

def sf(savefolder,figName,Dpi):
    First_part="C:\\Users\\Keiin\\OneDrive\\Desktop\\figures\\01\\Diatoms\\"
    Second_part=savefolder+"\\"
    Figure_name=str(figName)
    Last_part=".png"
    savefig(First_part+Second_part+Figure_name+Last_part,dpi=Dpi)
