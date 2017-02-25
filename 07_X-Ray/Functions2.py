import Functions1 as fn1
import numpy as np
import csv

###############################################################################

# CONSTANTS

path = "C:\\Users\\Jacob\\Google Drive\\School\\Physics Degree\\Current Classes\\PHY 334 Advanced Lat 1\\X-Ray\\Excel Files\\Part2\\"

d = 2.8201E-10 # m, d_{2,0,0} of NaCl

I = 1 # mA, source current
dt = 1 # s, dwell time
dB = 1, # degree, degree step for goniometer
B_min,B_max = 2.5,30 # degrees, range of goniometer
W74_alpha,W74_beta = 59.31,67.23 # keV, Tungsten
Mo_alpha,Mo_beta = 17.48,19.61

h = 4.135667662E-18 # keV s, planck's constant
c =  299792458 # m, speed of light in vacuum

###############################################################################




###############################################################################

# CLASSES

###############################################################################

class class1:
    def __init__(self,index,CSV,legend,U,source,cap,plot):
        self.index = index
        self.CSV = CSV
        self.legend = legend
        self.U = U
        self.source = source
        self.cap = cap
        self.plot = plot # [ subplot , legend , color ]
    #    index,CSV,         legend,        U,  source,cap,  plot
I1 =  class1(0,'Mo.csv',    [1,'Mo'],      35, 'Mo', 'no',  [221,'r'])
I2 =  class1(1,'MoZr.csv',  [1,'Mo - Zr'], 35, 'Mo', 'yes', [222,'orange'])
I3 = class1(2,'W74.csv',    [1,'W74'],     35, 'W74','no',  [223,'g'])
I4 = class1(3,'W74Zr.csv',  [1,'W74 - Zr'],35, 'W74','yes', [224,'b'])
II1 =  class1(4,'Mo.csv',   [1,'U = 35'],  35, 'Mo', 'no',  [111,'r'])
II2 = class1(5,'MoU30.csv', [1,'U = 30'],  30, 'Mo', 'no',  [111,'b'])
#II2 = class1(5,'MoU25.csv',[1,'U = 25'],  25, 'Mo', 'no',  [111,'g'])
#II3 = class1(6,'MoU20.csv',[1,'U = 20'],  20, 'Mo', 'no',  [111,'b'])
#parts = [I1,I2,I3,I4,II1,II2,II3]
parts = [I1,I2,I3,I4]


###############################################################################

# MISC FUNCTIONS

###############################################################################

def OpenCSV(TITLE,LL,UL): # Imports Table from csv file in python file folder
    # LL,UL give the interval of the table rows to import
    DATA = []
    with open(TITLE, newline='', encoding='utf-8') as d:
        reader = csv.reader(d)
        for row in reader:
            DATA.append(row)
    return DATA[LL:UL]
        
def XY(data,x_i,y_i):
    Xd,Xr,Y = [],[],[]
    for row in data:
        Xd.append(float(row[x_i])) # degrees
        Xr.append(float(row[x_i]) * (np.pi/180)) # radians
        Y.append(float(row[y_i])) # units?
    Xd = np.array(Xd).astype(float)
    Xr = np.array(Xr).astype(float)
    Y = np.array(Y).astype(int)
    return Xd,Xr,Y
    



###############################################################################

# LAB FUNCTIONS

###############################################################################

def Bragg(n,E):
    rad = np.arcsin( (n*h*c)/(2*d*E) )
    theta = rad * (180/np.pi)
    return theta