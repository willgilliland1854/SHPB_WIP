import numpy as np
import math
import pandas as pd
import os
import matplotlib.pyplot as plt 

print('')
print('Data Analysis...INITIATED')

#%% User-Defined Functions
def smooth(x,N):
    return np.convolve(x, np.ones(N), 'valid')/N

def intersectCheck(x1,y1,x2,y2):
    """
    intersectCheck takes array inputs and checks of two line segments intersect
    if line segments intersect returns 1, if they do not intersect returns 0"""
    u1 = np.array([(x1[1]-x1[0]), (y1[1]-y1[0])])
    v1 = np.array([(x2[0]-x1[0]), (y2[0]-y1[0])])
    w1 = np.array([(x2[1]-x1[0]), (y2[1]-y1[0])])
    u2 = np.array([(x2[1]-x2[0]), (y2[1]-y2[0])])
    v2 = np.array([(x1[0]-x2[0]), (y1[0]-y2[0])])
    w2 = np.array([(x1[1]-x2[0]), (y1[1]-y2[0])])
    if (np.cross(u1,v1)*np.cross(u1,w1) < 0)&(np.cross(u2,v2)*np.cross(u2,w2) < 0):
        intersect = 1
    else:   
        intersect = 0
    return intersect

def intersections2(x1,y1,x2,y2):
    check = intersectCheck(x1,y1,x2,y2)
    if check == 1:
        m1 = (y1[1] - y1[0])/(x1[1]-x1[0])
        m2 = (y2[1] - y2[0])/(x2[1]-x2[0])
        b1 = y1[0] - m1*x1[0]
        b2 = y2[0] - m1*x2[0]
        x_int = np.array([(b2 - b1)/(m1 - m2)])
        y_int = np.array([(m1*x_int + b1)])
        return x_int, y_int

def find_first(vec, item):
    """return the index of the first occurence of item in vec"""
    for i in range(len(vec)):
        if vec[i] >= item:
            return i
    return -1
    

#%% Code to automatically check for data sets that haven't been run


#%% 1. Read in Experiment Variables and Raw Data
# filename = input('Input the epxeriment log filename (Example: 20180423_None_1_10psi_8in_OTK)     ')

filename = '20180421_Calib_05_12psi_8in_OTK'

#os.chdir('/Volumes/kingstedtlab/SHPB - Data and Analysis/Raw Data and Experiment Log Forms')
os.chdir('K:/SHPB - Data and Analysis/Calibration Data') # Change directory

bar_mat1    = 1               # Bar Material
D_bar       = 0.75 # Bar diameter (m)
R_bar       = D_bar/2                                                                                                                           # Bar radius (m)
L_spec      = pd.read_excel(''.join([filename,'_LOG.xlsx']),'Sheet1',usecols = 'G',skiprows = 19,nrows = 1,header = None).values * 10**-3       # Specimen length (m)
D_spec      = pd.read_excel(''.join([filename,'_LOG.xlsx']),'Sheet1',usecols = 'O',skiprows = 19,nrows = 1,header = None).values * 10**-3       # Specimen diameter (m)
L_striker   = 8       # Striker bar length (m) 
gg_psi      = 12               # Gas gun pressure psi
V_ex        = 12               # Excitation Voltage
Gain = 5.0 / (20.833* V_ex/1000)   # Gain
test        = pd.read_excel(''.join([filename,'_LOG.xlsx']),'Sheet1',usecols = 'M',skiprows = 35,nrows = 1,header = None).values
ps          = pd.read_excel(''.join([filename,'_LOG.xlsx']),'Sheet1',usecols = 'E',skiprows = 25,nrows = 1,header = None).values                # Pulse shaping used? Yes-1, No-0

# Read in Raw Data
Data        = pd.read_excel(''.join([filename,'_DATA.xlsx']),usecols = 'A:E',skiprows = 16,skipfoot = 125017,header = None).values


#%% 2. Fixed SHPB Bar Geometry Values
if bar_mat1 == 1 and D_bar >= 0.019:                    # Aluminum Bars 3/4"
    c_0         = np.arange(4800, 5101)                 # Calibrated bar wave speed
    L_inc       = 98.0*0.0254                           # Incident bar length (m)
    dx_inc      = 48.75 * 0.025                         # Distance form strain gauges to specimen interface (m)
    L_trans     = 78.0 * 0.0254                         # Transmited bar length (m)
    dx_trans    = (78.0 - (38.0 + 15.0/16.0) ) * 0.0254 # Distance from specimen interface to second set of strain gauges
    v_bar       = 0.33                                  # Poisson's Ratio of the bar
    E_bar       = 71.7 *10**9                           # Elastic Modulus (GPa)
    GF          = 2.125                                 # Strain Gauge Gauge Factor
    A_bar       = math.pi*(D_bar/2)**2                  # Bar Cross-section (m^2)
    bridge      = 'full'                                # Type of Wheatstone Bridge Used
    
    # Allocate imported array to column variable names
    SHPB_time_raw       = Data[:,0]
    SHPB_inc_raw        = Data[:,1]
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw    = -SHPB_inc_raw
    SHPB_trans_raw      = Data[:,3]
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]
    
elif bar_mat1 == 1 and D_bar < 0.019:    # Aluminum Bars 1/2"
    c_0         = 5050
    L_inc       = 1
    dx_inc      = 1
    L_trans     = 1
    dx_trans    = 1
    v_bar       = 0.33
    GF          = 1
    bridge      = 'quarter'
    
    # Allocate imported array to column variable names
    SHPB_time_raw = Data[:,0]
    SHPB_inc_raw1 = Data[:,1]   # sign is negative so the compressive pulse is positive
    SHPB_inc_raw2 = Data[:,3]
    SHPB_trans_raw1 = Data[:,5]
    SHPB_trans_raw2 = Data[:,7]
    
    SHPB_inc_raw = np.add(SHPB_inc_raw1 + SHPB_inc_raw2)/2
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw   = -SHPB_inc_raw
    SHPB_trans_raw = np.add(SHPB_trans_raw1,SHPB_trans_raw2)/2
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]

elif bar_mat1 == 2 and D_bar >= 0.019:      # Steel Bars 3/4"
    c_0         = np.arange(4800, 5101)
    L_inc       = 96.0*0.0254
    dx_inc      = 47.91666667*0.0254
    L_trans     = 96.0*0.0254
    dx_trans    = 48.875*0.0254
    v_bar       = 0.3
    E_bar       = 190.0*10**9               # Value retrieved from Makeitfrom.com
    GF          = 2.130
    A_bar       = math.pi*(D_bar/2)**2      # Bar Cross-section (m^2)
    
    # Allocate imported array to column variable names
    SHPB_time_raw       = Data[:,0]
    SHPB_inc_raw        = Data[:,1]
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw    = -SHPB_inc_raw
    SHPB_trans_raw      = Data[:,3]
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw  = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]
    
elif bar_mat1 == 2 and D_bar < 0.019:    # Steel bars 1/2"
    c_0         = 5050
    L_inc       = 1
    dx_inc      = 1
    L_trans     = 1
    dx_trans    = 1
    v_bar       = 0.33
    GF          = 1
    bridge      = 'quarter'
    
    # Allocate imported array to column variable names
    SHPB_time_raw = Data[:,0]
    SHPB_inc_raw1 = Data[:,1]
    SHPB_inc_raw2 = Data[:,3]
    SHPB_trans_raw1 = Data[:,5]
    SHPB_trans_raw2 = Data[:,7]
    
    SHPB_inc_raw = np.add(SHPB_inc_raw1 + SHPB_inc_raw2)/2
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw   = -SHPB_inc_raw
    SHPB_trans_raw = np.add(SHPB_trans_raw1,SHPB_trans_raw2)/2
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]
#%%

y_SHPB = smooth(SHPB_inc_raw, 100)
x_SHPB = smooth(SHPB_time_raw, 100)

# if ps == 1: # a pule shaped signal - longer duration
#     x_SHPB_1 = SHPB_time_0[round(0.03 * len(SHPB_time_0)):round(0.5 * len(SHPB_time_0))]
#     x_SHPB = x_SHPB_1[49:int(len(x_SHPB_1)-50)]
#     #x_SHPB   = x_SHPB_2[0:58651]
#     y_SHPB_1 = SHPB_inc_raw[round(0.03 * len(SHPB_time_0)):round(0.5 * len(SHPB_time_0))]
#     y_SHPB = y_SHPB_1[49:int(len(y_SHPB_1)-50)]
#     #y_SHPB   = y_SHPB_2[0:58651]
# else: 
#     x_SHPB_1 = SHPB_time_0[round(0.03 * len(SHPB_time_0)):round(0.35 * len(SHPB_time_0))]
#     #x_SHPB_2 = x_SHPB_1[51:,]
#     x_SHPB = x_SHPB_1[49:int(len(x_SHPB_1)-50)]
#     #x_SHPB   = x_SHPB_2[0:58651]
#     y_SHPB_1 = SHPB_inc_raw[round(0.03 * len(SHPB_time_0)):round(0.35 * len(SHPB_time_0))]
#     #y_SHPB_2 = y_SHPB_1[51:,]
#     y_SHPB = y_SHPB_1[49:int(len(y_SHPB_1)-50)]
#     #y_SHPB   = y_SHPB_2[0:58651]
    
# y_s = smooth(y_SHPB_1,100)

#%% Normalize the SHPB data traces prior to taking the Hough Transform
x_norm = 2*(x_SHPB - min(x_SHPB)) / (max(x_SHPB - min(x_SHPB))) - 1
y_norm = 2*(y_SHPB - min(y_SHPB)) / (max(y_SHPB - min(y_SHPB))) - 1
# y_s_norm = 2*(y_s - min(y_s)) / (max(y_s - min(y_s))) - 1

plt.figure(figsize = [6.4, 4.8])
plt.plot(x_norm, y_norm)
plt.title('Normalized Incident Bar Signal')
plt.xlabel('Time [sec]'), plt.ylabel('Incident Bar Voltage [V]')
plt.grid()

down_sampling = round(len(x_norm)/300.0)
# Downsample x_norm and y_norm
# x_norm_ds = x_norm[::down_sampling]
# y_norm_ds = y_s_norm[::down_sampling]

# #%%
# plt.plot(SHPB_time_raw, SHPB_inc_raw)


#%%
filename = '20180918_IN718_BD_C1_23_50psi_8in_NK'

#os.chdir('/Volumes/kingstedtlab/SHPB - Data and Analysis/Raw Data and Experiment Log Forms')
os.chdir('K:/SHPB - Data and Analysis/Raw Data and Experiment Log Forms') # Change directory
#%%

bar_mat1    = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'C',skiprows = 21,nrows = 1,header = None).values                # Bar Material
D_bar       = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'M',skiprows = 21,nrows = 1,header = None).values * 0.0254       # Bar diameter (m)
R_bar       = D_bar/2                                                                                                                           # Bar radius (m)
L_spec      = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'G',skiprows = 19,nrows = 1,header = None).values * 10**-3       # Specimen length (m)
D_spec      = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'O',skiprows = 19,nrows = 1,header = None).values * 10**-3       # Specimen diameter (m)
L_striker   = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'E',skiprows = 31,nrows = 1,header = None).values * 0.0254       # Striker bar length (m) 
gg_psi      = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'N',skiprows = 31,nrows = 1,header = None).values                # Gas gun pressure psi
V_ex        = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'G',skiprows = 35,nrows = 1,header = None).values                # Excitation Voltage
Gain = 5.0 / (pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'M',skiprows = 35,nrows = 1,header = None).values * V_ex/1000)   # Gain
test        = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'M',skiprows = 35,nrows = 1,header = None).values
ps          = pd.read_excel(''.join([filename,'_Log.xlsx']),'Sheet1',usecols = 'E',skiprows = 25,nrows = 1,header = None).values                # Pulse shaping used? Yes-1, No-0

# Read in Raw Data
Data        = pd.read_excel(''.join([filename,'_DATA.xlsx']),usecols = 'A:E',skiprows = 16,skipfoot = 125017,header = None).values


#%% 2. Fixed SHPB Bar Geometry Values
if bar_mat1 == 1 and D_bar >= 0.019:                    # Aluminum Bars 3/4"
    c_0         = np.arange(4800, 5101)                 # Calibrated bar wave speed
    L_inc       = 98.0*0.0254                           # Incident bar length (m)
    dx_inc      = 48.75 * 0.025                         # Distance form strain gauges to specimen interface (m)
    L_trans     = 78.0 * 0.0254                         # Transmited bar length (m)
    dx_trans    = (78.0 - (38.0 + 15.0/16.0) ) * 0.0254 # Distance from specimen interface to second set of strain gauges
    v_bar       = 0.33                                  # Poisson's Ratio of the bar
    E_bar       = 71.7 *10**9                           # Elastic Modulus (GPa)
    GF          = 2.125                                 # Strain Gauge Gauge Factor
    A_bar       = math.pi*(D_bar/2)**2                  # Bar Cross-section (m^2)
    bridge      = 'full'                                # Type of Wheatstone Bridge Used
    
    # Allocate imported array to column variable names
    SHPB_time_raw       = Data[:,0]
    SHPB_inc_raw        = Data[:,1]
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw    = -SHPB_inc_raw
    SHPB_trans_raw      = Data[:,3]
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]
    
elif bar_mat1 == 1 and D_bar < 0.019:    # Aluminum Bars 1/2"
    c_0         = 5050
    L_inc       = 1
    dx_inc      = 1
    L_trans     = 1
    dx_trans    = 1
    v_bar       = 0.33
    GF          = 1
    bridge      = 'quarter'
    
    # Allocate imported array to column variable names
    SHPB_time_raw = Data[:,0]
    SHPB_inc_raw1 = Data[:,1]   # sign is negative so the compressive pulse is positive
    SHPB_inc_raw2 = Data[:,3]
    SHPB_trans_raw1 = Data[:,5]
    SHPB_trans_raw2 = Data[:,7]
    
    SHPB_inc_raw = np.add(SHPB_inc_raw1 + SHPB_inc_raw2)/2
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw   = -SHPB_inc_raw
    SHPB_trans_raw = np.add(SHPB_trans_raw1,SHPB_trans_raw2)/2
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]

elif bar_mat1 == 2 and D_bar >= 0.019:      # Steel Bars 3/4"
    c_0         = np.arange(4800, 5101)
    L_inc       = 96.0*0.0254
    dx_inc      = 47.91666667*0.0254
    L_trans     = 96.0*0.0254
    dx_trans    = 48.875*0.0254
    v_bar       = 0.3
    E_bar       = 190.0*10**9               # Value retrieved from Makeitfrom.com
    GF          = 2.130
    A_bar       = math.pi*(D_bar/2)**2      # Bar Cross-section (m^2)
    
    # Allocate imported array to column variable names
    SHPB_time_raw       = Data[:,0]
    SHPB_inc_raw        = Data[:,1]
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw    = -SHPB_inc_raw
    SHPB_trans_raw      = Data[:,3]
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw  = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]
    
elif bar_mat1 == 2 and D_bar < 0.019:    # Steel bars 1/2"
    c_0         = 5050
    L_inc       = 1
    dx_inc      = 1
    L_trans     = 1
    dx_trans    = 1
    v_bar       = 0.33
    GF          = 1
    bridge      = 'quarter'
    
    # Allocate imported array to column variable names
    SHPB_time_raw = Data[:,0]
    SHPB_inc_raw1 = Data[:,1]
    SHPB_inc_raw2 = Data[:,3]
    SHPB_trans_raw1 = Data[:,5]
    SHPB_trans_raw2 = Data[:,7]
    
    SHPB_inc_raw = np.add(SHPB_inc_raw1 + SHPB_inc_raw2)/2
    if np.mean(SHPB_inc_raw[: int(len(SHPB_inc_raw)/2)]) < 0: # flips pulse sign so compression is positive
        SHPB_inc_raw   = -SHPB_inc_raw
    SHPB_trans_raw = np.add(SHPB_trans_raw1,SHPB_trans_raw2)/2
    if np.mean(SHPB_trans_raw[int(len(SHPB_inc_raw)/2):round(0.75*len(SHPB_inc_raw))]) < 0: # flips pulse sign so compression is positive
        SHPB_trans_raw = -SHPB_trans_raw
        
    # Shift the time array to start at zero
    SHPB_time_0     = SHPB_time_raw[0:] - SHPB_time_raw[0]
    dt              = SHPB_time_0[1]

#%% 3.0 Complete the Hough Transform of the SHPB Incident Wave  
    
y_SHPB2 = smooth(SHPB_inc_raw, 100)
x_SHPB2 = smooth(SHPB_time_raw, 100)

x_norm2 = 2*(x_SHPB2 - min(x_SHPB2)) / (max(x_SHPB2 - min(x_SHPB2))) - 1
y_norm2 = 2*(y_SHPB2 - min(y_SHPB2)) / (max(y_SHPB2 - min(y_SHPB2))) - 1
# y_s_norm = 2*(y_s - min(y_s)) / (max(y_s - min(y_s))) - 1

plt.figure(figsize = [6.4, 4.8])
plt.plot(x_norm, y_norm, color = 'blue', label = 'Calibration Data w/ Pulse Shaper')
plt.plot(x_norm2, y_norm2, color = 'red', label = 'IN718 Data w/o Pulse Shaper')
plt.title('Normalized Incident Bar Signal')
plt.xlabel('Time [sec]'), plt.ylabel('Incident Bar Voltage [V]')
plt.legend(loc = 'best')
plt.grid()
    
    


