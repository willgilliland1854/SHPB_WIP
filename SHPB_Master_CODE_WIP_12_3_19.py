# -*- coding: utf-8 -*-
"""
@author: William Gilliland
Last edited: 12/2/19
"""
"""
 Nate Briggs and Owen Kingstedt

Last Modified: 4/23/2018
Version Number: 1

DESCRIPTION: The following code ingests a SHPB signal record, finds the
beginning and end points of the incident signal, performs the dispersion
correction, and calculates the material stress-strain behavior.

Late determination:
The windowed signals might show late determination. Check the results and
adjust the windowed signals by shifting the transmitted signal as needed
for each case.

CODE OUTLINE:
  0. Check to see if the file has been run
  1. Read in Experiment Variables and Raw Data
  2. Fixed SHPB Bar Geometry Values 
  3. Read in the Raw Experimental Data
  4. Complete the Hough Transform of the SHPB Incident Wave
  5. Calculate stress-strain behavior wihtout dispersion correction
  6. Perform the dispersion correction of the SHPB data
  7. Calculate the dispersion corrected stress-strain behavior
  8. Save the work space
  9. Save plots
      a. Save plot of raw data
      b. Save plot of Hough Space
      c. Save plot of the signal with intersection points
      d. Save plot of stress-strain behavior
      e. Save plot of strain-rate vs. time
      f. Save plot of the force balance
      g. Save plot of the frequency content
      h. Save plot of the corrected stress-strain behavior
      i. Save plot of strain-rate vs. time
      j. Save plot of the force balance
      k. Save plot comparing the corrected solution to the stress-strain
      behavior
"""
import numpy as np
import math
import pandas as pd
import os
import matplotlib as mpl
import matplotlib.pyplot as plt 
from scipy.optimize import fsolve

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

def intersectFinder(x1,y1,x2,y2):
    """
    Finds intersection of lines defined by points (x1[0], y1[0]) to (x1[1], y1[1])
    and (x2[0], y2[0]) and (x2[1], y2[1])"""   
    m1 = (y1[0]-y1[1])/(x1[0]-x1[1])
    m2 = (y2[0]-y2[1])/(x2[0]-x2[1])
    b1 = y1[0] - m1*x1[0]
    b2 = y2[0] - m2*x2[0]
    x_int = (b2-b1)/(m1-m2)
    y_int = m1*x_int + b1   
    return x_int, y_int


def poly_x_poly(x1,y1,x2,y2):
    """
    Uses a sweep-line to determine intersection of two poly lines that are both
    functions of x"""
    hits = 0
    # Sweep line loop
    for i in range (0,len(x1)-1):
    # Find the lower and upper indicies of x2 that are inside [x1[i], x1[i+1]] 
    # This interval defines the area of swept line
        lowerBound = np.min(np.where(x2 >= x1[i]))
        upperBound = np.max(np.where(x2 <= x1[i+1]))
        if lowerBound <= 0:
            interval = [lowerBound, upperBound+1]
        elif upperBound == len(x2)-1:
            interval = [lowerBound-1, upperBound]
        else:
            interval = [lowerBound-1, upperBound+1]       
        # Loop to check for intersections inside the interval defining swept line
        for j in range(interval[0], interval[1]):
            # Check points inside interval for intersections
            ax = np.array([x1[i], x1[i+1]])
            ay = np.array([y1[i], y1[i+1]])
            bx = np.array([x2[j], x2[j+1]])
            by = np.array([y2[j], y2[j+1]])
            int_check = intersectCheck(ax,ay,bx,by)             
            # If lines do intersect, determine point of intersection
            if int_check == 1:
                ints = intersectFinder(ax,ay,bx,by)
                if hits == 0:
                    intersections = np.array(ints)
                else:
                    intersections = np.vstack([intersections, ints])                    
                hits = hits + 1
    return intersections


def find_first(vec, item):
    """return the index of the first occurence of item in vec"""
    for i in range(len(vec)):
        if vec[i] >= item:
            return i
    return -1
    

def getPhiFcn(c0_loop,A,B,C,D,E,F,akomega):
    fcn = lambda ck: (ck/c0_loop-(A+(B/(C*(akomega/ck)**4+D*(akomega/ck)**3+E*(akomega/ck)**2+F*(akomega/ck)**(1.5)+1))))
    return fcn

#%% Code to automatically check for data sets that haven't been run


#%% 1. Read in Experiment Variables and Raw Data
# filename = input('Input the epxeriment log filename (Example: 20180423_None_1_10psi_8in_OTK)     ')

filename = '20190131_IN718_BD_C2_2_HT_720C_24h_37psi_8in_TL'

#os.chdir('/Volumes/kingstedtlab/SHPB - Data and Analysis/Raw Data and Experiment Log Forms')
os.chdir('C:/Kingstedt') # Change directory
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
    
    ###########################################################################
    ######### if something is broken look at the slices!!!!!!!!!!!!! ##########
    ###########################################################################
    
    if ps == 1: # a pule shaped signal - longer duration
        x_SHPB_1 = SHPB_time_0[round(0.03 * len(SHPB_time_0)):round(0.5 * len(SHPB_time_0))]
        x_SHPB = x_SHPB_1[49:int(len(x_SHPB_1)-50)]
        #x_SHPB   = x_SHPB_2[0:58651]
        y_SHPB_1 = SHPB_inc_raw[round(0.03 * len(SHPB_time_0)):round(0.5 * len(SHPB_time_0))]
        y_SHPB = y_SHPB_1[49:int(len(y_SHPB_1)-50)]
        #y_SHPB   = y_SHPB_2[0:58651]
    else: 
        x_SHPB_1 = SHPB_time_0[round(0.03 * len(SHPB_time_0)):round(0.35 * len(SHPB_time_0))]
        #x_SHPB_2 = x_SHPB_1[51:,]
        x_SHPB = x_SHPB_1[49:int(len(x_SHPB_1)-50)]
        #x_SHPB   = x_SHPB_2[0:58651]
        y_SHPB_1 = SHPB_inc_raw[round(0.03 * len(SHPB_time_0)):round(0.35 * len(SHPB_time_0))]
        #y_SHPB_2 = y_SHPB_1[51:,]
        y_SHPB = y_SHPB_1[49:int(len(y_SHPB_1)-50)]
        #y_SHPB   = y_SHPB_2[0:58651]
        
    y_s = smooth(y_SHPB_1,100) # Smooth the data trace to limit noise effects 
    
    
#%% Normalize the SHPB data traces prior to taking the Hough Transform
    x_norm = (x_SHPB - np.min(x_SHPB)) / (np.max(x_SHPB - np.min(x_SHPB)))
    y_norm = (y_SHPB - np.min(y_SHPB)) / (np.max(y_SHPB - np.min(y_SHPB)))
    y_s_norm = (y_s - np.min(y_s)) / (np.max(y_s - np.min(y_s)))
    
#    plt.plot(x_norm, y_s_norm)
#    plt.title('Python Smoothed Incident Wave')
#    plt.xlabel('time norm'), plt.ylabel('y_s norm')
    
    down_sampling = round(len(x_norm)/300.0)
    # Downsample x_norm and y_norm
    x_norm_ds = x_norm[::down_sampling]
    y_norm_ds = y_s_norm[::down_sampling]
        
    theta_steps = 180
    theta = np.arange(0,math.pi,math.pi/theta_steps)
    d_theta = theta[1]-theta[0]
    deg = np.array(theta)*180/math.pi
    
    # Initialize variables
    C = np.zeros((len(x_norm_ds),1))
    theta_B = np.zeros((len(theta),len(x_norm_ds)))
    theta_A = np.zeros((len(theta),len(x_norm_ds)))
    A = np.zeros((len(theta),len(x_norm_ds)))
    A_x = np.zeros((len(theta),len(x_norm_ds)))
    A_y = np.zeros((len(theta),len(x_norm_ds)))
    B = np.zeros((len(theta),len(x_norm_ds)))
    theta_step = np.zeros((len(theta)))
    



#%% 3.2 The following loop conducts the Hough transform 
k = np.arange(0,int(len(x_norm_ds)))
for j in range(0,int(len(x_norm_ds))):
    for i in range(0,int(len(theta))):
        x = x_norm_ds[k[j]]
        y = y_norm_ds[k[j]]
        C[j] = math.sqrt(x**2 + y**2)
        theta_1 = math.degrees(math.atan(y/x))
        theta_2 = 90 - theta_1
        theta_step[i] = deg[i] # theta vector in Hough space
        theta_B[i,j] = theta_2 + theta_step[i]
        theta_A[i,j] = 90 - theta_B[i,j]
        A[i,j] = C[j]*math.sin(math.radians(theta_A[i,j])) # r vector in Hough space
        B[i,j] = C[j]*math.sin(math.radians(theta_B[i,j]))
        A_x[i,j] = -A[i,j]*math.sin(math.radians(theta_step[i])) # used for creating a movie of the Hough Transform
        A_y[i,j] = A[i,j]*math.cos(math.radians(theta_step[i]))  # used for creating a movie of the Hough Transform
    # end for i
# end for j

R = np.absolute(A)

#plt.plot(R,theta)
#plt.title('Python Hough Transform')
#plt.xlabel('r'), plt.ylabel('theta')
    
del A, A_x, A_y, B, C, theta_B, theta_A, theta_1, theta_2, x, y # clear unnecessary large variables
 
#%% 3.3 Determine the line intersections in Hough Space

r_intersect = np.empty(0) # initialize vectors
theta_intersect = np.empty(0)

print('')
print('Hough Space Calculation....In progress')

k = np.shape(R) # total number of lines to calculate intersections over (same as the number of data points in the transform)
n = k[0]
k = k[1]

for i in range(0,len(np.arange(1,k))):
    for j in range((i+1),len((np.arange(2,k)))):
        if i != j:
            # loop steps forwards for each line,i=2 is compared to 3:end,
            # following loop is i=3 is compared to 4:end, etc.
#            x1 = R[:,i]
#            x2 = R[:,j]
#            r_int, theta_int = segmentIntersect(x1,deg,x2,deg)
#            r_intersect     = np.append(r_intersect, r_int)
#            theta_intersect = np.append(theta_intersect, theta_int)
            for z in range(0,(len(np.arange(0,n))-1)): # Check if line segments intersect using cross product
                x1          = [R[z,i], R[(z+1),i]]
                y1          = [deg[z], deg[(z+1)]]
                x2          = [R[z,j], R[(z+1),j]]
                y2          = [deg[z], deg[(z+1)]]
                intersect   = intersectCheck(x1,y1,x2,y2)
                if intersect == 1:
                    m_i             = (deg[z+1] - deg[z])/(R[z+1,i] - R[z,i])
                    m_j             = (deg[z+1] - deg[z])/(R[z+1,j] - R[z,j])
                    b_i             = deg[z] - m_i*R[z,i]
                    b_j             = deg[z] - m_j*R[z,j]
                    r_int           = np.array([(b_j - b_i)/(m_i - m_j)])
                    theta_int       = np.array([m_i*r_int + b_i])
                    r_intersect     = np.append(r_intersect, r_int)
                    theta_intersect = np.append(theta_intersect, theta_int)                
                else: # doesn't intersect
                    continue
        elif i == j: # prevent self comparison
            continue 
    PercentComplete = ((sum(range(2,int(len(x_norm_ds))))) - sum(range(2,(int(len(x_norm_ds)))-i)))/(sum(range(2,int(len(x_norm_ds)))))*100
        # PercentComplete = (int(len(x_norm_ds) - 2) - (int(len(x_norm_ds)-2) - i)) / int(len(x_norm_ds) - 2) * 100        
    print('    Hough Space Calculation....%.2f%%'%PercentComplete )
print('Hough Space Calculation....COMPLETE')        

del r_int, theta_int, intersect, m_i, m_j, b_i, b_j

# Calculate the bivariate histogram, needed to know the max number of intersections
r_bins = round(len(x_norm_ds)/2);
theta_bins = round((len(deg))/1)
N1 = np.histogram2d(r_intersect.ravel(), theta_intersect.ravel(),bins=(r_bins, theta_bins))
N = N1[0]

#%% Save data from previous section
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
np.save('r_intersect2.npy', r_intersect)
np.save('theta_intersect2.npy', theta_intersect)
#np.savetxt('r_bins.txt', r_bins)
#np.savetxt('theta_bins',theta_bins) 

#%% Load r_intersect and theta_intersect data, create histogram
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
r_intersect     = np.load('r_intersect2.npy')
theta_intersect = np.load('theta_intersect2.npy')

r_bins = round(len(x_norm_ds)/2);
theta_bins = round((len(deg))/1)
N1 = np.histogram2d(r_intersect, theta_intersect,bins=(r_bins, theta_bins))
#N1 = plt.hist2d(r_intersect, theta_intersect,bins=(r_bins, theta_bins))
N = N1[0]
#N = N[:, 1:len(np.transpose(N))-1]
#N = N1[:,0:int(len(np.transpose(N1)))]


#%% 3.4 Feature extraction from Hough space

if ps == 1:
    # There are two slopes of interest, the rising slope and the falling slope.
    # Thus there are only two points needed in the Hough Transform. Initial
    # code is limited to just the rising slope.
    pass1 = 0
    i = 0
#    r = np.empty
#    c = np.empty
    r = np.zeros([1],dtype = int)
    c = np.zeros([1],dtype = int)
    while pass1 < 1:
        val = np.max(N)
        [r1,c1] = np.where(N == val)
        # Loop addresses the case when the same number of intersections occurs 
        # in multiple bin locations in the hough space
        if len(r1) == 1:
            #if i == 0:
            r = np.insert(r,i,r1+1)
            c = np.insert(c,i,c1+1)
            #r[0] = r1
            #c[0] = c1
                
        else:
            r = np.insert(r,i,r1[0]+1)
            c = np.insert(c,i,c1[0]+1)
#        else:
#            if i == 0:
           # r[i] = r1[0]
           # c[i] = c1[0]
#            else:
#                r = np.insert(r,i, r1[0])
#                c = np.insert(c,i, c1[0])

        # Loop places limits on only searching for cases of a positive
        # slope. This search is more restrictive than the pulse shaped case
        if c[i] < (50/(180/theta_bins)):
            N[r[i]-1, c[i]-1] = 0
        elif c[i] > (90/(180/theta_bins)):
            N[r[i]-1, c[i]-1] = 0
        elif r[i] > math.floor(len(x_norm_ds)/2): # this may set a falling slope equal to zero
            N[r[i]-1, c[i]-1] = 0
        else: 
            pass1 = 1
        i = i + 1
        
    r = r[0:int(len(r)-1)]
    c = c[0:int(len(c)-1)]
    
    # Leading edge position
    r_interest_rise     = r[-1] * max(r_intersect)/r_bins
    theta_interest_rise = c[-1] * max(theta)/theta_bins * 180/math.pi
    theta_interest_rise_adjust = theta_interest_rise - 90
    
    pass2 = 0
    N1 = np.histogram2d(r_intersect, theta_intersect,bins=(r_bins, theta_bins))
    N = N1[0]
    #N = N[:, 1:len(np.transpose(N))-1]
    val2 = np.zeros(1)
    i = 0
#    r = np.empty
#    c = np.empty
    r = np.zeros([1],dtype = int)
    c = np.zeros([1],dtype = int)
    while pass2 < 1:
        val2 = np.insert(val2,i,np.max(N))
        r1,c1 = np.where(N == val2[i])
        # Loop addresses the case when the same number of intersections occurs 
        # in multiple bin locations in the hough space
        if len(r1) == 1:
            #if i == 0:
            r = np.insert(r,i,r1+1)
            c = np.insert(c,i,c1+1)
#                r[0] = r1
#                c[0] = c1
        else:
            #else:
            r = np.insert(r,i,r1[0]+1)
            c = np.insert(c,i,c1[0]+1)
#            r[i] = r1[0]
#            c[i] = c1[0]
#        else:
#            if i == 0:
#                r[i] = r1[i]
#                c[i] = c1[i]
#            else:
#                r = np.insert(r,i, r1[0])
#                c = np.insert(c,i, c1[0])
        # Loop places limits on only searching for cases of a positive slope
        if c[i] < 90 / (180/theta_bins):
            N[r[i]-1, c[i]-1] = 0
        elif r[i] < math.floor(r_bins/2):
            N[r[i]-1, c[i]-1] = 0
        elif r[i] > r_bins - 3:
            N[r[i]-1, c[i]-1] = 0
        else:
            pass2 = 1;
        i = i + 1
        
    r = r[0:int(len(r)-1)]
    c = c[0:int(len(c)-1)]
    
    # Falling edge position
    r_interest_fall     = r[-1] * max(r_intersect)/r_bins
    theta_interest_fall = c[-1] * max(theta)/theta_bins * 180/math.pi
    theta_interest_fall_adjust = theta_interest_fall - 90

# Case WITHOUT pulse shaping       
else:
    # There are two slopes of interest, the rising slope and the falling slope.
    # Thus there are only two points needed in the Hough Transform. Initial
    # code is limited to just the rising slope.
    pass1 = 0
    i = 0
    r = 0
    c = 0    
    while pass1 < 1:
        val = np.max(N)
        [r1,c1] = np.where(N == val)
        # Loop addresses the case when the same number of intersections occurs 
        # in multiple bin locations in the hough space
        if len(r1) == 1:
            # if i == 0:
                # r[0] = r1
                # c[0] = c1
            #else:
            r = np.insert(r,i,r1)
            c = np.insert(c,i,c1)
        else:
            #if i == 0:
            #    r[i] = r1[i]
            #    c[i] = c1[i]
            #else:
            r = np.insert(r,i, r1[0])
            c = np.insert(c,i, c1[0])

        # Loop places limits on only searching for cases of a positive
        # slope. This search is more restrictive than the pulse shaped case
        if c[i] < 80/(180/theta_bins):
            # N[r[i], c[i]] = 0
            np.insert(N, ([r[i],c[i]]), 0)
        elif c[i] > 90/(180/theta_bins):
            # N[r[i], c[i]] = 0
            np.insert(N, ([r[i],c[i]]), 0)
        elif r[i] > math.floor(r_bins/2): # this may set a falling slope equal to zero
            # N[r[i], c[i]] = 0
            np.insert(N, ([r[i],c[i]]), 0)
        else: 
            pass1 = 1
        i = i + 1

    # Leading edge position
    r_interest_rise     = r[-1] * max(r_intersect)/r_bins
    theta_interest_rise = c[-1] * max(theta)/theta_bins * 180/math.pi
    theta_interest_rise_adjust = theta_interest_rise - 90
    
    pass2 = 0
    N = N1[0]
    N = N[:, 1:len(np.transpose(N))-1]
    i = 0
    r = np.zeros([1],dtype = int)
    c = np.zeros([1],dtype = int)
    while pass2 < 1:
        val[i] = np.insert(val,i,np.max(N))
        [r1,c1] = np.where(N == val[i])
        # Loop addresses the case when the same number of intersections occurs 
        # in multiple bin locations in the hough space
        if len(r1) == 1:
            #if i == 0:
            np.insert(r,i,r1)
            np.insert(c,i,c1)
              #  r[0] = r1
              #  c[0] = c1
            #else:
        else:
#            if i == 0:
#                r[i] = r1[i]
#                c[i] = c1[i]
#            else:
            r = np.insert(r,i, r1[0])
            c = np.insert(c,i, c1[0])
        # Loop places limits on only searching for cases of a positive slope
        if c[i] < 90 / (180/theta_bins):
            #N[r[i], c[i]] = 0
            np.insert(N, [r[i],c[i]], 0)
        elif r[i] < math.floor(r_bins/2):
            #N[r[i], c[i]] = 0
            np.insert(N, [r[i],c[i]], 0)
        elif r[i] > r_bins - 3:
            #N[r[i], c[i]] = 0
            np.insert(N, [r[i],c[i]], 0)
        else:
            pass2 = 1;
        i = i + 1
    
    # Falling edge position
    r_interest_fall     = r[-1] * max(r_intersect)/r_bins
    theta_interest_fall = c[-1] * max(theta)/theta_bins * 180/math.pi
    theta_interest_fall_adjust = theta_interest_fall - 90
#%%
# The value for theta starts at 12-noon on a clock (i.e., vertical axis). 
# Subtract 90 degrees so that the angular measurements are from the 
# horizontal.
x_h_rise = r_interest_rise*(math.cos(math.radians(theta_interest_rise_adjust)))
y_h_rise = r_interest_rise*(math.sin(math.radians(theta_interest_rise_adjust)))        
slope_perp_rise = -x_h_rise/y_h_rise   

x_h_fall = r_interest_fall*(math.cos(math.radians(theta_interest_fall_adjust)))    
y_h_fall = r_interest_fall*(math.sin(math.radians(theta_interest_fall_adjust)))    
slope_perp_fall = -x_h_fall/y_h_fall
#%%
# Determine slopes of possible solutions
slope_Hough = np.vstack([slope_perp_rise, slope_perp_fall])
num_slopes = len(slope_Hough)     
      
# Down size to only account for possible slopes
x_h = np.vstack([x_h_rise, x_h_fall])
y_h = np.vstack([y_h_rise, y_h_fall])        
        
y_Hough = np.zeros([len(x_norm_ds), (num_slopes+1)])
for j in range(0, (num_slopes)):
    for i in range(0, (len(x_norm_ds-1))):
        y_Hough[i,j] = y_h[j] + slope_Hough[j]*(x_norm_ds[i] - x_h[j])
           
#%%
# Recreate the signal in x-y space
x_intersect = np.zeros([0],dtype = int)
y_intersect = np.zeros([0],dtype = int)

intersectAssign = 0
for j in range(0, (num_slopes)):
    intersects = poly_x_poly(x_norm_ds, y_Hough[:,j], x_norm, y_s_norm)
    if intersectAssign == 0: 
        all_intersects = intersects
        intersectAssign = 1
    else:
        all_intersects = np.vstack([all_intersects, intersects])
    
# zero_amp_intersects = poly_x_poly(x_norm_ds, y_Hough[:,2], x_norm, y_s_norm) 

x_intersect = all_intersects[:,0]
y_intersect = all_intersects[:,1]              
           
#%%
intersections = np.transpose(np.vstack([x_intersect, y_intersect]))
intersections = intersections[intersections[:,0].argsort()]

#%% Save intersection data
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
np.save('x_intersect.npy',x_intersect)
np.save('y_intersect.npy',y_intersect)
np.save('intersections.npy',intersections)

#%% Load recreated x-y signal
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
x_intersect = np.load('x_intersect.npy')
y_intersect = np.load('y_intersect.npy')
intersections = np.load('intersections.npy')

#%% 3.5 Calculate Intersctions to find the windows to use for data analysis

# Incident Signal Window initialization
inc_window_length       = np.zeros([len(c_0)-1],dtype = float)
inc_window_start        = np.zeros([len(c_0)-1],dtype = float)
inc_window_end          = np.zeros([len(c_0)-1],dtype = float)  
inc_window_start_index  = np.zeros([len(c_0)-1],dtype = int)
inc_window_end_index    = np.zeros([len(c_0)-1],dtype = int)
inc_pts                 = np.zeros([len(c_0)-1],dtype = float) 

# Reflected Signal Window initialization
ref_window_start        = np.zeros([len(c_0)-1],dtype = float)
ref_window_end          = np.zeros([len(c_0)-1],dtype = float)  
ref_window_start_index  = np.zeros([len(c_0)-1],dtype = int)
ref_window_end_index    = np.zeros([len(c_0)-1],dtype = int)
ref_pts                 = np.zeros([len(c_0)-1],dtype = float) 

# Transmitted Signal Window initialization
trans_window_start        = np.zeros([len(c_0)-1],dtype = float)
trans_window_end          = np.zeros([len(c_0)-1],dtype = float)   
trans_window_start_index  = np.zeros([len(c_0)-1],dtype = int)
trans_window_end_index    = np.zeros([len(c_0)-1],dtype = int)
trans_pts                 = np.zeros([len(c_0)-1],dtype = float)

# Incident Signal Window
for m in range(0, (len(c_0)-1)):
    intersection_index1 = find_first(x_norm_ds, intersections[0,0])
    intersection_index2 = find_first(x_norm_ds, intersections[-1,0])        
    
    x_unnorm_ds = (x_norm_ds)*(max(x_SHPB) - min(x_SHPB)) + min(x_SHPB)

    inc_window_length[m] = ((2*L_striker)/c_0[m])   
    if ps == 1:
        inc_window_start[m]         = x_unnorm_ds[intersection_index1] - (inc_window_length[m]*0.2)
        inc_window_end[m]           = x_unnorm_ds[intersection_index2] + (inc_window_length[m]*0.2)
        inc_window_start_index[m]   = find_first(SHPB_time_0, inc_window_start[m])
        inc_window_end_index[m]     = find_first(SHPB_time_0, inc_window_end[m])
        inc_pts[m]                  = (inc_window_end_index[m] - inc_window_start_index[m] + 1)
    else:
        inc_window_start[m]         = x_unnorm_ds[intersection_index1] - (inc_window_length[m]*0.2)
        inc_window_end[m]           = x_unnorm_ds[intersection_index2] + (inc_window_length[m]*0.2)
        inc_window_start_index[m]   = find_first(SHPB_time_0, inc_window_start[m])
        inc_window_end_index[m]     = find_first(SHPB_time_0, inc_window_end[m])
        inc_pts[m]                  = (inc_window_end_index[m] - inc_window_start_index[m] + 1)

# Reflected Signal Window
    ref_window_start[m]         = inc_window_start[m]+(2*dx_inc)/c_0[m]
    ref_window_end[m]           = inc_window_end[m]+(2*dx_inc)/c_0[m]
    ref_window_start_index[m]   = find_first(SHPB_time_0, ref_window_start[m])
    ref_window_end_index[m]     = find_first(SHPB_time_0, ref_window_end[m])
    ref_pts[m]                  = ref_window_end_index[m] - ref_window_start_index[m] + 1
    
# Transmitted Signal Window    
    trans_window_start[m]         = inc_window_start[m] + (dx_inc + dx_trans)/c_0[m]
    trans_window_end[m]           = inc_window_end[m] + (dx_inc + dx_trans)/c_0[m]
    trans_window_start_index[m]   = find_first(SHPB_time_0, trans_window_start[m])
    trans_window_end_index[m]     = find_first(SHPB_time_0, trans_window_end[m])
    trans_pts[m]                  = trans_window_end_index[m] - trans_window_start_index[m] + 1
    
    if m % 20 == 0:
        percent_signal_window = m/(len(c_0)-1)*100
        print('Window for data analysis... %.2f %%'%percent_signal_window)
print(' ') 
print('Window for data analysis... COMPLETE')
print(' ') 
    
#%% 4. Calculate stress-strain behavior without dispersion correction 
    
print('Non-dispersion Corrected Data...Initiated')

pts = np.vstack([max(inc_pts), max(ref_pts), max(trans_pts)])
pts = int(max(pts))

# Initialize Variables
inc_time = np.zeros([pts, len(c_0)], dtype = float)
inc_volt = np.zeros([pts, len(c_0)], dtype = float)
Vr_inc = np.zeros([pts, len(c_0)], dtype = float)

ref_time = np.zeros([pts, len(c_0)], dtype = float)
ref_volt = np.zeros([pts, len(c_0)], dtype = float)
Vr_ref = np.zeros([pts, len(c_0)], dtype = float)

trans_time = np.zeros([pts, len(c_0)], dtype = float)
trans_volt = np.zeros([pts, len(c_0)], dtype = float)
Vr_trans = np.zeros([pts, len(c_0)], dtype = float)

strain_bar_inc = np.zeros([pts, len(c_0)], dtype = float)
strain_bar_ref = np.zeros([pts, len(c_0)], dtype = float)
strain_bar_trans = np.zeros([pts, len(c_0)], dtype = float)

P1 = np.zeros([pts, len(c_0)], dtype = float)
P2 = np.zeros([pts, len(c_0)], dtype = float)

strain_rate_spec = np.zeros([pts, len(c_0)], dtype = float)
stress_spec = np.zeros([pts, len(c_0)], dtype = float)
strain_spec = np.zeros([pts, len(c_0)], dtype = float)

for m in range(0, len(c_0)-1):
    i_time = SHPB_time_raw[(inc_window_start_index[m]):(inc_window_end_index[m])]
    i_volt = SHPB_inc_raw[(inc_window_start_index[m]):(inc_window_end_index[m])]
    i_time_len = len(i_time)
    
    if i_time_len < pts:
        inc_time[0:i_time_len,m] = i_time
        inc_volt[0:i_time_len,m] = i_volt
    else:
        inc_time[:,m] = i_time
        inc_volt[:,m] = i_volt
    
    Vr_inc[:,m] = inc_volt[:,m] / (V_ex*Gain)
    
    r_time =  SHPB_time_raw[(ref_window_start_index[m]):(ref_window_end_index[m])]
    r_volt = -SHPB_inc_raw[(ref_window_start_index[m]):(ref_window_end_index[m])] # negative due to wave being tensile
    r_time_len = len(r_time)
    
    if r_time_len < pts:
        ref_time[0:r_time_len,m] = r_time
        ref_volt[0:r_time_len,m] = r_volt
    else:
        ref_time[:,m] = r_time
        ref_volt[:,m] = r_volt
    
    Vr_ref[:,m] = ref_volt[:,m] / (V_ex*Gain)
    
    t_time = SHPB_time_raw[(trans_window_start_index[m]):(trans_window_end_index[m])]
    t_volt = SHPB_trans_raw[(trans_window_start_index[m]):(trans_window_end_index[m])]
    t_time_len = len(t_time)
    
    if t_time_len < pts:
        trans_time[0:t_time_len,m] = t_time
        trans_volt[0:t_time_len,m] = t_volt
    else:
        trans_time[:,m] = t_time
        trans_volt[:,m] = t_volt
        
    Vr_trans[:,m] = trans_volt[:,m] / (V_ex*Gain)
    
    if m % 20 == 0:
        percent_Vr = m/(len(c_0)-1)*100
        print('Creating voltage and time arrays... %.2f %%'%percent_Vr)
print('Creating voltage and time arrays... COMPLETE')
print(' ') 
#%%

# Initialize stress/strain arrays
strain_bar_inc = np.zeros([len(Vr_inc[:,0]), len(c_0)])
strain_bar_ref = np.zeros([len(Vr_ref[:,0]), len(c_0)])
strain_bar_trans = np.zeros([len(Vr_trans[:,0]), len(c_0)])
strain_rate_spec = np.zeros([len(Vr_trans[:,0]), len(c_0)])
stress_spec = np.zeros([len(Vr_trans[:,0]), len(c_0)])

# Calculate the engineering strain values from the SHPB signals
for m in range(0, len(c_0)-1):
    percent_stress_strain = m/(len(c_0)-1)*100
    print('Calculate Stress/Strain Values... %.2f %%'%percent_stress_strain)
    
    if (bar_mat1 == 1) & (D_bar >= 0.019): # Aluminum 3/4"
        for i in range(0, len(Vr_inc[:,m])-1):
            strain_bar_inc[i,m]   = (2*Vr_inc[i,m])   / (GF*((v_bar + 1) - Vr_inc[i,m]*(v_bar - 1)))
            strain_bar_ref[i,m]   = (2*Vr_ref[i,m])   / (GF*((v_bar + 1) - Vr_ref[i,m]*(v_bar - 1)))
            strain_bar_trans[i,m] = (2*Vr_trans[i,m]) / (GF*((v_bar + 1) - Vr_trans[i,m]*(v_bar - 1)))
    elif (bar_mat1 == 2) & (D_bar >= 0.019): # Steel 3/4"
        for i in range(0, len(Vr_inc[:,m])-1):
            strain_bar_inc[i,m]   = (2*Vr_inc[i,m])   / (GF*((v_bar + 1) - Vr_inc[i,m]*(v_bar - 1)))
            strain_bar_ref[i,m]   = (2*Vr_ref[i,m])   / (GF*((v_bar + 1) - Vr_ref[i,m]*(v_bar - 1)))
            strain_bar_trans[i,m] = (2*Vr_trans[i,m]) / (GF*((v_bar + 1) - Vr_trans[i,m]*(v_bar - 1)))
    elif (bar_mat1 == 1) & (D_bar < 0.019): # Aluminum 1/2"
        for i in range(0, len(Vr_inc[:,m])-1):
            strain_bar_inc[i,m]   = (4*Vr_inc[i,m])   / (GF*((v_bar + 1) - Vr_inc[i,m]*(v_bar - 1)))
            strain_bar_ref[i,m]   = (4*Vr_ref[i,m])   / (GF*((v_bar + 1) - Vr_ref[i,m]*(v_bar - 1)))
            strain_bar_trans[i,m] = (4*Vr_trans[i,m]) / (GF*((v_bar + 1) - Vr_trans[i,m]*(v_bar - 1)))
    elif (bar_mat1 == 2) & (D_bar < 0.019): # Steel 1/2"
        for i in range(0, len(Vr_inc[:,m])-1):
            strain_bar_inc[i,m]   = (4*Vr_inc[i,m])   / (GF*((v_bar + 1) - Vr_inc[i,m]*(v_bar - 1)))
            strain_bar_ref[i,m]   = (4*Vr_ref[i,m])   / (GF*((v_bar + 1) - Vr_ref[i,m]*(v_bar - 1)))
            strain_bar_trans[i,m] = (4*Vr_trans[i,m]) / (GF*((v_bar + 1) - Vr_trans[i,m]*(v_bar - 1)))

    # The following assumes 1-D wave propagation, calculate the forces of the 
    # incident-specimen and the transmitted-specimen interfaces
    P1[:,m] = E_bar * A_bar * (strain_bar_inc[:,m] - strain_bar_ref[:,m]) # Force at the incident-specimen interface
    P2[:,m] = E_bar * A_bar * (strain_bar_trans[:,m]) # Force at the transmitted-specimen interface
    
    # Calculate the Nominal Specimen Strain Rate
    strain_rate_spec[:,m] = (2*c_0[m]*strain_bar_ref[:,m])/L_spec
    
    # Calculate the Nominal Specimen Stress
    stress_spec[:,m] = P2[:,m] / (math.pi*(D_spec/2)**2)
    
    # Calculate the Specimen Strain
    for i in range(0,len(strain_rate_spec[:,m])-2):
        strain_spec[i,m] = np.trapz(ref_time[0:i,m], strain_rate_spec[0:i,m])
        
#    percent_stress_strain = m/(len(c_0)-1)*100
#    print('Calculate Stress/Strain Values... %.2f %%'%percent_stress_strain)

print(' ')    
print('Non-Dispersion Correct Data... COMPLETE')

#%% Save strain values
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
np.save('strain_bar_inc.npy',strain_bar_inc)
np.save('strain_bar_ref.npy',strain_bar_ref)
np.save('strain_bar_trans.npy',strain_bar_trans)
np.save('P1.npy',P1)
np.save('P2.npy',P2)
np.save('strain_rate_spec.npy',strain_rate_spec)
np.save('stress_spec.npy',stress_spec)
np.save('strain_spec.npy',strain_spec)

#%% Load strain values
"""ONLY USED DURING WRITING. DELETE WHEN FINSIHED"""
strain_bar_inc = np.load('strain_bar_inc.npy')
strain_bar_ref = np.load('strain_bar_ref.npy')
strain_bar_trans = np.load('strain_bar_trans.npy')
P1 = np.load('P1.npy')
P2 = np.load('P2.npy')
strain_rate_spec = np.load('strain_rate_spec.npy')
stress_spec = np.load('stress_spec.npy')
strain_spec = np.load('strain_spec.npy')

    
#%% 6. Perform the dispersion correction of the SHPB data
print('Dispersion Correction...Initiated')

# Constants of the Gong approximation to the Pochhammer-Chree Solution
if bar_mat1 == 1: # Aluminum Bars
    A =  0.5710
    B =  0.4301
    C =  16.764
    D =  19.679
    E = -5.854
    F =  2.153
elif bar_mat1 == 2: # Steel Bars (Poisson's ratio of 0.3)
    A =  0.5746
    B =  0.4268
    C =  19.530
    D =  19.702
    E = -7.005
    F =  2.389

kmax = 100; # Sets the number of frequencies to use in dispersion correction

#%% 6.1 Dispersion Correct the Incident signal

inc_volt_dc_loop = np.zeros([len(inc_volt), len(c_0)-1]) # Initialize dc voltage array

for m in range(0, len(c_0)-1):   
    inc_len      = len(inc_volt[:,m])
    pad          = np.zeros(inc_len)
    inc_volt_pad = np.concatenate((pad,inc_volt[:,m],pad))
    N_inc        = len(inc_volt_pad)
    n_inc        = np.arange(1, N_inc+1)
    Ti_inc       = 1/(N_inc*dt)
    omega0_inc   = 2*math.pi*Ti_inc
    k            = np.arange(1,kmax+1)  # Vectorized kmax
    ck0          = 0.5*c_0[m]       # Initial guess for zero solver

    # Generating the ck/c0 -- a/gamma dispersion curve (solving for phi)
    # Vectorised a*k*omega/(2*pi)
    akomega = np.transpose(R_bar*k*Ti_inc)
    
    # FSOLVE is computationaly expensive.  With the assumtion that ck for
    # frequencies beyond 1250/a are constant, limit loop to minidx and loop
    # backwards, adjusting initial value parameter of FSOLVE
    minidx = int(np.min([np.ceil(((1250/R_bar)*1000)/Ti_inc), len(akomega)-1]))
    
    c0_loop = c_0[m]
    # Get end value, set initial search start point
    lowck = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[minidx]),ck0)

    # Initialize ck array with low ck values
    ck = np.ones(kmax)*lowck
    
    # Loop fsolve solver backwards, updating ck value
    for i in range(minidx,0,-1):
        ck[i] = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[i]),ck0)
        ck0 = ck[i]
    
    # Limit values to c0
    ck[ck > c_0[m]] = c_0[m]
    
    # Solve for the phase change
    phi_inc = k * omega0_inc * (dx_inc/ck - dx_inc/c_0[m])
    
    # Forward Fourier Transform
    M1 = np.tile((omega0_inc*n_inc*dt), (kmax, 1))
    cos_term = np.transpose(np.cos(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    sin_term = np.transpose(np.sin(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    
    Ak = 2* Ti_inc * np.sum(np.apply_along_axis(np.multiply,0, inc_volt_pad, cos_term)*dt, axis = 1)
    Bk = 2* Ti_inc * np.sum(np.apply_along_axis(np.multiply,0, inc_volt_pad, sin_term)*dt, axis = 1)

    # Shifted & Inverse Fourier Transform
    M2 = np.tile((omega0_inc*k*dt),(N_inc, 1))
    term1 = np.apply_along_axis(np.multiply,0,n_inc, np.transpose(M2))
    I1 = np.apply_along_axis(np.subtract,0, term1, phi_inc)
    
    I1_sin_term = np.apply_along_axis(np.multiply,0,Bk,np.transpose(np.sin(I1)))
    I1_cos_term = np.apply_along_axis(np.multiply,0,Ak,np.transpose(np.cos(I1)))
    
    inc_volt_dc = Ti_inc * np.sum(inc_volt_pad*dt) + np.sum((I1_sin_term + I1_cos_term), axis = 1)
    inc_volt_dc = inc_volt_dc[inc_len:-inc_len]
    inc_volt_dc_loop[:,m] = inc_volt_dc
    
#    inc_volt_dc_loop[:,m] = inc_volt_dc
    if m % 6 == 0:
        percent_Inc_Fourier = m/(len(c_0)-1)*100
        print('Preforming Fourier Transform of Incident Signal... %.1f %%'%percent_Inc_Fourier)
print('Preforming Fourier Transform of Incident Signal... COMPLETE')


#%% 6.2 Dispersion Correct the Reflected Signal

ref_volt_dc_loop = np.zeros([len(ref_volt), len(c_0)-1]) # Initialize dc voltage array

for m in range(0, len(c_0)-1):   
    ref_len      = len(ref_volt[:,m])
    pad          = np.zeros(ref_len)
    ref_volt_pad = np.concatenate((pad,ref_volt[:,m],pad))
    N_ref        = len(ref_volt_pad)
    n_ref        = np.arange(1, N_ref+1)
    Ti_ref       = 1/(N_ref*dt)
    omega0_ref   = 2*math.pi*Ti_ref
    k            = np.arange(1,kmax+1)  # Vectorized kmax
    ck0          = 0.5*c_0[m]       # Initial guess for zero solver

    # Generating the ck/c0 -- a/gamma dispersion curve (solving for phi)
    # Vectorised a*k*omega/(2*pi)
    akomega = np.transpose(R_bar*k*Ti_ref)
    
    # FSOLVE is computationaly expensive.  With the assumtion that ck for
    # frequencies beyond 1250/a are constant, limit loop to minidx and loop
    # backwards, adjusting initial value parameter of FSOLVE
    minidx = int(np.min([np.ceil(((1250/R_bar)*1000)/Ti_ref), len(akomega)-1]))
    
    c0_loop = c_0[m]
    # Get end value, set initial search start point
    lowck = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[minidx]),ck0)

    # Initialize ck array with low ck values
    ck = np.ones(kmax)*lowck
    
    # Loop fsolve solver backwards, updating ck value
    for i in range(minidx,0,-1):
        ck[i] = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[i]),ck0)
        ck0 = ck[i]
    
    # Limit values to c0
    ck[ck > c_0[m]] = c_0[m]
    
    # Solve for the phase change
    phi_ref = k * omega0_ref * (dx_inc/ck - dx_inc/c_0[m])
    
    # Forward Fourier Transform
    M1 = np.tile((omega0_ref*n_ref*dt), (kmax, 1))
    cos_term = np.transpose(np.cos(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    sin_term = np.transpose(np.sin(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    
    Ak = 2* Ti_ref * np.sum(np.apply_along_axis(np.multiply,0, ref_volt_pad, cos_term)*dt, axis = 1)
    Bk = 2* Ti_ref * np.sum(np.apply_along_axis(np.multiply,0, ref_volt_pad, sin_term)*dt, axis = 1)

    # Shifted & Inverse Fourier Transform
    M2 = np.tile((omega0_ref*k*dt),(N_ref, 1))
    term1 = np.apply_along_axis(np.multiply,0,n_ref, np.transpose(M2))
    I1 = np.apply_along_axis(np.subtract,0, term1, phi_ref)
    
    I1_sin_term = np.apply_along_axis(np.multiply,0,Bk,np.transpose(np.sin(I1)))
    I1_cos_term = np.apply_along_axis(np.multiply,0,Ak,np.transpose(np.cos(I1)))
    
    ref_volt_dc = Ti_ref * np.sum(ref_volt_pad*dt) + np.sum((I1_sin_term + I1_cos_term), axis = 1)
    ref_volt_dc = ref_volt_dc[ref_len:-ref_len]
    ref_volt_dc_loop[:,m] = ref_volt_dc
    
#    inc_volt_dc_loop[:,m] = inc_volt_dc
    if m % 6 == 0:
        percent_ref_Fourier = m/(len(c_0)-1)*100
        print('Preforming Fourier Transform of Reflected Signal... %.1f %%'%percent_ref_Fourier)
print('Preforming Fourier Transform of Reflected Signal... COMPLETE')


#%% 6.3 Dispersion Correct the Transmitted Signal

trans_volt_dc_loop = np.zeros([len(trans_volt), len(c_0)-1]) # Initialize dc voltage array

for m in range(0, len(c_0)-1):   
    trans_len      = len(trans_volt[:,m])
    pad            = np.zeros(trans_len)
    trans_volt_pad = np.concatenate((pad,trans_volt[:,m],pad))
    N_trans        = len(trans_volt_pad)
    n_trans        = np.arange(1, N_ref+1)
    Ti_trans       = 1/(N_trans*dt)
    omega0_trans   = 2*math.pi*Ti_trans
    k              = np.arange(1,kmax+1)  # Vectorized kmax
    ck0            = 0.5*c_0[m]       # Initial guess for zero solver

    # Generating the ck/c0 -- a/gamma dispersion curve (solving for phi)
    # Vectorised a*k*omega/(2*pi)
    akomega = np.transpose(R_bar*k*Ti_trans)
    
    # FSOLVE is computationaly expensive.  With the assumtion that ck for
    # frequencies beyond 1250/a are constant, limit loop to minidx and loop
    # backwards, adjusting initial value parameter of FSOLVE
    minidx = int(np.min([np.ceil(((1250/R_bar)*1000)/Ti_trans), len(akomega)-1]))
    
    c0_loop = c_0[m]
    # Get end value, set initial search start point
    lowck = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[minidx]),ck0)

    # Initialize ck array with low ck values
    ck = np.ones(kmax)*lowck
    
    # Loop fsolve solver backwards, updating ck value
    for i in range(minidx,0,-1):
        ck[i] = fsolve(getPhiFcn(c0_loop,A,B,C,D,E,F,akomega[i]),ck0)
        ck0 = ck[i]
    
    # Limit values to c0
    ck[ck > c_0[m]] = c_0[m]
    
    # Solve for the phase change
    phi_trans = k * omega0_trans * (dx_inc/ck - dx_inc/c_0[m])
    
    # Forward Fourier Transform
    M1 = np.tile((omega0_trans*n_trans*dt), (kmax, 1))
    cos_term = np.transpose(np.cos(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    sin_term = np.transpose(np.sin(np.apply_along_axis(np.multiply,0, k, np.transpose(M1))))
    
    Ak = 2* Ti_trans * np.sum(np.apply_along_axis(np.multiply,0, trans_volt_pad, cos_term)*dt, axis = 1)
    Bk = 2* Ti_trans * np.sum(np.apply_along_axis(np.multiply,0, trans_volt_pad, sin_term)*dt, axis = 1)

    # Shifted & Inverse Fourier Transform
    M2 = np.tile((omega0_trans*k*dt),(N_trans, 1))
    term1 = np.apply_along_axis(np.multiply,0,n_ref, np.transpose(M2))
    I1 = np.apply_along_axis(np.subtract,0, term1, phi_trans)
    
    I1_sin_term = np.apply_along_axis(np.multiply,0,Bk,np.transpose(np.sin(I1)))
    I1_cos_term = np.apply_along_axis(np.multiply,0,Ak,np.transpose(np.cos(I1)))
    
    trans_volt_dc = Ti_trans * np.sum(trans_volt_pad*dt) + np.sum((I1_sin_term + I1_cos_term), axis = 1)
    trans_volt_dc = ref_volt_dc[trans_len:-trans_len]
    trans_volt_dc_loop[:,m] = trans_volt_dc
    
#    inc_volt_dc_loop[:,m] = inc_volt_dc
    if m % 6 == 0:
        percent_trans_Fourier = m/(len(c_0)-1)*100
        print('Preforming Fourier Transform of Transmitted Signal... %.1f %%'%percent_trans_Fourier)
print('Preforming Fourier Transform of Transmitted Signal... COMPLETE')

#%% 7. Calculate the Dispersion Corrected Stress-Strain Behavior

strain_bar_inc_dc   = np.zeros([max(inc_pts), len(c_0)])
strain_bar_ref_dc   = np.zeros([max(inc_pts), len(c_0)])
strain_bar_trans_dc = np.zeros([max(inc_pts), len(c_0)])

P1_dc = np.zeros([max(inc_pts), len(c_0)])
P2_dc = np.zeros([max(inc_pts), len(c_0)])

strain_rate_spec_dc = np.zeros([max(inc_pts), len(c_0)])
stress_spec_dc      = np.zeros([max(inc_pts), len(c_0)])
strain_spec_dc      = np.zeros([max(inc_pts), len(c_0)])

Vr_inc_dc   = inc_volt_dc_loop / (V_ex*Gain)
Vr_ref_dc   = ref_volt_dc_loop / (V_ex*Gain)
# Vr_trans_ds = trans_volt_dc_loop / (V_ec*Gain)

# Calculate the engineering strain values form the SHPB signals
for m in range(0,len(c_0)-1): # Aluminum 3/4"
    if (bar_mat1 == 1) & (D_bar >= 0.019):
        for i in range(0, len(Vr_inc_dc)-1):
            strain_bar_inc_dc[i,m] = (2*Vr_inc_dc[i,m])/(GF* ( (v_bar + 1) - Vr_inc_dc[i,m] * (v_bar - 1)) )
            strain_bar_ref_dc[i,m] = (2*Vr_ref_dc[i,m])/(GF* ( (v_bar + 1) - Vr_ref_dc[i,m] * (v_bar - 1)) )
            #strain_bar_trans_dc[i,m] = (2*Vr_trans_dc[i,m])/(GF* ( (v_bar + 1) - Vr_trans_dc[i,m] * (v_bar - 1)) )
    
    elif (bar_mat1 == 2) & (D_bar >= 0.019): # Steel 3/4"
        for i in range(0, len(Vr_inc_dc)-1):
            strain_bar_inc_dc[i,m] = (2*Vr_inc_dc[i,m])/(GF* ( (v_bar + 1) - Vr_inc_dc[i,m] * (v_bar - 1)) )
            strain_bar_ref_dc[i,m] = (2*Vr_ref_dc[i,m])/(GF* ( (v_bar + 1) - Vr_ref_dc[i,m] * (v_bar - 1)) )
            #strain_bar_trans_dc[i,m] = (2*Vr_trans_dc[i,m])/(GF* ( (v_bar + 1) - Vr_trans_dc[i,m] * (v_bar - 1)) )
    
    elif (bar_mat1 == 1) & (D_bar < 0.019): # Aluminum 1/2"
        for i in range(0, len(Vr_inc_dc)-1):
            strain_bar_inc_dc[i,m] = (4*Vr_inc_dc[i,m])/(GF* ( 1 + 2*Vr_inc_dc) )
            strain_bar_ref_dc[i,m] = (4*Vr_ref_dc[i,m])/(GF* ( 1 + 2*Vr_ref_dc) )
            #strain_bar_trans_dc[i,m] = (4*Vr_trans_dc[i,m])/(GF* ( (v_bar + 1) - Vr_trans_dc[i,m] * (v_bar - 1)) )
     
    elif (bar_mat1 == 2) & (D_bar < 0.019): # Steel 1/2"
        for i in range(0, len(Vr_inc_dc)-1):
            strain_bar_inc_dc[i,m] = (4*Vr_inc_dc[i,m])/(GF* ( (v_bar + 1) - Vr_inc_dc[i,m] * (v_bar - 1)) )
            strain_bar_ref_dc[i,m] = (4*Vr_ref_dc[i,m])/(GF* ( (v_bar + 1) - Vr_ref_dc[i,m] * (v_bar - 1)) )
            #strain_bar_trans_dc[i,m] = (4*Vr_trans_dc[i,m])/(GF* ( (v_bar + 1) - Vr_trans_dc[i,m] * (v_bar - 1)) )
            
    if m % 6 == 0:
        percent_corrected = m/(len(c_0)-1)*100
        print('Calculating Dispersion Corrected Stress/Strain... %.1f %%'%percent_corrected)
        
    # The following assumes 1-D wave propagation, calculate the forces of the
    # incident-specimen and the transmitted-specimen interface
    P1_dc[:,m] = E_bar * A_bar * (strain_bar_inc_dc[:,m] - strain_bar_ref_dc[:,m]) # Force at the incident-specimen interface
    P2_dc[:,m] = E_bar * A_bar * (strain_bar_trans_dc[:,m]) # Froces at the transmitted-specimen interface

    # Calculate the Nominal Specimen Strain Rate
    strain_rate_spec_dc[:,m] = (2*c_0[m] * strain_bar_ref_dc[:,m])/L_spec
    
    # Calculate the Nominal specimen stress
    stress_spec_dc[:,m] = P2_dc[:,m]/(math.pi*(D_spec/2)**2)
    
    # Calculate the Specimen Stress
    '''
    Could speed this section up by using the pervious segment's area and adding 
    on the current section's area.
    '''
    for i in range(0,len(strain_rate_spec_dc[:,m])-1):
        strain_spec_dc[i,m] = np.trapz(ref_time[0:(i+1), m], strain_rate_spec_dc[0:(i+1), m])
        if i % 5000 == 0:
            percent_i = i/(len(strain_rate_spec_dc[:,m])-1)*100
            print('Integrating Strain Rate Spec... %.1f %%'%percent_i)
    
#    if m % 6 == 0:
#        percent_corrected = m/(len(c_0)-1)*100
#        print('Calculating Dispersion Corrected Stress/Strain... %.1f %%'%percent_corrected)
print('Calculating Dispersion Corrected Stress/Strain... COMPLETE')

#%% 8. Save Plots
print('Saving Plots...')

os.chdir('C:\Kingstedt\SHPB - Data and Analysis\Processed Data')
#os.mkdir(filename)
os.chdir('C:/Kingstedt/SHPB - Data and Analysis/Processed Data/{}'.format(filename))

#%%
# 8.1 Raw experimental data
plt.figure(figsize = (8.5,6.5))
h1 = plt.plot(SHPB_time_0, SHPB_trans_raw,'b',linewidth=0.5, label='Trans. Bar')
h2 = plt.plot(SHPB_time_0, SHPB_inc_raw,'r',linewidth=0.5, label='Inc. Bar')
plt.title('Raw Data',fontsize=14)
plt.xlabel('Time [s]',fontsize=12)
plt.ylabel('Strain Gauge Output [V]',fontsize=12)
plt.xlim((0, np.max(SHPB_time_0)))
plt.xticks(np.arange(0,0.0011,0.0001))
plt.ticklabel_format(axis = 'x', style = 'scientific',scilimits=(-2,2), useMathText = True)
plt.grid(True)
plt.legend(loc='best',fontsize=12, shadow=True)
plt.show

plt.savefig(('{}_Raw_Data.png'.format(filename)), format = 'png')

#%%
# 8.2 Incident Signal in Hough Space
plt.figure(figsize = (7,7))
for i in range(0,len(R[0,:])):
    if np.mod(i,10) == 0:
        plt.plot(R[:,i], deg, 'r.-',linewidth = 2)
    else:
        plt.plot(R[:,i], deg, 'b',linewidth = 0.5)
plt.grid(True)
plt.title('Hough Space Lines',fontsize=14)
plt.xlabel('Distance - r',fontsize=12)
plt.ylabel('Theta [deg]',fontsize=12)
plt.xlim((-0.03, math.sqrt(2)+0.03))
plt.ylim((-3,183))
plt.show

plt.savefig(('{}_Hough_Lines.png'.format(filename)), format = 'png')

#%%
# Hough Space Line Intersections
plt.figure(figsize=(7,6), constrained_layout=True)
fig3 = plt.hist2d(r_intersect, theta_intersect,bins=(len(x_norm_ds), len(deg)),
                  cmap=plt.cm.gist_yarg,)
cbar = plt.colorbar(fig3[3])
cbar.ax.set_ylabel('Counts', fontsize=12)
plt.title('Hough Space Line Intersections',fontsize=14)
plt.xlabel('Distance - r', fontsize=12)
plt.ylabel('Theta [deg]', fontsize=12)
plt.xlim((-0.03, math.sqrt(2)+0.03))
plt.ylim((-3,183))
plt.show

plt.savefig(('{}_Hough_Intersections.png'.format(filename)), format = 'png')

#%%
# Intersections in Normalized Space
plt.figure(figsize=(8.5,6.5))
plt.plot(x_intersect,y_intersect, 'ko')
plt.plot(x_norm_ds, y_norm_ds, 'r.-')
plt.plot(x_norm_ds, y_Hough[:,0], 'b')
plt.plot(x_norm_ds, y_Hough[:,1], 'b')
for i in range(0,len(all_intersects[:,0])):
    plt.plot(intersections[i,0], intersections[i,1],'kX', markersize = 7.5)
plt.grid(True)
plt.title('Intersections in Normalized Space',fontsize=14)
plt.xlabel('Time [s]',fontsize=12)
plt.ylabel('Normalized Strain Gauge Output [V]',fontsize=12)
plt.xlim((-0.03, 1.03))
plt.ylim((-0.03,1.03))
plt.show

plt.savefig(('{}_Norm_Intersections.png'.format(filename)), format = 'png')

#%%
# Raw Data with Hough Transform Predicted Windows

plt.figure(figsize=(8.5,6.5))
plt.plot(x_SHPB,y_SHPB, 'b', linewidth = 0.5)
plt.plot(SHPB_time_0, SHPB_inc_raw,'b', linewidth = 0.5)
plt.plot(SHPB_time_0, SHPB_trans_raw,'r',linewidth = 0.5)
#plot(inc_window_start[c_0_index,1])



plt.show

#%% Non-dispersion corrected stress-strain
plt.figure(figsize=(8.5,6.5))
plt.plot(strain_spec[:,0], stress_spec[:,0])
plt.grid(True)
plt.title('Non-Dispersion Corrected Stress-Strain w/o Window')
plt.xlabel('Strain')
plt.ylabel('Stress [MPa]')
plt.show

plt.savefig(('{}_Stress-Strain_no_dc.png'.format(filename)), format = 'png')

#%% Dispersion Corrected Incident Voltage
plt.figure(figsize=(8.5,6.5))
plt.plot(inc_time[:,0], inc_volt_dc_loop[:,0])
plt.grid(True)
plt.title('Dispersion Corrected Incident Bar Voltage w/o Window')
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.show

plt.savefig(('{}_Dispersion_Corrected_Inc_Volt.png'.format(filename)), format = 'png')




