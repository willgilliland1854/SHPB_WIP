import numpy as np
 
def cosd(theta):
    return np.cos(np.deg2rad(theta))

def sind(theta):
    return np.sin(np.deg2rad(theta))

def atand(theta):
   return np.rad2deg(np.tan(theta))

def houghTransform(x,y,theta = np.arange(0,180)):
    return x*sind(theta) + y*cosd(theta)
 
class Hough:
    
    def getAccumulator(self, x, y, theta = np.arange(0,180)):
        """
        Conducts the Hough transform for the data set given by x, y and creates
        an Hough accumulator array for that data
    
        Parameters
        ----------
        x : array or list
            the x values of data being transformed
            
        y : array or list
            the y values of data being transformed
        
        theta : array or list
            the angles theta of which the hough transform should iterate
        
        """
        # Allocate Space 
        self.accumulator = np.zeros((180, 142))
        for i in range(len(x)):
            # Rounding multiplying rho values for easy indexing for this
            # example
            self.rho = (x[i]*sind(theta) + y[i]*cosd(theta))*100
            for j in range(len(theta)):
                # Using np.round to bin values of rho before being added to the
                # accumulator to greatly increase the speed. 
                self.accumulator[j, int(np.round(self.rho[j]))] += 1
     
        return self
    
#%% A way to use histograms
bin_num = 200
rho_bins = np.linspace(-1.4, 1.4, bin_num);
rho = np.zeros(180)
for j in range(len(x_norm_ds)-1):
    for i in range(180):
        rho[i] = x_norm_ds[j]*sind(i) + y_norm_ds[j]*cosd(i)
    if j == 0:
        rho_total = rho
        theta_total = np.arange(0,180)
    else:
        rho_total = np.hstack([rho_total, rho])
        theta_total = np.hstack([theta_total, np.arange(0,180)])
        

accumulator = np.histogram2d(theta_total, rho_total, bins = [np.arange(0,180), rho_bins])

# finding the maxima in accumulator:
i,j = np.unravel_index(accumulator[0].argmax(), (accumulator[0].shape))
accumulator[1][i] # angle
accumulator[2][j] # radius
        