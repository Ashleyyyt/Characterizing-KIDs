import numpy as np
import matplotlib.pyplot as plt

# define normalized 2D gaussian
def gaus2d(x, y, mean_x, mean_y, sigma_x, sigma_y, amplitude):
    first_frac = ((x-mean_x)**2)/(2*sigma_x*2)
    second_frac = ((y-mean_y)**2)/(2*sigma_y*2)
    return amplitude*np.exp(-(first_frac+second_frac))

def calculate_Beam_Filling_Factor(length=100, width=20, FWHM=20, points=5000):

    #Define x and y
    x = np.linspace(0, length, points)
    y = np.linspace(0, length, points)
    x, y = np.meshgrid(x, y) # get 2D variables instead of 1D

    #mean: square box of (longest length)/2 
    mean = [length/2, length/2]

    #sigma = FWHM/(sqrt(8log2))
    divisor = np.sqrt(8*np.log(2))
    sigma = [FWHM/divisor, FWHM/divisor]
    z = gaus2d(x, y, mean[0], mean[1], sigma[0], sigma[1], 1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z)

    ax.set_xlabel('x / mm')
    ax.set_ylabel('y / mm')
    ax.set_zlabel('Intensity')

    #Integral = sum over all points
    Gauss_integral = np.sum(z)

    #Bar dimensions centered at length/2
    sub = width/2
    first_cuttoff = int(((mean[0]-sub)/length)*points)
    second_cuttoff = int(((mean[0]+sub)/length)*points)

    #Set outside hotbar dimensions = 0
    z[0:first_cuttoff,:] = 0
    z[second_cuttoff:,:] = 0

    #Integrate over all points 
    hotbar_integral = np.sum(z)

    #Calculate BFF
    Beam_Filling_Factor = hotbar_integral/Gauss_integral
    
    return Beam_Filling_Factor
