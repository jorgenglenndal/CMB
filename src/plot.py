import matplotlib.pyplot as plt
import numpy as np

m           = 1.0
s           = 1.0
km          = 1e3 * m
Mpc         = 3.08567758e22 * m
H0_over_h   = 100# * km/s/Mpc
loaded_file = np.loadtxt('../cosmology.txt')
#print(loaded_file[:,0])


#plt.plot(loaded_file[:,0],loaded_file[:,7]+loaded_file[:,8])
#plt.plot(loaded_file[:,0],loaded_file[:,4]+loaded_file[:,5])
#plt.plot(loaded_file[:,0],loaded_file[:,6])
#plt.show()

#plt.semilogy(loaded_file[:,0],loaded_file[:,1]/3.08567758e22)
#plt.ylim(1e0,5*1e4)
#plt.xlim(-12,0)
#plt.show()

#plt.semilogy(loaded_file[:,0],loaded_file[:,2]/(100*1000)*3.08567758e22)
#plt.ylim(1e-1,1e+3)
#plt.xlim(-12,0)
#plt.show()

#plt.plot(loaded_file[:,0],loaded_file[:,1]*loaded_file[:,2]/2.99792458e8)
#plt.xlim(-14.5,0)
#plt.ylim(0.8,3.2)
#plt.show()


loaded_file2 = np.loadtxt('../data/supernovadata.txt')
#print(loaded_file2)
z_data = loaded_file2[:,0]
d_L_data = loaded_file2[:,1] # [Gpc]
error_data = loaded_file2[:,2] # [Gpc]

#plt.errorbar(z_data,d_L_data,error_data)
#plt.show()

def convert_x_to_z(x):
    return np.exp(-x)-1

#print(convert_x_to_z(loaded_file[:,0]))
#plt.plot(convert_x_to_z(loaded_file[:,0]),loaded_file[:,11]/(3.08567758e22*1e3))
#plt.xlim(0,1.4)
#plt.ylim(-1,10)
#plt.show()


'''
chi squared part
'''

loaded_fitting_file = np.loadtxt('../results.txt',skiprows=200)
#print(loaded_fitting_file[:,0])
argmin = np.argmin(loaded_fitting_file[:,0])
#print(loaded_fitting_file(argmin))
chi2_min,h,OmegaM,OmegaK = loaded_fitting_file[argmin]
all_OmegaM = loaded_fitting_file[:,2]
all_h = loaded_fitting_file[:,1]
one_sigma_OmegaM = all_OmegaM[loaded_fitting_file[:,0] < chi2_min + 3.53]
all_OmegaL = 1 -(loaded_fitting_file[:,3]+loaded_fitting_file[:,2]) # assumes Omega_(Radiation) = 0
one_sigma_OmegaL = all_OmegaL[loaded_fitting_file[:,0] < chi2_min + 3.53]

two_sigma_OmegaM = all_OmegaM[loaded_fitting_file[:,0] < chi2_min + 8.02]
two_sigma_OmegaL = all_OmegaL[loaded_fitting_file[:,0] < chi2_min + 8.02]

#plt.scatter(two_sigma_OmegaM,two_sigma_OmegaL)
#plt.scatter(one_sigma_OmegaM,one_sigma_OmegaL)
#plt.show()

'''
The standard deviation
'''
h_avg = np.average(loaded_fitting_file[:,1])
h_squared_sum = 0
for i in loaded_fitting_file[:,1]:
    h_squared_sum += (i-h_avg)**2
sigma_h = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*h_squared_sum)


'''
Not needed
H0_avg = h_avg*H0_over_h
H0_squared_sum = 0
for i in loaded_fitting_file[:,1]*H0_over_h:
    H0_squared_sum += (i-H0_avg)**2
sigma_H0 = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*H0_squared_sum)
'''


OmegaM_avg = np.average(loaded_fitting_file[:,2])
OmegaM_squared_sum = 0
for i in loaded_fitting_file[:,2]:
    OmegaM_squared_sum += (i-OmegaM_avg)**2
sigma_OmegaM = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaM_squared_sum)

OmegaK_avg = np.average(loaded_fitting_file[:,3])
OmegaK_squared_sum = 0
for i in loaded_fitting_file[:,3]:
    OmegaK_squared_sum += (i-OmegaK_avg)**2
sigma_OmegaK = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaK_squared_sum)


OmegaL_avg = np.average(all_OmegaL)
OmegaL_squared_sum = 0
for i in all_OmegaL:
    OmegaL_squared_sum += (i-OmegaL_avg)**2
sigma_OmegaL = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaL_squared_sum)
#print(OmegaM,sigma_OmegaM)

#plt.hist(all_OmegaL,density=True,bins=35)
#OmegaL_linspace = np.linspace(np.amin(all_OmegaL),np.amax(all_OmegaL),800)
#plt.plot(OmegaL_linspace,1/(sigma_OmegaL*np.sqrt(2*np.pi))*np.exp(-(OmegaL_linspace-OmegaL_avg)**2/(2*sigma_OmegaL**2)))
#plt.show()

plt.hist(all_h,density=True,bins=35) # easy to go from h to H0...
h_linspace = np.linspace(np.amin(all_h),np.amax(all_h),800)
plt.plot(h_linspace,1/(sigma_h*np.sqrt(2*np.pi))*np.exp(-(h_linspace-h_avg)**2/(2*sigma_h**2)))
plt.show()


'''
Not needed
plt.hist(all_h*H0_over_h,density=True,bins=35) # easy to go from h to H0...
H0_linspace = np.linspace(np.amin(all_h*H0_over_h),np.amax(all_h*H0_over_h),800)
plt.plot(H0_linspace,1/(sigma_H0*np.sqrt(2*np.pi))*np.exp(-(H0_linspace-H0_avg)**2/(2*sigma_H0**2)))
plt.show()
'''


plt.hist(all_h*H0_over_h,density=True,bins=35) # easy to go from h to H0...
#H0_linspace = np.linspace(np.amin(all_h*H0_over_h),np.amax(all_h*H0_over_h),800)
plt.plot(h_linspace*H0_over_h,1/(sigma_h*H0_over_h*np.sqrt(2*np.pi))*np.exp(-(h_linspace*H0_over_h-h_avg*H0_over_h)**2/(2*(sigma_h*H0_over_h)**2)))
plt.show()


# How to analyze the resulting chains:
# * Load the chains and skip the first few hundred samples (the burnin of the chains). E.g. loadtxt(file,skiprows=200) in python
# * Find the minimum chi2 and the corresponding best-fit parameters (you can use np.argmin to get index of the minvalue in python)
# * Select all samples of OmegaM and OmegaLambda (computed from OmegaM and OmegaK) that satisfy chi2 < chi2_min + 3.53 
#   (e.g. OmegaM[chi2 < chi2min + 3.53] in python)
# * Scatterplotting these gives you the 1sigma (68.4%) confidence region
# * Find the standard deviation of the samples to get the 1sigma confidence region of the parameters (assuming the posterior is a gaussian)
# * Make and plot a histogram of the samples for the different parameters (OmegaM, OmegaK, OmegaLambda, H0)
# * You can also compute the mean and standard deviation of the chain values and use this to overplot a gaussian with the same mean and variance for comparison.
#


#  0  fp << x                  << " ";
#  1  fp << eta_of_x(x)        << " ";
#  2  fp << Hp_of_x(x)         << " ";
#  3  fp << dHpdx_of_x(x)      << " ";
#  4  fp << get_OmegaB(x)      << " ";
#  5  fp << get_OmegaCDM(x)    << " ";
#  6  fp << get_OmegaLambda(x) << " ";
#  7  fp << get_OmegaR(x)      << " ";
#  8  fp << get_OmegaNu(x)     << " ";
#  9  fp << get_OmegaK(x)      << " ";

