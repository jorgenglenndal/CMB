import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker


m           = 1.0
s           = 1.0
km          = 1e3 * m
Mpc         = 3.08567758e22 * m
H0_over_h   = 100# * km/s/Mpc  [not SI]
c = 2.99792458e8 * m/s


#   0  fp << x                                        << " ";
#   1  fp << eta_of_x(x)                              << " ";
#   2  fp << Hp_of_x(x)                               << " ";
#   3  fp << dHpdx_of_x(x)                            << " ";
#   4  fp << get_OmegaB(x)                            << " ";
#   5  fp << get_OmegaCDM(x)                          << " ";
#   6  fp << get_OmegaLambda(x)                       << " ";
#   7  fp << get_OmegaR(x)                            << " ";
#   8  fp << get_OmegaNu(x)                           << " ";
#   9  fp << get_OmegaK(x)                            << " ";
#  10  fp << t_of_x(x)                                << " ";
#  11  fp << get_luminosity_distance_of_x(x)          << " ";
#  12  fp << ddHpddx_of_x(x)                          << " ";


loaded_file = np.loadtxt('cosmology.txt')
#print(loaded_file[:,0])
x = loaded_file[:,0]
eta_x = loaded_file[:,1]
Hp_x = loaded_file[:,2]
dHpdx_x = loaded_file[:,3]
OmegaB_x = loaded_file[:,4]
OmegaCDM_x = loaded_file[:,5]
OmegaL_x = loaded_file[:,6]
OmegaR_x = loaded_file[:,7]
OmegaNu_x = loaded_file[:,8]
OmegaK_x = loaded_file[:,9]
t_x = loaded_file[:,10]
luminosity_distance_x = loaded_file[:,11]
ddHpddx_of_x = loaded_file[:,12]


'''
Demonstrates that the code works
'''
#fig,ax = plt.subplots(3,1,figsize=(8, 16))
#ax[0].plot(x,1/Hp_x*dHpdx_x)
#ax[0].plot(x,np.ones(len(x))*-1,'--',label='Radiation Dominated')
#ax[0].plot(x,np.ones(len(x))*-1/2,'--',label='Matter Dominated')
#ax[0].plot(x,np.ones(len(x))*1,'--',label='Dark Energy Dominated')
#ax[0].legend()
#ax[0].set_title('$\\frac{1}{\mathcal{H}(x)} \\cdot   \\frac{d\mathcal{H}(x)}{dx}$')
#ax[1].plot(x,1/Hp_x*ddHpddx_of_x)
#ax[1].set_title('$\\frac{1}{\mathcal{H}(x)} \\cdot   \\frac{d^2\mathcal{H}(x)}{dx^2}$')
#ax[1].plot(x[(-20 < x) & (x< -4)],np.ones(len(x[(-20 < x) & (x < -4)]))*1,'--',label='Radiation Dominated')
#ax[1].plot(x[(-5 < x) & (x< 5)],np.ones(len(x[(-5 < x) & (x < 5)]))*1,'--',label='Dark Energy Dominated')
#ax[1].plot(x,np.ones(len(x))*1/4,'--',label='Matter Dominated')
#ax[1].legend()
#ax[2].plot(x,eta_x*Hp_x/c)
#ax[2].set_title('$\\frac{\eta (x) \mathcal{H}(x)}{c}$')
#ax[2].set_ylim(0.8,3.2)
#ax[2].set_xlim(-14.5,0)
#ax[2].plot(x,np.ones(len(x))*1,'--',label='Radiation Dominated')
#ax[2].legend()
#ax[2].set_xlabel('x')
#plt.tight_layout()
##plt.savefig('demonstrate_that_code_works.pdf')
#plt.show()

'''
The Omegas
'''
#fig,ax = plt.subplots()
#ax.set_title('Evolution of $\Omega$')
#ax.plot(x,OmegaR_x + OmegaNu_x,label='$\Omega_R$ Radiation')
#ax.plot(x,OmegaB_x+OmegaCDM_x,label='$\Omega_M$ Matter')
#ax.plot(x,OmegaL_x,label='$\Omega_\Lambda$ Dark Energy')
#ax.set_xlabel('x')
#plt.legend()
##plt.savefig('Omegas.pdf')
#plt.show()

#plt.title('$t(x)\ \ \ (yr) $')
#plt.semilogy(x[(-18 < x) & (x<=0)],t_x[(-18 < x) & (x<=0)]/(60*60*24*356))
#plt.xlabel('x')
##plt.savefig('t_x.pdf')
#plt.show()


#fig,ax = plt.subplots()
#plt.title('$\\frac{\eta (x)}{c}\ \ \ (\mathrm{yr})$')
#plt.semilogy(x,eta_x/(c*60*60*24*365))
#plt.xlim(-18,0)
#plt.xlabel('x')
##plt.savefig('eta_x.pdf')
#plt.show()

#plt.title('$\mathcal{H}(x)\ \ \ \\left(\\frac{100\ km/s)}{Mpc}\\right)$')
#plt.semilogy(x,Hp_x/(100*km)*Mpc)
#plt.ylim(1e-1,1e+3)
#plt.xlim(-12,0)
#plt.xlabel('x')
##plt.savefig('Hp_x.pdf')
#plt.show()


loaded_file2 = np.loadtxt('data/supernovadata.txt')
#print(loaded_file2)
z_data = loaded_file2[:,0]
d_L_data = loaded_file2[:,1] # [Gpc]
error_data = loaded_file2[:,2] # [Gpc]



def convert_x_to_z(x):
    return np.exp(-x)-1

#fig,ax = plt.subplots()
#fig.suptitle('$\\frac{d_L}{z}\ \ \ (\mathrm{Gpc})$',fontsize=16)
#ax.semilogx(convert_x_to_z(x),luminosity_distance_x/(convert_x_to_z(x)*(Mpc*1e3)),label='Theoretical $d_L$')
#ax.errorbar(z_data,d_L_data/z_data,error_data/z_data,fmt='none',ecolor='r',capsize=5,label='Data')
#ax.set_xlim(0.005,2.26)
#ax.set_ylim(3.5,8)
#ax.set_xlabel('z')
#ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
#ax.legend()
##plt.savefig('dataplot.pdf')
#plt.show()


'''
chi squared part
'''
loaded_fitting_file = np.loadtxt('results.txt',skiprows=200)
#print(loaded_fitting_file[:,0])
argmin = np.argmin(loaded_fitting_file[:,0])
#print(loaded_fitting_file(argmin))
chi2_min,h_fit,OmegaM_fit,OmegaK_fit = loaded_fitting_file[argmin]
chains_chi2 = loaded_fitting_file[:,0]
chains_OmegaM = loaded_fitting_file[:,2]
chains_OmegaK = loaded_fitting_file[:,3]
chains_h = loaded_fitting_file[:,1]
one_sigma_OmegaM = chains_OmegaM[chains_chi2 < chi2_min + 3.53]
chains_OmegaL = 1 -(chains_OmegaK+chains_OmegaM) # assumes Omega_(Radiation) = 0
one_sigma_OmegaL = chains_OmegaL[chains_chi2 < chi2_min + 3.53]

two_sigma_OmegaM = chains_OmegaM[chains_chi2 < chi2_min + 8.02]
two_sigma_OmegaL = chains_OmegaL[chains_chi2 < chi2_min + 8.02]

#plt.scatter(two_sigma_OmegaM,two_sigma_OmegaL,label='$2\sigma$ constraint')
#plt.scatter(one_sigma_OmegaM,one_sigma_OmegaL,label='$1\sigma$ constraint')
#plt.legend()
#plt.xlabel('$\Omega_M$')
#plt.ylabel('$\Omega_\Lambda$')
##plt.savefig('constraints.pdf')
#plt.show()

'''
The standard deviation
'''
h_avg = np.average(chains_h)
h_squared_sum = 0
for i in chains_h:
    h_squared_sum += (i-h_avg)**2
sigma_h = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*h_squared_sum)


OmegaM_avg = np.average(chains_OmegaM)
OmegaM_squared_sum = 0
for i in chains_OmegaM:
    OmegaM_squared_sum += (i-OmegaM_avg)**2
sigma_OmegaM = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaM_squared_sum)


OmegaK_avg = np.average(chains_OmegaK)
OmegaK_squared_sum = 0
for i in chains_OmegaK:
    OmegaK_squared_sum += (i-OmegaK_avg)**2
sigma_OmegaK = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaK_squared_sum)


OmegaL_avg = np.average(chains_OmegaL)
OmegaL_squared_sum = 0
for i in chains_OmegaL:
    OmegaL_squared_sum += (i-OmegaL_avg)**2
sigma_OmegaL = np.sqrt(1/(len(loaded_fitting_file[:,0])-1)*OmegaL_squared_sum)


#OmegaL_linspace = np.linspace(np.amin(chains_OmegaL),np.amax(chains_OmegaL),800)
#plt.hist(chains_OmegaL,density=True,bins=35)
#plt.plot(OmegaL_linspace,1/(sigma_OmegaL*np.sqrt(2*np.pi))*np.exp(-(OmegaL_linspace-OmegaL_avg)**2/(2*sigma_OmegaL**2)))
#plt.show()


h_linspace = np.linspace(np.amin(chains_h),np.amax(chains_h),800)
#plt.title('Probability density $p(H_0) \ \ \ (100\ \mathrm{km}/\mathrm{s}/\mathrm{Mpc})^{-1}$')
#plt.hist(chains_h,density=True,bins=35,alpha=0.6)
#plt.plot(h_linspace,1/(sigma_h*np.sqrt(2*np.pi))*np.exp(-(h_linspace-h_avg)**2/(2*sigma_h**2)),linewidth=3,color='r')
#plt.xlabel('$H_0$  (100 km/s/Mpc)')
##plt.savefig('prob_dist_H0.pdf')
#plt.show()


'''
milsestone 2
'''

rec_loaded = np.loadtxt('recombination.txt')

x_rec = rec_loaded[:,0]
Xe_x = rec_loaded[:,1]
ne_x = rec_loaded[:,2]
tau_x = rec_loaded[:,3]
dtaudx_x = rec_loaded[:,4]
ddtau_ddx_x = rec_loaded[:,5]
g_tilde_x = rec_loaded[:,6]
dg_tildedx_x = rec_loaded[:,7]
ddg_tildeddx_x = rec_loaded[:,8]

plt.semilogy(x_rec,Xe_x)
plt.xlim(-12,0)
plt.show()


plt.semilogy(x_rec,tau_x)
plt.semilogy(x_rec,-dtaudx_x)
plt.semilogy(x_rec,ddtau_ddx_x)
plt.xlim(-12,0)
plt.ylim(10**-8,10**8)
plt.show()

plt.plot(x_rec,g_tilde_x)
plt.xlim(-12,0)
plt.show()

plt.plot(x_rec,dg_tildedx_x)
plt.xlim(-12,0)
plt.show()

plt.plot(x_rec,ddg_tildeddx_x)
plt.xlim(-12,0)
plt.show()

g_tilde_test = np.trapz(g_tilde_x,x_rec,dx=(x_rec[1]-x_rec[0]))
#print(g_tilde_test)

g_tilde_max = np.amax(g_tilde_x)
#print(g_tilde_max)
x_decoupling = x_rec[np.argmax(g_tilde_x)]
#print(x_decoupling)

#fp << x                    << " ";
#fp <<std::setprecision(15) <<Xe_of_x(x)           << " ";
#fp << ne_of_x(x)           << " ";
#fp << tau_of_x(x)          << " ";
#fp << dtaudx_of_x(x)       << " ";
#fp << ddtauddx_of_x(x)     << " ";
#fp << g_tilde_of_x(x)      << " ";
#fp << dgdx_tilde_of_x(x)   << " ";
#fp << ddgddx_tilde_of_x(x) << " ";
#fp << "\n";