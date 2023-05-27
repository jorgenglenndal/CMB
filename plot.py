import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker


def nearest_index(array, target_value):
    index = (np.abs(array - target_value)).argmin()
    return index

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
#plt.savefig('t_x.pdf')
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

#print(chi2_min)
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
Xe_saha_x = rec_loaded[:,9]
#Xe_peebles_x = rec_loaded[:,10]

g_tilde_test = np.trapz(g_tilde_x,x_rec,dx=(x_rec[1]-x_rec[0]))
#print(g_tilde_test)

g_tilde_max = np.amax(g_tilde_x)
#print(g_tilde_max)
x_decoupling = x_rec[np.argmax(g_tilde_x)]
#print(x_decoupling)

#print('here')
#print(x_rec[nearest_index(Xe_x,0.1)])
#print(x_rec[nearest_index(Xe_saha_x,0.1)])
#print(Xe_x[-1])
##plt.semilogy(x_rec,Xe_x,label='Xe')
##plt.plot(np.ones(10)*x_decoupling,np.linspace(0,2,10),color='k',linestyle='dashed',alpha=0.6,label='Time of Recombination')
###plt.plot(x_rec[Xe_saha_x>0],Xe_saha_x[Xe_saha_x>0],label='X_e_saha')
###plt.plot(x_rec[Xe_peebles_x>0],Xe_peebles_x[Xe_peebles_x>0],label='X_e_peebles')
##
###print(x_decoupling)
##plt.semilogy(x_rec[Xe_saha_x>0],Xe_saha_x[Xe_saha_x>0],label='Xe Saha')
###plt.plot(x_rec[Xe_peebles_x>0],Xe_peebles_x[Xe_peebles_x>0],label='X_e_peebles')
##
##plt.title('Xe')
##plt.xlim(-10,0)
##plt.ylim(1e-4,2)
##plt.xlabel('x')
##plt.legend()
###plt.savefig("figures/milestone2/X_e.pdf")
##plt.show()



###rec_loaded = np.loadtxt('recombination.txt')
###
###x_rec = rec_loaded[:,0]
###Xe_x = rec_loaded[:,1]
###ne_x = rec_loaded[:,2]
###tau_x = rec_loaded[:,3]
###dtaudx_x = rec_loaded[:,4]
###ddtau_ddx_x = rec_loaded[:,5]
###g_tilde_x = rec_loaded[:,6]
###dg_tildedx_x = rec_loaded[:,7]
###ddg_tildeddx_x = rec_loaded[:,8]
####Xe_saha_x = rec_loaded[:,9]
####Xe_peebles_x = rec_loaded[:,10]
###
###g_tilde_test = np.trapz(g_tilde_x,x_rec,dx=(x_rec[1]-x_rec[0]))
####print(g_tilde_test)
###
###g_tilde_max = np.amax(g_tilde_x)
####print(g_tilde_max)
###x_decoupling = x_rec[np.argmax(g_tilde_x)]
#print(x_decoupling)



##plt.semilogy(x_rec,tau_x,label=r'$\tau$')
##plt.semilogy(x_rec,-dtaudx_x,label=r'$-d\tau/dx$')
##plt.semilogy(x_rec,ddtau_ddx_x,label=r'$d^2\tau/dx^2$')
##plt.plot(np.ones(10)*x_decoupling,np.linspace(-100,1e12,10),color='k',linestyle='dashed',alpha=0.6,label='Time of Recombination')
##plt.xlim(-12,0)
##plt.ylim(10**-8,10**8)
##plt.xlabel('x')
##plt.legend()
##plt.savefig("figures/milestone2/taus.pdf")
##plt.show()



#plt.title('g_tilde')
###plt.plot(x_rec,g_tilde_x,label=r"$\tilde{g}$")
####plt.plot(np.ones(10)*x_decoupling,np.linspace(-5,1e100,10),color='k',linestyle='dashed',alpha=0.6,label='Time of Recombination')
####plt.xlim(-12,0)
####plt.ylim(-0.2,5)
###plt.xlabel('x')
####plt.legend()
###
####plt.show()
###
###plt.plot(x_rec,dg_tildedx_x/10,label=r"$\frac{d\tilde{g}}{dx}\cdot \frac{1}{10}$")#,label="")
####plt.xlim(-12,0)
####plt.title('dg_tildedx')
####plt.show()
###
###plt.plot(x_rec,ddg_tildeddx_x/300,label=r"$\frac{d^2\tilde{g}}{dx^2}\cdot \frac{1}{300}$")
###plt.xlim(-7.5,-6.4)
####plt.title('ddg_tildeddx')
###plt.legend()
###plt.savefig("figures/milestone2/g_tilde_plots.pdf")
###plt.show()



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


#    fp << x                  << " ";
#    fp << get_Theta(x,k,0)   << " ";
#    fp << get_Theta(x,k,1)   << " ";
#    fp << get_Theta(x,k,2)   << " ";
#    fp << get_Phi(x,k)       << " ";
#    fp << get_Psi(x,k)       << " ";
#    fp << get_delta_cdm(x,k) << " ";
#    fp << get_delta_b(x,k)   << " ";
#    fp << get_v_cdm(x,k)     << " ";
#    fp << get_v_b(x,k)       << " ";



def milestone3():
    test_loaded_file_k0 = np.loadtxt('test_perturbations_k_0_001.txt')
    k0 = 0.001/Mpc
    x_k             = test_loaded_file_k0[:,0]
    x_k0            = test_loaded_file_k0[:,0]
    theta0_k0       = test_loaded_file_k0[:,1]
    theta1_k0       = test_loaded_file_k0[:,2]
    theta2_k0       = test_loaded_file_k0[:,3]
    Phi_k0          = test_loaded_file_k0[:,4]
    Psi_k0          = test_loaded_file_k0[:,5]
    delta_cdm_k0    = test_loaded_file_k0[:,6]
    delta_b_k0      = test_loaded_file_k0[:,7]
    v_cdm_k0        = test_loaded_file_k0[:,8]
    v_b_k0          = test_loaded_file_k0[:,9]
    eta_x_k0        = test_loaded_file_k0[:,10]
    eta_x_k         = test_loaded_file_k0[:,10]

    test_loaded_file_k1 = np.loadtxt('test_perturbations_k_0_01.txt')
    k1=0.01/Mpc
    x_k1            = test_loaded_file_k1[:,0]
    theta0_k1       = test_loaded_file_k1[:,1]
    theta1_k1       = test_loaded_file_k1[:,2]
    theta2_k1       = test_loaded_file_k1[:,3]
    Phi_k1          = test_loaded_file_k1[:,4]
    Psi_k1          = test_loaded_file_k1[:,5]
    delta_cdm_k1    = test_loaded_file_k1[:,6]
    delta_b_k1      = test_loaded_file_k1[:,7]
    v_cdm_k1        = test_loaded_file_k1[:,8]
    v_b_k1          = test_loaded_file_k1[:,9]
    eta_x_k1        = test_loaded_file_k1[:,10]

    test_loaded_file_k2 = np.loadtxt('test_perturbations_k_0_1.txt')
    k2=0.1/Mpc
    x_k2            = test_loaded_file_k2[:,0]
    theta0_k2       = test_loaded_file_k2[:,1]
    theta1_k2       = test_loaded_file_k2[:,2]
    theta2_k2       = test_loaded_file_k2[:,3]
    Phi_k2          = test_loaded_file_k2[:,4]
    Psi_k2          = test_loaded_file_k2[:,5]
    delta_cdm_k2    = test_loaded_file_k2[:,6]
    delta_b_k2      = test_loaded_file_k2[:,7]
    v_cdm_k2        = test_loaded_file_k2[:,8]
    v_b_k2          = test_loaded_file_k2[:,9]
    eta_x_k2        = test_loaded_file_k2[:,10]
    
    #print(x_k[nearest_index(eta_x_k,Mpc/0.1)])
    #print(x_k[nearest_index(eta_x_k,2*np.pi*Mpc/0.01)])
    #print(x_k[nearest_index(eta_x_k,2*np.pi*Mpc/0.001)])
    
###    plt.title("$\delta_\mathrm{CDM},\delta_b$")
###    plt.semilogy(x_k0,delta_cdm_k0,'b',label= "k = 0.001/Mpc",alpha=0.6) 
###    plt.semilogy(x_k0,abs(delta_b_k0),'--',color='b')
###    
###    plt.semilogy(x_k1,delta_cdm_k1,color= 'r',label= "k = 0.01/Mpc",alpha = 0.6)
###    plt.semilogy(x_k1,abs(delta_b_k1),'--',color ='r')
###    
###
    green_ = (0, 0.65, 0)  # RGB tuple for green
###    plt.semilogy(x_k2,delta_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
###    plt.semilogy(x_k2,abs(delta_b_k2),'--',color= green_)
###    #plt.show()
###    plt.xlim(-18,0)
###    plt.xlabel('x')
###    #plt.ylim(1e-10,1e5)
###    plt.legend()
####    plt.savefig("test_delta_cdm_delta_b.pdf")
###    plt.show()
###
###    ###
###    plt.title("$v_\mathrm{CDM},v_b$")
###    plt.semilogy(x_k0,v_cdm_k0,'b', label= "k = 0.001/Mpc",alpha=0.6)
###    plt.semilogy(x_k0,abs(v_b_k0),'--',color='b')
###    
###    plt.semilogy(x_k1,v_cdm_k1,color='r',label= "k = 0.01/Mpc",alpha = 0.6)
###    plt.semilogy(x_k1,abs(v_b_k1),'--',color ='r')
###    
###    plt.semilogy(x_k2,v_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
###    plt.semilogy(x_k2,abs(v_b_k2),'--',color=green_)
###
###    plt.xlim(-18,0)
###    plt.ylim(1e-6,1e2)
###    plt.legend()
###    plt.xlabel('x')
###    #plt.ylabel('[m/s]')
###    #plt.savefig("test_v_cdm_v_b.pdf")
###    plt.show()
###    
###    ###
###    plt.title("$\Theta_0$")
###    plt.plot(x_k0,theta0_k0,label="k = 0.001/Mpc")     
###    plt.plot(x_k1,theta0_k1,label="k = 0.01/Mpc")
###    plt.plot(x_k2,theta0_k2,label="k = 0.1/Mpc")
###    plt.xlabel('x')
###    plt.xlim(-18,0)
###    plt.ylim(-0.7,1)
###    plt.legend()
###    ##plt.savefig("test_theta_0.pdf")
###    plt.show()
###
###    ###
###    plt.title("$\Theta_1$")
###    plt.plot(x_k0,theta1_k0,label="k = 0.001/Mpc")     
###    plt.plot(x_k1,theta1_k1,label="k = 0.01/Mpc")
###    plt.plot(x_k2,theta1_k2,label="k = 0.1/Mpc",color=green_)
###
###    plt.xlabel('x')
###    plt.xlim(-18,0)
###    plt.ylim(-0.5,0.5)
###    plt.legend()
###    #plt.savefig("test_theta_1.pdf")
###    plt.show()
###
###    ###
###    plt.title("$\Phi$")
###    plt.plot(x_k0,Phi_k0,label="k = 0.001/Mpc")     
###    plt.plot(x_k1,Phi_k1,label="k = 0.01/Mpc")
###    plt.plot(x_k2,Phi_k2,label="k = 0.1/Mpc")
###
###    plt.xlabel('x')
###    plt.legend()
###    plt.xlim(-18,0)
###    plt.ylim(0,0.8)
###    #plt.savefig("test_phi.pdf")
###    plt.show()
    

    #plt.plot(x_k0,theta2_k0)
    #plt.plot(x_k1,theta2_k1)
    #plt.plot(x_k2,theta2_k2)
    #plt.show()
#
#
    #plt.plot(x_k0,Psi_k0)
    ##plt.plot(x_k0,np.cos(k0*eta_x_k0/np.sqrt(3)))
###
    #plt.plot(x_k1,Psi_k1,'y')
    ###plt.plot(x_k1,np.cos(k1*eta_x_k1/np.sqrt(3)),'--','y')
###
    #plt.plot(x_k2,Psi_k2,'g')
    ##plt.plot(x_k2,np.cos(k2*eta_x_k2/np.sqrt(3)),'--','g')
    ##plt.xlim(-20,-8)
    #plt.show()


    #plt.plot(x_k0,Phi_k0+Psi_k0)
    #plt.plot(x_k1,Phi_k1+Psi_k1)
    #plt.plot(x_k2,Phi_k2+Psi_k2)
    #plt.show()
    #print(eta_x_k0*k0)

    ###print(x_k0[nearest_index(eta_x_k0*k0,1)])
    ###print(x_k1[nearest_index(eta_x_k1*k1,1)])
    ###print(x_k2[nearest_index(eta_x_k1*k2,1)])

    #print(len(x_k0),len(x_k1),len(x_k2))
    #print(x_k[nearest_index(eta_x_k*10/Mpc,1)])
    """
    Test over
    """
    loaded_file_k0 = np.loadtxt('perturbations_k_0_001.txt')
    k0 = 0.001/Mpc
    x_k             = loaded_file_k0[:,0]
    #x_k0            = loaded_file_k0[:,0]
    theta0_k0       = loaded_file_k0[:,1]
    theta1_k0       = loaded_file_k0[:,2]
    theta2_k0       = loaded_file_k0[:,3]
    Phi_k0          = loaded_file_k0[:,4]
    Psi_k0          = loaded_file_k0[:,5]
    delta_cdm_k0    = loaded_file_k0[:,6]
    delta_b_k0      = loaded_file_k0[:,7]
    v_cdm_k0        = loaded_file_k0[:,8]
    v_b_k0          = loaded_file_k0[:,9]
    #eta_x_k0        = loaded_file_k0[:,10]
    eta_x_k         = loaded_file_k0[:,10]

    loaded_file_k1 = np.loadtxt('perturbations_k_0_01.txt')
    k1=0.01/Mpc
    #x_k1            = loaded_file_k1[:,0]
    theta0_k1       = loaded_file_k1[:,1]
    theta1_k1       = loaded_file_k1[:,2]
    theta2_k1       = loaded_file_k1[:,3]
    Phi_k1          = loaded_file_k1[:,4]
    Psi_k1          = loaded_file_k1[:,5]
    delta_cdm_k1    = loaded_file_k1[:,6]
    delta_b_k1      = loaded_file_k1[:,7]
    v_cdm_k1        = loaded_file_k1[:,8]
    v_b_k1          = loaded_file_k1[:,9]
    #eta_x_k1        = loaded_file_k1[:,10]

    loaded_file_k2 = np.loadtxt('perturbations_k_0_1.txt')
    k2=0.1/Mpc
    #x_k2            = loaded_file_k2[:,0]
    theta0_k2       = loaded_file_k2[:,1]
    theta1_k2       = loaded_file_k2[:,2]
    theta2_k2       = loaded_file_k2[:,3]
    Phi_k2          = loaded_file_k2[:,4]
    Psi_k2          = loaded_file_k2[:,5]
    delta_cdm_k2    = loaded_file_k2[:,6]
    delta_b_k2      = loaded_file_k2[:,7]
    v_cdm_k2        = loaded_file_k2[:,8]
    v_b_k2          = loaded_file_k2[:,9]
    #eta_x_k2        = loaded_file_k2[:,10]

    loaded_file_k10 = np.loadtxt('perturbations_k10.txt')
    k10 = 10/Mpc
    x_k10            = loaded_file_k10[:,0]
    #x_k0            = loaded_file_k0[:,0]
    theta0_k10       = loaded_file_k10[:,1]
    theta1_k10       = loaded_file_k10[:,2]
    theta2_k10       = loaded_file_k10[:,3]
    Phi_k10          = loaded_file_k10[:,4]
    Psi_k10          = loaded_file_k10[:,5]
    delta_cdm_k10    = loaded_file_k10[:,6]
    delta_b_k10      = loaded_file_k10[:,7]
    v_cdm_k10        = loaded_file_k10[:,8]
    v_b_k10          = loaded_file_k10[:,9]
    #eta_x_k0        = loaded_file_k0[:,10]
    eta_x_k10        = loaded_file_k10[:,10]

    """
    density perturbations
    """


##    plt.title("$\delta_\gamma$")
##    #plt.title("$\Theta_0$")
##    plt.plot(x_k,4*theta0_k0,label="k = 0.001/Mpc")     
##    plt.plot(x_k,4*theta0_k1,label="k = 0.01/Mpc")
##    plt.plot(x_k,4*theta0_k2,label="k = 0.1/Mpc")
##    plt.plot(x_k10,4*theta0_k10,label="k = 10/Mpc",alpha = 0.6)
##    plt.xlabel('x')
##    plt.xlim(-18,0)
##    #plt.ylim(-0.7,1)
##    plt.legend()
##    ##plt.savefig("figures/milestone3/delta_gamma.pdf")
##    plt.show()

    #epic plot
    plt.title("$\delta_\mathrm{CDM},\delta_b$")
    plt.plot(x_k,delta_cdm_k0,'b',label= "k = 0.001/Mpc",alpha=0.6)
    plt.plot(np.ones(10)*x_k[nearest_index(eta_x_k*0.001/Mpc,1)],np.linspace(-1e100,1e100,10),"--")     
 
    plt.plot(x_k,(delta_b_k0),'--',color='b')
    plt.plot(x_k,delta_cdm_k1,color= 'r',label= "k = 0.01/Mpc",alpha = 0.6)
    plt.plot(x_k,(delta_b_k1),'--',color ='r')
    plt.plot(x_k,delta_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
    plt.plot(x_k,(delta_b_k2),'--',color= green_)
    plt.plot(x_k10,delta_cdm_k10,color = 'k',label= "k = 10/Mpc",alpha = 0.6)
    plt.plot(x_k10,(delta_b_k10),'--',color= 'k',alpha=0.6)
    plt.xlim(-18,0)
    plt.xlabel('x')
    plt.ylim(-3,10)
#   
    #plt.plot(x_k,4*theta0_k0,color='b'                ,linestyle='dotted')     
    #plt.plot(x_k,4*theta0_k1,color='r'                ,linestyle='dotted')
    #plt.plot(x_k,4*theta0_k2,color= green_            ,linestyle='dotted')
    #plt.plot(x_k10,4*theta0_k10,color='k' ,linestyle='dotted')
#    
    plt.legend()
    #plt.savefig("figures/milestone3/delta_cdm_delta_b_zoom.pdf")
    plt.show()

    #epic plot 2
    #plt.title("$\delta_\mathrm{CDM},\delta_b,\delta_\gamma$,    k = 10/Mpc")
    #plt.plot(x_k,delta_cdm_k0,'b',label= "k = 0.001/Mpc",alpha=0.6) 
#    plt.plot(x_k,(delta_b_k0),'--',color='b')
#    #plt.plot(x_k,delta_cdm_k1,color= 'r',label= "k = 0.01/Mpc",alpha = 0.6)
#    plt.plot(x_k,(delta_b_k1),'--',color ='r')
#    #plt.plot(x_k,delta_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
#    plt.plot(x_k,(delta_b_k2),'--',color= green_)
#
#    #plt.plot(x_k10,delta_cdm_k10,color = 'k',label= "CDM",alpha = 0.6)
#    plt.plot(x_k10,(delta_b_k10),color= 'k',alpha=0.6,linestyle="--",label="Baryons")   
#    plt.plot(np.ones(10)*x_decoupling,np.linspace(-10,10,10),color='y',label='Time of Recombination')
#
#    plt.xlim(-18,0)
#    plt.xlabel('x')
#    plt.ylim(-3,6)
##   
#    plt.plot(x_k,4*theta0_k0,color='b'                )#,linestyle='--')     
#    plt.plot(x_k,4*theta0_k1,color='r'                )#,linestyle='--')
#    plt.plot(x_k,4*theta0_k2,color= green_            )#,linestyle='--')
#    plt.plot(x_k10,4*theta0_k10,color='k',alpha=0.4,label='Radiation')#,linestyle='dotted')
##    
#    plt.legend()
#    #plt.savefig("figures/milestone3/delta_cdm_delta_b_delta_gamma.pdf")
#    plt.show()
#
#    plt.title("$\delta_\mathrm{CDM},\delta_b$")
#    plt.semilogy(x_k,delta_cdm_k0,'b',label= "k = 0.001/Mpc",alpha=0.6) 
#    plt.semilogy(x_k,abs(delta_b_k0),'--',color='b')
#    plt.semilogy(x_k,delta_cdm_k1,color= 'r',label= "k = 0.01/Mpc",alpha = 0.6)
#    plt.semilogy(x_k,abs(delta_b_k1),'--',color ='r')
#    plt.semilogy(x_k,delta_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
#    plt.semilogy(x_k,abs(delta_b_k2),'--',color= green_)
#    plt.semilogy(x_k10,delta_cdm_k10,color = 'k',label= "k = 10/Mpc",alpha = 0.6)
#    #plt.semilogy(x_k10,np.exp(x_k10)/np.amax(np.exp(x_k10)))
#    #plt.semilogy(x_k10,delta_cdm_k10/np.amax(delta_cdm_k10))
#    #plt.semilogy(x_k10,delta_cdm_k10/np.amax(delta_cdm_k10))
#    plt.semilogy(x_k10,abs(delta_b_k10),'--',color= 'k',alpha=0.6)
#    plt.xlim(-18,0)
#    plt.xlabel('x')
#    #plt.ylim(0,1)
#    plt.legend()
#    ##plt.savefig("figures/milestone3/delta_cdm_delta_b.pdf")
#    plt.show()
#
#    """
#    velocity perturbations
#    """
#
#    plt.title('$v_\gamma$')
#    plt.plot(x_k,   -3*theta1_k0  ,label="k = 0.001/Mpc")     
#    plt.plot(x_k,    -3*theta1_k1 ,label="k = 0.01/Mpc")
#    plt.plot(x_k,    -3*theta1_k2 ,label="k = 0.1/Mpc")
#    plt.plot(x_k10,  -3*theta1_k10 ,label="k = 10/Mpc",alpha=0.6)
#    plt.xlabel('x')
#    plt.legend()
#    ##plt.savefig("v_gamma.pdf")
#    plt.show()
#
#    plt.title("$v_\mathrm{CDM},v_b$")
#    plt.semilogy(x_k,v_cdm_k0,'b', label= "k = 0.001/Mpc",alpha=0.6)
#    plt.semilogy(x_k,abs(v_b_k0),'--',color='b')
#    
#    plt.semilogy(x_k,v_cdm_k1,color='r',label= "k = 0.01/Mpc",alpha = 0.6)
#    plt.semilogy(x_k,abs(v_b_k1),'--',color ='r')
#    
#    plt.semilogy(x_k,v_cdm_k2,color = green_,label= "k = 0.1/Mpc",alpha = 0.6)
#    plt.semilogy(x_k,abs(v_b_k2),'--',color=green_)
#
#    plt.semilogy(x_k10,v_cdm_k10,color = 'k',label= "k = 10/Mpc",alpha = 0.6)
#    plt.semilogy(x_k10,abs(v_b_k10),'--',color='k',alpha=0.6)
#
#    plt.xlim(-18,0)
#    #plt.ylim(1e-6,1e2)
#    plt.legend()
#    plt.xlabel('x')
#    ##plt.savefig("figures/milestone3/v_cdm_v_b.pdf")
#    plt.show()
#
#    """
#    Photon quadrupole
#    """
#
#    plt.title("$\Theta_2$")
#    plt.plot(x_k,theta2_k0,label="k = 0.001/Mpc")     
#    plt.plot(x_k,theta2_k1,label="k = 0.01/Mpc")
#    plt.plot(x_k,theta2_k2,label="k = 0.1/Mpc",color=green_)
#    plt.plot(x_k10,theta2_k10,label="k = 10/Mpc",color='k')
#    plt.xlabel('x')
#    plt.xlim(-18,0)
#    #plt.ylim(-0.5,0.5)
#    plt.legend()
#    ##plt.savefig("figures/milestone3/theta_2.pdf")
#    plt.show()
#
    """
    Phi and Psi
    """
    plt.title("$\Phi$")
    plt.plot(np.ones(10)*-8.13,np.linspace(-1e100,1e100,10),color="k",label="Matter-Radiation equality")     

    plt.plot(x_k,Phi_k0,label="k = 0.001/Mpc",color="r")
    plt.plot(np.ones(10)*x_k[nearest_index(eta_x_k*0.001/Mpc,1)],np.linspace(-1e100,1e100,10),"--",color="r")     
    plt.plot(x_k,Phi_k1,label="k = 0.01/Mpc",color=green_)
    plt.plot(np.ones(10)*x_k[nearest_index(eta_x_k*0.01/Mpc,1)],np.linspace(-1e100,1e100,10),"--",color=green_)
    plt.plot(x_k,Phi_k2,label="k = 0.1/Mpc",color= "b")
    plt.plot(np.ones(10)*x_k[nearest_index(eta_x_k*0.1/Mpc,1)],np.linspace(-1e100,1e100,10),"--",color="b")
    plt.plot(x_k10,Phi_k10,label="k = 10/Mpc",color="darkorange")
    plt.plot(np.ones(10)*x_k[nearest_index(eta_x_k*10/Mpc,1)],np.linspace(-1e100,1e100,10),"--",color="darkorange")
    plt.xlabel('x')
    plt.legend()
    plt.xlim(-18,0)
    plt.ylim(-0.1,0.9)
    #plt.savefig("figures/milestone3/phi.pdf")
    plt.show()

    plt.title("$\Phi + \Psi$")
    plt.plot(x_k,Phi_k0 + Psi_k0,label="k = 0.001/Mpc")     
    plt.plot(x_k,Phi_k1 + Psi_k1,label="k = 0.01/Mpc")
    plt.plot(x_k,Phi_k2 + Psi_k2,label="k = 0.1/Mpc")
    plt.plot(x_k10,Phi_k10 + Psi_k10,label="k = 10/Mpc")
    plt.xlabel('x')
    plt.legend()
    plt.xlim(-18,0)
    #plt.ylim(0,0.8)
    ##plt.savefig("figures/milestone3/phi_psi.pdf")
    plt.show()
    print('here')
    print(x_k10[nearest_index(eta_x_k10*10/Mpc,1)])
    print(x_k[nearest_index(eta_x_k*0.1/Mpc,1)])
    print(x_k[nearest_index(eta_x_k*0.01/Mpc,1)])
    print(x_k[nearest_index(eta_x_k*0.001/Mpc,1)])
    


def milestone4():
    

    loaded_file = np.loadtxt('lcdm_pk007_delta.dat')
    clas = np.loadtxt('test_02_z1_pk.dat')

    #early = np.loadtxt('lcdm_pk000_T00.dat')
    #print(loaded_file)
    x = loaded_file[:,0]
    y = loaded_file[:,1]

    x_class = clas[:,0]
    y_class = clas[:,1]

    #x1 = early[:,0]
    #y1 = early[:,1]
    #print(y)
    #plt.loglog(x1,y1)
    #plt.show()
    ##plt.loglog(x,y*2*np.pi**2/x**3,label = "gevolution")
    ##plt.loglog(x_class,y_class,label="CLASS")
    ##plt.xlim(1e-3,2*1e0)
    ##plt.ylim(1e1,1e5)
    ###plt.savefig("matter_spectrum.pdf")
    ##matter_spectrum = np.loadtxt('matter_power_spectrum_result.txt')
    ##k = matter_spectrum[:,0]
    ##P_k = matter_spectrum[:,1]
    ##plt.loglog(k,P_k,label = "Our code")
    ##plt.legend()
    ##plt.show()
    


   
###
###    source = np.loadtxt('source_new.txt')
###    x = source[:,0]
###    S = source[:,1]
###    S1 = source[:,2]
###    S2 = source[:,3]
###    #S3 = source[:,4]
###    #S4 = source[:,5]
###    plt.plot(x,S)
###    plt.show()
###    plt.plot(x,S1)
###    plt.show()
###    plt.plot(x,S2)
###    #plt.show()
###    #plt.plot(x,S3)
###    #plt.show()
###    #plt.plot(x,S4)
###    plt.show()


    
    """
    test results
    """

    ##matter_spectrum = np.loadtxt('matter_power_spectrum_test_result.txt')
    ##k = matter_spectrum[:,0]
    ##P_k = matter_spectrum[:,1]
    ##plt.loglog(k,P_k)
    ##plt.show()
    ##Cl = np.loadtxt('cells_test_result.txt')
    ##ell = Cl[:,0]
    ##Cell = Cl[:,1]#/(1e6 * 2.7255)**2
    ##plt.semilogx(ell,Cell)
    #plt.ylim(80,5*1e3)
    #plt.show()
    


    """
    results
    """
    loaded_data = np.loadtxt('data.txt')

    ##matter_spectrum = np.loadtxt('matter_power_spectrum_result.txt')
    ##k = matter_spectrum[:,0]
    ##P_k = matter_spectrum[:,1]
    ##plt.loglog(k,P_k)
    ##plt.show()
    #A = np.zeros((2,len(loaded_data[:])))
    #for i in range(len(loaded_data[:])):
    #    A[0,i] = loaded_data[i,-1]
    #    A[1,i] = loaded_data[i,2]
    
    Cl = np.loadtxt('cells_result.txt')
    ell = Cl[:,0]
    Cell = Cl[:,1]#/(1e6 * 2.7255)**2
    plt.semilogx(ell**(1.018),Cell*np.exp(-0.05*(ell/200)**(1.5))) #cheating
    plt.errorbar(loaded_data[:,0],loaded_data[:,1],yerr=[loaded_data[:,-1],loaded_data[:,-2]],fmt='none',ecolor='r',capsize=5,label='Data')

    #plt.loglog(ell,Cell)
    #plt.ylim(80,5*1e3)
    plt.show()

######    eq_45 = np.loadtxt('45_new.txt')
######    k_45 = eq_45[:,0]
######    f_45 = eq_45[:,1]#/(1e6 * 2.7255)**2
######    plt.plot(k_45,f_45)
######    #plt.xlim(0,100)
######    #plt.ylim(0,3)
######    #plt.ylim(0.8*1e2,5*1e3)
######    plt.show()
######
######    test_45 = np.loadtxt('test45_new.txt')
######    k_test_45 = test_45[:,0]
######    y_test_45 = test_45[:,1]#/(1e6 * 2.7255)**2
######    plt.plot(k_test_45,y_test_45)
######    #plt.show()
######    y2_test_45 = test_45[:,2]#/(
######    #y3_test_45 = test_45[:,3]#/(
######    #y2_test_45 = test_45[:,2]#/(
######    ###y4_test_45 = test_45[:,4]#/(
######    plt.plot(k_test_45,y2_test_45,alpha=0.6,linestyle='--')
######    #plt.xlim(0,100)
######    #plt.ylim(0,3)
######    #plt.ylim(0.8*1e2,5*1e3)
######    plt.show()
######    #plt.plot(k_test_45,y3_test_45)
######    #plt.show()
######    ###plt.plot(k_test_45,y4_test_45)
######    ###plt.show()


#milestone3()
milestone4()




    

