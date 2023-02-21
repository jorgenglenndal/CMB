import matplotlib.pyplot as plt
import numpy as np
loaded_file = np.loadtxt('../cosmology.txt')
#print(loaded_file[:,0])


plt.plot(loaded_file[:,0],loaded_file[:,7]+loaded_file[:,8])
plt.plot(loaded_file[:,0],loaded_file[:,4]+loaded_file[:,5])
plt.plot(loaded_file[:,0],loaded_file[:,6])
plt.show()

plt.semilogy(loaded_file[:,0],loaded_file[:,1]/3.08567758e22)
plt.ylim(1e0,5*1e4)
plt.xlim(-12,0)
plt.show()

plt.semilogy(loaded_file[:,0],loaded_file[:,2]/(100*1000)*3.08567758e22)
plt.ylim(1e-1,1e+3)
plt.xlim(-12,0)
plt.show()

plt.plot(loaded_file[:,0],loaded_file[:,1]*loaded_file[:,2]/2.99792458e8)
plt.xlim(-14.5,0)
plt.ylim(0.8,3.2)
plt.show()


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