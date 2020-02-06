from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import matplotlib.cm
from mpl_toolkits.mplot3d import axes3d
from DiffusionLandscape import Landscape

class ActivatorInhibitorLandscape(Landscape):
    """docstring for ActivatorInhibitorLandscape."""

    def set_activator_params(self, D = .005, rho = .01, mu = .01, sigma = 0, kappa = 0):
        self.D = D
        self.rho = rho
        self.mu = mu
        self.sigma = sigma
        self.kappa = kappa

    def set_inhibitor_params(self, D = .2, rho = .02, mu = .02, sigma = 0):
        self.D = D
        self.rho = rho
        self.mu = mu
        self.sigma = sigma

    def activator_half_step(self, inhibitor_concentration, TimeVector):
        self.concentration = self.concentration + (self.rho*self.concentration**2/((1 + self.kappa*self.concentration**2)*inhibitor_concentration) - self.mu*self.concentration + self.sigma)*TimeVector.delta_t/2

    def inhibitor_half_step(self, activator_concentration, TimeVector):
        self.concentration = self.concentration + (self.rho*activator_concentration**2 - self.mu*self.concentration + self.sigma)*TimeVector.delta_t/2

    def diffusion_half_step(self, anti_concentration, TimeVector, half_step_method):
        self.concentration = np.dot(self.diff_operator, self.concentration)
        self.concentration = self.half_step_method(anti_concentration, TimeVector)
        self.concentration = np.dot(self.diff_operator_inv, self.concentration)

    def activator_inhibitor_step(self, TimeVector):
        pass


if __name__ == '__main__':

    import TimeVector as TV
    import DiffusionLandscape as DL

    t_i, t_f, delta_t = 0, 1, 0.01
    start, stop, delta = -5, 5, 0.25

    TimeVector = TV.TimeVector(t_i, t_f, delta_t)

    ActivatorLandscape = ActivatorInhibitorLandscape(start = start, stop = stop, delta = delta, initial_conditions = 'random')
    InhibitorLandscape = ActivatorInhibitorLandscape(start = start, stop = stop, delta = delta, initial_conditions = 'random')

    Activator_concentrations = np.zeros((TimeVector.N_time_steps, len(ActivatorLandscape.axis), len(ActivatorLandscape.axis)))
    Inhibitor_concentrations = np.zeros((TimeVector.N_time_steps, len(InhibitorLandscape.axis), len(InhibitorLandscape.axis)))

    Activator_concentrations[0,:,:] = ActivatorLandscape.concentration
    Inhibitor_concentrations[0,:,:] = InhibitorLandscape.concentration

    ActivatorLandscape.set_activator_params()
    InhibitorLandscape.set_inhibitor_params()

    ActivatorLandscape.create_diff_operators(TimeVector)
    InhibitorLandscape.create_diff_operators(TimeVector)


    for time_step in range(TimeVector.N_time_steps):

        ActivatorLandscape.diffusion_half_step(InhibitorLandscape.concentration, TimeVector, activator_half_step)
        InhibitorLandscape.diffusion_half_step(ActivatorLandscape.concentration, TimeVector, inhibitor_half_step)

        ActivatorLandscape.diffusion_half_step(InhibitorLandscape.concentration, TimeVector, activator_half_step)
        InhibitorLandscape.diffusion_half_step(ActivatorLandscape.concentration, TimeVector, inhibitor_half_step)

        Activator_concentrations[time_step,:,:] = ActivatorLandscape.concentration
        Inhibitor_concentrations[time_step,:,:] = InhibitorLandscape.concentration

    #setup figure
    fig = plt.figure()
    ax1=fig.add_subplot(1,2,1)

    #set up list of images for animation
    ims=[]
    for time in xrange(Activator_concentrations.shape[0]):
        im = ax1.imshow(Activator_concentrations[time,:,:])
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_titile('Activator-Inhibitor System')
        ax2.legend()

        ims.append([im])

    #run animation
    ani = anim.ArtistAnimation(fig,ims, interval=50,blit=False)
    ani.save('ActivatorInhibitor.mp4')
    plt.show()

"""

    ''' initialize plot '''

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    im = ax.imshow(Ua[0,:,:], cmap = matplotlib.cm.gray)
    plt.show(block=False)


    ''' calculate diffusion '''

    for n in range(1,len(t)):

        for k in range(len(y)):

            Ua[n,:,k] = np.dot(Cya,Ua[n-1,:,k])
            Uh[n,:,k] = np.dot(Cyh,Uh[n-1,:,k])

        for i in range(len(x)):
            for j in range(len(y)):

                a = Ua[n-1,i,j]
                h = Uh[n-1,i,j]

                Ua[n,i,j] = Ua[n,i,j]+(rho_a*a**2/((1+kappa_a*a**2)*h)-mu_a*a+sigma_a)*dt/2
                Uh[n,i,j] = Uh[n,i,j]+(rho_h*a**2-mu_h*h+sigma_h)*dt/2

        for j in range(len(x)):

            Ua[n,j,:] = np.dot(Cxa_inv,Ua[n,j,:])
            Uh[n,j,:] = np.dot(Cxh_inv,Uh[n,j,:])

        for j in range(len(x)):

            Ua[n,j,:] = np.dot(Cxa,Ua[n,j,:])
            Uh[n,j,:] = np.dot(Cxh,Uh[n,j,:])

        for i in range(len(x)):
            for j in range(len(y)):

                a = Ua[n,i,j]
                h = Uh[n,i,j]

                Ua[n,i,j] = Ua[n,i,j]+(rho_a*a**2/((1+kappa_a*a**2)*h)-mu_a*a+sigma_a)*dt/2
                Uh[n,i,j] = Uh[n,i,j]+(rho_h*a**2-mu_h*h+sigma_h)*dt/2

        for k in range(len(y)):

            Ua[n,:,k] = np.dot(Cya_inv,Ua[n,:,k])
            Uh[n,:,k] = np.dot(Cyh_inv,Uh[n,:,k])

        im.set_array(Ua[n,:,:])
        fig.canvas.draw()


    rs = 1
    cs = 1

    X,Y = np.meshgrid(x,y)

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(221, projection='3d')
    ax.plot_surface(X, Y, Ua[0,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig1.add_subplot(222, projection='3d')
    ax.plot_surface(X, Y, Ua[1000,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig1.add_subplot(223, projection='3d')
    ax.plot_surface(X, Y, Ua[2000,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig1.add_subplot(224, projection='3d')
    ax.plot_surface(X, Y, Ua[-1,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')

    fig2 = plt.figure(2)
    ax = fig2.add_subplot(221, projection='3d')
    ax.plot_surface(X, Y, Uh[0,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig2.add_subplot(222, projection='3d')
    ax.plot_surface(X, Y, Uh[1000,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig2.add_subplot(223, projection='3d')
    ax.plot_surface(X, Y, Uh[2000,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig2.add_subplot(224, projection='3d')
    ax.plot_surface(X, Y, Uh[-1,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    '''
    fig2 = plt.figure(2)
    ax = fig2.add_subplot(221, projection='3d')
    ax.plot_surface(X, Y, Ua[-1,:,:], rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    ax = fig2.add_subplot(222, projection='3d')
    ax.plot_surface(X, Y, fa_fcn(X,Y,10), rstride=rs, cstride=cs)
    plt.xlabel('x')
    plt.ylabel('y')
    '''
    plt.show()
"""
