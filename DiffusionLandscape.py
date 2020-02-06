from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import matplotlib.cm

class Landscape(object):
    """docstring for ."""

    def __init__(self, start = -1, stop = 1, delta = .2, initial_conditions='random', bounday_conditions='continuous', D = 1):
        self.start = start
        self.stop = stop
        self.delta = delta
        self.D = D
        self.axis = np.arange(start, stop, delta)
        self.bounday_conditions = bounday_conditions

        if initial_conditions == 'random':
            self.concentration = np.random.rand(len(self.axis), len(self.axis))
        elif initial_conditions == 'zero':
            self.concentration = np.zeros((len(self.axis), len(self.axis)))
        else:
            self.concentration = initial_conditions(self.axis[:, None], self.axis[None,:])

    def create_diff_operators(self, TimeVector):
        A_center = -2 * np.diag(np.ones(len(self.axis)), 0)
        A_upper = np.diag(np.ones(len(self.axis) - 1), 1)
        A_lower = np.diag(np.ones(len(self.axis) - 1), -1)
        A = A_center + A_upper + A_lower

        if self.bounday_conditions == 'unspecified':
            A = 1/self.delta**2*A
        elif self.bounday_conditions == 'continuous':
            A[0,-1] = 1
            A[-1,0] = 1
            A = 1/self.delta**2*A

        s = self.D*TimeVector.delta_t/2
        I = np.identity(len(self.axis))

        self.diff_operator = I + s*A
        self.diff_operator_inv = inv(I - s*A)


    def diffusion_step(self):
        self.concentration = np.dot(np.dot(self.diff_operator,np.dot(self.diff_operator_inv,np.dot(self.concentration,self.diff_operator))),self.diff_operator_inv)



if __name__ == '__main__':

    import matplotlib.animation as anim
    import time

    import TimeVector as TV

    t_i, t_f, delta_t = 0, 1, 0.01

    TimeVector = TV.TimeVector(t_i, t_f, delta_t)

    start, stop, delta = -5, 5, 0.25
    zero_axis = int((stop - start)/(2*delta))


    def gaussian(x_value, y_value, A = 1, mu_x = 0, mu_y = 0, sigma = 1):
    	return (A/(2*np.pi*sigma**2))*np.exp(-(((x_value - mu_x)**2 + (y_value - mu_y)**2)/(2*sigma**2)))

    def gaussian_analytical(x_value, y_value, t, A = 1, D = 1):
        return (A/(2*np.pi*D*t)*np.exp(-(x_value**2 + y_value**2)/(2*D*t)))


    Landscape = Landscape(start = start, stop = stop, delta = delta, initial_conditions = gaussian, bounday_conditions = 'dirichlet', D = 1)

    Landscape.create_diff_operators(TimeVector)
    #diff_operator = Landscape.diff_operator
    #diff_operator_inv = Landscape.diff_operator_inv

    concentrations = Landscape.concentration

    for time_step in range(TimeVector.N_time_steps):
        print('time_step = ' + str(time_step))
        Landscape.diffusion_step()
        concentrations = np.dstack((concentrations, Landscape.concentration))


    #setup figure
    fig = plt.figure()
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)

    #set up list of images for animation
    ims=[]
    for time in range(concentrations.shape[2]):
        im = ax1.imshow(concentrations[:,:,time])
     #   im2 = ax2.imshow(image2[time,:,:])
        im2, = ax2.plot(concentrations[zero_axis,:,time], color='red', label='numerical solution')
        im3, = ax2.plot(gaussian_analytical(Landscape.axis, 0, 1 + time*TimeVector.delta_t), color='blue', label='analytical solution')
        ax2.set_xlabel('x')
        ax2.set_ylabel('concentration')
        ax2.legend()
     #   def setvisible(self,vis):
     #       for c in self.collections: c.set_visible(vis)
     #    im2.set_visible = types.MethodType(setvisible,im2,None)
     #    im2.axes = plt.gca()

        ims.append([im, im2, im3])

    #run animation
    ani = anim.ArtistAnimation(fig,ims, interval=50,blit=False)
    ani.save('DiffusionLandscape.mp4')
    plt.show()
