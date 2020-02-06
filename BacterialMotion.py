from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import matplotlib.cm

'''Guided Bacteria Random Walk'''

class Bacteria(object):
    """docstring for Bacteria."""

    def __init__(self, N_bacteria = 20, v = 1, t_tumble = .1, p_tumble = .5, sensitivity = 50, len_memory = 1):
        self.N_bacteria = N_bacteria
        self.v = v
        self.t_tumble = t_tumble
        self.p_tumble = p_tumble
        self.sensitivity = sensitivity
        self.len_memory = len_memory

    def release_bacteria(self, Landscape, TimeVector, x_i, y_i):
        self.x_positions = np.full((self.N_bacteria, 1), x_i)
        self.y_positions = np.full((self.N_bacteria, 1), y_i)
        self.len_shift_registers = int(self.len_memory/TimeVector.delta_t)
        self.shift_registers = np.full((self.N_bacteria, self.len_shift_registers), Landscape.concentration[x_i, y_i])

    def update_shift_registers(self, Landscape):
        Energy = np.zeros((self.N_bacteria))

        for bacteria in range(self.N_bacteria):
            x = self.x_positions[bacteria,-1]
            y = self.y_positions[bacteria,-1]
            Energy[bacteria] = Landscape.concentration[int(x/Landscape.delta), int(y/Landscape.delta)]

        updated_shift_registers = np.zeros((self.N_bacteria, self.len_shift_registers))
        updated_shift_registers[:,0] = Energy
        updated_shift_registers[:,1:] = self.shift_registers[:,:-1]
        self.shift_registers = updated_shift_registers

    def update_positions(self, TimeVector):
        t_half  = 1 + self.sensitivity*(self.shift_registers[:,0] - self.shift_registers[:,-1])          # define half_life as function of E
        tau_run = np.where(t_half > 0.2, t_half/np.log(2), 0.2/np.log(2))
        p_tumble = 1 - np.exp(np.divide(TimeVector.delta_t, tau_run))
        alphas = 2 * np.pi * np.random.rand(self.N_bacteria,1)
        r = np.random.rand(self.N_bacteria, 1)

        updated_x_positions = np.zeros((self.N_bacteria, 1))
        updated_y_positions = np.zeros((self.N_bacteria, 1))

        for bacteria in range(self.N_bacteria):
            if r[bacteria] < p_tumble[bacteria]:
                updated_x_positions[bacteria], updated_y_positions[bacteria] = self.x_positions[bacteria,-1], self.y_positions[bacteria,-1]
                alphas[bacteria] = 2 * np.pi * np.random.random()
            else:
                updated_x_positions[bacteria] = self.x_positions[bacteria,-1] + np.cos(alphas[bacteria]) * self.v * TimeVector.delta_t
                updated_y_positions[bacteria] = self.y_positions[bacteria,-1] + np.sin(alphas[bacteria]) * self.v * TimeVector.delta_t
        self.x_positions = updated_x_positions
        self.y_positions = updated_y_positions

    def bacteria_step(self, Landscape, TimeVector):
        self.update_shift_registers(Landscape)
        self.update_positions(TimeVector)

if __name__ == '__main__':

    import DiffusionLandscape
    import TimeVector as TV
    import matplotlib.animation as anim

    def calc_msd_origin(x_positions, y_positions):
        msd_origin = np.average(x_positions**2 + y_positions**2)
        return msd_origin

    def calc_msd_release(x_positions, y_positions, x_i, y_i):
        msd_release = np.average((x_positions - x_i)**2 + (y_positions - y_i)**2)
        return msd_release

    t_i, t_f, delta_t = 0, 10, 0.1
    start, stop, delta = -5, 5, .1
    x_i, y_i = 1, 1

    TimeVector = TV.TimeVector(t_i, t_f, delta_t)
    Landscape = DiffusionLandscape.Landscape(start = start,
                                            stop = stop,
                                            delta = delta,
                                            initial_conditions = 'zero',
                                            bounday_conditions = 'dirichlet',
                                            D = 1)
    Bacteria = Bacteria(N_bacteria=100)
    Bacteria.release_bacteria(Landscape, TimeVector, x_i, y_i)

    log_x_positions = np.zeros((Bacteria.N_bacteria, TimeVector.N_time_steps))
    log_y_positions = np.zeros((Bacteria.N_bacteria, TimeVector.N_time_steps))

    log_msd_origin = np.zeros(TimeVector.N_time_steps)
    log_msd_release = np.zeros(TimeVector.N_time_steps)

    for time_step in range(TimeVector.N_time_steps):
        print('time_step = ' + str(time_step))
        Bacteria.bacteria_step(Landscape, TimeVector)

        log_x_positions[:, time_step] = Bacteria.x_positions.reshape(Bacteria.N_bacteria,)
        log_y_positions[:, time_step] = Bacteria.y_positions.reshape(Bacteria.N_bacteria,)

        log_msd_origin[time_step] = calc_msd_origin(Bacteria.x_positions, Bacteria.y_positions)
        log_msd_release[time_step] = calc_msd_release(Bacteria.x_positions, Bacteria.y_positions, x_i, y_i)


    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(1,2,1)
    ax2 = fig1.add_subplot(1,2,2)

    #set up list of images for animation
    ims=[]
    for time in xrange(TimeVector.N_time_steps):
        im1 = ax1.imshow(Landscape.concentration, extent = (start, stop, start, stop), origin = 'lower', cmap = matplotlib.cm.gray)
        im1, = ax1.plot(log_x_positions[:, time], log_y_positions[:, time], 'bo')

        im2, = ax2.plot(log_msd_origin[:time], color='red', label='MSD origin')
        im3, = ax2.plot(log_msd_release[:time], color='blue', label='MSD release')
        im4, = ax2.plot(np.fromfunction(lambda t: 2*(Bacteria.v*TimeVector.delta_t)**2/(2*TimeVector.delta_t)*TimeVector.delta_t*t, (time,)), color='gray', label='MSD theoretical')
        ax2.set_xlabel('time')
        ax2.set_ylabel('MSD')
        ax2.legend()

        ims.append([im1, im2, im3, im4])

    #run animation
    ani = anim.ArtistAnimation(fig1,ims, interval=50,blit=False)
    #ani.save('BacterialMotion.mp4')
    plt.show()


    fig2 = plt.figure(2)
    fig2_ax1 = fig2.add_subplot(1,2,1)
    fig2_ax2 = fig2.add_subplot(1,2,2)
    im = fig2_ax1.imshow(Landscape.concentration, extent = (start, stop, start, stop), origin = 'lower', cmap = matplotlib.cm.gray)

    for bacteria in range(Bacteria.N_bacteria):
        fig2_ax1.plot(log_x_positions[bacteria], log_y_positions[bacteria])

    im2 = fig2_ax2.plot(log_msd_origin)
    im2 = fig2_ax2.plot(log_msd_release)

    plt.show()
