import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as cm

import TimeVector as TV
import DiffusionLandscape
import BacterialMotion

def calc_msd_origin(x_positions, y_positions):
    msd_origin = np.average(x_positions**2 + y_positions**2)
    return msd_origin

def calc_msd_release(x_positions, y_positions, x_i, y_i):
    msd_release = np.average((x_positions - x_i)**2 + (y_positions - y_i)**2)
    return msd_release

def gaussian(x_value, y_value, A = 1e3, mu_x = -1, mu_y = -1, sigma = 2):
	return (A/(2*np.pi*sigma**2))*np.exp(-(((x_value - mu_x)**2 + (y_value - mu_y)**2)/(2*sigma**2)))

def f(x, y):
    return 200 - np.sqrt(x**2 + y**2)

t_i, t_f, delta_t = 0, 20, 0.1
start, stop, delta = -10, 10, .1
zero_axis = int((stop - start)/(2*delta))
x_i, y_i = 1, 1

TimeVector = TV.TimeVector(t_i, t_f, delta_t)
Landscape = DiffusionLandscape.Landscape(start = start,
                                        stop = stop,
                                        delta = delta,
                                        initial_conditions = f,
                                        bounday_conditions = 'dirichlet',
                                        D = 10)

Landscape.create_diff_operators(TimeVector)

Bacteria = BacterialMotion.Bacteria(N_bacteria=1000, sensitivity = 1)
Bacteria.release_bacteria(Landscape, TimeVector, x_i, y_i)

log_concentrations = np.zeros((len(Landscape.axis), len(Landscape.axis), TimeVector.N_time_steps))

log_x_positions = np.zeros((Bacteria.N_bacteria, TimeVector.N_time_steps))
log_y_positions = np.zeros((Bacteria.N_bacteria, TimeVector.N_time_steps))

log_msd_origin = np.zeros(TimeVector.N_time_steps)
log_msd_release = np.zeros(TimeVector.N_time_steps)

for time_step in range(TimeVector.N_time_steps):
    print('time_step = ' + str(time_step))
    Landscape.diffusion_step(Landscape.diff_operator, Landscape.diff_operator_inv)
    Bacteria.bacteria_step(Landscape, TimeVector)

    log_concentrations[:,:,time_step] = Landscape.concentration

    log_x_positions[:, time_step] = Bacteria.x_positions.reshape(Bacteria.N_bacteria,)
    log_y_positions[:, time_step] = Bacteria.y_positions.reshape(Bacteria.N_bacteria,)

    log_msd_origin[time_step] = calc_msd_origin(Bacteria.x_positions, Bacteria.y_positions)
    log_msd_release[time_step] = calc_msd_release(Bacteria.x_positions, Bacteria.y_positions, x_i, y_i)


fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,3,1)
ax2 = fig1.add_subplot(1,3,2)
ax3 = fig1.add_subplot(1,3,3)

#set up list of images for animation
ims=[]
for time in xrange(TimeVector.N_time_steps):
    im1 = ax1.imshow(log_concentrations[:,:,time], extent = (start, stop, start, stop), origin = 'lower', cmap = cm.gray)
    im1, = ax1.plot(log_x_positions[:, time], log_y_positions[:, time], 'bo')

    im2, = ax2.plot(log_msd_origin[:time], color='red', label='MSD origin')
    im3, = ax2.plot(log_msd_release[:time], color='blue', label='MSD release')
    im4, = ax2.plot(np.fromfunction(lambda t: 2*(Bacteria.v*TimeVector.delta_t)**2/(2*TimeVector.delta_t)*TimeVector.delta_t*t, (time,)), color='gray', label='MSD theoretical')
    ax2.set_xlabel('time')
    ax2.set_ylabel('MSD')
    ax2.legend()

    im5, = ax3.plot(log_concentrations[zero_axis,:,time], color='red', label='numerical solution')
    ax3.legend()

    ims.append([im1, im2, im3, im4, im5])

    #run animation
ani = anim.ArtistAnimation(fig1, ims, interval=50, blit=False)
ani.save('BacterialMotion.mp4')
plt.show()
