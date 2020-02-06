import numpy as np

class TimeVector(object):
    """docstring for TimeVector."""

    def __init__(self, t_i = 0, t_f = 10, delta_t = 0.1):
        self.t_i = t_i
        self.t_f = t_f
        self.delta_t = delta_t
        self.N_time_steps = int((t_f - t_i)/delta_t)
