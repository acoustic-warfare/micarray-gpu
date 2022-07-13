# Global Configuration
f_sampling = 1000           # sampling frequency in Hz
t_start = 0                 # start time of simulation 
t_end = 1                   # end time of simulation
away_distance = 700         # distance between the array and sources

# Antenna array setup variables
r_a1 = [-0.08, 0, 0]        # coordinate position of origin of array1
r_a2 = [0.08, 0, 0]         # coordinate position of origin of array2
r_a3 = [-0.24, 0, 0]        # coordinate position of origin of array3
r_a4 = [0.24, 0, 0]         # coordinate position of origin of array4
rows = 8                    # number of rows
columns = 8                 # number of columns
distance = 20 * 10**(-3)    # distance between elements (m)

# Filter variables
filter_order = 200          # filter order
scale_factor = 10000        # scale factor, adjusting filter width
f_bands_N = 45                                  # number of frequency bands
bandwidth = [100, f_sampling/2-f_sampling/100]  # bandwidth of incoming audio signal

# Beamforming resolution
x_res = 10                          # resolution in x
y_res = 10                          # resolution in y



verbose = True

# Backend
backend = "cpu"  # Or gpu

# Index of which GPU device to use
gpu_device = 0