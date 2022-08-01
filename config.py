# --- EMULATING DATA variables ---
f_sampling = 15625          # sampling frequency in Hz
samples = 300
t_start = 0                 # start time of simulation 
t_end = samples/f_sampling                   # end time of simulation

# Source variables
away_distance = 700         # distance between the array and sources
# source1
f_start1 = 4000              # lowest frequency that the source emitts
f_end1 = 4000                # highest frequency that the source emitts
f_res1 = 20                 # resolution of frequency
theta_deg1 = 40              # theta angel of source placement, relative to origin of array
phi_deg1 = 20               # phi angel of source placement, relative to origin of array
t_start1 = 0                # start time of emission
t_end1 = t_end              # end time of emission

# source1
f_start2 = 3000              # lowest frequency that the source emitts
f_end2 = 3000                # highest frequency that the source emitts
f_res2 = 20                 # resolution of frequency
theta_deg2 = -30            # theta angel of source placement, relative to origin of array
phi_deg2 = 0               # phi angel of source placement, relative to origin of array
t_start2 = 0                # start time of emission
t_end2 = t_end                  # end time of emission


# --- ANTENNA ARRAY setup variables ---
r_a1 = [0, 0, 0]        # coordinate position of origin of array1
r_a2 = [0.08, 0, 0]         # coordinate position of origin of array2
r_a3 = [-0.24, 0, 0]        # coordinate position of origin of array3
r_a4 = [0.24, 0, 0]         # coordinate position of origin of array4
rows = 8                    # number of rows
columns = 8                 # number of columns
elements = rows*columns     # number of elements
distance = 20 * 10**(-3)    # distance between elements (m)


# --- FILTER variables ---
filter_order = 200          # filter order
scale_factor = 1000        # scale factor, adjusting filter width
f_bands_N = 45              # number of frequency bands


# --- OTHER variables ---
# Modes for adaptive weights
modes = 7

# Beamforming resolution
x_res = 10                          # resolution in x
y_res = 10                          # resolution in y

# Listening values
theta_listen = 0
phi_listen = 45


# --- RECORDED OR EMULATED DATA ---
#   choose if the program should work with recorded or emulated data
audio_signals = 'recorded'              # 'emulated' or 'recorded'
active_arrays = 1                       # number of active arrays

#  Emulated data
sources = 1                             # number of sources to emulate data from

# Recorded data
#filename = 'test2907.bin'  # filename of recorded data
filename = 'studion_A2_sound_0deg.bin'
path = '/home/batman/github/micarray-gpu/sample_data'
#path = '/home/batman/BatSignal/data/studion2607/'
#path = '/home/batman/BatSignal/'
initial_values = 40000                  # how many values to wait for until 


verbose = True

# Backend
backend = "cpu"  # Or gpu

# Index of which GPU device to use
gpu_device = 0