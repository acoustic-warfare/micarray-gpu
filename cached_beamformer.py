import multiprocessing
import numpy as np
import config
from scipy import signal
from scipy.io import wavfile
import signal as interupt
import queue

import time

propagation_speed = 340  # Speed of sound in air


class Source:
    """Audio source"""

    def __init__(self, f_start, f_end, f_res, theta_deg, phi_deg, rho, t_start, t_end):
        self.theta = theta_deg*np.pi/180
        self.phi = phi_deg*np.pi/180
        self.frequency = np.linspace(f_start, f_end, f_res)
        self.t_start = t_start
        self.t_end = t_end
        self.rho = rho


class Array:
    """Each array"""

    def __init__(self, r_a, element_distance, row_elements, column_elements):
        self.row_elements = row_elements
        self.column_elements = column_elements
        self.uni_distance = element_distance
        self.elements = row_elements * column_elements
        self.r_prime = np.zeros((3, self.elements))

        # place all microphone elements at the right position
        element_index = 0
        for i in range(row_elements):
            for j in range(column_elements):
                self.r_prime[0, element_index] = i * self.uni_distance + r_a[0]
                self.r_prime[1, element_index] = j * self.uni_distance + r_a[1]
                element_index += 1

        # center matrix in origin (0,0)
        self.r_prime[0, :] = self.r_prime[0, :] - \
            self.row_elements*self.uni_distance/2 + self.uni_distance/2
        self.r_prime[1, :] = self.r_prime[1, :] - \
            self.column_elements*self.uni_distance/2 + self.uni_distance/2


def calculate_filter_coefficients(f_sampling, frequency_bands, scale_factor, n_bands, filter_order):
    f_coefficients = np.zeros((n_bands, filter_order))
    for freq_ind in range(n_bands):
        nu_0 = 2*frequency_bands[freq_ind]/f_sampling
        cut_off = [nu_0 - nu_0/scale_factor, nu_0 + nu_0/scale_factor]
        b = signal.firwin(filter_order, cut_off, window="hamming",
                          pass_zero=False)  # filter coefficients
        f_coefficients[freq_ind, :] = b
    return f_coefficients


def weight_index(frequency):
    # calculates what mode to use, depending on the wavelength of the signal
    uni_distance = config.distance              # distance between elements
    # relative wavelength to distance between microphone elements
    wavelength_rel = frequency*uni_distance/propagation_speed

    if wavelength_rel > 0.1581:
        mode = 1
    elif (wavelength_rel <= 0.156) and (wavelength_rel > 0.0986):
        mode = 3
    elif (wavelength_rel <= 0.0986) and (wavelength_rel > 0.085):
        mode = 5
    elif (wavelength_rel <= 0.085) and (wavelength_rel > 0.07):
        mode = 6
    else:
        mode = 7
    return mode


def adaptive_array_config_matrix(matrix_array):
    # Creates the weight matrix
    row_elements = matrix_array.row_elements
    column_elements = matrix_array.row_elements

    # weight_matrix = np.zeros((7, row_elements*column_elements))
    weight_matrix = np.zeros((7, row_elements*column_elements))

    for mode in range(1, config.modes+1):
        # weight = np.zeros((1, row_elements*column_elements))
        weight = np.zeros((1, row_elements*column_elements))
        row_lim = np.ceil(row_elements/mode)
        column_lim = np.ceil(column_elements/mode)
        for i in range(int(row_lim)):
            for j in range(int(column_lim)):
                # this calculation could be wrong thanks to matlab and python index :))
                element_index = (mode*i*row_elements + mode*j)
                weight[0, element_index] = 1
        weight_matrix[mode-1, :] = weight
    return weight_matrix


def load_calibration_weights(array, elements, f_bands):
    # placeholder function, to be completed later
    # function should load calibration weights form file
    # returns matrix with calibration weightsfor all microphones, at all calibration frequencies
    weights = np.ones((f_bands, elements))
    return weights


def r_vec(theta, phi):
    return np.array([(np.sin(theta)*np.cos(phi)),
                     np.sin(theta)*np.sin(phi), np.cos(theta)])


class FIFO:

    def __init__(self, shape):
        self.buffer = np.empty(shape)

    def put(self, data):
        self.buffer = np.concatenate((self.buffer, data))

    def get(self, size):
        data = self.buffer[:size]

        self.buffer = self.buffer[size:]
        return data

    def peek(self, size):
        return self.buffer[:size]

    def getvalue(self):
        # peek with no copy
        return self.buffer

    def __len__(self):
        return len(self.buffer)

class RTCachedBeamformer(object):

    """Real-Time Cached Beamformer."""

    elements = config.rows*config.columns

    # Buffer is empty of data, but prepared to fill with each element
    receive_buffer = FIFO((0, elements))
    output_buffer = FIFO((0, 1))
    
    running = True

    def __init__(self, arrays, theta, phi, window: int):
        self.array_matrices = arrays
        self.theta = theta
        self.phi = phi
        self.window = window

        # Setup listening direction
        self.x_factor = np.sin(self.theta) * np.cos(self.phi)
        self.y_factor = np.sin(self.theta) * np.sin(self.phi)

        # load adaptive weights, calibration weights and filter coefficients
        self.frequency_bands = np.linspace(
            config.bandwidth[0], config.bandwidth[1], config.f_bands_N)

        # load adaptive weights, calibration weights and filter coefficients
        self.adaptive_weights_matrix = adaptive_array_config_matrix(
            self.array_matrices[0])

        self.calibration_weights = load_calibration_weights(
            0, config.rows*config.columns, config.f_bands_N)

        self.filter_coefficients = calculate_filter_coefficients(
            config.f_sampling, self.frequency_bands, config.scale_factor, config.f_bands_N, config.filter_order)
        # self.filter_coefficients = calculate_filter_coefficients(config.f_sampling, self.frequency_bands)
        self.freq_jobs = [freq_ind for freq_ind in range(
            len(self.filter_coefficients[:, 0]))]  # Jobs to do in parallel

        self.manager = multiprocessing.Manager()

        # Gather results
        self.return_queue = self.manager.dict()

        # Signal Catcher for Ctrl-C to stop the program gracefully
        interupt.signal(interupt.SIGINT, self.exit_gracefully)
        interupt.signal(interupt.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, *args, **kwargs):
        """Break loop on `SIGINT` or `SIGTERM`"""
        self.running = False


    def kernel(self, freq_ind, multi_window, r_prime, return_queue=None):
        b = self.filter_coefficients[freq_ind, :]
        frequency = self.frequency_bands[freq_ind]   # center frequency
        k = 2*np.pi*frequency/propagation_speed   # the narrowband frequency
        ny = frequency/config.f_sampling               # normalized frequency
        weights = self.adaptive_weights_matrix[weight_index(
            frequency)-1, :]  # Wieghts for microphones

        # Filter all 64 signals using vectorization
        filtered_signals = self.calibration_weights[freq_ind,
                                                    ] * signal.lfilter(b, 1.0, multi_window[:, ])
        phase_shift_value = -k * \
            (r_prime[0, ] * self.x_factor + r_prime[1, ]*self.y_factor)

        x_length = len(filtered_signals)
        y = np.zeros_like(filtered_signals)

        cossed = np.cos(phase_shift_value)
        sinned_ny = np.sin(phase_shift_value) / 2*np.pi*ny

        for i in range(1, x_length-1):
            y[i] = cossed * filtered_signals[i] + \
                sinned_ny * (filtered_signals[i+1]/2 - filtered_signals[i-1]/2)

        mic_data = np.sum(weights[:, ] * y, axis=(1)).reshape((x_length, 1))

        norm_coeff = 1/sum(weights)
        return_queue[freq_ind] = mic_data * norm_coeff

    def build_cache(self, frames):
        for frame in frames:
            self.get(frame)

    def old_get(self, frame):
        """Append to buffer and once buffer reaches window size, it will remove first item to only keep the newest"""
        self.buffer.append(frame)
        if len(self.buffer) > self.window_size:
            self.buffer.pop(0)
        return np.asarray(self.buffer)

    def oldput(self, audio_data):
        self.buffer.extend(audio_data)
    def put(self, audio_data):
        self.buffer = np.concatenate((self.buffer, audio_data))

    def get(self, size):
        audio_data = self.buffer[:size]

        self.buffer = self.buffer[size:]
        return audio_data

    def peek(self, size):
        return self.buffer[:size]

    def getvalues(self):
        # Peek without deletion
        return self.buffer

    def __len__(self):
        return len(self.buffer)

    def loop(self):
        while self.running:
            last = time.time()
            # Get the data in realtime
            audio_data = self.receive_buffer.get(self.window)
            if len(audio_data) == 0:
                print("Queue is empty...")
                time.sleep(0.5)
                continue

            for array in self.array_matrices:
                r_prime = array.r_prime
                jobs = []
                for freq_ind in self.freq_jobs:
                    p = multiprocessing.Process(
                        target=self.kernel, args=(freq_ind, audio_data, r_prime), kwargs=dict(return_queue=self.return_queue))

                    jobs.append(p)

                for p in jobs:
                    p.start()

                for proc in jobs:
                    proc.join()

            output = np.sum(self.return_queue.values(), axis=(0))
            self.output_buffer.put(output)
            crunch_time = time.time() - last
            # print(f"Current crunch-time: {round(crunch_time, 3)}", end="\r")
            rate = (config.window_size / crunch_time / config.f_sampling)**-1
            print(f"Real-Time Crunch Rate: {round(rate, 3)}s", end="\r")
            # break
        print()
        result = self.output_buffer.getvalue()
        result = result.reshape(len(result))
        result /= np.max(result)
        print(result, result.shape)
        print(max(result))
        # result = result.reshape()
        wavfile.write("result.wav", config.f_sampling, result.astype(np.float32))

if __name__ == "__main__":
    arrays = [Array(config.r_a1, config.distance, config.rows, config.columns)]

    beamformer = RTCachedBeamformer(arrays, config.theta_listen, config.phi_listen, config.window_size)

    filename = "array_audio_signals1.npy"
    array_audio_signals = np.load(filename, allow_pickle=True)
    print("Loading from Memory: " + filename)

    # Generate an arbitrary amount of cache and process it
    seconds = 1
    for _ in range(seconds):
        beamformer.receive_buffer.put(array_audio_signals[0])
    
    # Process it until it is interupted
    beamformer.loop()
