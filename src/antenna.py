#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# File: micarray-gpu/src/antenna.py
# Author: Irreq
# Date: 11/07-2022

import signal
import config

# Load Library
if config.backend == "cpu":
    import numpy as xp
    cubridge_lib = None
    if config.verbose:
        print("Will use CPU backend")
elif config.backend == "gpu":
    import ctypes
    import cupy as xp
    xp.cuda.runtime.setDevice(config.gpu_device)
    print("Will use GPU backend")
    try:
        cubridge_lib = ctypes.cdll.LoadLibrary(
            config.path + "/build/cubridge.so")
        # # B. Specify function signatures
        # bridge_fcn = cubridge_lib.pythonCudaBridgeWrapper
        # bridge_fcn.argtypes = [ctypes.c_int]

        # # bridge_fcn(int(input("Size: ")))
        # bridge_fcn(9)
        print("Successfully loaded Cuda APIs")
    except OSError:
        raise OSError("Have you compiled the Cuda Libs? run make")
else:
    exit("No backend has been defined")

arguments = {
    "rows": int,
    "columns": int,
    "distance": float,
    "r_a": list,
}


def test():
    import time
    print(time.time())
    time.sleep(1)


class Antenna(object):
    """Main Antenna Class
    
    Variables can be defined at runtime, but variables defined in `config.py`
    have higher priority.
    
    """
    running = True

    def __init__(self, rows=8, columns=8, idx=None, distance=2e-2, r_a=[0.1, 0.0, 0.0], args=None, callback=None):
        self.rows = rows
        self.cols = columns
        self.id = idx
        self.distance = distance
        self.r_a = r_a

        # Sanity check arguments
        if args is not None:
            self._sanity_check(args)

        # Print time if no callback is defined.
        if callback is None:
            self.callback = test
        else:
            self.callback = callback

        # Not configurable
        self._elements = self.rows * self.columns
        self._r_prime = xp.zeros((3, self._elements))

        self.array = self.setup()

        # Signal Catcher for Ctrl-C to stop the program gracefully
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, *args, **kwargs):
        """Break loop on `SIGINT` or `SIGTERM`"""
        self.running = False

    def _sanity_check(self, args):
        for tag in arguments:
            tmp = getattr(args, tag)
            try:
                assert isinstance(tmp, arguments.get(
                    tag)), 'Argument of wrong type!'
                self.__dict__[tag] = tmp
            except AssertionError:
                print(f"{tag} should be: {arguments.get(tag).__name__}")
                exit(1)

    def _get_audio_signal(self):
        return self.audio_signal

    def setup(self):
        element_index = 0

        for i in range(self.rows):
            for k in range(self.columns):
                self._r_prime[0][element_index] = i * \
                    self.distance + self.r_a[0]
                self._r_prime[1][element_index] = k * \
                    self.distance + self.r_a[1]
                element_index += 1

        return self._r_prime

    def loop(self):
        while self.running:
            self.callback()
