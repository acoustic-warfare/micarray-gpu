#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

import config
config.path = os.path.dirname(os.path.realpath(__file__))

from src.antenna import Antenna

if __name__ == "__main__":
    at = Antenna(args=config)
    at.loop()

