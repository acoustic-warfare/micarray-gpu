import config

from src.antenna import Antenna

if __name__ == "__main__":
    at = Antenna(args=config)
    at.loop()
