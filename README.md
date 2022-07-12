# micarray-gpu

A GPU Accelerated version of Acoustic Warfare 


## Requirements
---------------
### Hardware

A Nvidia GPU with CUDA support.

## Software

**NVIDIA Drivers**

It is preffered to install display drivers using the distribution's native package management tool, i.e `apt`. If not installed already, NVIDIA display drivers can be installed from [NVIDIA Download Center](https://www.nvidia.com/Download/index.aspx?lang=en-us).

**1. CUDA** Runtime and Toolkit

To install CUDA Toolkit `CUDA 11.7` from Nvidia on Ubuntu 20.04:

    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
    sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
    wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb

    sudo dpkg -i cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
    sudo cp /var/cuda-repo-ubuntu2004-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/

    sudo apt-get update
    sudo apt-get -y install cuda

**Conda**

Setup Conda Installation (Miniconda)

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    bash Miniconda3-latest-Linux-x86_64.sh

    conda init

Disable auto-init

    conda config --set auto_activate_base false

Recreate an environment inside micarray-gpu

    conda env create -f environment.yml

Activate and enter the newly created environment

    conda activate aw
    ...
    (aw) $ 