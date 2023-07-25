# Transmon_system_simulator

## Table of Contents

1. [About](#about)
2. [User_Guide](#user_guide)

## About <a name="about"></a>

A simulator aim at simulating and solving Transmon based qubit-coupler-resonator system.

This program is direct by Dr. Fei Yan, BAQIS. The original code structe is adoppted from Mr. Pan Shi, Wuhan University & BAQIS, and contributed by Mr. Jiheng Duan and Dr. Ruixia Wang.

The original code by Mr. Shi is shown in `~\old_ver\`. The revised version using `enr_` construction from `qutip` is shown in `~\reconstructed_enr_ver\`, which is the version that we currently develop.

## User_Guide <a name = "user_guide"></a>

### Requirement

You are required to install all these python packages in the system for running this simulator.
>* `Python`, version >= 3.10
>* `qutip`, version >= 4.6.3
>* `matplotlib`
>* `pickle`
>* `jupyter lab` or `jupyter notebook`, optional

### Download

Just downloading the entire folder of desired version.

### Usage
> We only consider the `reconstructed_enr_ver`.

There are three `.py` stores functions for setup system, plotting, and waveform library and one `.ipynb` files for debugging.

The detailed tutorial are shown in `debug.ipynb`, including required imports, define a system, what inside a system, scan zz-coupling, scan energy level, and driven system dynamics.

### Features

Our simulator includes the following features:

* Support multiple qubits
* Multiple initial state simulation
* Internal plot function
* t.b.c
