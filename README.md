# QuSim, a transmon system simulator

<img src="QuSim_logo.png" alt="Alt text" width="330" height="300">

## Table of Contents

1. [About](#about)
2. [User_Guide](#user_guide)

## About <a name="about"></a>

A simulator aim at simulating and solving Transmon based qubit-coupler-resonator system.

This program is directed by Dr. Fei Yan, BAQIS. The original code structe is adoppted from Mr. Pan Shi, Wuhan University & BAQIS, and contributed by Mr. Jiheng Duan and Dr. Ruixia Wang, and Dr. Zhen Chen.

## User_Guide <a name = "user_guide"></a>

### Install (use pip)

```
pip install qusim
```

### Requirement

You are required to install all these python packages in the system for running this simulator. See `requirements.txt`.
>* `Python`, version >= 3.10
>* `qutip`, version >= 4.6.3
>>* `qutip 4.6.3` is strongly recommanded. As we use `qutip 4.7.5`, the `propagator` solver will crash your terminal.
>>* Nocie that `qutip 4.6.3` requires `scipy 1.10.1` exactly. Otherwise it will raise up error in function `fast_crs_matrix`. (https://stackoverflow.com/questions/76560698/python-3-10-qutip-attributeerror-cant-set-attribute-format)
>* `matplotlib`
>* `pickle`, any version
>* `jupyter lab` or `jupyter notebook`, optional

### Download

Just downloading the entire folder. 

There is also a release.

### Usage
> We have two version: the `transmon_system.py` and `arb_qubit_system.py`

The detailed tutorial are shown in folder `tutorial/~`, including required imports, define a system, what inside a system, scan zz-coupling, scan energy level, driven system dynamics, single qubit gate calibration, and DRAG calibration.

### Features

Our simulator includes the following features:

* Support multiple qubits
* Multiple initial state simulation
* Internal plot function
* DRAG
* Generalized for any qubit
* Good data saving and loading functions
* Tunable coupler helper
* Aviod crossing helper
