# QuSim, a transmon system simulator

<img src="https://github.com/RunawayFancy/QuSim/blob/main/QuSim_logo.png" alt="Alt text" width="330" height="300">

## Table of Contents

1. [About](#about)
2. [User_Guide](#user_guide)
3. [Citation](#citation)
4. [Contribution](#contribution)

## About <a name="about"></a>

A simulator aim at simulating and solving Transmon based qubit-coupler-resonator system.

This program is directed by Dr. Fei Yan, BAQIS. The original code structe is adoppted from Mr. Pan Shi, Wuhan University & BAQIS, and contributed by Mr. Jiheng Duan and Dr. Ruixia Wang, and Dr. Zhen Chen.

Look up the PyPI page: [https://pypi.org/project/qusim/](https://pypi.org/project/qusim/)

## User_Guide <a name = "user_guide"></a>

### Install (use pip)

```
pip install qusim
```

### Requirement

You are required to install all these python packages in the system for running this simulator. See `requirements.txt`.
>* `Python`, version >= 3.10
>> Installing virtual environment using conda `conda -n env_name python=python_ver -ipython`
>* `qutip`, version >= 4.6.3
<!-- >>* `qutip 4.6.3` is strongly recommanded. As we use `qutip 4.7.5`, the `propagator` solver will crash your terminal. -->
>>* Nocie that `qutip 4.6.3` requires `scipy 1.10.1` exactly. Otherwise it will raise up error in function `fast_crs_matrix`. (https://stackoverflow.com/questions/76560698/python-3-10-qutip-attributeerror-cant-set-attribute-format)
>* `matplotlib==3.8.3`
>* `mpmath==1.3.0`,
>* `pillow==10.2.0`,
>* `scipy==1.12.0`,
>* `sympy==1.12`,
>* `tqdm==4.66.2`,
>* `h5py==3.11.0`,
>* `pyqtgraph==0.13.7`,
>* `PyQt5==5.15.10`,
>* `vispy==0.14.3`,
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

### Todo list
> *Stop lazying*
>>                    ---------Jiheng Duan

* Pulse sequence compiler
* A better gate set paramter management system
* Bentchmakring tool box
* State/process/gate set tomography tool box


## Citation <a name="citation"></a>

Thank you for taking the time to try our package out. If you found it useful in your research, please cite us as follows:

```
@software{jhd2024qusim,
  author = {Jiheng Duan, Pan Shi, Zhen Chen and Fei Yan},
  title  = {QuSim: A object-oriented software toolkit for qubits simulation and control},
  url    = {https://github.com/RunawayFancy/QuSim},
  version = {1.0.3},
  year   = {2024},
}
```

## Contribution & Contact <a name="contribution"></a>

This package is open-source and welcomes all types of contributions. Whether you want to report a bug, request a feature, or submit a pull request, we encourage you to get involved. We’re also open to deeper collaborations to make the tool as useful as possible for everyone. If you’d like to discuss further, please email jiheng.duan@rochester.edu to schedule a meeting.