from setuptools import setup, find_packages

setup(
    name='qusim',
    version='1.1.0',
    packages=find_packages(),
    install_requires=[
        'matplotlib==3.8.3',
        'mpmath==1.3.0',
        'numpy==1.26.4',
        'pillow==10.2.0',
        'qutip==4.7.5',
        'scipy==1.12.0',
        'sympy==1.12',
        'tqdm==4.66.2',
        'h5py==3.11.0',
        'pyqtgraph==0.13.7',
        'PyQt5==5.15.10',
        'vispy==0.14.3',
    ],
    author='Jiheng Duan',
    author_email='jiheng.duan@rochester.edu',
    description='A simulator aim at simulating and solving Transmon based qubit-coupler-resonator system',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/RunawayFancy/QuSim',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
    ],
)
