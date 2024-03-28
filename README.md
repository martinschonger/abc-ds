# ABC-DS: obstacle Avoidance with Barrier-Certified polynomial Dynamical Systems

Accompanying code for the paper **Learning Barrier-Certified Polynomial Dynamical Systems for Obstacle Avoidance with Robots** by
Martin Schonger<sup>1*</sup>,
Hugo T. M. Kussaba<sup>1*</sup>,
Lingyun Chen<sup>1</sup>,
Luis Figueredo<sup>2</sup>,
Abdalla Swikir<sup>1</sup>, 
Aude Billard<sup>3</sup>,
and Sami Haddadin<sup>1</sup>, which has been accepted for presentation at ICRA 2024.

<sup>1</sup>Munich Institute of Robotics and Machine Intelligence (MIRMI), Technical University of Munich (TUM), Germany. Abdalla Swikir is also with the Department of Electrical and Electronic Engineering, Omar Al-Mukhtar University (OMU), Albaida, Libya.\
<sup>2</sup>School of Computer Science, University of Nottingham, UK. Luis Figueredo is also an Associated Fellow at the MIRMI, TUM.\
<sup>3</sup>Learning Algorithms and Systems Laboratory, EPFL, Switzerland.\
<sup>*</sup>These authors contributed equally to the paper.


### Setup
Install MATLAB (tested with R2023a).

Install the required MathWorks toolboxes:
Control System Toolbox,
Robust Control Toolbox,
Optimization Toolbox,
Signal Processing Toolbox,
Symbolic Math Toolbox,
Statistics and Machine Learning Toolbox.

Install the required third party tools:
[YALMIP](https://yalmip.github.io/),
[PENLAB](https://web.mat.bham.ac.uk/kocvara/penlab/)+[PENBMI](http://www.penopt.com/penbmi.html),
[MOSEK](https://www.mosek.com/),
[GUROBI](https://www.gurobi.com/),
(Optional for plotting: [crameri colormaps](https://de.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)).

> **Note**
> Make sure that the non-toolbox paths are before/on top of the toolbox paths.

Run:
```bash
git clone https://github.com/martinschonger/abc-ds.git
cd abc-ds
git submodule init
git submodule update
```

### Usage
Open the `abc-ds` folder in MATLAB.

Configure the desired experiments in `main2.m` and run this script.

Check the `output` folder for results and logs.

(Optionally, recreate the plots from the paper with `generate_plots.m`, and the animations from the video with `generate_plots_video.m`.)


### Contact
martin.schonger@tum.de


This software was created as part of Martin Schonger's master's thesis in Computer Science at the Technical University of Munich's (TUM) School of Computation, Information and Technology (CIT).


Copyright © 2023 Martin Schonger  
This software is licensed under the GPLv3.
