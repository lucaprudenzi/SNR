# SNR
The aim of this project is computing the amplitude, the frequency and the signal-to-noise ratio of gravitational waves measured by one of the gravitational interferometers, Ligo Hanford, Ligo Livingston or Virgo, considering their orientation. At the current development state, the difference between the various interferometers is not in their sensibiity curve (it is considered the same for all), but in their different orientation in space: for this reason every interferometer observes the source under different angles, so their pattern functions F+ and Fx are not alike. 

Given the position of the source in a geocentric system, with the X-axis (the x-arm of an hypotetic interferometer) oriented towards the Greenwitch meridian and the Z-axis (the normal to the plane of this hypotetic interferometer) towards the North Pole, the first step is computing under which angles the source is seen by an interferometer in a certain position on earth surface, so with respect to the x-y plane tangent to the Earth surface and rotated with respect to the initial system. An Asymptote program shows the situation in 3D space and computes for the interferometer system the angles Theta and Phi (that identify the source point on the sky) and the angle Psi. The latter describes how the x-axis of the interferometer, projected on the binary rotational plane,is oriented respect to the semi-major axis of the binary that lies on the binary rotational plane.

The second part of the project, written in python in snr_module.py, computes the frequency, the amplitude and the SNR measured by the interferometer that sees the sources under the angles computed in the previous part. The gravitational signal is computed up to the isco frequency using linearized equations.
<p align="center">
<img src="/images/asy.png">
<img src="/images/gw150914.png">
<img width="460" height="300" src="/images/terminal.png">
</p>

## Requirements
- Python3
	- matplotlib
	- tabulate
	- numpy
- Asymptote

## How to run the project

- The only thing to modify is a text file, geoparams.txt. The parameters written in the file are:
	- mass1
	- mass2
	- luminosity distance in Mpc
	- redshift
	- theta (angle between Z-axis and soure position) in geocentric system
	- phi (angle between X-axis and source projected position on xy plane) in geocentric system
	- psi (angle between projected X-axis on binary plane and semimajor axis of the binary) in geocentric system
	- iota (angle between the normal to the interferometer arms and the orbital angular momentum)
	- the interferometer for which you want to compute the angles above: 1 for Hanford, 2 for Livingston, 3 for Virgo

Once written the desidered parameters, you have to run a bash script that executes the Asymptote module and then the python module.
To execute the script:
```
bash snr.sh
```
The Asymptote script is executed, then after having closed the Asymptote window, the python script is executed and the result SNR is printed.

- If you only want to see how the different planes are oriented, you can execute directly the Asymptote program
```
asy -V spherelite.asy
```
The three systems of axes follow the convention of X-axis: blue, Y-axis: green, Z-axis: red. At the center of the sphere there are two systems: the detector system and the geocentric system (they are approximatly at the same point if we consider the distance of the gravitational source). On the sphere surface there is the binary, with the X-axis towards the semimajor axis of the binary.
The asymptote program returns theta, phi and psi for the selected detector.
