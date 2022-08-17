# Interface Recombination Calculator

This interface recombination calculator is programmed as a Matlab App using the extended Shockley-Read-Hall interface recombination model,  with Girisch and Aberle’s iterative algorithm.

To use, downlowad: srv_appXX.mlapp and place in Matlab console folder. Run by executin gsrv_appXX.mlapp.

## Use

A tutorial of the theoretical basics and use of this app is given in https://youtu.be/bHikP9CuTh4

A second tutorial has been created to support the use of models to interpret KP measurements of CPD in the dark and under illumination. https://youtu.be/S1aIW3m5tQY 

## Method

Shockley-Read-Hall formalism for defect-mediated recombination is detailed in [1–3], and it is extended to an arbitrary trap level density function in [4]. 

To calculate surface recombination velocity (SRV) we use the algorithm below as described in [4,5].

The silicon space charge density is calculated after the work in [6,7].

This app was developed in support of the publication "Charge fluctuations at the Si–SiO2 interface and its effect on surface recombination in solar cells", Ruy S Bonilla, Isabel Al-Dhahir, Mingzhe Yu, Phillip Hamer, Pietro Altermatt. Solar Energy Materials and Solar Cells, 2020.

Additionally, the Kelvin probe part is fully described on "Modelling of Kelvin probe surface voltage and photovoltage in dielectric-semiconductor interfaces" Ruy Sebastian Bonilla, Published 5 August 2022, Materials Research Express, Volume 9, Number 8.


![image](https://user-images.githubusercontent.com/53188769/84266444-6038f180-ab1c-11ea-85d3-5735829bd662.png)

## Preview

An example calculation using this Matlab App is shown below:

![image](https://user-images.githubusercontent.com/53188769/84266721-db020c80-ab1c-11ea-9b9b-9161eae96c97.png)

## References

1]	C.T. Sah, R.N. Noyce, W. Shockley, Carrier Generation and Recombination in P-N Junctions and P-N Junction Characteristics, Proceedings of the Institute of Radio Engineers. 45 (1957) 1228–1243. doi:10.1109/JRPROC.1957.278528.

[2]	W. Shockley, W.T. Read, Statistics of the Recombinations of Holes and Electrons, Physical Review. 87 (1952) 835–842.

[3]	R.N. Hall, Electron-Hole Recombination in Germanium, Physical Review. 87 (1952) 387.

[4]	R.B.M. Girisch, R.P. Mertens, R.F. Dekeersmaecker, Determination of Si-SiO2 Interface Recombination Parameters using a Gate-Controlled Point-Junction Diode under Illumination, Ieee Transactions on Electron Devices. 35 (1988) 203–222.

[5]	A.G. Aberle, S. Glunz, W. Warta, Impact of illumination level and oxide parameters on Shockley-Read-Hall recombination at the Si-SiO2 interface, Journal of Applied Physics. 71 (1992) 4422–4431. doi:10.1063/1.350782.

[6]	A.S. Grove, D.J. Fitzgerald, Surface effects on p-n junctions: Characteristics of surface space-charge regions under non-equilibrium conditions, Solid-State Electronics. 9 (1966) 783–806. doi:10.1016/0038-1101(66)90118-3.

[7]	C.E. Young, Extended curves of the space charge, electric field, and free carrier concentration at the surface of a semiconductor, and curves of the electrostatic potential inside a semiconductor, Journal of Applied Physics. 32 (1961) 329–332. doi:10.1063/1.1736007.


