# Cahn-Hilliard-Navier-Stokes
This repository contains the code for the mathematical model presented in the paper: 'A Bi-Component Model to Assess the Rheology of Soft Cellular Aggregates Probed Using the Micropipette Aspiration Technique'

_Abstract_. The micro-pipette aspiration technique is a classical experiment used to characterize the physical properties of inert fluids and biological soft materials such as cellular aggregates. The physical parameters of the fluid, as viscosity and interfacial tension, are obtained by studying how the fluid enters the pipette when the suction pressure is increased and how it relaxes when the suction pressure is put to zero. A mathematical model representative of the experiment is needed to extrapolate the physical parameters of the fluid-like matter; however, for biological materials as cells or cell aggregates these models are always based on strong starting hypotheses that impact the significance of the identified parameters. In this article, starting from the bi-constituent nature of the cell aggregate, we derive a general mathematical model based of a Cahn-Hilliard-Navier-Stokes set of equations. The model is applied to describe quantitatively the aspiration-retraction dynamics of a cell-aggregate into and out of a pipette. We demonstrate the predictive capability of the model and highlight the impact of the assumptions made on the identified parameters by studying two cases: one with a non-wetting condition between the cells and the wall of the pipette (classical assumption in the literature) and the second one, which is more realistic, with a partial wetting condition (contact angle θs = 150°). Furthermore, our results provide a purely physical explanation to the asymmetry between the aspiration and retraction responses which is alternative to the proposed hypothesis of an mechano-responsive alteration of the surface tension of the cell aggregate.

## Archive organization
* _Experimental data_ contains the excel file with the data about the aspiration and retraction phase of the micropipette aspiration experiment;
* _Images_ contains all the images reported in the [paper](https://www.i2m.u-bordeaux.fr/L-institut-UMR-5295/pages_perso/Sciume-Giuseppe); 
* _Fenics_ contains the codes used for non-wetting and partially-wetting cases including the final set of parameters identified.

## Version
The numerical simulations have been performed using the legacy FEniCS environnement which has been installed following the instructions available in the [FEniCS website](https://fenicsproject.org/download/)

## Experiment
This numerical study is based on the experimental results published by [Guevorkian _et al._(2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.218101). In this paper the authors propose the use of the MPA technique to study the surface tension and mechanical properties of cell aggregates. Due to its viscoelastic properties when the cell aggregate is aspirated, it initially responds as an elastic solid, then as a viscous fluid at times larger than a certain characteristic time of about few tens of min. The micropipette has been chemically treated to limit the adhesion of cells to the internal wall of the pipette, leading to the assumption of a 180° contact angle. During the experiment, the aspiration of the aggregate in the pipette of radius Rp was monitored by tracking the position L(t) of the end of the tongue (Figure 1.a). In the example presented in Figure 1.a, Rp = 35 µm while the initial radius of the cells aggregate (just before aspiration) is R(0) = 175 µm. A suction pressure ∆p = −1180 Pa has been applied within the pipette during 180 minutes to aspirate the aggregate and monitor the advancement of the tongue; then the pressure ∆p has been set to zero which leads the retraction of the aggregate. During the aspiration phase after a fast initial deformation (elastic phase) the advancement of the front reaches a quite constant velocity as shown in Figure 1.b. The behavior is similar in the initial stage of the retraction phase; then in the final stage of the retraction the behavior slightly evolves (see data points in Figure 1.b after 225 min). 

## Mathematical model
xxxxxx xx  xxx xxx

## Geometry and boundary conditions
xxxxxx xx  xxx xxx

## References

Guevorkian, MJ Colbert, M Durth, S Dufour, F Brochard-Wyart (2010) Aspiration of biological viscoelastic drops. Phys. Rev. Lett. 104, 218101.



