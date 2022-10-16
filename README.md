# Diffusive Behaviors of Polymer Chains under Forced Extension
## Abstract: 
Previous work has shown that a single locus, or a specific DNA segment, on a free chromosome undergoes a subdiffusive Brownian diffusion in vitro and in vivo. As most biological functions involve multiple such loci to interact with each other, as a start, we focus on the diffusive behaviors of two loci under the condition of forced extension. In this study, the polymer chains were modeled with the Rouse model and programmatic simulation was carried to solve the fractional Langevin equation representing the forced extension conditions. Preliminary results indicated the extension of the chains starts from both ends where the forces are applied and propagates towards the middle of the chain. The mean square distance (MSD) between the nth locus and the center of the chain undergoes the transition from a short term diffusive behavior to a mid-term sub-diffusive behavior, then to an accelerated diffusive behavior when the chain extension reaches the nth locus. The center-of-mass of the entire does not seem to be perturbed by the opposite but equal pulling force and remains to diffuse with an alpha of 1 as expected.
Key Words: Rouse Model, Brownian motion, polymer, mean square distance

## Research Result
[Simulation Code](python_test/bm_test/test_3d.py)<br> 
[Experiment Report](data/dynamics/force_extension_dynamics/dynamics_of_force_extension_of_dna.md)<br> 
[Highlight Graph of the Nth Momomer MSD wrt the Center of the Chain](https://docs.google.com/spreadsheets/d/e/2PACX-1vQ0MQw2Pa8abpHk2KiH6BZaIXhbKDfTusS5cA4SvKvIuDEg80QPQF26xqr-rOOgEevIqeUlIlV-2yPD/pubchart?oid=397119336&format=interactive)

## Screenshots
[Chain Forced Extension](resources/PolymerChainsForceExtension.png)<br> 
[Simple Brownian Motion Ensemble MSD](resources/simple_brownian_motion_simulation.png)<br> 
[Brownian Motion in Harmonic Field](resources/BrownianMotionInHarmonicField.png)

## Reading List
- [Soft Matter Physics, Doi](https://kupdf.net/download/cgxnqsoftmatterphysics_59b0bb83dc0d609e1e568edb_pdf)
- [MIT: Stochastic Differential Equations](https://ocw.mit.edu/courses/mathematics/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/lecture-notes/MIT18_S096F13_lecnote21.pdf)
- [Nonequilibrium Statistical Mechanics](https://chz276.ust.hk/public/Cloud::siqin/References/Robert%20Zwanzig%20Nonequilibrium%20Statistical%20Mechanics.pdf): Solving first order linear differential equations
- [Papoulis Pillai: Probability Random Variables and Stochastic Processes-4th Edition](http://ce.sharif.edu/courses/97-98/1/ce181-1/resources/root/Text_Books_References/Papoulis_Pillai_Probability_RandomVariables_and_Stochastic_Processes-4th_Edition_2002.pdf)
- [Papoulis Pillai Slides](https://www.mhhe.com/engcs/electrical/papoulis/)
- [Brownian dynamics simulations of bead-rod and bead-spring chains: numerical algorithms and coarse-graining issues](http://www.polyhub.org/pub/Documentation/CoarseGrainedAndMultiscaleSimulation/bamin_09.pdf)
- [Twisting and stretching single DNA molecules](https://www.sciencedirect.com/science/article/pii/S0079610700000183): This has [coth](https://www.wolframalpha.com/input/?i=derivative+of+y%3Dln%281%2Fx+*+sinh%28x%29%29)
