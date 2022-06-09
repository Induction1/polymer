# Dynamics of Force Extension of DNA

## Experiment Abstract
- Simulations have been done to study DNA chains under forced extension. In each experiment, an ensemble of
100 random-walked chains were subject to the opposite pulling forces of magnitude p at both ends. The ensemble averaged 
MSD was calculated and recorded at each step of delta_t. Various DNA lengths and extension forces were experimented 
while keeping the Hooke at 1000.
- The final result reveals that
    - The MSD of the chains are unperturbed in short time scale of t < 1e-2, and then transition 
into an alpha=1 extension phase in long time scale until the chain reaches the maximum equilibrium length.
    - Varying pulling force or the chain length does not change the long-term alpha, even though the short-long
    term transition duration can vary.

## Preliminary Studies
- Preliminary experiments have been done to find the optimal delta_t, so that the result converges within reasonable amount 
of time. It has been found that for delta_t > 5e-4, the Euler method does not converge. The result does converge for 
delta_t <= 1e-4. At the end, delta_t = 5e-5 was chosen for the experiment, so that a simulation of n = 100 completes with
in around 25 min on my machine.
- Preliminary experiments on the number of chains used for the ensemble average in each experiment have also been conducted. 
The results of 50, 100, 200 chains were compared. After comparison, 100-chain experiments were select for reasonably low
noise and simulation time. 

## Experiment Parameters
- n: The chain length
- k: Hooke constant
- p: Pulling force
- delta_t: The simulation step time interval

| Experiment | n | k | p | delta_t | # of steps | # of chains | Data Directory | Comment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 50 | 1000 | 100 | 5e-5 | 500000 | 100 | n_50_k_1000_p_100_delta_t_5e-5 | |
| 2 | 100 | 1000 | 100 | 5e-5 | 500000 | 100 | n_100_k_1000_p_100_delta_t_5e-5 | |
| 3 | 100 | 1000 | 200 | 5e-5 | 500000 | 100 | n_100_k_1000_p_200_delta_t_5e-5 | |
| 4 | 50 | 1000 | 100 | 5e-5 | 500000 | 500 | n_50_k_1000_p_100_delta_t_5e-5 | Center of Mass MSD and 0th, 13th and 25th monomer MSD |
| 5 | 100 | 1000 | 100 | 5e-5 | 500000 | 500 | n_100_k_1000_p_100_delta_t_5e-5 | Center of Mass MSD and 0th, 25th and 50th monomer MSD |

## Experiment Results
- The raw data and the down-sampled data are in data/dynamics/force_extension_dynamics directory
- [Fig 1: Simulation Diagram](https://docs.google.com/spreadsheets/d/e/2PACX-1vQ0MQw2Pa8abpHk2KiH6BZaIXhbKDfTusS5cA4SvKvIuDEg80QPQF26xqr-rOOgEevIqeUlIlV-2yPD/pubchart?oid=1729019559&format=interactive)
- [Fig 2: COM and nth Monomer MSD without Pulling Force](https://docs.google.com/spreadsheets/d/e/2PACX-1vQ0MQw2Pa8abpHk2KiH6BZaIXhbKDfTusS5cA4SvKvIuDEg80QPQF26xqr-rOOgEevIqeUlIlV-2yPD/pubchart?oid=156674255&format=interactive)
    * Without pulling force, the COM (center-of-mass) MSD obeys diffuses with alpha=1
    * Monomer MSD shows an alpha of 1 when time scale < 0.001 and an alpha of about 0.5 at time scale > 0.01.
    * The power law of COM MSD and long term monomer MSD are consistent with literature.
    * The 0th monomer has a slightly greater MDS than inner monomer.
    * There is no noticeable MSD difference between the 25th and 50th monomer.
- [Fig 3: COM and nth Monomer MSD](https://docs.google.com/spreadsheets/d/e/2PACX-1vQ0MQw2Pa8abpHk2KiH6BZaIXhbKDfTusS5cA4SvKvIuDEg80QPQF26xqr-rOOgEevIqeUlIlV-2yPD/pubchart?oid=397119336&format=interactive)
    * Compared with the no-pulling force result, when there is a pulling force in x-axis, it appears the chain starts extending from both ends of the chain and slowly propagate toward the center.
    * The nth monomer MSD is essentially unperturbed until the extension propagated to the nth monomer. The closer the monomer to the center the later it is affected by the pulling effect.
    * As the pulling forces are opposite in direction and equals in magnitude, the COM MSD is unaffected by the forces.
## Appendix
### Simulation Environment
- OS: MacOS 11.2.3 (20D91)
- Processor: 3.3 GHz 6-Core Intel Core i5
- Memory: 128 GB 2667 MHz DDR4

### Data Post-Processing and Considerations
- Each raw data output file contains 500k records in its entirety. In order to visualize the data, down-sampling is 
performed in a logarithmic fashion so that each decade has roughly the same number of data points.
- The final results are show in double-log MSD vs Time plot.

