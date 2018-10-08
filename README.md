# Event Driven simulations on Coupled Oscillators
A discrete phase coupled oscillators model was [introduced](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.145701 "Universality of Synchrony: Critical Behavior in a Discrete Model of Stochastic Phase-Coupled Oscillators")
by Wood et al as an example of a non-equilibrium system which presents
equilibrium-like properties near a phase transition.

Here we present some code to perform a dynamical simulation of this model on
arbitrary connection graphs between oscillators.

## How it works
The connection graph between oscillators is defined during compile time.
The current Makefile will automatically generate a ring-lattice topology with
optional random rewiring (see the Watts Strogatz model[[1]](https://www.nature.com/articles/30918 "Small World Networks")).

If you just type `make` in the root folder, a basic simulation file will be
generated for a network with a small number of nodes. In order to generate
simulation files for a custom network, run `make N=X K=Y p=Z seed=W` where `N`
is the number of nodes, `K` is half the number of neighbors per node in a ring
lattice, `p` is the probability of rewiring edges as explained in [[1]](https://www.nature.com/articles/30918 "Small World Networks")
and `seed` is the random number generator seed used when producing the random
networks. Note that `X, Y, W` are integers while `p` is a real number.

*IMPORTANT NOTE: when passing the probability `p`, Z must be given with EXACTLY
6 decimal places*

## Compiling and Running
Generate a simulation file with `make N=X K=Y p=Z gseed=W`.
Typing `make N=1000 K=15 p=0.001000 gseed=W` will generate a file named
`sim-01000-0015-0_001000-gseed_42`, which can be run with the commands:

    `./sim-01000-0015-0_001000-gseed_42 -t`

For the time evolution of a single realization of the simulation

    `./sim-01000-0015-0_001000-gseed_42 -b`

For a simulation containing many **batches** of simulations and their respective
order parameters, frequency, etc.

To see more options of command line parameters run

    `./sim-01000-0015-0_001000-gseed_42 --help`
