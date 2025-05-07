# Anticoherent States Generation

In this repository, you will find the parameters which are necessary to generate anticoherent spin states.

## Description of the protocol

Starting from a coherent state pointing in the $x$ direction, our protocol consists of a sequence of $n_{C}$ cycles, where each cycle (except the first) applies a rotation around the $y$-axis followed by a squeezing operation along $z$. The corresponding operations are described by the operators

$$ R_{y}(\theta) = e^{-iJ_{y}\theta}, \qquad S_{z}(\eta) = e^{-iJ_{z}^{2}\eta} $$

where $\theta$ and $\eta$ are the amplitudes of the rotation and the squeezing. Note that in the first cycle only the squeezing operation is applied. The final state of the system after $n_C$ cycles is then given by

$$ |\psi_{n_C}\rangle = \left(\prod_{i=1}^{n_C}S_{z}(\eta_i)\,R_{y}(\theta_i)\right)|\psi_0\rangle. $$

with $\theta_1=0$.

## Data file structure

In the Data folder, you will find the files storing the parameters $\theta_i$ and $\eta_i$ ($i=1,2,\dots,n_C$) which maximise the anticoherence measures of spin states obtained via numerical optimisation. There is two types of file structure :

- Analytical parameters : For each number of qubits $N=2j$ (first column) are given the squeezing parameters $\eta_2$ and $\eta_3$ (second and third columns) as explained in our paper.
- Other files : For a given number of cycles (first column), we give the maximal anticoherence measure we obtained (second column) followed by the $n_C-1$ rotation parameters $\theta$ and then the $n_C$ squeezing parameters $\eta$. The files preceded by "OptimisedSqueezing" have the same structure and are obtained from an optimisation on the total squeezing duration.

## Julia code

In the Julia folder, you will find the functions necessary to optimise the protocol in order to generate the anticoherence measure. An example of use is given in ACStatesGeneration.jl.
