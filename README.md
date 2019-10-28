# Characterization of Conditional Von Neumann Entropy

This repository contains all the codes used to generate results in the paper
'Characterization of Conditional Von Neumann Entropy' written in `MATLAB` The pdf of the paper is available in the repository. 
It also contains a few quantum functions useful for 2-qubit systems written in `sage`. 

Functions from the following libraries have been used to write the functions in this repository. 
- cvx
- cvxquad
- quantinf

The repository primarily contains codes to generate analytical and numerical witnesses to separate a given state from the CVENN class. 

Clone the repository:


## To generate analytical witness 
Open matlab in the same directory as all the files are present:

To use functions in analytical witness file:

```aw = analytical_witness```

Define the matrix `rho_s` as the state you want to separate from the CVENN class, or use the predefined function to 
generate a Werner like state `alpha|phi^+> <phi+| + (1-alpha)I/d^2`. Please ensure the state you generate lies outside the CVENN class. Refer to the `alpha vs dimensions graph` in the paper and pick an `alpha` such that the state lies outside the class. Eg: `alpha = 0.8, d = 5`

```rho_s = aw.werner_like_state(alpha, dimensions)```

Generate a witness for the state:

```W_a = aw.witness(rho_s, dimensions)```

Check the state against the witness:

```check_witness(W_a, rho_s)```

Test any other state against the witness to see if states inside CVENN give a positive value. 

## To generate a numerical witness
Note that unlike the analytical witness which quickly finds solutions for large dimensions, the numerical witness takes a 
long time to converge for greater than a 2-qubit system. To find a witness for a 2-qubit or 2-qutrit system, do the following. Note that the 3-qutrit system will take 8-9 minutes to find the closest state. 

We would have to use the `cvx` library to converge to the closest CVENN state. To setup cvx, open MATLAB in the cvx folder of this repository and run `cvx_setup` in MATLAB. 

Now open MATLAB in the root of the repository. To add all the functions in `cvxquad` to the MATLAB path, run:

```addpath(cvxquad)```

Once again, either define `rho_s`, the state you wish to separate from CVENN yourself, or generate a Werner like state. 

```rho_s = aw.werner_like_state(alpha, dimensions) % where dimensions = 2 or 3``` 

Find the numerical witness:

```W_n = numerical_witness(rho_s, dimensions)```

Check the state against the witness

```check_witness(W_n, rho_s)```

Check other states against this witness

```check_witness(W_n, x)```

## Further work:

The distance from the witness is an example of an entanglement measure of the state. Comparing how the analytical and numerical witness perform as such a measure could be an interesting project. Note that the numerical witness finds distance from the closest CVENN state, whereas the analytical witness finds distance from any state at the boundary. (Read paper for more details).



