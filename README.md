# The-two-dimensional-lattice-basis-reduction
This library implements the algorithm(s) proposed in the paper “Algorithms for the Shortest Vector Problem in 2-dimensional Lattices, Revisited”.

## Requirements

- gcc
- gmp
- Ubuntu

## How to use

All programs can be executed independently. If the lattice basis data are too large, they can be loaded from files. The reduced basis produced by the algorithms can also be written to files. The entire project is compiled and executed on Ubuntu. Data generation can be achieved by generating data files.

## Example

Compilation：

```
gcc CrossEUC.c -o alg3 -lgmp
```

Run：

```
./alg3
```

Input lattice:
```
23476 21505
18355 16814
```
Output lattice:
```text
1 35
-33 34
```
