# Determinants by Geometric Memorization

This project presents a didactic and visual method for computing determinants of square matrices.
Instead of relying on Laplace expansion, the determinant is constructed through a recursive process
guided by pivots and minors, emphasizing geometric and structural understanding.

The method was developed with educational clarity in mind and is especially suited for interactive
visualization and learning environments.

---

## Core Idea

The determinant is interpreted as a recursive construction driven by pivots.
At each step:

- A pivot element is selected
- The corresponding minor matrix is formed
- The process continues until reaching order one

This approach highlights the internal structure of the determinant and avoids the formal complexity
of classical expansions.

> Pivots may be chosen freely, except for zero, which cannot serve as a valid pivot.

---

## Interactive Visualization

The main component of this project is an interactive HTML page that allows the user to:

- Visualize square matrices of arbitrary order
- Select pivots dynamically
- Observe how minors are generated
- Follow the recursive construction of the determinant

The visualization favors geometric memorization rather than symbolic repetition.

---

## Documentation

The project includes:

- A PDF article in Portuguese describing the method
- A PDF article in English presenting the same ideas for international readers
- A MATLAB script implementing the recursive determinant computation

All documents are available in the `docs/` and `matlab/` folders.

---

## MATLAB Script

The MATLAB script demonstrates how the method can be implemented programmatically.
It follows the same pivot-based recursive logic used in the visual model.

Updated English documentation was added to the original script to clarify the method and
its educational purpose.

---

## Project Structure
/ ├── index.html ├── docs/ │   ├── Determinants_Geometric_Memorization_EN.pdf │   └── Determinantes_Menorização_24062020.pdf └── matlab/ └── detp.m

---

## Educational Perspective

This project is intended for students, educators, and enthusiasts interested in alternative
approaches to linear algebra.

The emphasis is not on speed, but on understanding.

---

## Credits

Concept, method, and visualization developed collaboratively.
The project values clarity, intuition, and respect for the learner.
