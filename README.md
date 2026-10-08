# QuantumFold: Quantum VQE for Protein Folding & Binding-Site Affinity

An end-to-end quantum computing pipeline addressing Use Case 01: Protein structure and binding-site prediction built for the IBM Qiskit Fall Fest.

## Project Overview
Simulating biological macromolecular folding exceeds classical tractability due to Levinthal's paradox and exponential conformational search spaces. 

QuantumFold maps the discrete 2D Hydrophobic-Polar lattice problem to an Ising Hamiltonian and computes the lowest-energy folded geometry using the Variational Quantum Eigensolver. Candidate small-molecule ligand affinities are subsequently screened using quantum state overlap fidelity against the folded pocket.

## Scientific Methodology
1. **Lattice Hamiltonian Mapping:** Contact interactions between non-bonded hydrophobic residues lower free energy, translated into a 2-qubit Ising Hamiltonian:
   $$H = 0.5 \cdot I - 0.5 \cdot Z_0 - 0.5 \cdot Z_1 + 1.0 \cdot (Z_0 \otimes Z_1)$$
2. **Variational Ansatz:** Employs parameterized RealAmplitudes circuits with full CNOT entanglement.
3. **Classical Optimization:** Classical feedback loop optimized using COBYLA with convergence tracking.
4. **Docking Evaluation:** Ligand poses modeled as parameterized phase rotations evaluated against target pocket state vectors via quantum fidelity:
   $$F(\psi_{pocket}, \psi_{ligand}) = \vert{}\langle \psi_{pocket} \vert{} \psi_{ligand} \rangle\vert{}^2$$
5. **Noise Verification:** Benchmarked under depolarizing error models on the Qiskit Aer simulator.

## Key Results
* **Exact Energy Ground Truth:** -0.5000
* **VQE Minimum Energy:** -0.5000 with zero approximation error
* **Optimal Conformation:** State |01> successfully achieves maximum contact energy while avoiding steric overlap
* **Docking Screening:** Successfully ranked 4 distinct candidate poses based on overlap affinity

## Project Architecture
* `protein_folding.ipynb`: End-to-end research notebook containing the Hamiltonian formulation, VQE solver, lattice decoder, and Aer noise simulation.
* `app.py`: Interactive Streamlit dashboard for real-time parameter tuning and visual folding inspection.
* `requirements.txt`: Frozen dependencies for full reproducibility.

## Quickstart
Clone the repository and install requirements:
```bash
git clone https://github.com/JPSAHOO123/QuantumFold.git
cd QuantumFold
pip install -r requirements.txt