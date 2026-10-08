import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp, Statevector
from scipy.optimize import minimize

st.set_page_config(page_title="QuantumFold", layout="wide")

st.title("QuantumFold: Protein Structure & Binding-Site Prediction")
st.caption("IBM Qiskit Fall Fest Hackathon | Use Case 01")

st.sidebar.header("Configuration")
optimizer_choice = st.sidebar.selectbox("Optimizer", ["COBYLA", "SLSQP"])
iterations = st.sidebar.slider("Max Iterations", 20, 120, 60)

hamiltonian = SparsePauliOp.from_list([
    ("II", 0.5),
    ("ZI", -0.5),
    ("IZ", -0.5),
    ("ZZ", 1.0)
])

ansatz = RealAmplitudes(num_qubits=2, reps=2, entanglement="full")

if st.button("Run Quantum VQE Optimization"):
    loss_history = []

    def cost_func(params):
        bound_circuit = ansatz.assign_parameters(params)
        sv = Statevector(bound_circuit)
        energy = sv.expectation_value(hamiltonian).real
        loss_history.append(energy)
        return energy

    init_params = np.random.uniform(0, 2 * np.pi, ansatz.num_parameters)
    res = minimize(cost_func, init_params, method=optimizer_choice, options={"maxiter": iterations})

    final_sv = Statevector(ansatz.assign_parameters(res.x))
    probs = final_sv.probabilities_dict()
    best_bit = max(probs, key=probs.get)

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("VQE Energy Convergence")
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.plot(loss_history, color="#006699", linewidth=2)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Energy Eigenvalue")
        ax.grid(True, linestyle="--", alpha=0.5)
        st.pyplot(fig)
        st.success(f"Final Ground Energy: {res.fun:.4f}")

    with col2:
        st.subheader("Optimal Folded Conformation")
        b0, b1 = int(best_bit[-1]), int(best_bit[-2])
        coords = [
            (0, 0),
            (1, 0),
            (1, 1 if b0 == 0 else -1),
            (0 if b1 == 1 else 2, 1 if b0 == 0 else -1)
        ]
        xs, ys = zip(*coords)
        labels = ["H0", "P1", "P2", "H3"]
        colors = ["#e63946", "#457b9d", "#457b9d", "#e63946"]

        fig2, ax2 = plt.subplots(figsize=(4, 4))
        ax2.plot(xs, ys, "-o", color="gray", linewidth=2)
        for i, (x, y) in enumerate(coords):
            ax2.scatter(x, y, color=colors[i], s=400, zorder=3)
            ax2.text(x, y, labels[i], color="white", weight="bold", ha="center", va="center")
        ax2.set_xlim(-1, 3)
        ax2.set_ylim(-2, 2)
        ax2.set_aspect("equal")
        ax2.grid(True, linestyle="--", alpha=0.5)
        st.pyplot(fig2)
        st.info(f"Predicted Ground State: |{best_bit}> ({probs[best_bit]:.1%} confidence)")

    st.markdown("---")
    st.subheader("Ligand Binding-Site Docking Screening")
    
    poses = [
        ("Pose A (Direct Alignment)", [0.0, 0.0]),
        ("Pose B (45° Tilt)", [np.pi / 4, np.pi / 2]),
        ("Pose C (Orthogonal)", [np.pi / 2, np.pi]),
        ("Pose D (Inverted)", [np.pi, np.pi])
    ]
    
    dock_data = []
    for name, angles in poses:
        test_qc = ansatz.assign_parameters(res.x).copy()
        for q_idx, ang in enumerate(angles):
            test_qc.rz(ang, q_idx)
        cand_sv = Statevector(test_qc)
        fidelity = np.abs(np.vdot(final_sv.data, cand_sv.data)) ** 2
        dock_data.append((name, -10.0 * fidelity))

    p_names, p_scores = zip(*dock_data)
    fig3, ax3 = plt.subplots(figsize=(8, 3))
    ax3.bar(p_names, p_scores, color="#2a9d8f")
    ax3.set_ylabel("Affinity Score (kcal/mol)")
    ax3.grid(axis="y", linestyle="--", alpha=0.5)
    st.pyplot(fig3)