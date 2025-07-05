"""
Unit and regression test for the photonic_circuit_solver package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import numpy as np

import random as rand

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, state_fidelity

import photonic_circuit_solver


def test_photonic_circuit_solver_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "photonic_circuit_solver" in sys.modules

def test_basic():
    indices = []
    for i in range(1,10):
        state = photonic_circuit_solver.Stabilizer(i)
        if state.signvector.all() == np.zeros(i).all():
            pass
        else:
            indices.append(i)
    assert len(indices)==0
def test_graph_states():
    state = photonic_circuit_solver.Stabilizer(edgelist = [[i,(i+1)%5] for i in range(5)])
    operations = photonic_circuit_solver.circuit_solver(state)
    phi = Statevector.from_label('0'*7)
    psi = Statevector.from_label('0'*7)
    ref = QuantumCircuit(7)
    for i in range(5):
        ref.h(i)
    for i in range(5):
        ref.cz(i,(i+1)%5)
    phi = phi.evolve(ref)
    qc = QuantumCircuit(7)
    for operation in operations:
        if operation[0]=='H':
            qc.h(operation[1])
        elif operation[0] == 'Emission':
            qc.cx(operation[1],operation[2])
        elif operation[0]=='Z':
            qc.z(operation[1])
        elif operation[0]=='X':
            qc.x(operation[1])
        elif operation[0]=='Y':
            qc.y(operation[1])
        elif operation[0]=='S':
            qc.s(operation[1])
        elif operation[0]=='CNOT':
            qc.cx(operation[1],operation[2])
        elif operation[0]=='Measure':
            qc.cx(operation[1],operation[2])
            qc.h(operation[1])
    psi = psi.evolve(qc)
    assert np.isclose([state_fidelity(psi,phi)],[1])


