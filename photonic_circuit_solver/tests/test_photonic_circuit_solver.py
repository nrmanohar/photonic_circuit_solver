"""
Unit and regression test for the photonic_circuit_solver package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import numpy as np

import random as rand

import math

#import sympy

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, state_fidelity

import photonic_circuit_solver

from photonic_circuit_solver import *


def test_photonic_circuit_solver_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "photonic_circuit_solver" in sys.modules

def test_supplemental_basic():
    edges = [[0, 1], [0, 2], [1, 2], [1, 3], [2, 3]]
    graph = Graph(edges)
    L = graph.draw()
    state2 = graph.state()
    state = Stabilizer(edgelist = edges)
    n_emitters = num_emitters(state)
    expected_state = tensor(Qubit('0'*n_emitters), graph.state())
    protocol = circuit_solver(state)
    circ = qiskit_circuit_solver(state)
    g1 = Rz(0,math.pi)
    g2 = Rx(0,math.pi)
    g3 = Rx(0,math.pi)
    print(circ)

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

def test_graph_states_extended():
    candidate = [[0,1],[1,2],[2,3],[3,0],[0,2]]
    state = photonic_circuit_solver.Stabilizer(edgelist = candidate)
    operations = photonic_circuit_solver.circuit_solver(state)
    phi = Statevector.from_label('0'*6)
    psi = Statevector.from_label('0'*6)
    ref = QuantumCircuit(6)
    for i in range(4):
        ref.h(i)
    for edge in candidate:
        ref.cz(edge[0],edge[1])
    phi = phi.evolve(ref)
    qc = QuantumCircuit(6)
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

def test_stabilizers():
    indices = []
    for n in range(2,100):
        state = Stabilizer(edgelist = [[i,(i+1)%n] for i in range(n)])
        stabs = state.stabilizers()
        for j in range(n):
            compare = ''
            for k in range(n):
                if k==j:
                    compare+='X'
                elif k==(j+1)%n or k==(j-1)%n:
                    compare+='Z'
                else:
                    compare+='I'
            if compare!=stabs[j]:
                indices.append(n)
                break
    assert len(indices)==0

def test_qiskit_methods():
    candidate = [[0,1],[1,2],[2,3],[3,0],[0,2]]
    state = Stabilizer(edgelist = candidate)
    q1 = qiskit_circuit_solver(state)
    q2 = qiskit_circuit_solver_alternate(state)
    q3 = qiskit_circuit_solver(state, True)
    q4 = qiskit_circuit_solver_alternate(state, True)
    assert q1 == q2
    assert q3 == q4


    candidate = [[i,(i+1)%6] for i in range(6)]
    candidate.append([0,3])
    candidate.append([2,5])
    for i in range(6):
        candidate.append([i,i+6])
    state = Stabilizer(edgelist = candidate)
    q1 = qiskit_circuit_solver(state)
    q2 = qiskit_circuit_solver_alternate(state)
    q3 = qiskit_circuit_solver(state, True)
    q4 = qiskit_circuit_solver_alternate(state, True)
    assert q1 == q2
    assert q3 == q4


def test_compare_qiskit():
    cap = 100
    depth = 150
    for k in range(cap):
        print(k)
        state = Stabilizer(9)
        phi = Statevector.from_label('0'*9)
        qc = QuantumCircuit(9)
        gates = ['H','S','CNOT','X','Z','Y']
        for j in range(depth):
            q1 = rand.randint(0,8)
            gate = gates[rand.randint(0,5)]
            q2 = q1
            if gate == 'CNOT':
                while q2==q1:
                    q2 = rand.randint(0,8)
            state.clifford(gate,q1,q2)
            if gate == 'H':
                qc.h(q1)
            elif gate == 'S':
                qc.s(q1)
            elif gate == 'X':
                qc.x(q1)
            elif gate == 'Z':
                qc.z(q1)
            elif gate == 'Y':
                qc.y(q1)
            elif gate == 'CNOT':
                qc.cx(q1,q2)
        for i in range(1,9):
            state.row_add(0,i)
        stabs = state.stabilizers()
        phi = phi.evolve(qc)
        for stab in stabs:
            ref = phi.copy()
            pc = QuantumCircuit(9)
            phase = 1
            if stab[0] == '-':
                phase = -1
                stab = stab.lstrip('-')
            for i in range(9):
                if stab[i]=='X':
                    pc.x(i)
                elif stab[i]=='Z':
                    pc.z(i)
                elif stab[i]=='Y':
                    pc.y(i)
            ref = ref.evolve(pc)
            print(stab)
            print(qc)
            assert ref == phase*phi

def test_stabilizer_methods():
    cap = 100
    depth = 100
    for k in range(cap):
        state = Stabilizer(9)
        gates = ['H','S','CNOT','X','Z','Y']
        for j in range(depth):
            q1 = rand.randint(0,8)
            gate = gates[rand.randint(0,5)]
            q2 = q1
            if gate == 'CNOT':
                while q2==q1:
                    q2 = rand.randint(0,8)
            state.clifford(gate,q1,q2)
        print(state)
        print(state.num_qubits)
        val = 0
        try:
            print(state.__stab)
            val = 10
        except:
            pass
        assert val == 0
        state.report()
        state2 = state.clone()
        state2.stabilizer_measurement()
        state3 = Stabilizer(2)
        state3.new_stab(2,'XX,-YY')
        state3.new_stab(2,'XX,ZZ')
        state3.measurement('ZI',[0])
        rref(state3)
        assert state3.stabilizers() == ['ZI','IZ']
        state3.measurement('XI,XX',[0,1])
        rref(state3)
        state3.flip()
        assert state3.stabilizers() == ['-IX','XI']

def test_rref():
    state = Stabilizer(edgelist = [[0,1],[1,2],[2,3],[3,0],[0,2]])
    rref(state)
    h = heightfunction(state)
    for i in range(len(h)):
        assert partial_height(state, i)==h[i]
    prot = circuit_solver(state)
    assert prot is not None
    assert emitter_cnot(state) == [num_emitters(state),num_cnots(state)]

def test_partially_circuit_building():
    state = Stabilizer(edgelist = [[0,1],[1,2],[2,3],[3,0],[0,2]])
    qc = state.circuit_builder()
    state.qiskit_stabilizers()
    qc2 = state.stabilizer_measurement()
    qc3 = state.build_and_measure()        
    state.draw_circuit('text')

def test_display():
    edges = [[0, 1], [0, 2], [1, 2], [1, 3], [2, 3]]
    graph = Graph(edges)
    stabilizer = Stabilizer(edgelist = edges)
    n_emitters = num_emitters(stabilizer)

    expected_state = tensor(Qubit('0'*n_emitters), graph.state())
    assert expected_state == expected_state