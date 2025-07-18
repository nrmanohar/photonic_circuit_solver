'Adding test and utility functions built by Tarek Razzaz. Almost all functions require additional optional packages beyond just NumPy'

try:
    from sympy import I, pi, sqrt, sin, cos, exp, sign, simplify, latex, Expr, Function, Add, Mul
    from sympy.physics.quantum import qapply
    from sympy.physics.quantum.qubit import Qubit, QubitBra, qubit_to_matrix
    from sympy.physics.quantum.gate import IdentityGate, X, Y, Z, H, S, T, CNOT, CPHASE
except:
    pass
try:
    import networkx as nx
except:
    pass
try:
    import matplotlib.pyplot as plt
except:
    pass
try:
    from IPython.display import clear_output
except:
    pass

def dephase(state):
    '''
    Removes the global phase from a SymPy statevector by dividing by the phase of the first component (Requires SymPy)

    :param state: The SymPy statevector to be dephased
    :type state: Expr

    :return: state without the global phase
    :rtype: Expr
    '''
    try:
        res = qubit_to_matrix(state)
        glob_phase = 1
        for elem in res:
            if elem != 0:
                glob_phase = sign(elem)
                break

        return state / glob_phase
    except:
        return state


def apply(state, operation = None, dp = False):
    '''
    Simplifies an expression describing a SymPy statevector

    :param state: The SymPy statevector to be simplified
    :type state: Expr

    :param operation: An operation (usually a quantum gate) to apply to the state, defaults to None
    :type operation: Expr

    :param dp: A boolean determining if the statevector should be dephased, defaults to False
    :type dp: bool

    :return: The simplified version of the SymPy statevector
    :rtype: Expr

    '''
    if operation is None:
        if dp:
            pre_simplified = dephase(qapply((state).doit()))
            try:
                return simplify(pre_simplified).expand()
            except:
                return pre_simplified.evalf()
        else:
            pre_simplified = qapply((state).doit())
            try:
                return simplify(pre_simplified).expand()
            except:
                return pre_simplified.evalf()
    else:
        if dp:
            pre_simplified = dephase(qapply((operation * state).doit()))
            try:
                return simplify(pre_simplified).expand()
            except:
                return pre_simplified.evalf()
        else:
            pre_simplified = qapply((operation * state).doit())
            try:
                return simplify(pre_simplified).expand()
            except:
                return pre_simplified.evalf()


def tensor(state1, state2):
    '''
    Computes the tensor product of two SymPy statevectors (not fully tested) (requires SymPy)

    :param state1: The first SymPy statevector in the tensor product
    :type state1: Expr

    :param state2: The second SymPy statevector in the tensor product
    :type state2: Expr

    :return: The SymPy statevector representing the tensor product of state1 and state2
    :rtype: Expr
    '''

    term1 = apply(state1)
    term2 = apply(state2)

    if isinstance(term1, Qubit):
        if isinstance(term2, Qubit):
            bits = term1.args + term2.args
            string = ''
            for bit in bits:
                string += str(bit)
            return Qubit(string)
        elif isinstance(term2, Add):
            result = tensor(term1, term2.args[0])
            for i in range(1, len(term2.args)):
                result += tensor(term1, term2.args[i])
            return result
        elif isinstance(term2, Mul) and isinstance(term2.args[-1], Qubit):
            result = term2.args[0]
            for i in range(1, len(term2.args) - 1):
                result *= term2.args[i]
            return result * tensor(term1, term2.args[-1])

    elif isinstance(term1, Add):
        result = tensor(term1.args[0], term2)
        for i in range(1, len(term1.args)):
            result += tensor(term1.args[i], term2)
        return result

    elif isinstance(term1, Mul) and isinstance(term1.args[-1], Qubit):
        result = term1.args[0]
        for i in range(1, len(term1.args) - 1):
            result *= term1.args[i]
        return result * tensor(term1.args[-1], term2)

    else:
        raise TypeError("Unexpected input")


def Rz(target: int, theta: float):
    '''
    Apply an RZ gate via SymPy (requires SymPy)

    :param theta: Rotation angle
    :type theta: float

    :param target: Qubit to apply the RZ gate to
    :type target: int
    '''
    return cos(theta/2)*IdentityGate(target) - I*sin(theta/2)*Z(target)

def Rx(target: int, theta: float):
    '''
    Apply an RZ gate via SymPy (requires SymPy)

    :param theta: Rotation angle
    :type theta: float

    :param target: Qubit to apply the RZ gate to
    :type target: int
    '''
    return cos(theta/2)*IdentityGate(target) - I*sin(theta/2)*X(target)

def Ry(target: int, theta: float):
    '''
    Apply an RZ gate via SymPy (requires SymPy)

    :param theta: Rotation angle
    :type theta: float

    :param target: Qubit to apply the RZ gate to
    :type target: int
    '''
    return cos(theta/2)*IdentityGate(target) - I*sin(theta/2)*Y(target)

class Graph:
    """
    This is a class that encodes graphs, and contains a few convenient functions (most useful methods require optional packages)

    :param edges: A list of edges describing the graph
    :type edges: list
    """

    def __init__(self, edges: list[tuple]) -> None:
        n_qubits = 0
        for edge in edges:
            if edge[0] > n_qubits:
                n_qubits = edge[0]
            if edge[1] > n_qubits:
                n_qubits = edge[1]
        n_qubits += 1

        qubits = range(n_qubits)
        for edge in edges:
            if edge[0] not in qubits:
                raise Exception("Unexpected qubit: " + str(edge[0]))
            if edge[1] not in qubits:
                raise Exception("Unexpected qubit: " + str(edge[1]))

        self.n_qubits = n_qubits
        self.edges = edges


    def draw(self) -> None:
        '''
        Method to generate an nx graph (requires NetworkX)
        
        '''
        G = nx.Graph()
        G.add_nodes_from([i for i in range(self.n_qubits)])
        G.add_edges_from(self.edges)
        return G


    def state(self, progress: bool = False):
        '''
        Uses SymPy to calculate the statevector associated with the graph (requires SymPy)

        :param progress: A boolean option to continuously print out the progress of calculating the state (useful for large graphs), defaults to False
        :type progress: bool

        :return: A SymPy expression representing the statevector
        :rtype: Expr
        '''
        state = Qubit('0'*self.n_qubits)

        photons = []
        for edge in self.edges:
            if edge[0] not in photons:
                photons.append(edge[0])
            if edge[1] not in photons:
                photons.append(edge[1])

        for i in range(len(photons)):
            if progress:
                clear_output(wait=True)
                print("Applying Hadamards: " + str(round(100*i/len(photons))) + "%")
            photon = photons[i]
            state = apply(state, H(photon))

        for i in range(len(self.edges)):
            edge = self.edges[i]
            if progress:
                clear_output(wait=True)
                print("Applying CZs: " + str(round(100*i/len(self.edges))) + "%")
            state = apply(state, CPHASE(edge[0], edge[1]))

        if progress:
            clear_output(wait=True)

        return state
