"Contains the functions needed to implement the emission circuit solving procedure"
import numpy as np
import math
from .stabilizer import *

def rref(state: Stabilizer):
    """
    Given a stabilizer state, implements the RREF gauge procedure (details in doi.org/10.1088/1367-2630/7/1/170)

    :param state: Stabilizer state on which to implement the RREF gauge procedure
    :type state: Stabilizer

    """    

    N=state.size
    K=N
    KU=0
    NL=0
    while NL<N-1 and KU<K-1:
        zeroitem = True
        oneitem = False
        twoitem = False
        r1=N
        r2=N
        for k in range(KU,K):
            if state.tab[k,NL]!=0 or state.tab[k,NL+N]!=0:
                r1 = k
                zeroitem = False
                oneitem = True
                break
        for k in range(r1,K):
            if state.tab[k,NL]!=0 or state.tab[k,NL+N]!=0:
                if state.tab[k,NL]!=state.tab[r1,NL] or state.tab[k,NL+N]!=state.tab[r1,NL+N]:
                    r2 = k
                    oneitem = False
                    twoitem = True
                    break
        if zeroitem:
            NL+=1
        elif oneitem:
            state.swap(KU,r1)
            for i in range(KU+1,K):
                if state.tab[i,NL]!=0 or state.tab[i,NL+N]!=0:
                    state.row_add(KU,i)
            KU+=1
            NL+=1
        elif twoitem:
            state.swap(KU,r1)
            state.swap(KU+1,r2)
            for i in range(KU+2,K):
                if state.tab[i,NL]!=0 or state.tab[i,NL+N]!=0:
                    if state.tab[i,NL]==state.tab[KU,NL] and state.tab[i,NL+N]==state.tab[KU,NL+N]:
                        state.row_add(KU,i)
                    elif state.tab[i,NL]==state.tab[KU+1,NL] and state.tab[i,NL+N]==state.tab[KU+1,NL+N]:
                        state.row_add(KU+1,i)
                    else:
                        state.row_add(KU,i)
                        state.row_add(KU+1,i)
            NL+=1
            KU+=2

def heightfunction(state: Stabilizer):
    """
    Given a stabilizer state, finds the height function (bipartite entanglement entropy) of the state (formula in doi.org/10.1038/s41534-022-00522-6, equation 1)

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :return: The outputs of the height function as a list
    :rtype: list

    """        

    rref(state)
    N=state.size
    leftmost = []
    for i in range(state.size):
        for j in range(state.size):
            if state.tab[i,j]!=0 or state.tab[i,j+N]:
                leftmost.append(j+1)
                break
    height = []
    for i in range(state.size+1):
        count = sum(j > i for j in leftmost)
        height.append(state.size-i-count)
    return height

def partial_height(state: Stabilizer,index: int):
    """
    Given a stabilizer state, finds the height function (bipartite entanglement entropy) of the state (formula in doi.org/10.1038/s41534-022-00522-6, equation 1) at a specified index

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :param index: The index at which to evaluate the height function
    :type index: int

    :return: The outputs of the height function at the index as an integer
    :rtype: int

    """    
    rref(state)
    N=state.size
    leftmost = []
    for i in range(state.size):
        for j in range(state.size):
            if state.tab[i,j]!=0 or state.tab[i,j+N]:
                leftmost.append(j+1)
                break
    return state.size-index-sum(j > index for j in leftmost)

def plot_height(state: Stabilizer):
    """
    Given a stabilizer state, finds the height function (bipartite entanglement entropy) of the state (formula in doi.org/10.1038/s41534-022-00522-6, equation 1) and plots it (need matplotlib installed)

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    """     

    try:
        import matplotlib.pyplot as plt
    except:
        raise ImportError('Matplotlib.pyplot not present')
    try:
        height = heightfunction(state)
        x_val = [i for i in range(state.size+1)]
        tickers = range(math.floor(min(height)), math.ceil(max(height))+1)
        plt.grid(color = 'blue', linewidth = 0.5)
        plt.plot(x_val,height,color='blue')
        plt.scatter(x_val,height,color='blue')
        plt.yticks(tickers)
        plt.title('Target Height Function')
        plt.xlabel('x')
        plt.ylabel('h(x)')
    except:
        pass


def num_emitters(state: Stabilizer):
    """
    Given a graph state with a certain fixed ordering, finds the minimal number of emitters required to emit that state

    :param state: Target photonic stabilizer state
    :type state: Stabilizer

    :return: Minimal number of photons needed to emit the input states
    :rtype: int

    """   
    height = heightfunction(state)
    emitters = max(height)
    return emitters

def circuit_solver(state: Stabilizer):
    """
    Given a stabilizer state, finds a circuit to generate it from a minimal set of emitters (details in doi.org/10.1038/s41534-022-00522-6)

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :return: Procedure to generate the state, of the form ['Action', q1, q2 (if applicable)]. 'Action' can be 'X', 'Z', 'H', 'S', 'CNOT', 'Measure', or 'Emission'. For 'CNOT', q1 is the control and q2 is the target. For emission and measure, q1 is the emitter and q2 is the photon.
    :rtype: list

    """   
    n_e = num_emitters(state)
    n_p = state.size
    N = n_e+n_p
    target_state = Stabilizer(N)
    for i in range(n_p):
        target_state.signvector[i] = state.signvector[i]
        for j in range(n_p):
            target_state.tab[i,j] = state.tab[i,j]
            target_state.tab[i,j+N] = state.tab[i,j+n_p]
    protocol = []
    for h in range(n_p,0,-1):
        photonindex = h-1
        d = partial_height(target_state,h)-partial_height(target_state,h-1)
        if d<0:
            rows = []
            for i in range(N):
                toggler = True
                for j in range(n_p):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
                if toggler:
                    rows.append(i)
            sums=[]
            for row in rows:
                sum=0
                for j in range(n_p,N):
                    if target_state.tab[row,j]!=0 or target_state.tab[row,j+N]!=0:
                        sum+=1
                sums.append(sum)
            row = rows[sums.index(min(sums))]
            for i in range(n_p,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    emit = i
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    break
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    protocol.append(['CNOT',i,emit])
                    target_state.clifford('CNOT',i,emit)
            if target_state.signvector[row]==1:
                target_state.clifford('X',emit)
                protocol.append(['X',emit])
            target_state.clifford('H',emit)
            target_state.clifford('CNOT',emit,photonindex)
            protocol.append(['Measure',emit,photonindex])
        for i in range(N):
            toggler = True
            if target_state.tab[i,photonindex]==0 and target_state.tab[i,photonindex+N]==0:
                toggler = False
            if toggler:
                for j in range(photonindex):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
            if toggler:
                row = i
                break
        emit = -1
        if target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==0:
            protocol.append(['H',photonindex])
            target_state.clifford('H',photonindex)
        elif target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==1:
            protocol.append(['S',photonindex])
            target_state.clifford('S',photonindex)
            target_state.clifford('Z',photonindex)
            protocol.append(['H',photonindex])
            target_state.clifford('H',photonindex)

        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',i])
                    target_state.clifford('S',i)
                    target_state.clifford('Z',i)
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                break
        if emit!= -1:
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        protocol.append(['S',i])
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        protocol.append(['H',i])
                        target_state.clifford('H',i)
                    protocol.append(['CNOT',i,emit])
                    target_state.clifford('CNOT',i,emit)
            if target_state.signvector[row]==1:
                target_state.clifford('X',photonindex)
                protocol.append(['X',photonindex])
            target_state.clifford('CNOT',emit,photonindex)
            protocol.append(['Emission',emit,photonindex])
        for i in range(N):
            if target_state.tab[i,photonindex+N]!=0 or target_state.tab[i,photonindex]!=0:
                if i!=row:
                    target_state.row_add(row,i)
    rref(target_state)
    sums = []
    for i in range(n_p,N):
        sum=0
        for j in range(n_p,N):
            if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                sum+=1
        sums.append(sum)
    
    if max(sums)==1:
        decoupled = True
    else:
        decoupled = False
    
    while not decoupled:
        minimum = min(i for i in sums if i > 1)
        row = n_p+sums.index(minimum)
        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',emit])
                    target_state.clifford('H',emit)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',emit])
                    target_state.clifford('S',emit)
                    target_state.clifford('Z',emit)
                    protocol.append(['H',emit])
                    target_state.clifford('H',emit)
                break
        for i in range(emit+1,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    protocol.append(['S',i])
                    target_state.clifford('S',i)
                    target_state.clifford('Z',i)
                    protocol.append(['H',i])
                    target_state.clifford('H',i)
                target_state.clifford('CNOT',i,emit)
                protocol.append(['CNOT',i,emit])
        for i in range(n_p,N):
            if target_state.tab[i,emit]!=0 or target_state.tab[i,emit+N]!=0:
                if i!= row:
                    target_state.row_add(row,i)
        sums = []
        for i in range(n_p,N):
            sum=0
            for j in range(n_p,N):
                if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                    sum+=1
            sums.append(sum)
        if max(sums)==1:
            decoupled = True
        else:
            decoupled = False
    for emitter in range(n_p,N):
        for row in range(N):
            toggler = True
            colcheck = True
            if target_state.tab[row,emitter]==0 and target_state.tab[row,emitter+N]==0:
                toggler = False
                colcheck = False
            if colcheck:
                for i in range(emitter):
                    if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                        toggler = False
                        break
            if toggler:
                target_state.swap(row,emitter)
                break
    for i in range(n_p,N):
        if target_state.tab[i,i]!=0:
            if target_state.tab[i,i]==1 and target_state.tab[i,i+N]==0:
                protocol.append(['H',i])
                target_state.clifford('H',i)
            elif target_state.tab[i,i]==1 and target_state.tab[i,i+N]==1:
                protocol.append(['S',i])
                target_state.clifford('S',i)
                target_state.clifford('Z',i)
                protocol.append(['H',i])
                target_state.clifford('H',i)
    
    for i in range(n_p,N):
        if target_state.signvector[i]==1:
            target_state.clifford('X',i)
            protocol.append(['X',i])
    checker_state = Stabilizer(N)
    if np.array_equal(checker_state.tab,target_state.tab) and np.array_equal(checker_state.signvector,target_state.signvector):
        protocol = protocol[::-1]
        return protocol
    else:
        print('Something went wrong')
        print(target_state.stabilizers())
        return None
     
def qiskit_circuit_solver(state: Stabilizer, simple: bool = False):
    """
    Given a stabilizer state, creates a qiskit circuit to generate it from a minimal set of emitters.

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :param simple: Uses one register for both photons and emitters instead of seperate registers. Defaults to false
    :type simple: bool (optional, defaults to False)

    :return: Qiskit circuit corresponding to the protocol
    :rtype: QuantumCircuit

    """
    try:
        from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    except:
        ImportError('Qiskit failed to import')  
    n_e = num_emitters(state)
    n_p = state.size
    N = n_e+n_p
    if simple:
        q = QuantumRegister(N, 'q')
        c = ClassicalRegister(N, 'c')
        qc = QuantumCircuit(q,c)
    else:
        p = QuantumRegister(n_p, 'p')
        e = QuantumRegister(n_e, 'e')
        c = ClassicalRegister(N, 'c')
        qc = QuantumCircuit(p,e,c) 
    target_state = Stabilizer(N)
    for i in range(n_p):
        target_state.signvector[i] = state.signvector[i]
        for j in range(n_p):
            target_state.tab[i,j] = state.tab[i,j]
            target_state.tab[i,j+N] = state.tab[i,j+n_p]
    for h in range(n_p,0,-1):
        photonindex = h-1
        d = partial_height(target_state,h)-partial_height(target_state,h-1)
        if d<0:
            rows = []
            for i in range(N):
                toggler = True
                for j in range(n_p):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
                if toggler:
                    rows.append(i)
            sums=[]
            for row in rows:
                sum=0
                for j in range(n_p,N):
                    if target_state.tab[row,j]!=0 or target_state.tab[row,j+N]!=0:
                        sum+=1
                sums.append(sum)
            row = rows[sums.index(min(sums))]
            for i in range(n_p,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    emit = i
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        qc.h(i)
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        qc.s(i)
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        qc.h(i)
                        target_state.clifford('H',i)
                    break
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        qc.h(i)
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        qc.s(i)
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        qc.h(i)
                        target_state.clifford('H',i)
                    qc.cx(i,emit)
                    target_state.clifford('CNOT',i,emit)
            if target_state.signvector[row]==1:
                target_state.clifford('X',emit)
                qc.x(emit)
            target_state.clifford('H',emit)
            target_state.clifford('CNOT',emit,photonindex)
            qc.reset(emit)
            with qc.if_test((c, 1)):
                qc.x(photonindex)
            qc.measure(emit,c[0])
        for i in range(N):
            toggler = True
            if target_state.tab[i,photonindex]==0 and target_state.tab[i,photonindex+N]==0:
                toggler = False
            if toggler:
                for j in range(photonindex):
                    if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                        toggler = False
                        break
            if toggler:
                row = i
                break
        emit = -1
        if target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==0:
            qc.h(photonindex)
            target_state.clifford('H',photonindex)
        elif target_state.tab[row,photonindex]==1 and target_state.tab[row,photonindex+N]==1:
            qc.s(photonindex)
            target_state.clifford('S',photonindex)
            target_state.clifford('Z',photonindex)
            qc.h(photonindex)
            target_state.clifford('H',photonindex)

        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    qc.h(i)
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    qc.s(i)
                    target_state.clifford('S',i)
                    target_state.clifford('Z',i)
                    qc.h(i)
                    target_state.clifford('H',i)
                break
        if emit!= -1:
            for i in range(emit+1,N):
                if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                    if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                        qc.h(i)
                        target_state.clifford('H',i)
                    elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                        qc.s(i)
                        target_state.clifford('S',i)
                        target_state.clifford('Z',i)
                        qc.h(i)
                        target_state.clifford('H',i)
                    qc.cx(i,emit)
                    target_state.clifford('CNOT',i,emit)
            if target_state.signvector[row]==1:
                target_state.clifford('X',photonindex)
                qc.x(photonindex)
            target_state.clifford('CNOT',emit,photonindex)
            qc.cx(emit, photonindex)
        for i in range(N):
            if target_state.tab[i,photonindex+N]!=0 or target_state.tab[i,photonindex]!=0:
                if i!=row:
                    target_state.row_add(row,i)
    rref(target_state)
    sums = []
    for i in range(n_p,N):
        sum=0
        for j in range(n_p,N):
            if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                sum+=1
        sums.append(sum)
    
    if max(sums)==1:
        decoupled = True
    else:
        decoupled = False
    
    while not decoupled:
        for i in range(2,N):
            if i in sums:
                minimum = i
                break
        row = n_p+sums.index(minimum)
        for i in range(n_p,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                emit = i
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    qc.h(emit)
                    target_state.clifford('H',emit)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    qc.s(emit)
                    target_state.clifford('S',emit)
                    target_state.clifford('Z',emit)
                    qc.h(emit)
                    target_state.clifford('H',emit)
                break
        for i in range(emit+1,N):
            if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                if target_state.tab[row,i]==1 and target_state.tab[row,i+N]==0:
                    qc.h(i)
                    target_state.clifford('H',i)
                elif target_state.tab[row,i]==1 and target_state.tab[row,i+N]==1:
                    qc.s(i)
                    target_state.clifford('S',i)
                    target_state.clifford('Z',i)
                    qc.h(i)
                    target_state.clifford('H',i)
                target_state.clifford('CNOT',i,emit)
                qc.cx(i,emit)
        for i in range(n_p,N):
            if target_state.tab[i,emit]!=0 or target_state.tab[i,emit+N]!=0:
                if i!= row:
                    target_state.row_add(row,i)
        sums = []
        for i in range(n_p,N):
            sum=0
            for j in range(n_p,N):
                if target_state.tab[i,j]!=0 or target_state.tab[i,j+N]!=0:
                    sum+=1
            sums.append(sum)
        if max(sums)==1:
            decoupled = True
        else:
            decoupled = False
    for emitter in range(n_p,N):
        for row in range(N):
            toggler = True
            colcheck = True
            if target_state.tab[row,emitter]==0 and target_state.tab[row,emitter+N]==0:
                toggler = False
                colcheck = False
            if colcheck:
                for i in range(emitter):
                    if target_state.tab[row,i]!=0 or target_state.tab[row,i+N]!=0:
                        toggler = False
                        break
            if toggler:
                target_state.swap(row,emitter)
                break
    for i in range(n_p,N):
        if target_state.tab[i,i]!=0:
            if target_state.tab[i,i]==1 and target_state.tab[i,i+N]==0:
                qc.h(i)
                target_state.clifford('H',i)
            elif target_state.tab[i,i]==1 and target_state.tab[i,i+N]==1:
                qc.s(i)
                target_state.clifford('S',i)
                target_state.clifford('Z',i)
                qc.h(i)
                target_state.clifford('H',i)
    
    for i in range(n_p,N):
        if target_state.signvector[i]==1:
            target_state.clifford('X',i)
            qc.x(i)
    checker_state = Stabilizer(N)
    if np.array_equal(checker_state.tab,target_state.tab) and np.array_equal(checker_state.signvector,target_state.signvector):
        qc.data = qc.data[::-1]
        return qc
    else:
        print('Something went wrong')
        print(target_state.stabilizers())
        return None

def qiskit_circuit_solver_alternate(state: Stabilizer, simple: bool = False):
    """
    Given a stabilizer state, creates a qiskit circuit to generate it from a minimal set of emitters. Uses circuit_solver to generate the protocol instead of doing it simultaneously.

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :param simple: Uses one register for both photons and emitters instead of separate ones
    :type simple: bool (optional, defaults to False)

    :return: Qiskit circuit corresponding to the protocol
    :rtype: QuantumCircuit

    """
    try:
        from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    except:
        ImportError('Qiskit failed to import')  
    photons = state.size
    emitters = num_emitters(state)
    if simple:
        q = QuantumRegister(photons+emitters, 'q')
        c = ClassicalRegister(photons+emitters, 'c')
        qc = QuantumCircuit(q,c) 
    else:
        p = QuantumRegister(photons, 'p')
        e = QuantumRegister(emitters,'e')
        c = ClassicalRegister(photons+emitters, 'c')
        qc = QuantumCircuit(p,e,c)
    operations = circuit_solver(state)
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
            qc.measure(operation[1],c[0])
            with qc.if_test((c, 1)):
                qc.x(operation[2])
            qc.reset(operation[1])
    return qc

def num_cnots(state: Stabilizer):
    """
    Given a graph state with a certain fixed ordering, finds the number of emitters-emitter CNOT gates required to emit that state using this implementation of the solving algorithm

    :param state: Target photonic stabilizer state
    :type state: Stabilizer

    :return: Number of emitter-emitter CNOT gates needed to emit the input state using this version of the algorithm
    :rtype: int

    """
    procedure = circuit_solver(state)
    cnots = 0
    for i in range(len(procedure)):
        if procedure[i][0] == 'CNOT':
            cnots+=1
    return cnots


def emitter_cnot(state: Stabilizer):
    """
    Given a stabilizer state, finds the number of emitters and the number of emitter-emitter CNOT gates needed to implement the circuit

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :return: A list of the form [emitters,cnots]
    :rtype: list

    """   
    emitter = num_emitters(state)
    procedure = circuit_solver(state)
    cnots = 0
    for i in range(len(procedure)):
        if procedure[i][0] == 'CNOT':
            cnots+=1
    data = [emitter,cnots]
    return data

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    state = Stabilizer(edgelist = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]])