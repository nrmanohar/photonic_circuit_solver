"Contains the Stabilizer class to represent stabilizer states"
import numpy as np
import math

class Stabilizer:
    '''
    This is a class that encodes the stabilizer state in terms of its stabilizers. If no input is given, it will initialize a two qubit bell state. If only the n is given, it will initialize the n qubit computational zero state
    
    :param n: Number of qubits
    :type n: int, Optional (only if providing edgelist)

    :param stabs: The stabilizers, either in a string or a list, in the format 'XX,-YY' or '[XX,-YY]'. Optional
    :type stabs: list or string, optional (initializes zero state by default)

    :param edgelist: A list of edges for a graph state. Optional
    :type edgelist: list

    :cvar size: The number of qubits, initial value: n
    :cvar __stabs: The stabilizers of the state, initial value: stabs (note, this is a dunder attribute, can't be directly called outside the class. Use the .stabilizers() method instead)
    :cvar tab: The tableau of the state
    :cvar signvector: The signvector of the state
    :cvar gauss: A nxn Gaussian matrix (used for empty_column calculations)
    
    '''
    def __init__(self, n: int = None, stabs: str | list = None, edgelist: list[list] = None):
        """Constructor method

        """
        
        if edgelist is None:    
            if n is None and stabs is None:
                n = 2
                stabs = 'XX,ZZ'
            
            elif n is not None and stabs is None:
                stabs = []
                for i in range(n):
                    str = ''
                    for j in range(n):
                        if i==j:
                            str = str+'Z'
                        else:
                            str = str+'I'
                    stabs.append(str)
            

            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
            list = self.tableau()
            self.tab = list[0]
            self.signvector = list[1]
            try:
                assert self.square()
            except:
                raise ValueError("Invalid input, number of qubits not equal to number of stabilizers")
            try:
                assert not self.empty_column()
            except:
                raise ValueError("Invalid input, free qubit (all stabilizers for some qubit is the identity)")
            try:
                assert self.commuter()
            except:
                raise ValueError("Invalid Inputs, Stabilizers do not commute")
            try:
                assert self.linear_independence()
            except:
                raise ValueError("Invalid Inputs, Stabilizers are not independant")
        else:
            self.graph_state(edgelist = edgelist)

    def __str__(self) -> str:
        """
        Returns the stabilizers as a string and outputs them

        """
        stabs = self.stabilizers()
        stabstring = stabs[0]
        for i in range(1,len(stabs)):
            stabstring+=', '+stabs[i]
        return stabstring

    def __repr__(self) -> str:
        """
        Returns the stabilizers as a string and outputs them

        """
        stabs = self.stabilizers()
        stabstring = stabs[0]
        for i in range(1,len(stabs)):
            stabstring+=', '+stabs[i]
        return stabstring

    def square(self) -> bool:
        toggler = True
        for i in range(len(self.__stab)):
            str = self.__stab[i]
            str = str.lstrip('-')
            if len(str)!=len(self.__stab):
                return False
        return toggler

    def commuter(self) -> bool:
        """
        Tests whether the stabilizers commute with each other

        :return: Whether or not they commute
        :rtype: boolean
        """
        for i in range(self.size):
            toggler=0
            for j in range(i+1,self.size):
                for k in range(self.size):
                    if self.tab[i,k]==self.tab[j,k] and self.tab[i,k+self.size]==self.tab[j,k+self.size]:
                        toggler = toggler
                    elif (self.tab[i,k+self.size]==0 and self.tab[i,k]==0) or (self.tab[j,k+self.size]==0 and self.tab[j,k]==0):
                        toggler = toggler
                    else:
                        toggler = toggler+1
                if toggler%2 != 0:
                    return False
        return True
    def num_qubits(self) -> int:
        """
        Returns the size of the stabilizer (the number of qubits)

        :return: The size of the stabilizer
        :rtype: int
        """
        return self.size
    
    def graph_state(self,edgelist: list[list] = [[0,1],[1,2],[2,3],[3,4],[4,0]]):
        """
        Generates a graph state based on inputed edgelist

        :param edgelist: The list of connections, defaults to [[0,1],[1,2],[2,3],[3,4],[4,0]]
        :type edgelist: Nested list

        """
        num=0
        for i in range(len(edgelist)):
            for j in range(len(edgelist[i])):
                if edgelist[i][j]>num:
                    num = edgelist[i][j]
        self.size = num+1
        tab = np.zeros(2*self.size*self.size)
        tab = tab.reshape(self.size,2*self.size)
        for i in range(self.size):
            tab[i,i]=1
        for i in range(len(edgelist)):
            q1 = edgelist[i][0]
            q2 = edgelist[i][1]
            tab[q1,q2+self.size]=1
            tab[q2,q1+self.size]=1
        sign = np.zeros(self.size)
        self.tab = tab
        self.signvector = sign
    
    def tableau(self) -> list:
        """
        Converts the stabilizers to a tableau and signvector

        :return: A list contained the tableau and the signvector
        :rtype: list
        """
        tab = np.zeros(2*self.size*self.size)
        tab = tab.reshape(self.size,2*self.size)
        sign = np.zeros(self.size)
        for i in range(len(self.__stab)):
            if self.__stab[i][0]=='-':
                sign[i]=1
                self.__stab[i]=self.__stab[i][1:]
            for j in range(len(self.__stab[i])):
                if self.__stab[i][j].upper()=='I':
                    pass
                elif self.__stab[i][j].upper()=='X':
                    tab[i,j]=1
                elif self.__stab[i][j].upper()=='Z':
                    tab[i,j+self.size]=1
                elif self.__stab[i][j].upper()=='Y':
                    tab[i,j]=1
                    tab[i,j+self.size]=1
                else:
                    raise ValueError("Invalid stabilizer (not a Pauli string)")
        return [tab,sign]
    def stabilizers(self)-> list[str]:
        """
        Returns a list of the stabilizers of the state, as per the tableau

        :return: A list of stabilizers
        :rtype: list  
        
        """
        self.__stab = []
        for i in range(self.size):
            str = ""
            if self.signvector[i]==1:
                str = str+"-"
            for j in range(self.size):
                if self.tab[i,j]==0 and self.tab[i,j+self.size]==0:
                    str = str+"I"
                elif self.tab[i,j]==1 and self.tab[i,j+self.size]==0:
                    str = str+"X"
                if self.tab[i,j]==0 and self.tab[i,j+self.size]==1:
                    str = str+"Z"
                if self.tab[i,j]==1 and self.tab[i,j+self.size]==1:
                    str = str+"Y"
            self.__stab.append(str)
        return self.__stab
    def new_stab(self,size: int =None,newstabs: str|list = None):
        """
        Resets the stabilizer and new tableau associated with it

        :param size: The size of the new state
        :type size: int (optional)

        :param newstabs: The new stabilizers
        :type newstabs: string or list
        """
        if size is None and newstabs is None:
            size = 2
            newstabs = 'XX,ZZ'
        
        if size is not None and newstabs is None:
            newstabs = []
            for i in range(size):
                str = ''
                for j in range(size):
                    if i==j:
                        str = str+'Z'
                    else:
                        str = str+'I'
                newstabs.append(str)
        self.size = size
        try:
            self.__stab = newstabs.split(',')
        except:
            self.__stab = newstabs
        list = self.tableau()
        self.tab = list[0]
        self.signvector = list[1]
        try:
            assert self.square()
        except:
            raise ValueError("Invalid input, number of qubits not equal to number of stabilizers")
        try:
            assert not self.empty_column()
        except:
            raise ValueError("Invalid input, free qubit (all stabilizers for some qubit is the identity)")
        try:
            assert self.commuter()
        except:
            raise ValueError("Invalid Inputs, Stabilizers do not commute")
        try:
            assert self.linear_independence()
        except:
            raise ValueError("Invalid Inputs, Stabilizers are not independant")

    def clifford(self,type: str,q1: int,q2: int=None):
        """
        Applies a clifford gate to the stabilizer

        :param type: The clifford gate to be operated, 'H', 'X', 'Y', 'Z', 'CNOT' or 'CX', 'SWAP', 'CZ', or 'S'
        :type type: string

        :param q1: The qubit to operate on, or the control qubit for entangling gates
        :type q1: int

        :param q2: The qubit to target, defaults to None
        :type q2: int
        """
        if type.lower() == 'h':
            self.tab[:,[q1,q1+self.size]] = self.tab[:,[q1+self.size,q1]]
            for i in range(self.size):
                if self.tab[i,q1]*self.tab[i,q1+self.size]==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'cnot' or type.lower() == 'cx':
            if q2 == None:
                raise ValueError('Second qubit not specified for cnot gate')
            elif q1 == q2:
                pass
            else:
                for i in range(self.size):
                    self.tab[i,q2] = (self.tab[i,q1]+self.tab[i,q2])%2
                    self.tab[i,self.size+q1] = (self.tab[i,q1+self.size]+self.tab[i,q2+self.size])%2
                    if self.tab[i,q1]==1 and self.tab[i,q2+self.size]==1:
                        if self.tab[i,q2]==self.tab[i,self.size+q1]:
                            self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'swap':
            if q2 == None:
                raise ValueError('Second qubit not specified for swap gate')
            elif q1 == q2:
                pass
            else:
                self.tab[:,[q1,q2]] = self.tab[:,[q2,q1]]
                self.tab[:,[q1+self.size,q2+self.size]] = self.tab[:,[q2+self.size,q1+self.size]]
        elif type.lower() == 'z':
            for i in range(self.size):
                if self.tab[i,q1]==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'x':
            for i in range(self.size):
                if self.tab[i,q1+self.size]==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'y':
            for i in range(self.size):
                if (self.tab[i,q1]+self.tab[i,q1+self.size])==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 's':
            for i in range(self.size):
                if self.tab[i,q1]==1:
                    self.signvector[i]=(self.signvector[i]+self.tab[i,q1+self.size])%2
                    self.tab[i,q1+self.size] = (self.tab[i,q1+self.size]+1)%2
        elif type.lower() == 'cz':
            if q2 == None:
                raise ValueError('Second qubit not specified for cnot gate')
            else:
                self.clifford('h',q2)
                self.clifford('cnot',q1,q2)
                self.clifford('h',q2)
        else:
            raise ValueError("Invalid gate inputted. Valid types are 'H' for Hadamard, 'S' for the phase gate, 'CNOT' for the Control Not, 'CZ' for the Control Z.")
    
    def row_commute(self, stab1: str, stab2: str) -> bool:
        """
        Checks if two stabilizers commute

        :param stab1: The first stabilizer
        :type stab1: string

        :param stab2: The second stabilizer
        :type stab2: string

        :return: Whether the stabilizers commute
        :rtype: boolean
        
        """
        stab1 = stab1.lstrip('-')
        stab2 = stab2.lstrip('-')
        if len(stab1)!=len(stab2):
            raise ValueError("Your stabilizers aren't of same length")
        toggler = 0
        for i in range(len(stab1)):
            if stab1[i] != 'I' and stab2[i] != 'I' and stab1[i] != stab2[i]:
                toggler+= 1
        if toggler%2 == 0:
            return True
        else:
            return False
    
    def measurement(self, stabilizers: list|str, outcomes:list = None):
        """
        Implements a measurement of a Pauli string with a specified outcome (Pauli string shouldn't have a sign on it, use the outcomes to denote sign. Ex, measuring -X is equivalent to measurement 'X', [1])

        If the state is measured by a Pauli that is a stabilizer of the state (modulo sign) then the measurement does nothing and is ignored. 
        Note, since the measurement stabilizer is already a member of the stabilizer group, there is only one valid measurement outcome. As of right now, the code doesn't verify that whether what you provided matches the expected outcome.

        :param stabilizers: The set of stabilizers that were measured, in order
        :type stabilizers: str

        :param outcomes: The outcomes of these measurements (defaults to 0 for all, can put 0 or 1)
        :type outcomes: list
        
        """
        try:
            stabilizers = stabilizers.split(',')
        except:
            stabilizers = list(stabilizers)

        for i in range(len(stabilizers)):
            if len(stabilizers[i])!= self.size:
                raise ValueError('Stabilizers are wrong, inaccurate size')

        if outcomes == None:
            outcomes = [0 for i in range(len(stabilizers))]
        for i in range(len(stabilizers)):
            stabs = self.stabilizers()
            for j in range(len(stabs)):
                if not self.row_commute(stabs[j],stabilizers[i]):
                    index = j
                    break
            if index is not None:
                for k in range(index+1,len(stabs)):
                    if not self.row_commute(stabs[k],stabilizers[i]):
                        self.row_add(index,k)
                if outcomes[i]==1:
                    self.signvector[index]=1
                for k in range(self.size):
                    if stabilizers[i][k] == 'I':
                        self.tab[index,k] = 0
                        self.tab[index,k+self.size] = 0
                    elif stabilizers[i][k]=='X':
                        self.tab[index,k] = 1
                        self.tab[index,k+self.size] = 0
                    elif stabilizers[i][k]=='Z':
                        self.tab[index,k] = 0
                        self.tab[index,k+self.size] = 1
                    elif stabilizers[i][k]=='Y':
                        self.tab[index,k] = 1
                        self.tab[index,k+self.size] = 1
            elif index is None:
                pass
            print(self.stabilizers())



    def report(self):
        """
        Prints the tableau and the signvector

        """
        print(self.tab)
        print(self.signvector)

    def gaussian(self):
        """
        Generates an array that contains information about where stabilizers are known
        
        """
        self.gauss = np.zeros(self.size*self.size)
        self.gauss = self.gauss.reshape(self.size,self.size)
        for i in range(self.size):
            for j in range(self.size):
                if self.tab[i,j]==1 or self.tab[i,j+self.size]==1:
                    self.gauss[i,j]=1
    
    def empty_column(self) -> bool:
        """
        Tests whether there are any empty stabilizers (free qubits)

        :return: Whether there is an empty column or not
        :rtype: boolean
        """
        self.gaussian()
        zed = self.gauss.sum(axis=0)
        empty = False
        for i in range(self.size):
            if zed[i]==0:
                empty = True
        return empty

    def linear_independence(self):
        """
        Checks if the generators are linearly independent

        """
        rank = np.linalg.matrix_rank(self.tab)
        rank = int(rank)
        if rank == self.size:
            return True
        else:
            return False

    def row_add(self,row1: int, row2: int):
        """
        Multiplies two stabilizers in the tableau together, and puts them into the second row

        """
        if row1==row2:
            pass
        elif row1>= self.size or row2>=self.size:
            pass
        else:
            phase_tracker = 0
            for i in range(self.size):
                if self.tab[row1,i]==0 and self.tab[row1,i+self.size]==0:
                    pass
                elif self.tab[row2,i]==0 and self.tab[row2,i+self.size]==0:
                    self.tab[row2,i]=self.tab[row1,i]
                    self.tab[row2,i+self.size]=self.tab[row1,i+self.size]
                elif self.tab[row1,i]==self.tab[row2,i] and self.tab[row1,i+self.size]==self.tab[row2,i+self.size]:
                    self.tab[row2,i] = 0
                    self.tab[row2,i+self.size] = 0

                else:
                    if self.tab[row1,i]==0 and self.tab[row1,i+self.size]==1:
                        if self.tab[row2,i]==1 and self.tab[row2,i+self.size]==0:
                            phase_tracker += 1
                        else:
                            phase_tracker -= 1
                    elif self.tab[row1,i]==1 and self.tab[row1,i+self.size]==0:
                        if self.tab[row2,i]==0 and self.tab[row2,i+self.size]==1:
                            phase_tracker -= 1
                        else:
                            phase_tracker +=1
                    else:
                        if self.tab[row2,i]==0 and self.tab[row2,i+self.size]==1:
                            phase_tracker += 1
                        else:
                            phase_tracker -= 1
                    self.tab[row2,i] = (self.tab[row2,i]+self.tab[row1,i])%2
                    self.tab[row2,i+self.size] = (self.tab[row2,i+self.size]+self.tab[row1,i+self.size])%2
            phase_tracker = (phase_tracker%4)/2
            self.signvector[row2] = (self.signvector[row1]+self.signvector[row2]+phase_tracker)%2
                

    def circuit_builder(self):
        """
        Uses reverse operations to build the stabilizer state (requires qiskit)

        :return: A Qiskit circuit that makes the stabilizer
        :rtype: QuantumCircuit        
        """
        reference = np.copy(self.tab)
        sign = np.copy(self.signvector)
        try:
            from qiskit import QuantumCircuit, transpile
        except:
            raise ImportError('Qiskit failed to Import')
        qc = QuantumCircuit(self.size)
        self.rref()
        sums = []
        for i in range(self.size):
            sum=0
            for j in range(self.size):
                if self.tab[i,j]!=0 or self.tab[i,j+self.size]!=0:
                    sum+=1
            sums.append(sum)
        if max(sums)==1:
            decoupled = True
        else:
            decoupled = False
        while not decoupled:
            minimum = min(i for i in sums if i > 1)
            row = sums.index(minimum)
            xcount = 0
            for i in range(self.size):
                xcount+=self.tab[row,i]

            if self.size-xcount>xcount:
                for i in range(self.size):
                    if self.tab[row,i]!=0 or self.tab[row,i+self.size]!=0:
                        if self.tab[row,i]==1 and self.tab[row,i+self.size]==0:
                            self.clifford('h',i)
                            qc.h(i)
                        elif self.tab[row,i]==1 and self.tab[row,i+self.size]==1:
                            self.clifford('s',i)
                            self.clifford('z',i)
                            qc.s(i)
                            self.clifford('h',i)
                            qc.h(i)
                        emit = i
                        break
                for i in range(emit+1,self.size):
                    if self.tab[row,i]!=0 or self.tab[row,i+self.size]!=0:
                        if self.tab[row,i]==1 and self.tab[row,i+self.size]==0:
                            self.clifford('h',i)
                            qc.h(i)
                        elif self.tab[row,i]==1 and self.tab[row,i+self.size]==1:
                            self.clifford('s',i)
                            self.clifford('z',i)
                            qc.s(i)
                            self.clifford('h',i)
                            qc.h(i)
                        self.clifford('cnot',i,emit)
                        qc.cx(i,emit)
            else:
                for i in range(self.size):
                    if self.tab[row,i]!=0 or self.tab[row,i+self.size]!=0:
                        if self.tab[row,i]==0 and self.tab[row,i+self.size]==1:
                            self.clifford('h',i)
                            qc.h(i)
                        elif self.tab[row,i]==1 and self.tab[row,i+self.size]==1:
                            self.clifford('s',i)
                            self.clifford('z',i)
                            qc.s(i)
                        emit = i
                        break
                for i in range(emit+1,self.size):
                    if self.tab[row,i]!=0 or self.tab[row,i+self.size]!=0:
                        if self.tab[row,i]==0 and self.tab[row,i+self.size]==1:
                            self.clifford('h',i)
                            qc.h(i)
                        elif self.tab[row,i]==1 and self.tab[row,i+self.size]==1:
                            self.clifford('s',i)
                            self.clifford('z',i)
                            qc.s(i)
                        self.clifford('cnot',emit,i)
                        qc.cx(emit,i)
            for i in range(self.size):
                if self.tab[i,emit]!=0 or self.tab[i,emit+self.size]!=0:
                    self.row_add(row,i)    
            sums = []
            for i in range(self.size):
                sum=0
                for j in range(self.size):
                    if self.tab[i,j]!=0 or self.tab[i,j+self.size]!=0:
                        sum+=1
                sums.append(sum)    
            if max(sums)==1:
                decoupled = True
            else:
                decoupled = False
        self.rref()
        for i in range(self.size):
            if self.tab[i,i]==1 and self.tab[i,i+self.size]==0:
                self.clifford('H',i)
                qc.h(i)
            elif self.tab[i,i]==1 and self.tab[i,i+self.size]==1:
                self.clifford('S',i)
                self.clifford('Z',i)
                self.clifford('H',i)
                qc.s(i)
                qc.h(i)
            if self.signvector[i]==1:
                self.clifford('X',i)
                qc.x(i)
        self.tab = np.copy(reference)
        self.signvector = np.copy(sign)
        qc.data = qc.data[::-1]
        basis_set = ['id', 'cx', 'h', 's', 'x', 'z', 'y', 'cz']
        return transpile(qc, basis_gates = basis_set)

    def draw_circuit(self, style = 'text', save = None):
        """
        Draws a circuit that can generate the given stabilizer state (requires qiskit, matplotlib and pylatexenc package)

        :param style: The type of output, 'mpl' for matplotlib, 'text' for ASCII drawing, 'latex_source' for raw latex output
        :type style: string, optional. Defaults to 'text'

        :param save: If you want to save the file to something (optional)
        :type save: string

        """
        if style == 'mpl':
            try:
                import matplotlib.pyplot as plt
                try:
                    circ = self.circuit_builder()
                    circ.draw(output = style, filename = save)
                    plt.show()
                except:
                    print("pylatexenc likely not installed")
            except:
                print("matplotlib package not installed")
        elif style == 'text':
            circ = self.circuit_builder()
            circ.draw(output = style, filename = save)
        elif style == 'latex_source':
            circ = self.circuit_builder()
            circ.draw(output = style, filename = save)

    
    def qiskit_stabilizers(self):
        """
        Asks Qiskit to return the stabilizers (requires qiskit)

        :return: A qiskit stabilizer state representation
        :rtype: StabilizerState (qiskit)

        """
        circ = self.circuit_builder()
        try:
            from qiskit.quantum_info import StabilizerState
        except:
            raise ImportError('Qiskit.quantum_info not installed')
        try:
            stab = StabilizerState(circ)
            return stab
        except:
            pass

    def stabilizer_measurement(self):
        """
        A circuit to measure the associated stabilizers of this state (requires qiskit)

        :return: A qiskit circuit for measureing stabilizer
        :rtype: QuantumCircuit

        """
        try:
            from qiskit import QuantumCircuit,ClassicalRegister
        except:
            raise ImportError('Qiskit not installed')
        qs = QuantumCircuit(2*self.size)
        bits = []
        for i in range(self.size):
            bits.append(i)
        reg = ClassicalRegister(self.size)
        qs.add_register(reg)
        stabs = self.stabilizers()
        for i in range(self.size):
            qs.h(self.size+i)
        for i in range(self.size):
            stabs[i]=stabs[i].lstrip('-')
            for j in range(self.size):
                if stabs[i][j]=='X':
                    qs.cx(self.size+i,j)
                elif stabs[i][j]=='Z':
                    qs.cz(self.size+i,j)
                elif stabs[i][j]=='Y':
                    qs.cy(self.size+i,j)
        
        for i in range(self.size):
            qs.h(self.size+i)
        for i in range(self.size):
            if self.signvector[i]==1:
                qs.x(self.size+i)
        for i in range(self.size):
            qs.measure(i+self.size,i)

        return qs
    
    def build_and_measure(self):
        """
        A circuit to implement the circuit and then to measure the associated stabilizers (requires qiskit).

        :return: A qiskit circuit for measureing stabilizer
        :rtype: QuantumCircuit

        """
        try:
            from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
        except:
            ImportError('Qiskit not imported')
        circ = self.circuit_builder()
        qs = QuantumCircuit(2*self.size)
        bits = []
        for i in range(self.size):
            bits.append(i)
        qs = qs.compose(circ,bits)
        reg = ClassicalRegister(self.size)
        qs.add_register(reg)
        qs.barrier()
        stabs = self.stabilizers()
        for i in range(self.size):
            qs.h(self.size+i)
        for i in range(self.size):
            stabs[i]=stabs[i].lstrip('-')
            for j in range(self.size):
                if stabs[i][j]=='X':
                    qs.cx(self.size+i,j)
                elif stabs[i][j]=='Z':
                    qs.cz(self.size+i,j)
                elif stabs[i][j]=='Y':
                    qs.cy(self.size+i,j)
        for i in range(self.size):
            qs.h(self.size+i)
        for i in range(self.size):
            if self.signvector[i]==1:
                qs.x(self.size+i)
        for i in range(self.size):
            qs.measure(i+self.size,i)
        
        return qs
    def swap(self, r1: int,r2: int):
        """
        Swaps two rows in the stabilizer

        :param r1: The first row
        :type type: int

        :param r2: The second row
        :type q1: int
        """
        self.tab[[r1, r2]] = self.tab[[r2, r1]]
        self.signvector[[r1, r2]] = self.signvector[[r2, r1]]
    def flip(self):
        """
        Flips the tableau over

        """
        self.tab = np.flip(self.tab,axis=0)
        self.signvector = np.flip(self.signvector,axis=0)
    def clone(self):
        """
        Generates a copy of the stabilizer state

        """
        newstab = self.stabilizers()
        int = self.size
        state = Stabilizer(n=int,stabs=newstab)
        return state
    def rref(self, KU:int = 0, NL:int = 0):
        """
        Implements the RREF gauge procedure (details in doi.org/10.1088/1367-2630/7/1/170)

        :param KU: Starting row
        :type KU: int

        :param NL: Starting column
        :type NL: int


        """    

        N=self.size
        K=N
        while NL<N-1 and KU<K-1:
            zeroitem = True
            oneitem = False
            twoitem = False
            r1=N
            r2=N
            for k in range(KU,K):
                if self.tab[k,NL]!=0 or self.tab[k,NL+N]!=0:
                    r1 = k
                    zeroitem = False
                    oneitem = True
                    break
            for k in range(r1,K):
                if self.tab[k,NL]!=0 or self.tab[k,NL+N]!=0:
                    if self.tab[k,NL]!=self.tab[r1,NL] or self.tab[k,NL+N]!=self.tab[r1,NL+N]:
                        r2 = k
                        oneitem = False
                        twoitem = True
                        break
            if zeroitem:
                NL+=1
            elif oneitem:
                self.swap(KU,r1)
                for i in range(KU+1,K):
                    if self.tab[i,NL]!=0 or self.tab[i,NL+N]!=0:
                        self.row_add(KU,i)
                KU+=1
                NL+=1
            elif twoitem:
                self.swap(KU,r1)
                self.swap(KU+1,r2)
                for i in range(KU+2,K):
                    if self.tab[i,NL]!=0 or self.tab[i,NL+N]!=0:
                        if self.tab[i,NL]==self.tab[KU,NL] and self.tab[i,NL+N]==self.tab[KU,NL+N]:
                            self.row_add(KU,i)
                        elif self.tab[i,NL]==self.tab[KU+1,NL] and self.tab[i,NL+N]==self.tab[KU+1,NL+N]:
                            self.row_add(KU+1,i)
                        else:
                            self.row_add(KU,i)
                            self.row_add(KU+1,i)
                NL+=1
                KU+=2
if __name__ == "__main__":
    # Do something if this file is invoked on its own
    state = Stabilizer(edgelist = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]])