"Contains the classes and function to manipulate stabilizer and graph states"
import numpy as np
import math

class Stabilizer:
    '''
    This is a class that encodes the stabilizer state in terms of its stabilizers. If no input is given, it will initialize a bell state. If only the n is given, it will initialize n qubits in the 0 state
    
    :param n: Number of qubits
    :type n: int, Optional (if providing edgelist)

    :param stabs: The stabilizers, either in a string or a list, in the format 'XX,-YY' or '[XX,-YY]' (case sensitive). Optional
    :type stabs: list or string, optional (initializes zero state by default)

    :param edgelist: A list of edges for a graph state. Optional
    :type edgelist: list

    :cvar size: The number of qubits, initial value: n
    :cvar __stabs: The stabilizers of the state, initial value: stabs (note, this is a dunder attribute, can't be directly called outside the class. Use the .stabilizers() method instead)
    :cvar tab: The tableau of the state
    :cvar signvector: The signvector of the state
    :cvar gauss: A nxn Gaussian matrix (used for empty_column calculations)
    
    '''
    def __init__(self, n = None, stabs = None, edgelist = None):
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

    def __str__(self):
        """
        Returns the stabilizers as a string and outputs them

        """
        stabs = self.stabilizers()
        stabstring = stabs[0]
        for i in range(1,len(stabs)):
            stabstring+=', '+stabs[i]
        return stabstring

    def __repr__(self):
        """
        Returns the stabilizers as a string and outputs them

        """
        stabs = self.stabilizers()
        stabstring = stabs[0]
        for i in range(1,len(stabs)):
            stabstring+=', '+stabs[i]
        return stabstring

    def square(self):
        toggler = True
        for i in range(len(self.__stab)):
            str = self.__stab[i]
            str = str.lstrip('-')
            if len(str)!=len(self.__stab):
                return False
        return toggler

    def commuter(self):
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
    def num_qubits(self):
        """
        Returns the size of the stabilizer (the number of qubits)

        :return: The size of the stabilizer
        :rtype: int
        """
        return self.size
    
    def graph_state(self,edgelist = [[0,1],[1,2],[2,3],[3,4],[4,0]]):
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
    
    def tableau(self):
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
                if self.__stab[i][j]=='I':
                    pass
                elif self.__stab[i][j]=='X':
                    tab[i,j]=1
                elif self.__stab[i][j]=='Z':
                    tab[i,j+self.size]=1
                elif self.__stab[i][j]=='Y':
                    tab[i,j]=1
                    tab[i,j+self.size]=1
                else:
                    raise ValueError("Invalid stabilizer (not a Pauli string)")
        return [tab,sign]
    def stabilizers(self):
        """
        Returns a list of the stabilizers of the state, as per the tableau

        :return: A list of operations to take a standard state to the given stabilizer state
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
    def new_stab(self,size=None,newstabs=None, ignore_commute = False):
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

    def clifford(self,type,q1,q2=None):
        """
        Applies a clifford gate to the stabilizer

        :param type: The clifford gate to be operated, 'H', 'X', 'Y', 'Z', 'CNOT', 'CZ', or 'S'
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
        elif type.lower() == 'cnot':
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
                self.clifford('h', q2)
        else:
            raise ValueError("Invalid gate inputted. Valid types are 'H' for Hadamard, 'S' for the phase gate, 'CNOT' for the Control Not, 'CZ' for the Control Z.")
    
    def row_commute(self, stab1, stab2):
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
    
    def measurement(self, stabilizers, outcomes=None):
        try:
            stabilizers = stabilizers.split(',')
        except:
            stabilizers = list(stabilizers)

        for i in range(len(stabilizers)):
            if len(stabilizers[i])!= self.size:
                raise ValueError('Stabilizers are wrong, inaccurate size')

        if outcomes == None:
            outcomes = [0 for i in range(len(stabilizers))]
        stabs = self.stabilizers()
        for i in range(len(stabilizers)):
            for j in range(len(stabs)):
                if not self.row_commute(stabs[j],stabilizers[i]):
                    index = j
                    break
            for k in range(index+1,len(stabs)):
                if not self.row_commute(stabs[k],stabilizers[i]):
                    self.row_add(index,k)

            stabs = self.stabilizers()
            if outcomes[i]==1:
                stabilizers[i] = '-'+stabilizers[i]
            stabs[index]=stabilizers[i]
            self.new_stab(self.size,stabs,True)

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
    
    def empty_column(self):
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

    def row_add(self,row1,row2):
        """
        Multiplies two stabilizers in the tableau together, specifying a new stabilizer, and puts them into the second row

        """
        if row1==row2:
            pass
        elif row1>= self.size or row2>=self.size:
            pass
        else:
            phase_tracker = 1
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
                            phase_tracker = phase_tracker*complex(1j)
                        else:
                            phase_tracker = phase_tracker*complex(-1j)
                    elif self.tab[row1,i]==1 and self.tab[row1,i+self.size]==0:
                        if self.tab[row2,i]==0 and self.tab[row2,i+self.size]==1:
                            phase_tracker = phase_tracker*complex(-1j)
                        else:
                            phase_tracker = phase_tracker*complex(1j)
                    else:
                        if self.tab[row2,i]==0 and self.tab[row2,i+self.size]==1:
                            phase_tracker = phase_tracker*complex(1j)
                        else:
                            phase_tracker = phase_tracker*complex(-1j)
                    self.tab[row2,i] = (self.tab[row2,i]+self.tab[row1,i])%2
                    self.tab[row2,i+self.size] = (self.tab[row2,i+self.size]+self.tab[row1,i+self.size])%2
        phase_tracker = (1-1*np.real(phase_tracker))/2
        self.signvector[row2] = (self.signvector[row1]+self.signvector[row2]+phase_tracker)%2
                

    def circuit_builder(self):
        """
        Uses reverse operations to build the stabilizer state

        :return: A Qiskit circuit that makes the stabilizer
        :rtype: QuantumCircuit        
        """
        reference = np.copy(self.tab)
        sign = np.copy(self.signvector)
        rev_operations = []

        broken = False

        for i in range(self.size):
            if self.tab[i,i]==0:
                if self.tab[i,i+self.size]==1:
                    rev_operations.append(['H',i])
                    self.clifford('H',i)
            if self.tab[i,i]==0:
                for j in range(i+1,self.size):
                    if self.tab[j,i]==1:
                        self.tab[[i,j]]=self.tab[[j,i]]
                        self.signvector[[i,j]]=self.signvector[[j,i]]
                        break
            if self.tab[i,i]==0:
                for j in range(i+1,self.size):
                    if self.tab[j,i+self.size]==1:
                        self.tab[[i,j]]=self.tab[[j,i]]
                        self.signvector[[i,j]]=self.signvector[[j,i]]
                        rev_operations.append(['H',i])
                        self.clifford('H',i)
                        break
            if self.tab[i,i]==0:
                for j in range(i):
                    if self.tab[j,i+self.size]==1:
                        self.row_add(j,i)
                        rev_operations.append(['H',i])
                        self.clifford('H',i)
                        break
            if self.tab[i,i]==0:
                broken = True
                break
            elif self.tab[i,i]==1:
                for j in range(self.size):
                    if self.tab[i,j]==1 and j!=i:
                        rev_operations.append(["CNOT",i,j])
                        self.clifford("CNOT",i,j)
                

        if broken:
            self.tab = np.copy(reference)
            self.signvector = np.copy(sign)
            print("Something went wrong in the building procedure. Check your stabilizers and maybe reformat them and try again")
            return None

        for i in range(self.size):
            if self.tab[i,i+self.size]==1:
                rev_operations.append(["S",i])
                self.clifford("S",i)

        for i in range(self.size):
            for j in range(self.size):
                if self.tab[i,j+self.size]==1:
                    rev_operations.append(["CZ",i,j])
                    self.clifford("CZ",i,j)

        for i in range(self.size):
            self.clifford('H',i)
            rev_operations.append(['H',i])
        
        for i in range(self.size):
            if self.signvector[i]==1:
                rev_operations.append(['X',i])
                self.clifford('X',i)
        
        self.tab = np.copy(reference)
        self.signvector = np.copy(sign)
        
        rev_operations.reverse()

        circuit = QuantumCircuit(self.size)

        
        for i in range(len(rev_operations)):
            if rev_operations[i][0]=='H':
                circuit.h(rev_operations[i][1])
            elif rev_operations[i][0]=='S':
                circuit.s(rev_operations[i][1])
                circuit.z(rev_operations[i][1])
            elif rev_operations[i][0]=='X':
                circuit.x(rev_operations[i][1])
            elif rev_operations[i][0]=='Y':
                circuit.y(rev_operations[i][1])
            elif rev_operations[i][0]=='Z':
                circuit.z(rev_operations[i][1])
            elif rev_operations[i][0]=='CNOT':
                circuit.cnot(rev_operations[i][1],rev_operations[i][2])
            elif rev_operations[i][0]=='CZ':
                circuit.cz(rev_operations[i][1],rev_operations[i][2])
        return circuit

    def draw_circuit(self, style = 'mpl', save = None):
        """
        Draws a circuit that can generate the given stabilizer state (requires matplotlib and pylatexenc package)

        :param style: The type of output, 'mpl' for matplotlib, 'text' for ASCII drawing, 'latex_source' for raw latex output
        :type style: String, optional. Defaults to 'mpl'

        :param save: If you want to save the file to something (optional)
        :type save: String

        """
        if style == 'mpl':
            try:
                import matplotlib
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
        Asks Qiskit to return the stabilizers

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
        A circuit to measure the associated stabilizers of this state

        :return: A qiskit circuit for measureing stabilizer
        :rtype: QuantumCircuit

        """
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
        A circuit to implement the circuit and then to measure the associated stabilizers.

        :return: A qiskit circuit for measureing stabilizer
        :rtype: QuantumCircuit

        """
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
    def swap(self, r1,r2):
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

def rref(state):
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

def heightfunction(state):
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

def partial_height(state,index):
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

def plot_height(state):
    """
    Given a stabilizer state, finds the height function (bipartite entanglement entropy) of the state (formula in doi.org/10.1038/s41534-022-00522-6, equation 1) and plots it (note, you need Matplotlib installed)

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
        plt.show()
    except:
        pass


def num_emitters(state):
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

def circuit_solver(state):
    """
    Given a stabilizer state, finds a circuit to generate it from a minimal set of emitters (details in doi.org/10.1038/s41534-022-00522-6)

    :param state: Stabilizer state from which to calculate the height function
    :type state: Stabilizer

    :return: Procedure to generate the state
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
        for i in range(2,N):
            if i in sums:
                minimum = i
                break
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

def num_cnots(state):
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


def emitter_cnot(state):
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