
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import QFT

import numpy as np

import utilities_multiples as ut_multiples
import operations as ops


def oracle_multiples(k:int, nqubits_input:int, approx_QFT:int=0, oracle:QuantumCircuit=None,
                     qubits_oracle:list=None, name:str=None, init_H:bool=False, classic_register:bool=None):
    '''
    Builder of the multiples oracle.

    Input:
        - k (int): Number whose multiples are to be computed.

        - nqubits_input (int): Number of qubits for the input register.
        These are the qubits in whose states the multiples of k are to be sought.

        - approx_QFT (int): Approximation degree of the QFTs for the additions.
        By default is 0 (no approximation).

        - oracle (QuantumCircuit): Oracle extra which is to be applied in combination with the rest.
        By default is None, so all the multiples would get a \pi-phase.

        - qubits_oracle (list): Qubits on which the oracle extra is to be applied.
        By default is None. If it is not provided, the oracle would be
        applied to all the qubits in the input register.

        - name (str): Name of the circuit/oracle. By default None. If none is provided then
        the number assigned is "Multiples of k".

        - init_H (bool): Bool to say wheter H gates are applied to the input qubits (True), in order
        to obtain the whole circuit, or just the oracle (False), to just get the oracle for marking.

    Output:
        - circuit (QuantumCircuit): QuantumCircuit which marks states representing natural numbers
        which are multiples of k from 0 to 2^{nqubits_input}-1.
        Note: In case oracle is provided, further changes may be done to the final states.
    '''

    # Number of qubits on the remainders register
    n_remainders = np.ceil(np.log2(k)).astype(int)


    if classic_register is None:
        classic_register = not init_H
    else:
        pass

    # Build circuit with the respective registers
    input_register = QuantumRegister(nqubits_input, 'input')
    remainders = QuantumRegister(n_remainders + 1, 'remainders')
    range_remainders = list(range(1, len(remainders)+1))
    ancilla = QuantumRegister(1, 'ancilla')

    if name:
        name = name
    else:
        name = ' Multiples of %d '%k

    if classic_register:
        classic_multiples = ClassicalRegister(nqubits_input)
        circuit = QuantumCircuit(input_register,
                                remainders,
                                ancilla,
                                classic_multiples,
                                name = name)
    else:
        circuit = QuantumCircuit(input_register,
                                remainders,
                                ancilla,
                                classic_multiples,
                                name = name)
    
    # Full superposition to the input register
    if init_H:
        circuit.h(input_register)
    else:
        pass

    # Equivalent to QFT when starting in 0
    circuit.h(remainders)

    # Obtain the remainders of the powers of 2
    list_remainders = ut_multiples.get_remainders_power_2(k, nqubits_input-1)

    # Dictionary of circuits for sum (in order not to construct repeteadly the circuits)
    circuits_sums = {}

    # For loop to build the circuits
    for i, remainder in enumerate(list_remainders):
        # Construct the circuits for the remainders (if not in the dictionary)
        if remainder not in circuits_sums:
            circuits_sums[remainder] = ops.c_phase_add_mod_K(target_register=range_remainders, control=0,
                                                        ancilla_qubit=len(range_remainders)+1, num_sum=remainder,
                                                        K=k, approx_QFT=approx_QFT, include_QFTs=False)
        else:
            pass

        # Append sum circuits to main circuit
        circuit.append(circuits_sums[remainder], [input_register[i]] + [remainders[j] for j in range(len(remainders))] + [ancilla[0]])

    # Inverse QFT to return to the amplitude space and have the remainders in amplitude and not in phase
    circuit.append(QFT(num_qubits=len(remainders), do_swaps=False, inverse=True, approximation_degree=approx_QFT), remainders)


    # X gates to use the |0...0> states to activate the gates
    circuit.x(remainders[:-1])

    # Apply the extra oracle
    if oracle:
        # Configure qubits on which to apply the oracle extra
        if not qubits_oracle:
            qubits_oracle = input_register
        else:
            pass
        circuit.append(oracle, qubits_oracle)
    else: #Otherwise, just mark the multiples
        circuit.append(ut_multiples.multi_control_z(len(remainders)-1), remainders[:-1])
    
    # Undo the X gates
    circuit.x(remainders[:-1])

    # Apply QFT to remainders register to sum on the phase space
    circuit.append(QFT(num_qubits=len(remainders), do_swaps=False, inverse=False, approximation_degree=approx_QFT), remainders)


    # For loop to build the circuits for undo the modulated sums
    # It uses the property of remainders ring that the inverse operation of +3 mod 7, is +4 mod 7
    # Generalization: + p mod k is the opposite to +(k-p) mod k
    for i, remainder_original in enumerate(list_remainders):
        # Calculate the inverse in the remainders ring
        inverse_remainder = k-remainder_original
        
        # Construct the circuits for the remainders (if not in the dictionary)
        if inverse_remainder not in circuits_sums:
            circuits_sums[inverse_remainder] = ops.c_phase_add_mod_K(target_register=range_remainders, control=0,
                                                                 ancilla_qubit=len(range_remainders)+1,
                                                                 num_sum=inverse_remainder, K=k,
                                                                 approx_QFT=approx_QFT, include_QFTs=False)
        else:
            pass

        # Append sum circuits to main circuit
        circuit.append(circuits_sums[inverse_remainder], [input_register[i]] + [remainders[j] for j in range(len(remainders))] + [ancilla[0]])


    # Inverse QFT to return to amplitude space
    circuit.append(QFT(num_qubits=len(remainders), do_swaps=False, inverse=True, approximation_degree=approx_QFT), remainders)

    # The oracle is built such that at this point the remainders register is at full superposition
    # Some states have phases, hence it should be studied whether this inverse QFT may be replaced with full Hs.

    return circuit

