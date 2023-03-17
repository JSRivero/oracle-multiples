'''
Script containing the functions needed to perform the basic operations
of adding and modulo adding (using Draper algorithm)
'''

from qiskit import QuantumCircuit
from qiskit.circuit.library import QFT

import utilities_multiples as ut_multiples

import numpy as np
import math



def get_angles_phase_addition(number_:int, nqubits:int) -> np.array:
    '''
    Calculates and returns an array of the angles required to add the value
    number to a register of nqubits qubits by using the Draper Algorithm.

    Input:
        - number (int): Value to be added.

        - nqubits (int): Number of qubits on which the addition is performed.
         This number must be at least the number of bits needed to number in binary,
         otherwise an error will be arisen.

    Outputs:
        - angles (np.array): Array with the angles needed to perform the Draper Addition.
    '''

    num_binary = ut_multiples.to_binary(number=number_, nbits=nqubits)
    angles = np.zeros(nqubits)

    angles = np.array([sum([
                    math.pow(2, nqubits-i-1-j) for j in range(nqubits-i-1, nqubits) if num_binary[j]=='1'])*np.pi
                     for i in range(nqubits)])

    return angles


def phase_add(circuit:QuantumCircuit, target_register:list, num_sum:int, inv:bool=False):
    '''
    This function performs the quantum addition using Drapper addition.

    This addition is realized using direct QFT before and inverse QFT after.
    #! IMPORTANT: These QFTs do need to NOT apply swap gates.

    Input:
        - circuit (QuantumCircuit): In which to perform the addition.

        - target_register (list): Quantum Register with the list of qubits where
            the addition is performed.

        - num_sum (int): Value to be added.

        - inv (bool): If value is False, performs an addition. If True,
            performs the opposite, a substraction. By default False.

    '''

    # Get number of qubits in the register
    nqubits = len(target_register)

    # If the Addition is to be done inversey, the angles must be negative
    if inv:
        angles = -1 * get_angles_phase_addition(num_sum, nqubits)
    else:
        angles = get_angles_phase_addition(num_sum, nqubits)

    # Perform the phase changes
    for i in range(nqubits):
        circuit.p(angles[i], target_register[i])


def c_phase_add(circuit:QuantumCircuit, control:int, target_register:list, num_sum:int, inv:bool=False):
    '''
    This function performs the quantum controlled addition using Drapper addition.

    This addition is realized using direct QFT before and inverse QFT after.
    #! IMPORTANT: These QFTs do need to NOT apply swap gates.

    Input:
        - circuit (QuantumCircuit): In which to perform the addition.

        - control (int): Control qubit of the circuit.

        - target_register (list): Quantum Register with the list of qubits where
            the addition is performed.

        - num_sum (int): Value to be added.

        - inv (bool): If value is False, performs an addition. If True,
            performs the opposite, a substraction.

    '''

    # Get number of qubits in the register
    nqubits = len(target_register)

    # If the Addition is to be done inversey, the angles must be negative
    if inv:
        angles = -1 * get_angles_phase_addition(num_sum, nqubits)
    else:
        angles = get_angles_phase_addition(num_sum, nqubits)

    # Perform the phase changes
    for i in range(nqubits):
        circuit.cp(theta=angles[i], control_qubit=control, target_qubit=target_register[i])


def c_phase_add_mod_K(target_register:list, control:int, ancilla_qubit:int, num_sum:int, K:int, approx_QFT:int=0, include_QFTs:bool=False, name:str=None):
    '''
    This function perforns the quantum controlled modulo K addition by using Draper addition.

    Input:
        - target_register (list): Quantum Register with the list of qubits where
            the addition is performed.

        - control (int): Qubit to be the control.

        - ancilla_qubit (int): Qubit to be the ancilla of the algorithm.

        - num_sum (int): Value to be added.

        - K (int): Number such that operation will be performed modulo K.

        - approx_QFT (int): Approximation degree of the QFT. By default is 0 (No approximation)

        - include_QFTs (bool): Apply the initial direct QFT and the final inverse QFT if True,
        not apply them if False. By default False.

        - name (str): Name of the circuit. If None is provided, the number would be
        " + num_sum mod N ". By default None.

    Output:
        - circuit (QuantumCircuit): Circuit which performs the modulated addition using Draper algorithm.

    '''

    nqubits = len(target_register)

    # Construction of the circuit
    # The number of qubits of the circuit is nqubits + 2, one for the control and another one for the ancilla.
    if name:# If name is provided give such name to the circuit
        circuit = QuantumCircuit(nqubits+2, name=name)
    else: # Otherwise, the name is just " + num_sum mod K"
        circuit = QuantumCircuit(nqubits+2, name = ' + %d mod %d'%(num_sum, K))

    # Define the QFTs to be used along the function to avoid constructing them every time needed
    directQFT = QFT(num_qubits=nqubits, do_swaps=False, approximation_degree=approx_QFT)
    inverseQFT = QFT(num_qubits=nqubits, do_swaps=False, approximation_degree=approx_QFT, inverse=True)

    # First Operation
    # Sum the num_sum value and substract the K value (mod number)

    if include_QFTs:
        #  Direct QFT
        circuit.append(directQFT, target_register)

    # Controlled addition
    c_phase_add(circuit=circuit, control=control, target_register=target_register, num_sum=num_sum, inv=False)
    # Substraction of K
    phase_add(circuit=circuit, target_register=target_register, num_sum=K, inv=True)

    #  Inverse QFT
    circuit.append(inverseQFT, target_register)   

    # Check whether it was overfilled (the substraction was not needed
    # because the value of the sum was not exceeding K)
    circuit.cx(target_register[nqubits-1], ancilla_qubit)

    # Second operation
    # Add K again and substract num_sum

    # Direct QFT (return to phase space)
    circuit.append(directQFT, target_register) 

    # Add K again if the ancilla qubit is 1 (the substraction overfilled)
    c_phase_add(circuit=circuit, target_register=target_register, control=ancilla_qubit, num_sum=K, inv=False)

    # Substract num_sum to return to initial value but keeping the mod K
    c_phase_add(circuit=circuit, target_register=target_register, control=control, num_sum=num_sum, inv=True)


    # Inverse QFT (return to phase space)
    circuit.append(inverseQFT, target_register) 

    # Return ancilla qubit to 0 if needed (only in case overfilled happened before)
    circuit.x(target_register[nqubits-1])
    circuit.cx(target_register[nqubits-1], ancilla_qubit)
    circuit.x(target_register[nqubits-1])

    # Third operation
    # Add num_sum finally

    # Add numsum again
    # Direct QFT (return to phase space)
    circuit.append(directQFT, target_register)

    # Add num_sum to return to initial value but keeping the mod K
    c_phase_add(circuit=circuit, target_register=target_register, control=control, num_sum=num_sum, inv=False)

    if include_QFTs:
        # Inverse QFT (return to phase space)
        circuit.append(inverseQFT, target_register)

    return circuit

