from qiskit import QuantumCircuit
import utilities_multiples as utmult


def oracle_less_than(number:int, nqubits:int, name=None) -> QuantumCircuit:

    '''
    This function builds a quantum circuit, an oracle, which marks with a pi-phase
    those states which represent numbers strictly smaler than the number given by parameter.

    The procedure is almost the same for all numbers, with the only exception of a difference
    if the first bit of the number in binary is 1 or 0.

    Input:
    number: integer (int) containing the objective number,
       or a string (str) with the binary representation of such number.
    nqubits: integer (int) number of qubits of the circuit.
       It must be larger than the number of digits of the binary representation of number.
    name: string (str), default None, name of the circuit.

    Output:
    circuit: QuantumCircuit which marks with fase pi the states which
    represent in binary the numbers strictly smaller than number.
    '''

    # Construction of the circuit
    if name:# If name is provided give such name to the circuit
        circuit = QuantumCircuit(nqubits, name=name)
    else: # Otherwise, the name is just " < number"
        circuit = QuantumCircuit(nqubits, name = ' < %d '%number)

    # Binary representation of the number
    num_binary = utmult.to_binary(number, nqubits)
    
    # Discard the 0s at the end, as they will not be used and save
    # unnecessary X gates
    num_binary = num_binary.rstrip('0')

    
    if num_binary[0] == '1':
        # If the first digit is 1
        # Mark all the states of the form |0q1...>
        circuit.x(nqubits-1)
        circuit.z(nqubits-1)
        circuit.x(nqubits-1)
    else:
        # If first digit is 0
        # Apply X gate to first qubit
        circuit.x(nqubits-1)
    
    # For loop on the remaining digits
    for position1, value in enumerate(num_binary[1:]):
        # Rename the position as it starts with 0 in the second bit and
        # we want it to be 1.
        position = position1 + 1

        if value == '0':
            # If the digit is 0
            # Just apply a X gate
            circuit.x(nqubits-position-1)
        else:
            # If the digit bi is 1
            # Apply a multi-controlled Z gate to mark states of the shape:
            # |bn...bi+1 0 qi-1...q1>
            # where bn,...,bi+1 are the first n-i bits of m, which is of the shape bn...bi+1 1 bi-1...b1
            # because we just checked that bi is 1.
            # Hence, the numbers of the form bn...bi+1 0 qi-1...q1 are smaller than m.
            circuit.x(nqubits-position-1)
            multi_z = utmult.multi_control_z(position + 1)
            circuit.append(multi_z.to_gate(), range(nqubits-1, nqubits-position-2, -1))
            circuit.x(nqubits-position-1)
    
    for position, value in enumerate(num_binary):
        # Apply X gates to qubits in position of bits with a 0 value
        if value == '0':
            circuit.x(nqubits-position-1)
        else:
            pass
    
    return circuit
