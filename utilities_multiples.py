from qiskit import QuantumCircuit
from collections import namedtuple
import qiskit
import numpy as np


def get_remainders_power_2(number, power):
    '''
    Remainders of powers of 2
    Given a number and a power, returns an array with the remainders of 
    dividing powers of 2 up to the given power by the provided number
    
    Input:
        - number (int): Integer number greater than 2
        - power (int): Integer, compute powers of 2 up to that exponent

    Ouput:
        - remainders (list): List of remainders of the power of 2.

    This functions conducts no checks on the input.
    '''

    remainders = [1, 2]

    for p in range(2, power + 1):
        check_num = remainders[-1]*2
        if check_num >= number:
            check_num = check_num - number
        else:
            pass
        if check_num == 1:
            complete_repetitions = power // len(remainders) + 1
            # rondas_parcial = power % len(remainders)
            # remainders = remainders * rondas_completas + remainders[:rondas_parcial+1]
            remainders = remainders * complete_repetitions
            remainders = remainders[:power+1]
            break
        else:
            pass

        remainders.append(check_num)

    return remainders


def to_binary(number:int, nbits:int = None):
    
    '''
    This fucntion transforms an integer to its binary form (string).
    If a determined number of bits is required (more than the needed ones),
    it can be passed as a parameter too, nbits, None by default.
    It is needed that the number of bits passed as a parameter is larger
    than the number of bits needed to write the number in binary. 

    Input:
        - number: integer (int).
        - nbits: integer (int), None by default

    Output:
        - binary: string (str) containing the number in its binary form.
        It writes 0s in front if nbits is larger than the number of bits needed
        to write the binary form.
    '''

    if nbits is None:
        return bin(number)[2:]
    else:
        binary = bin(number)[2:]
        if nbits < len(binary):
            print('Error, a larger nbits number is required.')
        else:
            return '0' * (nbits - len(binary)) + binary


def grover_diffuser(nqubits:int):
    
    '''
    Diffuser of grover

    Input: 
    nqubits (int): Number of qubits.

    Output:
    circuit: QuantumCircuit containing the diffuser.
    '''

    circuit=QuantumCircuit(nqubits, name='Diffuser')
    circuit.h(range(nqubits))
    circuit.x(range(nqubits))
    control_z=multi_control_z(nqubits)
    circuit.append(control_z.to_gate(), range(nqubits))
    circuit.x(range(nqubits))
    circuit.h(range(nqubits))
    return circuit


"""
n-qubit controlled gate
"""

def multi_control_z(nqubits):
    circuit = QuantumCircuit(nqubits,name=' CZ (%d)' %(nqubits))
    xgate = np.array([[1,0],[0,-1]])
    controls = list(range(nqubits-1))
    target = nqubits-1
    mc_gate(gate=xgate, circuit=circuit, controls=controls, targ=target)
    return circuit


def mc_gate(gate: np.ndarray, circuit: qiskit.QuantumCircuit, controls: list, targ: int):
    """
    Parameters
    ----------
    gate: 2 X 2 unitary gate
    circuit: qiskit.QuantumCircuit
    controls: list of control qubits
    targ: target qubit
    Returns
    -------
    """

    n_qubits = len(controls) + 1
    gate_circuit = qiskit.QuantumCircuit(n_qubits, name="T" + str(targ))
    gate_circuit.permutation = list(range(len(controls) + 1))
    _c1c2(gate, n_qubits, gate_circuit)
    _c1c2(gate, n_qubits, gate_circuit, step=-1)

    _c1c2(gate, n_qubits - 1, gate_circuit, False)
    _c1c2(gate, n_qubits - 1, gate_circuit, False, step=-1)

    circuit.compose(gate_circuit, controls + [targ], inplace=True)
    circuit.permutation = gate_circuit.permutation


def _c1c2(gate, n_qubits, qcirc, first=True, step=1):
    pairs = namedtuple("pairs", ["control", "target"])

    if step == 1:
        start = 0
        reverse = True
    else:
        start = 1
        reverse = False

    qubit_pairs = [
        pairs(control, target)
        for target in range(n_qubits)
        for control in range(start, target)
    ]

    qubit_pairs.sort(key=lambda e: e.control + e.target, reverse=reverse)

    for pair in qubit_pairs:
        exponent = pair.target - pair.control
        if pair.control == 0:
            exponent = exponent - 1
        param = 2**exponent
        signal = -1 if (pair.control == 0 and not first) else 1
        signal = step * signal
        if pair.target == n_qubits - 1 and first:
            csqgate = _gate_u(gate, param, signal)
            qcirc.compose(csqgate, qubits=[pair.control, pair.target], inplace=True)
        else:
            qcirc.crx(signal * np.pi / param, pair.control, pair.target)


def _gate_u(agate, coef, signal):
    param = 1 / np.abs(coef)

    values, vectors = np.linalg.eig(agate)
    gate = np.power(values[0] + 0j, param) * vectors[:, [0]] @ vectors[:, [0]].conj().T
    gate = (
        gate
        + np.power(values[1] + 0j, param) * vectors[:, [1]] @ vectors[:, [1]].conj().T
    )

    if signal < 0:
        gate = np.linalg.inv(gate)

    sqgate = qiskit.QuantumCircuit(1, name="U^1/" + str(coef))
    sqgate.unitary(gate, 0)  # pylint: disable=maybe-no-member
    csqgate = sqgate.control(1)

    return csqgate

