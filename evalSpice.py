import numpy as np


def evalSpice(filename):
    lines = readfile(filename)  # Function to parse the input file
    matA, matB, nodes_dict, aux_dict = nodal_equations(
        lines
    )  # Function to creating the matrix
    return solve(matA, matB, nodes_dict, aux_dict)  # Solving the matrix


def solve(matA, matB, nodes_dict, aux_dict):
    try:
        sol = np.linalg.solve(matA, matB)
    except:
        raise ValueError("Circuit error: no solution")
    del nodes_dict["GND"]
    for i in nodes_dict.keys():
        nodes_dict[i] = float(sol[nodes_dict[i] - 1])
    nodes_dict["GND"] = 0
    Iout = {}
    for i in aux_dict.keys():
        Iout[i] = float(sol[aux_dict[i][2] - 1])
        if (
            abs(Iout[i]) > 1e290
        ):  # To check if current is not blowing up in any shorted wire(R = 0)
            raise ValueError("Circuit error: no solution")
    return nodes_dict, Iout


def readfile(filename):
    try:
        with open(filename, "r") as f:
            content = f.read().split("\n")
    except:  # If file is not present
        raise FileNotFoundError("Please give the name of a valid SPICE file as input")
    for line in content:
        if line == "":
            continue
        i = line.split()
        if i[0] == ".circuit":
            start = content.index(line)
        if i[0] == ".end":
            end = content.index(line)
    try:  # Returns the list which actually contains the circuit data
        if start < end:
            return content[start + 1 : end]
        else:  # .circuit and .end in the wrong order
            raise ValueError("Malformed circuit file")
    except:  # If .circuit and .end is not present
        raise ValueError("Malformed circuit file")


def nodal_equations(lines):
    nodes = []
    nodes_dict = {}
    node_count = 0
    for line in lines:
        if line == "" or line[0] == "#":
            continue
        i = line.split()
        for j in range(1, 3):
            if i[j] == "GND" and i[j] not in nodes:
                nodes.append(i[j])
                nodes_dict[i[j]] = node_count  # Setting the ground node to have index 0
                node_count += 1
    for line in lines:
        if line == "":
            continue
        i = line.split()
        for j in range(1, 3):
            if i[j] not in nodes:
                nodes.append(i[j])
                nodes_dict[i[j]] = node_count  # Each node has an associated node index
                node_count += 1
    if "GND" not in nodes_dict.keys():  # Ground node must be present
        raise ValueError("Malformed circuit file")
    aux_dict = {}  # For auxillary equations due to current through to voltage sources
    aux_count = node_count
    for line in lines:
        if line == "":
            continue
        i = line.split()
        if i[0][0] == "V":
            aux_dict[i[0]] = (i[1], i[2], aux_count)
            aux_count += 1
    names = []
    for line in lines:
        if line == "":
            continue
        names.append(line.split()[0])
    if len(set(names)) != len(names):  # If any name is repeated
        raise ValueError("Malformed circuit file")
    matA, matB = create_mat(lines, aux_count, nodes_dict, aux_dict)
    matA = matA[1:, 1:]  # To eliminate redundant equations for the ground node
    matB = matB[1:]
    return matA, matB, nodes_dict, aux_dict


def create_mat(lines, aux_count, nodes_dict, aux_dict):
    matB = np.zeros(aux_count)
    matA = np.zeros((aux_count, aux_count))
    for j in nodes_dict.keys():
        for line in lines:
            if line == "":
                continue
            i = line.split()
            for k in range(1, 3):
                if i[k] == j:
                    if i[0][0] == "R" and float(i[3]) != 0:
                        if (
                            len(i) > 4 and i[4] != "#"
                        ):  # To detect any non commented text
                            raise ValueError("Malformed circuit file")
                        if k == 1:  # Populating the conductance matrix
                            matA[nodes_dict[j]][nodes_dict[i[1]]] += 1 / float(i[3])
                            matA[nodes_dict[j]][nodes_dict[i[2]]] -= 1 / float(i[3])
                        elif k == 2:
                            matA[nodes_dict[j]][nodes_dict[i[2]]] += 1 / float(i[3])
                            matA[nodes_dict[j]][nodes_dict[i[1]]] -= 1 / float(i[3])
                    elif i[0][0] == "R" and float(i[3]) == 0:
                        if (
                            len(i) > 4 and i[4] != "#"
                        ):  # To detect any non commented text
                            raise ValueError("Malformed circuit file")
                        if (
                            k == 1
                        ):  # Zero resistance is replaced by the smallest possible floating point value
                            matA[nodes_dict[j]][nodes_dict[i[1]]] += 1 / float("1e-308")
                            matA[nodes_dict[j]][nodes_dict[i[2]]] -= 1 / float("1e-308")
                        elif k == 2:
                            matA[nodes_dict[j]][nodes_dict[i[2]]] += 1 / float("1e-308")
                            matA[nodes_dict[j]][nodes_dict[i[1]]] -= 1 / float("1e-308")
                    elif i[0][0] == "I":
                        if len(i) > 5 and i[5] != "#":
                            raise ValueError("Malformed circuit file")
                        if k == 1:  # Populating the current vector
                            matB[nodes_dict[j]] = -float(i[4])
                        elif k == 2:
                            matB[nodes_dict[j]] = float(i[4])
                    elif i[0][0] == "V":
                        if len(i) > 5 and i[5] != "#":
                            raise ValueError("Malformed circuit file")
                        if (
                            k == 1
                        ):  # Populating the extended conductance matrix for the auxillary equations
                            matA[nodes_dict[j]][aux_dict[i[0]][2]] = 1
                        elif k == 2:
                            matA[nodes_dict[j]][aux_dict[i[0]][2]] = -1
                        matA[aux_dict[i[0]][2]][nodes_dict[i[1]]] = 1
                        matA[aux_dict[i[0]][2]][nodes_dict[i[2]]] = -1
                        matB[aux_dict[i[0]][2]] = float(i[4])
                    else:  # To detect invalid circuit elements
                        raise ValueError("Only V, I, R elements are permitted")
    return matA, matB
