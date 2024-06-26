# Running the code:
Run the evalSpice.py file which returns a dictionary containing the node voltages by giving it a .ckt file as input. The <mark> teest_evalSpice.py </mark> runs the source file on the test circuits in the <mark> testdata </mark> folder. The testdata folder is also expected to contain .exp files with the correct node voltages.

# Format of the .ckt and .exp files:
Consider the following .ckt file:  
.circuit  
V1   1 GND  dc 2  
R1   1   2     1  
R2   2 GND     1  
.end  
This represents a 2V dc source V1 connected between nodes 1 and GND, and resistors R1 and R2 (each of 1 ohm) connected between nodes 1 and 2, and 2 and GND respectively. It will have the corresponding .exp file as:  
({'1': 2.0, '2': 1.0, 'GND': 0.0}, {'V1': -1.0})  
which is a tuple of dictionaries, where the first dictionary gives the voltage at each node, and the second dictionary gives the current through each voltage source.
Note that current sources can also be included in the .ckt files.
# Working of the Code:

## Parsing the file  
The input file has been split into a list of lines using newline as a delimiter.
The actual circuit data(starting from .circuit upto .end) has been returned. 

## Creating the node dictionary and the auxillary dictionary  
Two dictionaries are used: one containing all the nodes(including GND) and one auxillary dictionary containing voltage sources. The nodes dictionary assigns an index to every node, with GND being index 0. The other indices are based on whichever node appears first while reading the file. The auxillary dictionary assigns each voltage source(uniquely identified by its name) a tuple which contains an auxillary index to be assigned to that source. The auxillary index is shifted by the number of nodes, because doing so will make it easier to construct the matrix. The tuple also contains the two terminal nodes of the voltage source. The required matrices are then constructed by calling the create_mat function inside the nodal_equations function. These matrices have the redundant equations for the ground node included within them, so they are trimmed accordingly to get the required equations. This has been done because it is easier to include all nodes while actually constructing the matrix.

## Creating the extended conductance matrix  
The extended conductance matrix has a size of (number of nodes + number of voltage sources). This is done so as to take into account the current through the voltage sources. The current source(s) between two nodes have been accounted for by adding (and subtracting, depending on the polarity) the current value in the indices corresponding to the terminal nodes in the extended current matrix. The bottom (number of voltage sources) entries of the extended current matrix contains the value of the corresponding voltage sources, in the order of their auxillary index, while the corresponding entries in the solution vector contain the current through the voltage sources. While going through the lines, if a voltage source is encountered, then the value of the voltage is added to its corresponding auxillary index in the extended current matrix, while the current through that source will be found in the soultion vector at that same auxillary index upon solving the equations. In addition, 1s and -1s have been added to the extended conductance matrix at the appropriate locations so as to get the equations for the voltage difference across the voltage source in terms of the node voltages. The resistors have been accounted for by adding (and subtracting, depending on which node the KCL is being written on) the conductances. While iterating through the lines, while considering the i^th and j^th indexed node(where the i^th and j^th node are connected), the conductance between them, if it exists, has been added to the (i, i) and (j, j) locations and subtracted from the (i, j) and (j, i) locations from the conductance matrix, in accordance with Kirchoff's Law.

## Solving the Matrix
The system of linear equations is solved using numpy.linalg.solve, and the indices associated with each of the nodes are replaced by the entry corresponding to the same index in the output vector. The remaining entries of the output vector are the currents through the voltage sources, which have been returned in a seperate dictionary. The current to voltage source correspondance is obtained through the auxillary dictionary. 

## Exception Handling and Edge Cases
The following exceptions and edge cases have been covered. The lines in the code where they have been covered are mentioned in the parenthesis.  
 If the numpy.linalg.solve return an error, it is because the circuit is invalid, and it has no solution. This covers any situations like two current sources in series or two voltage sources in parallel. (7)
 The input file is malformed if it does not have .circuit and .end, or has them in the wrong order. (35)
 The input file is also malformed if it does not have a ground node. (64)
 Any text which is not commented results in a malformed file. (97, 106, 115, 122)
 A file is malformed if the name of any element is repeated. This is checked by checking for any duplicates in the list of names. (80)
 Any elements other than V, I, R are not permitted. (132)
 Short circuit: If a resistor of 0 ohms is connected across two nodes, the two nodes will become shorted. If a zero resistance is given, it has been replaced by the smallest floating point number in python3, which is 1e-308(105). A short may or may not result in an invalid circuit, so if the current through any two nodes blows up (>1e+290), an invalid circuit error is raised(18). This has its limitations, and may fail if the value of the voltage source is very low, or the value of other resistors in the circuit is very low. Two extra test cases- one where the short is invalid, and one where it is valid- have been included. The first is test_short_v.ckt, which has two resistors and a battery in series connected to GND on oth sides, with one resistor being 0. This is a valid circuit, and the .exp file has also been included. The other one is test_short.ckt, which has a battery and zero ohm resistor in series, with both ends being connected to ground. This is an invalid circuit, and the corresponding error has been raised.
