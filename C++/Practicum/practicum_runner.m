%practicum_runner.m
matrix_neo= ([9,13,12,-23;-8,7,-6,5])
bend_threshold=6
cypher = magic_matrix(matrix_neo,bend_threshold)

access_granted = check_prime(cypher)
