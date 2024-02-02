'''
Date        : 2024-02-02
Author      : Gautam Singh
Roll Number : CS21BTECH11018
File        : CS5610_CS21BTECH11018_HW2_crt.py
Purpose     : Implement the Chinese Remainder Theorem to solve a simulataneous 
              system of congruences if the moduli are co-prime, else output -1.
'''

# Matrix multiplication algorithm using pure Python lists
def matrix_mul(A: list, B: list) -> list:
    # Dimension check
    if len(A[0]) != len(B):
        raise ArithmeticError(f'Incompatible dimensions for matrix multiplication, ({len(A)}, {len(A[0])}) * ({len(B)}, {len(B[0])})')
    # Perform the matrix multiplication itself
    zip_B = list(zip(*B))
    return [[sum(aij*bij for aij, bij in zip(A_rows, B_cols)) for B_cols in zip_B] for A_rows in A]
    

# Extended Euclid's algorithm
def euclid_ext(a: int, b: int) -> (int, int, int):
    # Accumulate result in M
    M = [[1, 0], [0, 1]]
    # Make a the smaller number
    is_swapped = False
    if a > b:
        a, b = b, a
        is_swapped = True
    while a > 0:
        # (a, b) -> (b%a, a)
        q = b // a
        M1 = [[-q, 1], [1, 0]]
        M = matrix_mul(M1, M)
        a, b = b%a, a
    x, y = M[1][0], M[1][1]
    # Swap x and y appropriately
    if is_swapped:
        x, y = y, x
    return x, y, b

# Constants
INPUT_FILE='testinput-crt.txt'

# Parse input from file
fh = open(INPUT_FILE, 'r')
testcases = [[int(x) for x in line.split(',')] for line in fh.readlines()]

for case in testcases:
    a, m, b, n = case
    # Apply Euclid's algorithm on m, n
    p, q, g = euclid_ext(m, n)
    if g != 1:
        print("-1")
    else:
        # We have x = a + mk1 = b + nk2
        # Thus, mk1 - nk2 = b - a
        # But mp + nq = 1.
        # Thus, k1 = p(b - a) and x = a + mp(b - a)
        x = a + m*p*(b - a)
        # Find x modulo mn
        mn = m*n
        x = (x%mn + mn)%mn
        print(x)
