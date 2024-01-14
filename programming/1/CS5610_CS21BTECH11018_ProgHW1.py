'''
Date        : 2024-01-14
Author      : Gautam Singh
Roll Number : CS21BTECH11018
File        : CS5610_CS21BTECH11018_HW1.py
Purpose     : Implement Euclid's algorithm to find gcd of two numbers a and b 
              and find integers x, y such that ax + by = gcd(a,b).
'''

# Imports
import csv

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
INPUT_FILE='testinput-gcd.csv'

# Read input using csv module
with open(INPUT_FILE, 'r') as fh:
    reader = csv.reader(fh)
    for row in reader:
        a, b = int(row[0]), int(row[1])
        # Apply Euclid's algorithm
        x, y, c = euclid_ext(a, b)
        print(f'x={x},y={y},c={c}')