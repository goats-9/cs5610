'''
Date        : 2024-03-03
Author      : Gautam Singh
Roll Number : CS21BTECH11018
File        : CS5610_CS21BTECH11018_ProgHW3.py
Purpose     : Program to implement the Tonelli-Shanks Algorithm to find the 
              smallest integer x such that x^2 = a in Zp, where p is a prime and
              0 < a < p.
'''

# Imports
import random

def modexp(
    a: int,
    b: int,
    n: int
) -> int:
    """
    Compute the value of (a^b) mod n using repeated squaring.
    """
    ret = 1
    a_pow_2 = a
    while b:
        if b&1:
            ret = (ret*a_pow_2)%n
        b >>= 1
        a_pow_2 = (a_pow_2*a_pow_2)%n
    return ret

def find_k(
    z: int,
    p: int
) -> int:
    """
    Returns the smallest k such that z^(2^k) mod p = 1.
    """
    k = 0
    # Maintain value of z^(2^k) here
    y = z
    while y != 1:
        k += 1
        y = (y*y)%p
    return k

def tonelli_shanks(
    a: int,
    p: int
) -> int:
    """
    Function that implements the Tonelli-Shanks algorithm, given 0 < a < p and p
    is a prime.
    """
    # Factorize p - 1 = 2^t * m
    m = p - 1
    t = 0
    while m % 2 == 0:
        t += 1
        m = m >> 1
    b = modexp(a, m, p)
    k = find_k(b, p)
    if k == t:
        return 0
    x = modexp(a, (m+1)>>1, p)
    r = 0
    while modexp(r, (p-1)>>1, p) != p-1:
        r = random.randint(1, p-1)
    s = modexp(r, m, p)
    S = modexp(s, 1<<(t-k), p)
    while k > 0:
        b = (b*S)%p
        x = (x*modexp(s, 1<<(t-k-1), p))%p
        k = find_k(b, p)
        S = modexp(s, 1<<(t-k), p)
    return x

# Constants
INPUT_FILE='inputSquareRoots.csv'

# Parse input from file
testcases = []
with open(INPUT_FILE, 'r', encoding='utf-8') as fh:
    testcases = [[int(x) for x in line.split(',')] for line in fh.readlines()]

# Solve for each testcase
for case in testcases:
    a, p = case
    # Call the Tonelli-Shanks algorithm
    ans = tonelli_shanks(a, p)
    # There are two solutions, if they exist: ans and p - ans.
    # Note that even if ans = 0 for no solution, 0 will be output.
    print(min(ans, p - ans))
