'''
Date        : 2024-03-08
Author      : Gautam Singh
Roll Number : CS21BTECH11018
File        : CS5610_CS21BTECH11018_ProgHW4.py
Purpose     : Program to implement Euclid's Extended Algorithm in Zp[x] to 
              compute the gcd of two given polynomials f(x) and g(x) and also
              compute polynomials u(x) and v(x) such that f(x)*u(x) + g(x)*v(x)
              = gcd(f(x), g(x)).
'''

# Imports
import itertools
import sys

sys.setrecursionlimit(1000)

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

class Polynomial:
    """
    Object-oriented implementation of a polynomial in Zp[x].
    """
    def __init__(
        self,
        p: int,
        coeff: tuple[int, ...] = tuple()
    ) -> None:
        self.p = p
        # Coefficients should lie in Zp
        self.coeff = [(c%p + p)%p for c in coeff]
    
    def degree(self):
        """ Return the degree of the polynomial. """
        return len(self.coeff) - 1

    def field_equal(self, other):
        """ Check if both polynomials belong to the same field Zp[x]. """
        assert self.p == other.p
    
    def derivative(self):
        """ Compute the derivative of the polynomial. """
        deg = self.degree()
        # Coefficients of the derived polynomial go here
        der_coeff = []
        for i, coeff in enumerate(self.coeff):
            # Constant term check
            if deg - i - 1 < 0:
                continue
            next_coeff = (deg - i - 1)*coeff%self.p
            # Leading zeros should not be appended
            if next_coeff != 0 or len(der_coeff) != 0:
                der_coeff.append(next_coeff)
        # Create a polynomial object and return
        return Polynomial(self.p, tuple(der_coeff))

    def square_free(self):
        """ Method to find the square-free part of the given polynomial. """

    def factorize(self):
        """ Factorize the polynomial using the Cantor-Zassenhaus algorithm. """

    def __add__(self, other):
        """ Add two polynomials. """
        if isinstance(other, int):
            other = Polynomial(self.p, (other,))
        self.field_equal(other)
        # Get coefficients of polynomial objects, with smallest degree first
        a = self.coeff[::-1]
        b = other.coeff[::-1]
        # c = a + b, with zeros padded for later coefficients
        c = [sum(t)%self.p for t in itertools.zip_longest(a, b, fillvalue=0)]
        c = c[::-1]
        cnt = len(c)
        for i, val in enumerate(c):
            if val != 0:
                cnt = i
                break
        c = c[cnt:]
        if len(c) == 0:
            c = [0]
        return Polynomial(self.p, tuple(c))

    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        return self + other

    def __mul__(self, other):
        """ Multiply two polynomials. """
        if isinstance(other, int):
            other = Polynomial(self.p, (other,))
        self.field_equal(other)
        res = Polynomial(self.p)
        # Perform naive multiplication in Zp[x]
        for _, bi in enumerate(other.coeff):
            res_add = self
            res_add.coeff = [c*bi%p for c in res_add.coeff]
            res += res_add
            # Shift the results for the next addition
            res.coeff.append(0)
        # The last shift is not required, if it exists
        if len(res.coeff):
            res.coeff.pop()
        return res

    def __rmul__(self, other):
        return self*other

    def __imul__(self, other):
        return self*other

    def __sub__(self, other):
        """ Subtract `other` from `self`. """
        if isinstance(other, int):
            other = Polynomial(self.p, (other,))
        self.field_equal(other)
        other1 = other*-1
        res = self + other1
        return res

    def __rsub__(self, other):
        return other + -1*self

    def __isub__(self, other):
        return self - other

    def __floordiv__(self, other):
        """ Return the quotient of `self` / `other`. """
        self.field_equal(other)
        quot = Polynomial(self.p, tuple())
        while self.degree() >= other.degree():
            c = self.coeff[0]*modexp(other.coeff[0], p-2, p)%p
            self.coeff = [sum(t)%p for t in itertools.zip_longest(self.coeff, [(p - c*oc)%p for oc in other.coeff], fillvalue=0)]
            cnt = len(self.coeff)
            for i, si in enumerate(self.coeff):
                if si != 0:
                    cnt = i
                    break
            self.coeff = self.coeff[cnt:]
            quot.coeff.append(c)
        return quot

    def __rfloordiv__(self, other):
        return other // self

    def __ifloordiv__(self, other):
        return self // other

    def __mod__(self, other):
        if isinstance(other, int):
            other = Polynomial(self.p, (other,))
        self.field_equal(other)
        while self.degree() >= other.degree():
            c = self.coeff[0]*modexp(other.coeff[0], p-2, p)
            self -= c*other
        return self

    def __imod__(self, other):
        return other % self

    def __rmod__(self, other):
        return self % other

    def __divmod__(self, other):
        """ Return both quotient and remainder on dividing `self` by `other`. """
        return self // other, self % other

    def __rdivmod__(self, other):
        return other // self, other % self

    def __idivmod__(self, other):
        return self // other, self % other

    def __eq__(self, other):
        """ Check if two polynomials are equal. """
        return self.p == other.p and self.coeff == other.coeff

    def __neq__(self, other):
        """ Check if two polynomials are not equal. """
        return not self.p == other.p

    def __str__(self, var: str = 'x') -> str:
        """ Print the polynomial, given a variable representation. """
        val = ''
        deg = self.degree()
        for i, ai in enumerate(self.coeff):
            if ai == 0:
                continue
            if deg == i:
                val += f'{ai}'
            elif deg == i + 1:
                if ai > 1:
                    val += f'{ai}*'
                val += var
            else:
                if ai > 1:
                    val += f'{ai}*'
                val += f'{var}^{deg-i}'
            val += ' + '
        return val[:-2]


def matrix_mul(A: list, B: list) -> list:
    """
    Function to implement matrix multiplication using pure Python lists.
    """
    # Dimension check
    if len(A[0]) != len(B):
        raise ArithmeticError(f'Incompatible dimensions for matrix multiplication, ({len(A)}, {len(A[0])}) * ({len(B)}, {len(B[0])})')
    # Perform the matrix multiplication itself
    zip_B = list(zip(*B))
    return [[sum(aij*bij for aij, bij in zip(A_rows, B_cols)) for B_cols in zip_B] for A_rows in A]

def euclid_ext_zp (
    p: int,
    f: Polynomial,
    g: Polynomial
) -> tuple[Polynomial, Polynomial, Polynomial]:
    """
    Function that implements Euclid's Extended Algorithm in Zp[x] to find
    the gcd of f(x) and g(x), as well as find polynomials u(x) and v(x) such
    that f(x)*u(x) + g(x)*v(x) = gcd(f(x), g(x)).
    """
    # Initialize M to be the identity matrix
    # Note that M@[a b].T = [b a%b].T
    zero_poly = Polynomial(p)
    M = [[Polynomial(p, (1,)), zero_poly], [zero_poly, Polynomial(p, (1,))]]
    swap_fl = False
    if f.degree() > g.degree():
        f, g = g, f
        swap_fl = True
    while f != zero_poly:
        q, r = divmod(g, f)
        M1 = [[-1*q, 1], [1, 0]]
        M = matrix_mul(M1, M)
        f, g = r, f
    if swap_fl:
        return g, M[1][1], M[1][0]
    else:
        return g, M[1][0], M[1][1]

# Constants
INPUT_FILE='input-polygcd2.csv'

# Parse input from file
p = 0
f = []
g = []

with open(INPUT_FILE, 'r', encoding='utf-8') as fh:
    L = fh.readlines()
    L = [l.split(',') for l in L]
    p = int(L[0][0])
    df = int(L[1][0])
    f = [int(x) for x in L[1][1:2+df]]
    dg = int(L[2][0])
    g = [int(x) for x in L[2][1:2+dg]]

# Input is of the form
# p
# deg(f) f
# deg(g) g

f = Polynomial(p, tuple(f))
g = Polynomial(p, tuple(g))

# Apply Euclid's Extended Algorithm
h, u, v = euclid_ext_zp(p, f, g)
h *= modexp(h.coeff[0], p-2, p)

# Print the outputs
print(f'GCD: {h}\nu: {u}\nv: {v}')
