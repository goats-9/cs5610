'''
Date        : 2024-04-01
Author      : Gautam Singh
Roll Number : CS21BTECH11018
File        : CS5610_CS21BTECH11018_ProgHW5.py
Purpose     : Program to implement the Cantor-Zassenhaus algorithm to factor 
              univariate polynomials in Zp[x].
'''

# Imports
import sys
import math
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

def modinv(
    a: int,
    n: int,
) -> int:
    """ Function to compute the modular inverse of a in Zn. Wrapper around the
    `modexp` function. """
    return modexp(a,n-2,n)

class Polynomial:
    """ Object-oriented implementation of a polynomial in Zp[x]. """
    def __init__(
        self,
        p: int,
        coeff: dict[int, int]
    ) -> None:
        # Check if p is prime
        self.p = p
        # Coefficients should lie in Zp
        self.coeff = {}
        for _key, _val in coeff.items():
            _val_p = (_val%self.p + self.p)%self.p
            if _val_p != 0:
                self.coeff[_key] = _val_p

    def degree(self) -> int:
        """ Return the degree of the polynomial. """
        if len(self.coeff) == 0:
            return -1
        return max(self.coeff.keys())

    def field_equal(self, other: 'Polynomial'):
        """ Check if both polynomials belong to the same field Zp. """
        assert self.p == other.p

    def pow_reduce(self, f, a: int):
        """ Return `(f**a) % self` using binary exponentiation. """
        _x_pow = f
        _res = Polynomial(self.p, {0:1})
        _a = a
        while _a:
            if _a&1:
                _res = (_res*_x_pow)%self
            _a >>= 1
            _x_pow = (_x_pow*_x_pow)%self
        return _res

    def gcd(self, g):
        """ Method to implement Euclid's Extended Algorithm in Zp[x] to find
        the gcd of f(x) and g(x). """
        _f = Polynomial(self.p, self.coeff)
        assert isinstance(g, Polynomial)
        self.field_equal(g)
        _g = Polynomial(g.p, g.coeff)
        if _f.degree() > _g.degree():
            _f, _g = _g, _f
        _zero_poly = Polynomial(self.p, {})
        while _f != _zero_poly:
            _f, _g = _g % _f, _f
        # Make the gcd monic
        return _g*modinv(_g.coeff[_g.degree()], p)

    def derivative(self):
        """ Compute the derivative of the polynomial. """
        # Coefficients of the derived polynomial go here
        _der_coeff = {}
        for _key, _val in self.coeff.items():
            # Constant term check
            if _key == 0:
                continue
            _der_coeff[_key - 1] = _key*_val
        # Create a polynomial object and return
        return Polynomial(self.p, _der_coeff)

    def pth_root(self):
        """ Method to find the pth root of `self`. """    
        _coeff_list = {}
        for _key, _val in enumerate(self.coeff.keys()):
            if _key%self.p:
                raise ArithmeticError(f'Cannot find the {self.p}th root of {self}.')
            _coeff_list[_key/self.p] = _val
        return Polynomial(self.p, _coeff_list)

    def square_free(self):
        """ Method to find the square-free part of the given polynomial. """
        _g = self
        _f = self
        _h = Polynomial(self.p, {})
        _zero_poly = Polynomial(self.p, {})
        while _g != _zero_poly:
            _h, _g = _g, _g.derivative()
        if _h.degree() <= 0:
            _drv = self.derivative()
            _drv_gcd = self.gcd(_drv)
            return self // _drv_gcd
        _f = _f // _h
        _h = _h.pth_root().square_free()
        return (_h*_f) // _f.gcd(_f.derivative())

    def distinct_degree_factors(self) -> dict[int, 'Polynomial']:
        """ Method to return a list of products of irreducible factors of
        `self`, split by degree. """
        # Maintain a list of product of irreducible factors of each degree
        _poly_dict = {}
        _i = 1
        _g = self
        while _g.degree() > 0:
            # Create the coefficients for x^(p^i) - x
            _xpi = self.pow_reduce(Polynomial(self.p, {1:1}), self.p**_i)
            _xpi -= Polynomial(self.p, {1:1})
            _f = _g.gcd(_xpi)
            if _f.degree() > 0:
                _poly_dict[_i] = _f
            _g = _g // _f
            _i = _i + 1
        return _poly_dict

    def irreducible_factors(self, deg: int) -> list:
        """ Method to return the set of irreducible factors of `self` with
        strictly positive degree `deg`. `self` must be square-free and have
        irreducible factors of degree `deg` only. """
        assert deg > 0, 'Degree of factors must be strictly positive.'
        if self.degree() == deg:
            return [self]
        _h1 = Polynomial(self.p, {0:1})
        _h2 = Polynomial(self.p, {0:1})
        while _h1 == 1 or _h2 == 1:
            # Pick random g of degree < self.degree()
            _random_coeffs = random.choices(population=range(self.p), k=self.degree())
            _g = Polynomial(self.p, dict(zip(list(range(self.degree())), _random_coeffs)))
            _g_pow = self.pow_reduce(_g, (self.p**deg - 1)//2) - 1
            _h1 = self.gcd(_g_pow)
            _h2 = self // _h1
        _l1 = _h1.irreducible_factors(deg)
        _l2 = _h2.irreducible_factors(deg)
        return _l1 + _l2

    def factorize(self):
        """ Factorize the polynomial using the Cantor-Zassenhaus algorithm. """
        # Phase 1: Finding the square-free part
        _sqfree_self = self.square_free()
        # Phase 2: Splitting the square-free part into products of irreducible
        # factors of the same degree.
        _irred_dict = _sqfree_self.distinct_degree_factors()
        # Phase 3: Find irreducible factors of the same degree
        _factor_list = []
        for _deg, _prod in _irred_dict.items():
            _factor_list_deg = _prod.irreducible_factors(_deg)
            _factor_list.extend(_factor_list_deg)
        return _factor_list

    def __add__(self: 'Polynomial', other: 'int | Polynomial') -> 'Polynomial':
        """ Add two polynomials. """
        _other = other
        if isinstance(other, int):
            _other = Polynomial(self.p, {0 : other})
        assert isinstance(_other, Polynomial)
        self.field_equal(_other)
        res = {}
        for key, val in _other.coeff.items():
            if key in res:
                res[key] += val
            else:
                res[key] = val
        for key, val in self.coeff.items():
            if key in res:
                res[key] += val
            else:
                res[key] = val
        return Polynomial(self.p, res)

    def __radd__(self, other):
        return self + other

    def __iadd__(self, other):
        return self + other

    def __mul__(self, other: 'int | Polynomial') -> 'Polynomial':
        """ Multiply two polynomials. """
        if isinstance(other, int):
            _other = Polynomial(self.p, {0: other})
        else:
            _other = Polynomial(other.p, other.coeff)
        assert isinstance(_other, Polynomial)
        self.field_equal(_other)
        _res = {}
        # Perform naive multiplication in Zp[x]
        for _key, _val in _other.coeff.items():
            for _key1, _val1 in self.coeff.items():
                if _key + _key1 in _res:
                    _res[_key + _key1] += _val * _val1
                else:
                    _res[_key + _key1] = _val * _val1
        return Polynomial(self.p, _res)

    def __rmul__(self, other):
        return self*other

    def __imul__(self, other):
        return self*other

    def __sub__(self, other: 'int | Polynomial') -> 'Polynomial':
        """ Subtract `other` from `self`. """
        _other = other
        if isinstance(other, int):
            _other = Polynomial(self.p, {0:other})
        assert isinstance(_other, Polynomial)
        self.field_equal(_other)
        other1 = _other*-1
        res = self + other1
        return res

    def __rsub__(self, other):
        return self - other

    def __isub__(self, other):
        return self - other

    def __floordiv__(self, other : 'int | Polynomial') -> 'Polynomial':
        """ Return the quotient of `self` / `other`. """
        _other = other
        if isinstance(other, int):
            _other = Polynomial(self.p, {0:other})
        assert isinstance(_other, Polynomial)
        self.field_equal(_other)
        # Set up the divisor and quotient
        _div = Polynomial(self.p, self.coeff)
        _quot = Polynomial(self.p, {})
        # Perform naive division
        while _div.degree() >= _other.degree():
            # Find the constant that needs to be multiplied
            _c = (_div.coeff[_div.degree()]*modexp(_other.coeff[_other.degree()], \
                                                   self.p - 2, self.p))%self.p
            # Add to quotient
            _quot_term = Polynomial(self.p, {_div.degree() - _other.degree() : _c})
            _quot += _quot_term
            _div -= _other*_quot_term
        return _quot

    def __ifloordiv__(self, other):
        return self // other

    def __mod__(self, other):
        """ Return the remainder on dividing `self` by `other`. """
        _other = other
        if isinstance(other, int):
            _other = Polynomial(self.p, {0:other})
        assert isinstance(_other, Polynomial)
        self.field_equal(_other)
        _quot = self // _other
        return self - (_other * _quot)

    def __imod__(self, other):
        return self % other

    def __divmod__(self, other):
        """ Return both quotient and remainder on dividing `self` by `other`. """
        return self // other, self % other

    def __idivmod__(self, other):
        return self // other, self % other

    def __pow__(self, other):
        """ Method to perform binary exponentiation of `self`, with nonnegative
        integer exponent `other`. """
        assert isinstance(other, int) and other >= 0, 'Exponent must be a nonnegative integer.'
        ans = 1
        self_pow = self
        while other:
            if other&1:
                ans = ans * self_pow
            other >>= 1
            self_pow = self_pow * self_pow
        return ans

    def __eq__(self, other):
        """ Check if two polynomials are equal. """
        _other = other
        if isinstance(other, int):
            _other = Polynomial(self.p, {0:other})
        return self.p == _other.p and self.coeff == _other.coeff

    def __neq__(self, other):
        """ Check if two polynomials are not equal. """
        return not self.p == other.p

    def __str__(self, var: str = 'x') -> str:
        """ Print the polynomial, given a variable representation. """
        # List of terms
        _terms_list = []
        for _key, _val in sorted(self.coeff.items(), reverse=True):
            if _key == 0:
                _terms_list.append(f'{_val}')
            elif _key == 1:
                if _val == 1:
                    _terms_list.append(f'{var}')
                else:
                    _terms_list.append(f'{_val}*{var}')
            else:
                if _val == 1:
                    _terms_list.append(f'{var}^{_key}')
                else:
                    _terms_list.append(f'{_val}*{var}^{_key}')
        if len(_terms_list) == 0:
            return '0'
        return ' + '.join(_terms_list)

    def __repr__(self) -> str:
        return self.__str__('x')

# Constants
INPUT_FILE='input-CZ.csv'

# Parse input from file
p = 0
L = []

with open(INPUT_FILE, 'r', encoding='utf-8') as fh:
    L = fh.readlines()
    L = [l.split(',') for l in L]
    p = int(L[0][0])
    L = L[1:]

for testcase in L:
    coeffs = [int(x) for x in testcase[1:]]
    coeffs = coeffs[::-1]
    f = Polynomial(p, dict(zip(list(range(len(coeffs))), coeffs)))
    f_factors = f.factorize()
    all_f_factors = []
    zero_poly = Polynomial(p, {})
    for poly in f_factors:
        cnt = 0
        while f.degree() > 0 and f%poly == zero_poly:
            f = f // poly
            cnt += 1
        if cnt == 1:
            all_f_factors.append(f'({poly})')
        else:
            all_f_factors.append(f'({poly})^{cnt}')
    print(' * '.join(all_f_factors))
