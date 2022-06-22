#!/usr/bin/env python3
from ctypes import sizeof
from re import X
import string
import matplotlib.pyplot as plt


"""
testing Curve:
# Field characteristic:
p = 5
# Curve coefficient:
a = 2
# Curve coefficient:
b = (-1)
# generatorPoint:
gp = (0, 3)
# skalar factor:
s = 5
"""


def inverse_mod(k, p):
    # simple approach => return pow(k, -1, p)
    """Returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


def add(point1, point2, p, a):
    """Returns the result of point1 + point2 with field characteristic p and curve coefficient a
    / points at infinity representation is None """

    if point1 is None:
        # O + point2 = point2
        return point2
    if point2 is None:
        # point1 + O = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + a) * inverse_mod(2 * y1, p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod((x1 - x2), p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % p,
              -y3 % p)

    return result


def double(point, p, a):
    return add(point, point, p, a)


def substract(P, Q, p, a):
    """P, Q:=(x,y), -Q = (x, -y mod p)
    R = P-Q = P+(-Q) = P + (x, -y mod p)"""
    (qx, qy) = Q
    return add(P, (qx, (-qy % p)), p, a)


#####################################################################################
# Naive-Method
def mult(k, P, p, a):
    counter = 0
    res = None
    while k > 0:
        res = add(res, P, p, a)
        counter += 1
        k = k-1
    return res, counter

#####################################################################################
# Binary-Methods (Double and Add)


def scalar_mult(s, point, p, a):
    """Returns s * point using the double-and-add algorithm over finite fields."""
    counter = 0
    result = None
    addend = point

    while s:  # bit-sequence bist MSB to LSB
        if s & 1:  # ‘&’ is a bitwise operator in Python => last bit 0 or 1 ?
            result = add(result, addend, p, a)
            counter += 1
        addend = double(addend, p, a)
        counter += 1
        s >>= 1  # shift right
    return result, counter


def scalar_mult_2(s, P, p, a):
    """Iterative algorithm, index increasing"""
    counter = 0
    bits = bin(s)[:1:-1]  # inversed bit-sequence bist LSB to MSB
    res = None
    temp = P
    for bit in bits:
        if bit == '1':
            res = add(res, temp, p, a)  # point add
            counter += 1
        temp = double(temp, p, a)
        counter += 1
    return res, counter

#################################################################################################
# m-ary Method


def int2base(x, base):
    """convert dezimal x to base repesentation"""
    digs = string.digits + string.ascii_letters
    if x < 0:
        sign = -1
    elif x == 0:
        return digs[0]
    else:
        sign = 1
    x *= sign
    digits = []
    while x:
        digits.append(digs[x % base])
        x = x // base
    if sign < 0:
        digits.append('-')
    digits.reverse()
    return ''.join(digits)


def scalar_mult_radix4(k, P, p, a):
    """m-ary Method with while main loop"""
    counter = 0
    r = 2
    m = (2 ** r)  # => base == radix
    # print("init radix4:", "k=", k, "P=", P, "r=", r, "m=", m)

    # Precomputation:
    table = [None, P]
    for i in range(2, m):
        table.append(add(table[i-1], P, p, a))
        counter += 1

    res = None

    # Main loop:
    radix_k = str(int2base(k, m))[::-1]  # LSB to MSB
    j = len(radix_k)-1
    while j >= 0:
        res = scalar_mult(m, res, p, a)[0]
        temp = table[int(radix_k[j])]
        res = add(res, temp, p, a)
        counter = (counter + m+1)
        j = j-1
    return res, counter


def scalar_mult_radix4_2(k, P, p, a):
    """m-ary Method with for in """
    counter = 0
    r = 2
    m = (2 ** r)  # => base == radix

    # Precomputation:
    table = [None, P]
    for i in range(2, m):
        table.append(add(table[i-1], P, p, a))
        counter += 1

    res = None

    # Main loop:
    radix_k = str(int2base(k, m))  # [::-1]  # LSB to MSB
    for kj in radix_k:  # for rd in radix_k:
        res = scalar_mult(m, res, p, a)[0]
        res = add(res, table[int(kj)], p, a)
        counter += 1
    return res, counter


#################################################################################################
#################################################################################################
# Methods with NAF representation
def naf(x):
    """convert int to NAF repesentation recusive"""
    if x == 0:
        return []
    z = 0 if x % 2 == 0 else 2 - (x % 4)
    return naf((x-z) // 2) + [z]


def scalar_mult_bin_naf_book(k, P, p, a):
    """Doubel and Add with NAF representation of k"""
    counter = 0
    klistNAF = naf(k)[::-1]
    # print(klistNAF)
    # print((klistNAF[::-1]))

    res = None
    # temp = P
    i = (len(klistNAF)-1)
    # for naf in kNAFs:
    while i >= 0:
        res = scalar_mult(2, res, p, a)[0]
        counter += 2
        if klistNAF[i] == 1:
            res = add(res, P, p, a)  # point add
            counter += 1
        if klistNAF[i] == -1:
            res = substract(res, P, p, a)
            counter += 1
        # temp = double(temp, p, a)
        i = i-1
    return res, counter

    #################################################################################################
# window Methods


def scalar_mult_sliding_window(k, P, p, a):
    """Sliding Windows Method with while main loop"""
    r = 2
    m = (2 ** r)  # => base == radix
    # print("init radix4:", "k=", k, "P=", P, "r=", r, "m=", m)

    # Precomputation:
    table = [None, P, double(P, p, a)]
    for i in range(1, 2**(r-1)):
        table.append(add(table[2*i-1], table[2], p, a))

    res = None

    # Main loop:
    radix_k = str(int2base(k, m))[::-1]  # LSB to MSB
    j = len(radix_k)-1
    while j >= 0:
        if radix_k[j] == 0:
            res = scalar_mult(2, res, p, a)[0]
            j = j-1
        else:
            hj = res = add(res, table[int(radix_k[j])], p, a)
        j = j-1
    return res


#################################################################################################
#################################################################################################
# Testung and assertion

def testLoop():
    for i in range(1, 70):
        print("Run", i, ":",
              mult(i, (0, 3), 5, 2)[0],
              scalar_mult(i, (0, 3), 5, 2)[0],
              scalar_mult_2(i, (0, 3), 5, 2)[0],
              scalar_mult_radix4(i, (0, 3), 5, 2)[0],
              scalar_mult_radix4_2(i, (0, 3), 5, 2)[0],
              scalar_mult_bin_naf_book(i, (0, 3), 5, 2)[0]
              )


testLoop()


def compare():
    klist = list(range(0, (700+1)))
    data1 = []
    data2 = []
    data3 = []
    for i in klist:
        r1, c1 = mult(i, (0, 3), 5, 2)
        data1.append(c1)
        r2, c2 = scalar_mult(i, (0, 3), 5, 2)
        data2.append(c2)
        r3, c3 = scalar_mult_bin_naf_book(i, (0, 3), 5, 2)
        data3.append(c3)

    plt.plot(klist, data1, label="naive mult")
    plt.plot(klist, data2, label="binary")
    plt.plot(klist, data3, label="m-ary")
    plt.title("count of add scaling k")
    plt.xlabel('k')
    plt.ylabel('number of operations')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)

    plt.savefig('plot.png')
    plt.show()


compare()


def logictest():
    tres = [None, (0, 3), (4, 4), (2, 4), (2, 1), (4, 1), (0, 2)]
    for i in range(0, 7):
        #print("Run:", i, "expected res:", tres[i])
        assert mult(i, (0, 3), 5, 2)[0].__eq__(tres[i]), "wrong result"
        assert scalar_mult(i, (0, 3), 5, 2)[0].__eq__(tres[i]), "wrong result"
        assert scalar_mult_2(i, (0, 3), 5, 2)[
            0].__eq__(tres[i]), "wrong result"
        assert scalar_mult_radix4(i, (0, 3), 5, 2)[0].__eq__(
            tres[i]), "wrong result"
        assert scalar_mult_radix4_2(
            i, (0, 3), 5, 2)[0].__eq__(tres[i]), "wrong result"


logictest()

# Operator-Desciption:
# <<	  Zero fill left shift	Shift left by pushing zeros in from the right and let the leftmost bits fall off
# >>	  Signed right shift	Shift right by pushing copies of the leftmost bit in from the left, and let the rightmost bits fall off
# [::-1]  Reverse String or bin
