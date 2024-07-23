from mpmath import iv

def my_arctan_floor(x, n):
    sum = iv.mpf(0)
    for i in range(0, (2 * n) + 2):
        x_term = x ** ((2 * i) + 1)
        numerator = iv.mpf((-1) ** i)
        denominator = iv.mpf((2 * i) + 1)
        term = x_term * (numerator / denominator)
        sum += term
    return sum

def my_arctan_ceiling(x, n):
    sum = iv.mpf(0)
    for i in range(0, (2 * n) + 1):
        x_term = x ** ((2 * i) + 1)
        numerator = iv.mpf((-1) ** i)
        denominator = iv.mpf((2 * i) + 1)
        term = x_term * (numerator / denominator)
        sum += term
    return sum

def my_pi_floor(n):
    return (iv.mpf(16) * my_arctan_floor(iv.mpf("1/5"), n)) - (iv.mpf(4) * my_arctan_ceiling(iv.mpf("1/239"), n))

def my_pi_ceiling(n):
    return (iv.mpf(16) * my_arctan_ceiling(iv.mpf("1/5"), n)) - (iv.mpf(4) * my_arctan_floor(iv.mpf("1/239"), n))

def my_atan_floor(x, n):
    first = (my_pi_floor(n)/iv.mpf(2))
    second = my_arctan_ceiling(iv.mpf(1)/x, n)
    return (first - second)

def my_atan_ceiling(x, n):
    first = (my_pi_ceiling(n)/iv.mpf(2))
    second = my_arctan_floor(iv.mpf(1)/x, n)
    return (first - second)

def my_atan(x):
    return iv.mpf([my_atan_floor(x, 10).a, my_atan_ceiling(x, 10).b])

def test_atan(x, n):
    return iv.mpf([my_atan_floor(x, n).a, my_atan_ceiling(x, n).b])
