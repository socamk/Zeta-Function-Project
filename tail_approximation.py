from mpmath import iv

def l(u):
    u = iv.mpf(u)
    return (iv.mpf("0.112") * iv.log(u)) + (iv.mpf("0.278") * iv.log(iv.log(u))) + iv.mpf("2.510")

def l1(u):
    u = iv.mpf(u)
    return (iv.mpf("0.059") * iv.log(u)) + iv.mpf("2.067")


def e1(y, tau):
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    c = y/iv.mpf("2")
    val = iv.mpf("1")/(c ** iv.mpf("2"))
    val = val + iv.mpf("1")/((y + tau) ** iv.mpf("2"))
    val = val * iv.mpf("0.006") * iv.mpf("4") * (iv.pi ** iv.mpf("2"))
    return val

def e2(x, y):
    x = iv.mpf(x)
    y = iv.mpf(y)
    c = y/iv.mpf("2")
    val = (iv.mpf("1") - (iv.mpf("2") * x))/(iv.mpf("2") * y)
    val = val * iv.log((iv.mpf("2") * y)/c)
    return val

def e3(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    c = y/iv.mpf("2")
    tau = iv.mpf(tau)
    val = ((iv.mpf("1") - x) ** iv.mpf("2"))/(iv.mpf("3") * (tau ** iv.mpf("3")))
    val = val + iv.mpf("1")/(y - c)
    val = val * (iv.mpf("1") - (iv.mpf("2") * x))
    return val

def e4(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    c = y/iv.mpf("2")
    tau = iv.mpf(tau)
    val = (iv.mpf("2") - (iv.mpf("4") * x))/(tau ** iv.mpf("2"))
    val = val + ((iv.mpf("2") - (iv.mpf("4") * x))/((y -c) ** iv.mpf("2")))
    input = iv.mpf("2") * y
    val = val * l(input)
    return val

def e5(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    val = iv.mpf("4") - (iv.mpf("8") * x)
    val = val/(tau ** iv.mpf("3"))
    input = iv.mpf("2") * y
    val = val * l1(input)
    return val

def e6(x, y):
    x = iv.mpf(x)
    y = iv.mpf(y)
    c = y/iv.mpf("2")
    val = ((iv.mpf("2") * y) - c)**iv.mpf("2")
    val = val * iv.log((iv.mpf("2") * y) - c)
    val = val - (((y - c) ** iv.mpf("2")) * iv.log(y - c))
    val = val/(iv.pi * ((y - c) ** iv.mpf("2")))
    val = val * ((iv.mpf("1") - (iv.mpf("2") * x))/(iv.mpf("2") * y))
    return val

def b(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    val = ((iv.mpf("1") - (iv.mpf("2") * x))/tau)
    val = val * e1(y, tau)
    val = val + e2(x, y)
    val = val + (e3(x, y, tau) * iv.log(y/(iv.mpf("2") * iv.pi)))
    val = val * (iv.mpf("1")/(iv.mpf("2") * iv.pi))
    val = val + ((e4(x, y, tau) + e5(x, y, tau))/iv.mpf("2"))
    return val

def B(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    val = (iv.mpf("1") - (iv.mpf("2") * x))/(tau)
    val = val * e1(y, tau)
    val = val + e2(x, y)
    val = val * (iv.mpf("1")/(iv.mpf("2") * iv.pi))
    val = val + ((e4(x, y, tau) + e5(x, y, tau))/iv.mpf("2"))
    val = val + e6(x, y)
    return val

def r(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    val = (iv.mpf("1") - (iv.mpf("2") * x))/(iv.mpf("2") * iv.pi * tau)
    val = val * iv.log(y/(iv.mpf("2") * iv.pi))
    val = val - b(x, y, tau)
    return val

def R(x, y, tau):
    x = iv.mpf(x)
    y = iv.mpf(y)
    tau = iv.mpf(tau)
    val = (iv.mpf("1") - (iv.mpf("2") * x))/(iv.mpf("2") * iv.pi * tau)
    val = val * iv.log(y/(iv.mpf("2") * iv.pi))
    val = val + B(x, y, tau)
    return val