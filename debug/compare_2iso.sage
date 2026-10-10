eA = 3
eB = 2
p = 2**eA * 3**eB - 1
Fp = GF(p)
R.<x> = PolynomialRing(Fp)
Fp2.<i> = GF(p**2, modulus=x**2 + 1)

def xDBL(A, X, Z):
    A24 = (A + 2) / 4
    add = X + Z
    sub = X - Z
    sq_add = add * add
    sq_sub = sub * sub
    four_XZ = sq_add - sq_sub
    next_X = sq_add * sq_sub
    next_Z = four_XZ * (sq_sub + A24 * four_XZ)
    return next_X, next_Z

def get_next_A_new(XK, ZK):
    XK2 = XK * XK
    ZK2 = ZK * ZK
    return (2 * ZK2 - 4 * XK2) / ZK2

def eval_new(XK, ZK, X, Z):
    t0 = (XK + ZK) * (X - Z)
    t1 = (XK - ZK) * (X + Z)
    return X * (t0 + t1), Z * (t0 - t1)

def on_curve(A, X, Z):
    if Z == 0:
        return True, "inf"
    xv = X / Z
    E = EllipticCurve(Fp2, [0, A, 0, 1, 0])
    try:
        P = E.lift_x(xv)
        return True, P.order()
    except Exception:
        return False, None

A = Fp2(6)
E0 = EllipticCurve(Fp2, [0, A, 0, 1, 0])
P = E0(Fp2(69 + 22*i), Fp2(70 + 49*i))
Q = E0(Fp2(25 + 21*i), Fp2(56 + 28*i))
Spt = P + Q
S = (Spt.x(), Fp2(1))

print("=== corrected chain ===")
for i in [2, 1, 0]:
    K = S
    for j in range(i):
        K = xDBL(A, K[0], K[1])
    ok, ordK = on_curve(A, K[0], K[1])
    K2 = xDBL(A, K[0], K[1])
    kx = K[0]/K[1] if K[1] != 0 else "inf"
    print("i=%s K on curve %s ord %s 2K inf %s Kx %s" % (i, ok, ordK, K2[1]==0, kx))
    A = get_next_A_new(K[0], K[1])
    S = eval_new(K[0], K[1], S[0], S[1])
    okS, ordS = on_curve(A, S[0], S[1])
    sx = S[0]/S[1] if S[1] != 0 else "inf"
    print("  A'=%s S on curve %s ord %s Sx %s" % (A, okS, ordS, sx))
