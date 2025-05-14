#
# Copyright (c) 2025 INRIA
#

import casadi as cs
import numpy as np


def LTL(M):
    n = M.shape[0]
    U = cs.SX(n, n)

    for i in range(n - 1, -1, -1):
        for j in range(n - 1, i, -1):
            s = sum(U[i, k] * U[j, k] for k in range(j + 1, n))
            U[i, j] = (M[i, j] - s) / U[j, j]

        s = sum(U[i, k] ** 2 for k in range(i + 1, n))
        U[i, i] = cs.sqrt(M[i, i] - s)

    return U


def backsolve(L, y):
    # return x such that L*x = y

    n = L.shape[0]
    x = cs.SX(n, y.shape[1])

    for i in range(n):
        s = sum(L[i, j] * x[j, :] for j in range(i))
        x[i, :] = (y[i, :] - s) / L[i, i]

    return x


def forwardsolve(L, y):
    # return x such that L^T*x = y

    n = L.shape[0]
    x = cs.SX(n, y.shape[1])

    for i in range(n - 1, -1, -1):
        s = sum(L[j, i] * x[j, :] for j in range(i + 1, n))
        x[i, :] = (y[i, :] - s) / L[i, i]

    return x


def LTL_solve(U, y):
    z = forwardsolve(U.T, y)
    x = backsolve(U.T, z)

    return x


def make_symmetric(M):
    m = M.shape[0]
    n = M.shape[1]
    assert m == n, "Matrix must be square"
    for i in range(m):
        for j in range(i):
            M[i, j] = M[j, i]

    return M


def get_num_operations(fun):
    return fun.n_instructions() - fun.nnz_in() - fun.nnz_out()


n = 4
A = cs.SX.sym("A", n, n)
b = cs.SX.sym("b", n, 1)
U = LTL(A)

U_fun = cs.Function("U_fun", [A], [U])
U_fun.generate(f"U_fun{n!s}_codegen.c", {"with_header": True})

A_inv_chol = make_symmetric(LTL_solve(U, cs.SX.eye(n)))
A_inv_chol_fun = cs.Function("A_inv_chol_fun", [A], [A_inv_chol])
A_inv_chol_fun.generate(f"A_inv_chol_fun{n!s}_codegen.c", {"with_header": True})

inv_A = make_symmetric(cs.inv(A))
inv_A_fun = cs.Function("inv_A_fun", [A], [inv_A])
inv_A_fun.generate(f"inv_A_fun{n!s}_codegen.c", {"with_header": True})

U_sym = cs.SX.sym("U", n, n)
x_sol = LTL_solve(U_sym, b)
x_sol_fun = cs.Function("x_sol_fun", [U_sym, b], [x_sol])
x_sol_fun.generate(f"x_sol_fun{n!s}_codegen.c", {"with_header": True})

A_val = np.random.rand(n, n)
A_val = A_val @ A_val.T
b_val = np.random.rand(n, 1)

U_val = U_fun(A_val)
print("U_val", U_val)

A_inv_chol_val = A_inv_chol_fun(A_val)
print("A_inv_chol_val", A_inv_chol_val)

inv_A_val = inv_A_fun(A_val)
print("inv_A_val", inv_A_val)
assert np.allclose(A_inv_chol_val, inv_A_val)

x_sol_val = x_sol_fun(U_val, b_val)
print("x_sol_val", x_sol_val)
print("A_inv@b_val", inv_A_val @ b_val)
assert np.allclose(x_sol_val, inv_A_val @ b_val)

print("Num operations Casadi inverse", get_num_operations(inv_A_fun))
print("Num operations inverse from LTL", get_num_operations(A_inv_chol_fun))
print("Num operations LTL", get_num_operations(U_fun))
print("Num operations LTL solve", get_num_operations(x_sol_fun))
