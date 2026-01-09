import sympy as sp
import math

def newton_solver(f0, var, x0_list, digits=8, taylor_order=7):

    # 原函数和导数
    df0 = sp.diff(f0, var)

    # 泰勒展开（在 x=0）
    f_taylor = sp.series(f0, var, 0, taylor_order).removeO()
    df_taylor = sp.diff(f_taylor, var)

    # 数值函数
    f = sp.lambdify(var, f0, "math")
    df = sp.lambdify(var, df0, "math")
    f_T = sp.lambdify(var, f_taylor, "math")
    df_T = sp.lambdify(var, df_taylor, "math")

    tol = 10 ** (-digits)
    eps_taylor = 1e-3   

    roots = []

    for idx, x0 in enumerate(x0_list, 1):
        print("\n" + "="*60)
        print(f"开始迭代第 {idx} 个初值: x0 = {x0}")
        x = float(x0)

        for k in range(1, 100):

            if abs(x) < eps_taylor:
                fx = f_T(x)
                dfx = df_T(x)
                mode = "Taylor"
            else:
                fx = f(x)
                dfx = df(x)
                mode = "Exact"

            x_new = x - fx / dfx

            print(
                f"迭代 {k:2d}: x = {x_new:.20f} | "
                f"mode = {mode}"
            )

            if abs(x_new - x) < tol:
                print("-"*60)
                print(f"收敛于第 {k} 次迭代")
                print(f"近似根 = {x_new:.{digits}f}")
                roots.append(x_new)
                break

            x = x_new
        else:
            print("未在最大迭代次数内收敛")

    return roots


if __name__ == "__main__":
    x = sp.symbols('x')
    f1 = 10*sp.exp(sp.sin(x)) - (x**2 - 5*x + 4)
    #f1 = sp.tan(x) - sp.tanh(x)
    roots = newton_solver(f1, x, [0.1], digits=8)