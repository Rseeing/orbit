import sympy as sp
import math

def newton_solver():
    var_name = input("请输入未知数符号（如 x, t）：").strip()
    expr_str = input("请输入方程 f(x)=0 的表达式(字符串格式)：").strip()
    x0_input = input("请输入初始估计值 x0（可以输入多个，用空格或逗号分隔）：")
    digits = int(input("请输入有效数字位数（如 8）："))

    # 支持多个初值
    x0_list = []
    for val in x0_input.replace(',', ' ').split():
        try:
            x0_list.append(float(val))
        except ValueError:
            print(f"无效输入：{val}，已跳过。")

    var = sp.symbols(var_name)   # 转符号变量
    expr = sp.sympify(expr_str)  # 字符串转符号表达式
    d_expr = sp.diff(expr, var)  # 求导

    f = sp.lambdify(var, expr, "math")  # 符号表达式转 Python 函数
    df = sp.lambdify(var, d_expr, "math")  

    tol = 10 ** (-digits)

    # 对每个初值迭代
    for idx, x0 in enumerate(x0_list, 1):
        print("\n" + "="*60)
        print(f"开始迭代第 {idx} 个初值: x0 = {x0}")
        x = x0
        for k in range(1, 100):
            fx = f(x)
            dfx = df(x)
            
            if abs(fx) < 1e-64 :  # #####################
                print("f(x) <1e-64，迭代终止")
                print(f"近似根 = {x_new:.{digits}f}")
                break
            
            x_new = x - fx / dfx
            print(f"迭代 {k:2d}: x = {x_new:.20f}, f(x) = {fx:.32f}, df(x) = {dfx:.15f}") 
            
            if abs(x_new - x) < tol or abs(fx/dfx) < tol:
                print("-"*60)
                print(f"收敛于第 {k} 次迭代")
                print(f"近似根 = {x_new:.{digits}f}")
                break
            
            x = x_new
        else:
            print("未在最大迭代次数内收敛")

if __name__ == "__main__":
    newton_solver()


#  10*exp(sin(x))-(x^2-5*x+4)
#  tan(x)-tanh(x)
#  E、e可能被识别为2.718281828459045
