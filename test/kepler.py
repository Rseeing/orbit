import numpy as np

def kepler_E(Me, e, tolerance=1e-8):
    
    # 检查
    if not 0 <= e < 1:
        raise ValueError(f"椭圆偏心率e必须在[0, 1)范围内，当前值: {e}")
    
    # 初值选择
    if Me <= np.pi:
        E = Me + e/2
    else:
        E = Me - e/2
        
    # 牛顿迭代
    while True:
        f = E - e * np.sin(E) - Me 
        f_D = 1 - e * np.cos(E)
        ratio = f / f_D 
        
        if abs(ratio) < tolerance:
            return E
        
        E = E - ratio
    

def kepler_F(Mh, e, tolerance=1e-8):
    
    # 检查
    if e <= 1:
        raise ValueError(f"双曲线偏心率e必须大于1，当前值: {e}")
    
    # 初值选择
    F = Mh
        
    # 牛顿迭代
    while True:
        f = e * np.sinh(F) - F - Mh
        f_D = e * np.cosh(F) - 1
        ratio = f / f_D
        
        if abs(ratio) < tolerance:
            return F
        
        F = F - ratio


        
######全局
def stumpff_S(z):
    if z > 0:
        sz = np.sqrt(z)
        return (sz - np.sin(sz)) / (sz * z)
    elif z < 0:
        sz = np.sqrt(-z)
        return (np.sinh(sz) - sz) / (sz * (-z))
    else:
        return 1/6
    
def stumpff_C(z):
    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        return 0.5

def kepler_chi(delta_t, r0, vr0, a, mu=1.0, tolerance=1e-12):
    
    #参数准备 + 初值
    sqrt_mu = np.sqrt(mu)
    if a != 0:
        alpha = 1.0 / a   #   椭圆轨道（a > 0）：α > 0 ； 双曲线轨道（a < 0）：α < 0 
    else:
        alpha = 0.0   #  抛物线轨道（a = ∞）：α = 0
    chi = sqrt_mu * delta_t * abs(alpha)

    while True:
        z = alpha * chi * chi
        C = stumpff_C(z)
        S = stumpff_S(z)
        
        f = r0 * vr0 / sqrt_mu * chi**2 * C + (1 - alpha * r0) * chi**3 * S + r0 * chi - sqrt_mu * delta_t
        f_D = r0 * vr0 / sqrt_mu * chi * (1 - alpha * chi**2 * S) + (1 - alpha * r0) * chi**2 * C + r0
        
        ratio = f / f_D
        chi = chi - ratio
        
        if abs(ratio) < tolerance:
            return chi
