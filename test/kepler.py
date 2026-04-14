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


