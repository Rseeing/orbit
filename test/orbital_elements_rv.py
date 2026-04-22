import numpy as np

# 地球引力常数 (m^3/s^2)
MU_EARTH = 3.986004418e5  # km³/s²

def orbital_elements(r, v, mu=MU_EARTH):
    
    r = np.array(r, dtype=float)
    v = np.array(v, dtype=float)
        
    # 1. 距离
    r_norm = np.linalg.norm(r)
    
    # 2. 速度大小
    v_norm = np.linalg.norm(v)
    
    # 3. 径向速度大小
    v_radial = np.dot(r, v) / r_norm
    
    # 4. 角动量矢量
    h_vec = np.cross(r, v)
    
    # 5. 比角动量的模
    h = np.linalg.norm(h_vec)
        
    # 6. 倾角
    inclination = np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0))
    
    # 7. N矢量（升交点矢量）
    K = np.array([0, 0, 1])
    N_vec = np.cross(K, h_vec)
    
    # 8. N的模
    N_norm = np.linalg.norm(N_vec)
    
    # 9. 升交点赤经
    if N_norm > 1e-12:
        Omega_RAAN = np.arccos(np.clip(N_vec[0] / N_norm, -1.0, 1.0))
        if N_vec[1] < 0:
           Omega_RAAN = 2 * np.pi - Omega_RAAN
    else:
        Omega_RAAN = 0.0
    
    # 10. 偏心率矢量
    e_vec = (1/mu) * ((v_norm**2 - mu/r_norm) * r - r_norm * v_radial * v)
    
    # 11. 偏心率
    e = np.linalg.norm(e_vec)
    
    # 12. 近地点幅角
    if N_norm > 1e-12 and e > 1e-12:
        omega_argp = np.arccos(np.clip(np.dot(N_vec, e_vec) / (N_norm * e), -1.0, 1.0))
        if e_vec[2] < 0:
            omega_argp = 2 * np.pi - omega_argp
    else:
        omega_argp = 0.0
    
    # 13. 真近点角
    if e > 1e-12:
        true_anomaly = np.arccos(np.clip(np.dot(e_vec, r) / (e * r_norm), -1.0, 1.0))
        if v_radial < 0:
            true_anomaly = 2 * np.pi - true_anomaly
    else:
        # 圆轨道
        if N_norm > 1e-12:
            true_anomaly = np.arccos(np.clip(np.dot(N_vec, r) / (N_norm * r_norm), -1.0, 1.0))
            if r[2] < 0:
                true_anomaly = 2 * np.pi - true_anomaly
        else:
            # 倾角为0
            true_anomaly = np.arccos(np.clip(r[0] / r_norm, -1.0, 1.0))
            if r[1] < 0:
                true_anomaly = 2 * np.pi - true_anomaly
                
        
    return {
        # 输入参数
        'position': r.tolist(),
        'velocity': v.tolist(),
        'mu': mu,
        
        # 1-3: 基本状态量
        'distance': float(r_norm),  # km
        'speed': float(v_norm),
        'radial_speed': float(v_radial),
        
        # 4-5: 角动量
        'angular_momentum_vector': h_vec.tolist(),
        'angular_momentum': float(h),
        
        # 6: 倾角
        'inclination': float(np.degrees(inclination)),
        'inclination_rad': float(inclination),
        
        # 7-8: 升交点矢量
        'node_vector': N_vec.tolist(),
        'node_magnitude': float(N_norm),
        
        # 9: 升交点赤经
        'RAAN': float(np.degrees(Omega_RAAN)),
        'RAAN_rad': float(Omega_RAAN),
        
        # 10-11: 偏心率
        'eccentricity_vector': e_vec.tolist(),
        'eccentricity': float(e),
        
        # 12: 近地点幅角
        'argument_of_perigee': float(np.degrees(omega_argp)),
        'argument_of_perigee_rad': float(omega_argp),
        
        # 13: 真近点角
        'true_anomaly': float(np.degrees(true_anomaly)),
        'true_anomaly_rad': float(true_anomaly),
        
    }

def rv(h, e, i, Omega, omega, theta, mu=MU_EARTH):
    
    # 位置向量（PQW近焦点坐标系）
    r_pqw = (h**2 / mu) / (1 + e * np.cos(theta)) * np.array([
    np.cos(theta),
    np.sin(theta),
    0
    ])
        
    # 速度向量（PQW近焦点坐标系）
    v_pqw = mu / h * np.array([
        -np.sin(theta),
        e + np.cos(theta),
        0
    ])
    
    # PQW近焦点坐标系 转 XYZ地心赤道坐标系
    # 旋转矩阵 
    c_Om, s_Om = np.cos(Omega), np.sin(Omega)
    c_i, s_i = np.cos(i), np.sin(i)
    c_om, s_om = np.cos(omega), np.sin(omega)
    
    Q_PX = np.array([
        [c_Om*c_om - s_Om*c_i*s_om, -c_Om*s_om - s_Om*c_i*c_om,  s_Om*s_i],
        [s_Om*c_om + c_Om*c_i*s_om, -s_Om*s_om + c_Om*c_i*c_om, -c_Om*s_i],
        [s_i*s_om,                   s_i*c_om,                   c_i     ]
    ])
    
    
    return Q_PX @ r_pqw, Q_PX @ v_pqw
        
def rv_deg(h, e, i_deg, Omega_deg, omega_deg, theta_deg, mu=MU_EARTH):
    ##  单位：度
    i = np.radians(i_deg)
    Omega = np.radians(Omega_deg)
    omega = np.radians(omega_deg)
    theta = np.radians(theta_deg)
    
    return rv(h, e, i, Omega, omega, theta, mu)


# 例4.1
if __name__ == "__main__":
    # 已知状态向量
    r = [-6045, -3490, 2500]
    v = [-3.457, 6.618, 2.533]
    
    elements = orbital_elements(r, v)
    
    print("轨道根数计算结果:")
    print("=" * 50)
    print(f"距离: {elements['distance']/1000:.3f} km")
    print(f"速度: {elements['speed']:.3f} m/s")
    print(f"径向速度: {elements['radial_speed']:.3f} m/s")
    print(f"角动量 h: {elements['angular_momentum']:.3e} m²/s")
    print(f"倾角 i: {elements['inclination']:.6f}°")
    print(f"升交点赤经 Ω: {elements['RAAN']:.6f}°")
    print(f"偏心率 e: {elements['eccentricity']:.8f}")
    print(f"近地点幅角 ω: {elements['argument_of_perigee']:.6f}°")
    print(f"真近点角 θ: {elements['true_anomaly']:.6f}°")
    print("=" * 50)
    
        
# 例4.2
if __name__ == "__main__":
    # 已知轨道根数
    h = 80000        
    e = 1.4           
    i = 30             
    Omega = 40        
    omega = 60         
    theta = 30         
    
    r, v = rv_deg(h, e, i, Omega, omega, theta)
    
    print("位置向量 r (km):")
    print(f"  x = {r[0]:.3f}")
    print(f"  y = {r[1]:.3f}")
    print(f"  z = {r[2]:.3f}")
    print(f"  距离 = {np.linalg.norm(r)/1000:.3f} km")
    
    print("\n速度向量 v (km/s):")
    print(f"  vx = {v[0]:.3f}")
    print(f"  vy = {v[1]:.3f}")
    print(f"  vz = {v[2]:.3f}")
    print(f"  速度 = {np.linalg.norm(v):.3f} m/s")
    
