import numpy as np
from pipeline.solver import KlausmeierSolver

def calculate_dispersion_relation(a, m, d1, d2, k_max=7.0):
    """
    funkcja obliczająca relację dyspersji lambda(k).
    """
    # wyznaczenie punktu równowagi (v* > 1)
    delta = a**2 - 4 * m**2
    if delta < 0:
        return None, None
    
    # wybieramy pierwiastek z plusem (stan v > 1)
    v_star = (a + np.sqrt(delta)) / (2 * m)
    u_star = m / v_star

    # macierz Jacobiego 
    J11 = -1 - v_star**2
    J12 = -2 * m          # uproszczone z -2uv
    J21 = v_star**2
    J22 = m               # uproszczone z 2uv - m

    # relacja dyspersji (z dyfuzją)
    k_values = np.linspace(0, k_max, 1000)
    re_lambda_max = []

    for k in k_values:
        # macierz M(k) = J - k^2 * D
        m11 = J11 - d1 * (k**2)
        m12 = J12
        m21 = J21
        m22 = J22 - d2 * (k**2)
        
        # obliczamy wartości własne macierzy
        # lambda^2 - Tr(M)*lambda + det(M) = 0
        tr_M = m11 + m22
        det_M = m11*m22 - m12*m21
        
        delta_lambda = tr_M**2 - 4*det_M
        
        # szukamy największej części rzeczywistej
        if delta_lambda >= 0:
            l1 = (tr_M + np.sqrt(delta_lambda)) / 2
            l2 = (tr_M - np.sqrt(delta_lambda)) / 2
            max_re = max(l1, l2)
        else:
            max_re = tr_M / 2
            
        re_lambda_max.append(max_re)

    return k_values, np.array(re_lambda_max)

def analyze_domain_size(N_values, dt, T_max, params): 
    mean_biomass_results = []
    max_biomass_results = []
    
    for N in N_values: 
        solver = KlausmeierSolver(N, dt, params)

        solver.U = np.ones((N, N)) * params['a'] 
        solver.V = np.random.uniform(0.5, 1.5, (N, N)) 
        
        for _ in range(int(T_max / dt)):
            solver.solve_step()

        v_center = solver.V[5:-5, 5:-5]
        if v_center.size > 0:
            mean_biomass_results.append(np.mean(v_center))
            max_biomass_results.append(np.max(v_center))
        else:
            mean_biomass_results.append(np.mean(solver.V))
            max_biomass_results.append(np.max(solver.V))
    
    return mean_biomass_results, max_biomass_results