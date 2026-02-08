import numpy as np

# algorytm Thomasa
def solve_vectorized(a, b, c, d):
    n = d.shape[0]
    
    # zapisujemy same przekątne bo poza nimi są zera
    ac = a[:, np.newaxis]        # dolna przekątma
    bc = b[:, np.newaxis].copy() # glowna przekątna
    cc = c[:, np.newaxis]        # gorna przekątna
    dc = d                       # prawa strona rownania
    
    # eliminacja
    for i in range(1, n):
        m = ac[i-1] / bc[i-1]
        bc[i] = bc[i] - m * cc[i-1]
        dc[i] = dc[i] - m * dc[i-1]
        
    # podstawianie wsteczne
    x = np.zeros_like(dc)
    x[-1] = dc[-1] / bc[-1]
    for i in range(n-2, -1, -1):
        x[i] = (dc[i] - cc[i] * x[i+1]) / bc[i]
    return x

class KlausmeierSolver:
    def __init__(self, N, DT, params):
        
        self.N = N
        self.DT = DT
        self.d1 = params['d1']
        self.d2 = params['d2']
        self.m = params['m']
        self.a = params['a']
        
        # macierze
        mu_u = (self.d1 * self.DT) / 2.0 
        self.main_u = (1 + 2 * mu_u) * np.ones(self.N) 
        self.main_u[0]=1; self.main_u[-1]=1
        self.off_u  = -mu_u * np.ones(self.N-1)       
        self.off_u[0]=0;  self.off_u[-1]=0

        mu_v = (self.d2 * self.DT) / 2.0
        self.main_v = (1 + 2 * mu_v) * np.ones(self.N) 
        self.main_v[0]=1; self.main_v[-1]=1
        self.off_v  = -mu_v * np.ones(self.N-1)         
        self.off_v[0]=0;  self.off_v[-1]=0

        self.U = np.zeros((self.N, self.N))
        self.V = np.zeros((self.N, self.N))

    def solve_step(self):
        u, v = self.U, self.V
        
        # limit
        u = np.clip(u, 0, 50); v = np.clip(v, 0, 50)
        
        uv2 = u * (v**2)
        du = self.a - u - uv2
        dv = uv2 - self.m * v
        u += du * self.DT
        v += dv * self.DT
        
        # Dirichlet
        u[0,:]=0; u[-1,:]=0
        u[:,0]=0; u[:,-1]=0
        v[0,:]=0; v[-1,:]=0
        v[:,0]=0; v[:,-1]=0

        # Y
        rhs_u = u + (self.d1 * self.DT / 2.0) * self.laplacian(u, 0)
        rhs_v = v + (self.d2 * self.DT / 2.0) * self.laplacian(v, 0)
        rhs_u[:,0]=0; rhs_u[:,-1]=0; rhs_v[:,0]=0; rhs_v[:,-1]=0
        u = solve_vectorized(self.off_u, self.main_u, self.off_u, rhs_u)
        v = solve_vectorized(self.off_v, self.main_v, self.off_v, rhs_v)

        # X
        rhs_u = u + (self.d1 * self.DT / 2.0) * self.laplacian(u, 1)
        rhs_v = v + (self.d2 * self.DT / 2.0) * self.laplacian(v, 1)
        rhs_u[0,:]=0; rhs_u[-1,:]=0; rhs_v[0,:]=0; rhs_v[-1,:]=0
        u = solve_vectorized(self.off_u, self.main_u, self.off_u, rhs_u.T).T
        v = solve_vectorized(self.off_v, self.main_v, self.off_v, rhs_v.T).T
        
        self.U, self.V = u, v

    def laplacian(self, arr, axis):
        res = np.zeros_like(arr)
        if axis == 0: 
            res[1:-1, :] = arr[:-2, :] - 2*arr[1:-1, :] + arr[2:, :]
        else:         
            res[:, 1:-1] = arr[:, :-2] - 2*arr[:, 1:-1] + arr[:, 2:]
        return res