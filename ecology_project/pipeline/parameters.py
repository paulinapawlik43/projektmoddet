import numpy as np

class KlausmeierParameters:
    """
    klasa przeliczająca parametry fizyczne modelu Klausmeiera na bezwymiarowe
    """
    def __init__(self, A, L, R, J, M, Dw, Dn):
        self.A = A   # opady [mm/h]
        self.L = L   # parowanie [1/h]
        self.R = R   # pobór wody [1/kg]
        self.J = J   # konwersja wody na biomasę [kg/mm]
        self.M = M   # śmiertelność [1/h]
        self.Dw = Dw # dyfuzja wody [m^2/h]
        self.Dn = Dn # dyfuzja biomasy [m^2/h]

    def get_scaling_factors(self):
        """ obliczanie stałych skalujących W0, N0, T0"""
        # L*T0 = 1
        T0 = 1.0 / self.L
        
        # R*T0*N0^2 = 1
        N0 = np.sqrt(self.L / self.R)
        
        # J*R*W0*T0*N0 = 1
        W0 = 1.0 / (self.J * self.R * T0 * N0)

        # przyjmujemy X0 = 1.0 
        # 1 jednostka siatki to 1 metr
        X0 = 1.0
        
        return W0, N0, T0, X0

    def get_dimensionless_params(self):
        W0, N0, T0, X0 = self.get_scaling_factors()
        
        a = (self.A * T0) / W0
        m = self.M * T0
        d1 = (self.Dw * T0) / (X0**2)
        d2 = (self.Dn * T0) / (X0**2)
        
        return {
            "a": a,
            "m": m,
            "d1": d1,
            "d2": d2
        }