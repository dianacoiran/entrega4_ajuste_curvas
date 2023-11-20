import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit

def df_dt(x, t, alpha, beta, gamma, delta):
    """Función del sistema en forma canónica"""
    dx = alpha * x[0] - beta * x[0] * x[1]
    dy = - gamma * x[1] + delta * x[0] * x[1]
    return np.array([dx, dy])

# Datos observados (simulados)
tf = 500
N = 50
t = np.linspace(0, tf, N)
conds_iniciales = np.array([10, 70])  # X0, presa. Y0, depredador
params_reales = [0.1, 0.02, 0.3, 0.01]
solucion_real = odeint(df_dt, conds_iniciales, t, args=tuple(params_reales))

# Añadir ruido a los datos simulados para imitar observaciones experimentales
np.random.seed(42)
datos_ruidosos = solucion_real + np.random.normal(0, 2, solucion_real.shape)

# Función de ajuste
def ajuste_lotka_volterra(t, alpha):
    params_fijos = [alpha, 0.02, 0.3, 0.01]  # Mantener beta, gamma y delta fijos
    solucion = odeint(df_dt, conds_iniciales, t, args=tuple(params_fijos))
    return solucion.flatten()

# Ajustar el modelo a los datos ruidosos para estimar el parámetro alpha
popt, pcov = curve_fit(ajuste_lotka_volterra, t, datos_ruidosos.flatten(), p0=[0.1])

plt.figure(figsize=(12, 8))

# Graficar resultados
plt.plot(t, datos_ruidosos[:, 0], 'o', label='Presa (Datos ruidosos)')
plt.plot(t, datos_ruidosos[:, 1], 'o', label='Depredador (Datos ruidosos)')
plt.plot(t, ajuste_lotka_volterra(t, *popt).reshape((N, 2))[:, 0], label='Presa (Ajuste)')
plt.plot(t, ajuste_lotka_volterra(t, *popt).reshape((N, 2))[:, 1], label='Depredador (Ajuste)')
plt.legend()
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Ajuste de Lotka-Volterra')
plt.show()

# Calcular y graficar residuales
residuales_presa = datos_ruidosos[:, 0] - ajuste_lotka_volterra(t, *popt).reshape((N, 2))[:, 0]
residuales_depredador = datos_ruidosos[:, 1] - ajuste_lotka_volterra(t, *popt).reshape((N, 2))[:, 1]

plt.subplot(1, 2, 1)
plt.plot(t, residuales_presa, 'o', label='Residuales Presa')
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.xlabel('Tiempo')
plt.ylabel('Residuales')
plt.title('Residuales de la Presa')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t, residuales_depredador, 'o', label='Residuales Depredador')
plt.axhline(0, color='black', linestyle='--', linewidth=1)
plt.xlabel('Tiempo')
plt.ylabel('Residuales')
plt.title('Residuales del Depredador')
plt.legend()

plt.tight_layout()
plt.show()


# Mostrar el valor estimado de alpha
#print(f"Valor estimado de alpha: {popt[0]}")
