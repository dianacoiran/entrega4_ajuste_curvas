import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

df = pd.read_csv("tabla_origen.csv")

## columna name se convierte a indice
df.set_index("name", inplace=True)

# Seleccionar las columnas relevantes
data = df[["G1_CD8_FP", "G2_CD8_FP", "G3_CD8_FP", "G1_CD8_FN", "G2_CD8_FN", "G3_CD8_FN"]]

## Realizar la prueba de t-Student para cada gen
p_value_list = []
t_statistic_list = []
name_list = []

for index, row in data.iterrows():
    CD8_FP = list(row[:3])  # columnas de datos positivos
    CD8_FN = list(row[3:])  # columnas de datos negativos

## termino sumando para evitar problemas con resultados igual a cero
    CD8_FP_median = row[:3].median() + 0.001
    CD8_FN_median = row[3:].median() + 0.001
    Fold_Change = CD8_FP_median / CD8_FN_median  # Calcular Fold Change

    t_statistic, p_value = ttest_ind(CD8_FP, CD8_FN)
    t_statistic_list.append(t_statistic)
    p_value_list.append(p_value)
    name_list.append(index)
    print(index, CD8_FP, CD8_FN, p_value)  # Imprimir resultados

## Crear un DataFrame con los resultados

result_df = pd.DataFrame({"gene": name_list, "t-statistic": t_statistic_list, "p-value": p_value_list})
print(result_df)

# Crear DataFrame con resultados
result_df = pd.DataFrame({"gene": name_list, "t-statistic": t_statistic_list, "p-value": p_value_list})

# Corrección de Bonferroni
reject, p_corrected, _, _ = multipletests(result_df['p-value'], method='bonferroni')

# Agregar columnas con los resultados corregidos al DataFrame
result_df['p-value_corrected'] = p_corrected
result_df['reject_null'] = reject

# Imprimir el DataFrame resultante con los resultados corregidos
#print(result_df)

# Graficar Volcano Plot
plt.figure(figsize=(10, 6))

# Puntos no significativos en amarillo
plt.scatter(result_df['t-statistic'][~result_df['reject_null']],
            -np.log10(result_df['p-value'][~result_df['reject_null']]),
            color='yellow', alpha=0.5, label='No significativo')

# Puntos significativos en rojo
plt.scatter(result_df['t-statistic'][result_df['reject_null']],
            -np.log10(result_df['p-value'][result_df['reject_null']]),
            color='red', label='Significativo')

# Línea para indicar umbral de significancia
plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, label='Umbral de significancia')

# Etiquetas y título
plt.xlabel('t-statistic')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot')
plt.legend()

# Anotar genes para los puntos significativos
for gene, x, y in zip(result_df[result_df['reject_null']]['gene'],
                      result_df[result_df['reject_null']]['t-statistic'],
                      -np.log10(result_df[result_df['reject_null']]['p-value'])):
    plt.annotate(gene, (x, y), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

# Mostrar el gráfico
plt.show()

