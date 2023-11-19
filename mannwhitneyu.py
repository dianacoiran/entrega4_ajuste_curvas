import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("tabla_origen.csv")

## columna name se convierte a indice
df.set_index("name", inplace=True)

## de la tabla original extraigo las columnas de interes
tabla_mod = df[["G1_CD8_FP","G2_CD8_FP","G3_CD8_FP",
                "G1_CD8_FN","G2_CD8_FN","G3_CD8_FN"]]

## listas para almacenar resultados
p_values_list = []
Fold_Change_list = []
name_list = []

## Iterar sobre filas del DataFrame
for index, row in tabla_mod.iterrows():
    CD8_FP = list(row[:3])  # columnas de datos positivos
    CD8_FN = list(row[3:])  # columnas de datos negativos

## termino sumando para evitar problemas con resultados igual a cero
    CD8_FP_median = row[:3].median() + 0.001
    CD8_FN_median = row[3:].median() + 0.001

    Fold_Change = CD8_FP_median / CD8_FN_median  # Calcular Fold Change

    _, p_value = stats.mannwhitneyu(CD8_FP, CD8_FN, alternative="two-sided")
    p_values_list.append(p_value)
    Fold_Change_list.append(Fold_Change)
    name_list.append(index)
    #print(index, CD8_FP, CD8_FN, p_value)  # Imprimir resultados

## Crear DataFrame con resultados
result_tabla_mod = pd.DataFrame({"gene": name_list, "FC": Fold_Change_list, "p-value": p_values_list})

# # Corrección de Bonferroni

reject, p_corrected, _, _ = multipletests(result_tabla_mod['p-value'], method='bonferroni')

# # Agregar columnas con los resultados corregidos al DataFrame
result_tabla_mod['p-value_corrected'] = p_corrected
result_tabla_mod['reject_null'] = reject

#  Imprimir el DataFrame resultante con los resultados corregidos
# #print(result_tabla_mod)

#### graficacion con bonferroni ###

result_tabla_mod['log2_FC'] = np.log2(result_tabla_mod['FC'])
result_tabla_mod['log10_p-value'] = -np.log10(result_tabla_mod['p-value_corrected'])

 ## Configurar el Volcano Plot
plt.figure(figsize=(10, 6))

# # Todos los puntos en amarillo (no significativos)
plt.scatter(result_tabla_mod['log2_FC'], result_tabla_mod['log10_p-value'], color='yellow', alpha=0.5, label='No significativo')

# # # Puntos significativos en rojo
significant_genes = result_tabla_mod[result_tabla_mod['reject_null']]
plt.scatter(significant_genes['log2_FC'], significant_genes['log10_p-value'], color='red', label='Significativo')

# # # Añadir líneas para indicar umbral de significancia
plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, label='Umbral de significancia')

# # Etiquetas y título
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot')
plt.legend()

# # # Añadir nombres de genes para genes significativos
for gene, x, y in zip(significant_genes['gene'], significant_genes['log2_FC'], significant_genes['log10_p-value']):
    plt.annotate(gene, (x, y), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8)

plt.show()


# ###### grafica sin bonferroni ####

##Calcular log2 del Fold Change y log10 del valor p
result_tabla_mod['log2_FC'] = np.log2(result_tabla_mod['FC'])
result_tabla_mod['log10_p-value'] = -np.log10(result_tabla_mod['p-value'])

plt.figure(figsize=(10, 6))
## Todos los puntos en amarillo (no significativos)
plt.scatter(result_tabla_mod['log2_FC'], result_tabla_mod['log10_p-value'], color='yellow', alpha=0.5, label='No significativo')
## Puntos significativos en rojo
significant_genes = result_tabla_mod[result_tabla_mod['p-value'] < 0.05]
plt.scatter(significant_genes['log2_FC'], -np.log10(significant_genes['p-value']), color='red', label='Significativo')

## Añadir líneas para indicar umbral de significancia
plt.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, label='Umbral de significancia')

## Etiquetas y título
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot')
plt.legend()

for gene, x, y in zip(significant_genes['gene'], significant_genes['log2_FC'], -np.log10(significant_genes['p-value'])):
    plt.annotate(gene, (x, y), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8)

plt.show()

#########################################################################




