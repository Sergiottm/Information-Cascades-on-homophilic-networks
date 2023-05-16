import networkx as nx
import ndlib.models.ModelConfig as mc
import ndlib.models.epidemics as ep
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Definimos los valores de grado y threshold que queremos probar
grados = np.linspace(1, 10, num=10, dtype=int)
thresholds = np.linspace(0.1, 1, num=10)



# Creamos una matriz de resultados para almacenar los valores de nodos infectados
resultados = np.zeros((len(grados), len(thresholds)))

# Iteramos sobre los valores de grado y threshold, ejecutando el modelo para cada combinación
for i, grado in enumerate(grados):
    for j, threshold in enumerate(thresholds):
        # Creamos la red de Erdos-Renyi
        g = nx.erdos_renyi_graph(1000, grado/1000)

        # Creamos el modelo de Independent Cascades
        model = ep.IndependentCascadesModel(g)

        # Configuramos los parámetros del modelo
        config = mc.Configuration()
        config.add_model_parameter('fraction_infected', 0.01)

        for e in g.edges():
            config.add_edge_configuration('threshold', e, threshold)

        model.set_initial_status(config)

        # Ejecutamos el modelo
        iterations = model.iteration_bunch(100)
        trends = model.build_trends(iterations)

        # Almacenamos el número total de nodos infectados
        total_infected = []
        for iteration in iterations:
            total_infected.append(iteration['node_count'][1])


        num_total_infectados = sum(total_infected)
        resultados[i, j] = num_total_infectados



gp, ti = [], []
for i in range(len(grados)):
    for j in range(1, len(thresholds)):
        if i > 0 and j > 0:
            gp.append(grados[i]/thresholds[j])
            ti.append(resultados[i, j])


# Creamos la figura y el eje tridimensional
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Creamos la superficie tridimensional
X, Y = np.meshgrid(grados, thresholds)
ax.plot_surface(X, Y, resultados, cmap='plasma')

# Configuramos los ejes y los títulos
ax.set_xlabel('Average Grade')
ax.set_ylabel('Threshold')
ax.set_zlabel('Infected nodes')


# Mostramos el gráfico
plt.savefig("ojo.png")

# Creamos la figura 2D
fig2 = plt.figure()
s = [10]*len(gp)
ax2 = fig2.add_subplot(111)
ax2.scatter(gp, ti, s=s)
ax2.set_xlabel('Grade/Threshold')
ax2.set_ylabel('Infected nodes')
ax2.set_title('Relationship between degree/threshold and infected nodes')

plt.savefig("critic.png")



