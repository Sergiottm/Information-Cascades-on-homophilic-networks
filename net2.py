import random
import networkx as nx
import ndlib.models.ModelConfig as mc
import ndlib.models.epidemics as ep
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from ndlib.viz.mpl.DiffusionTrend import DiffusionTrend

def calculate_taa(n, g, c, sa, steps):
    sb=sa

    for i in range(steps):
        # Escogemos un nodo al azar
        focal_node = random.choice(list(g.nodes()))

        # Seleccionamos un nodo candidato como vecino
        if random.random() < c:
            # Triadic closure
            neighbors = list(g.neighbors(focal_node))
            if len(neighbors) > 0:
                neighbor_node = random.choice(neighbors)
                while neighbor_node == focal_node: 
                    neighbor_node = random.choice(neighbors)

                # Buscamos un nodo candidato adicional que tenga un vecino común con el nodo vecino seleccionado y el nodo focal original
                additional_candidate_nodes = [n for n in list(g.neighbors(neighbor_node)) if n != focal_node and n not in neighbors and len(set(g.neighbors(n)).intersection(neighbors)) > 0]
                if additional_candidate_nodes:
                    candidate_node = random.choice(additional_candidate_nodes)
                else:
                    continue
                    
            else:
                continue
                    
        else:
            # Si no sale ninguno escogemos uno aleatorio
            candidate_node = random.choice(list(g.nodes()))
            neighbor_node=focal_node
            while candidate_node == focal_node: 
                candidate_node = random.choice(list(g.nodes()))


        # Aceptamos o rechazamos el nodo en función del bias
        if g.nodes[focal_node]["group"] == g.nodes[candidate_node]["group"]:
            # Mismo grupo
            if g.nodes[focal_node]["group"] == "A":
                accept_prob = sa
            else:
                accept_prob = sb
        else:
            # Diferente
            if g.nodes[focal_node]["group"] == "A":
                accept_prob = 1 - sa
            else:
                accept_prob = 1 - sb

        if random.random() < accept_prob:
            # Añadimos arista entre el nodo focal y el vecino
            if not g.has_edge(focal_node, candidate_node):
                # Quitamos una arista aleatoria del nodo focal
                # Quitamos una arista aleatoria del nodo focal, que no sea la que lo une con neighbor_node
                focal_edges = list(g.edges(focal_node))
                if len(focal_edges) > 1:
                    if g.has_edge(focal_node, neighbor_node):
                        neighbors_edge=list(g.edges(neighbor_node))
                        if len(neighbors_edge)>1:
                            focal_edges.remove((focal_node, neighbor_node))
                            if not g.has_edge(focal_node, neighbor_node):
                                g.add_edge(focal_node, candidate_node)

                        # Si aún quedan aristas, seleccionamos una aleatoria y la eliminamos
                    if len(focal_edges)>1:
                            removed_edge = random.choice(focal_edges)
                            re_0=list(g.edges(removed_edge[0]))
                            re_1=list(g.edges(removed_edge[1]))
                            if len(re_0)>1 and len(re_1)>1:
                                g.remove_edge(removed_edge[0], removed_edge[1])
                                if not g.has_edge(removed_edge[0], removed_edge[1]):
                                    g.add_edge(focal_node, candidate_node)

                        
                
              

    # Calculamos la homofilia y el bias de probabilidad
    aa_edges = 0
    bb_edges = 0
    ab_edges = 0

    for u, v in g.edges():
        if g.nodes[u]['group'] == 'A' and g.nodes[v]['group'] == 'A':
            aa_edges += 1
        elif g.nodes[u]['group'] == 'B' and g.nodes[v]['group'] == 'B':
            bb_edges += 1
        else:
            ab_edges += 1

    total_edges = g.number_of_edges()

    Paa =aa_edges/total_edges
    Pab = ab_edges/total_edges
    Pbb=bb_edges/total_edges

    Taa= 2*Paa/(2*Paa +Pab)
    pos = nx.spring_layout(g)
    return g, Taa

        



def run_simulation(g, num_iterations):
    a_list = [node for node, data in g.nodes(data=True) if data['group'] == 'A']
    infected_nodes = [random.choice(a_list)]

    # Crear un objeto de configuración del modelo IC
    model = ep.IndependentCascadesModel(g)

    

    # Configurar las probabilidades de adopción para todos los nodos
    config = mc.Configuration()
    config.add_model_initial_configuration("Infected", infected_nodes)

    # Creamos un diccionario con los valores de threshold para cada combinación de grupos
    group_thresholds = {('A', 'A'): 0.1,
                        ('A', 'B'): 0.5,
                        ('B', 'A'): 0.5,
                        ('B', 'B'): 0.5}

    # Recorremos todos los enlaces del grafo
    for e in g.edges():

        # Obtenemos los grupos de los nodos que conecta el enlace
        group1 = g.nodes[e[0]]['group']
        group2 = g.nodes[e[1]]['group']

        # Asignamos el threshold correspondiente en función de los grupos de los nodos
        if (group1, group2) in group_thresholds:
            threshold = group_thresholds[(group1, group2)]

        # Configuramos el threshold para el enlace
        config.add_edge_configuration("threshold", e, threshold)

    # Establecemos la configuración inicial del modelo
    model.set_initial_status(config)
        # Ejecutar la simulación
    iterations = model.iteration_bunch(num_iterations)
    trends = model.build_trends(iterations)

    total_infected = []
    a_infected = []
    b_infected = []
    a_infected_count = 0
    b_infected_count = 0
    a_inf = []
    b_inf = []

    for iteration in iterations:
        total_infected.append(iteration['node_count'][1])
        for node, status in iteration['status'].items():
            if status == 1 and g.nodes[node]["group"] == "A":
                a_infected_count += 1

            if status == 1 and g.nodes[node]["group"] == "B":
                b_infected_count += 1

        a_inf.append(a_infected_count)
        b_inf.append(b_infected_count)

    num_total_infectados = sum(total_infected)

    return total_infected, a_inf, b_inf




sa=0.75
steps=200000
p=0.005
c=0.75
n=1000




h = nx.erdos_renyi_graph(n, p)

node_labels = ["A"] * (n//2) + ["B"] * (n//2)
random.shuffle(node_labels)
for i, node in enumerate(h.nodes()):
    group = node_labels[i]
    h.nodes[node]["group"] = group

g, Taa = calculate_taa(n, h, c, sa, steps)
a_inf_list = []
total_infected_list = []
b_inf_list = []

for i in range(1, 1001):
    # Calcular la TAA

    # Crear las listas para guardar los resultados
    a_inf = []
    total_infected = []
    b_inf = []

    # Correr la simulación
    total_infected, a_inf, b_inf = run_simulation(g, 40)

    # Asignar los nombres de las listas basadas en el valor de i
    locals()[f"a_inf_{i}"] = a_inf
    locals()[f"total_infected_{i}"] = total_infected
    locals()[f"b_inf_{i}"] = b_inf

    # Agregar las listas a la lista correspondiente
    a_inf_list.append(a_inf)
    total_infected_list.append(total_infected)
    b_inf_list.append(b_inf)

# Calcular las listas promedio
a_inf_avg = [sum(x) / len(x) for x in zip(*a_inf_list)]
total_infected_avg = [sum(x) / len(x) for x in zip(*total_infected_list)]
b_inf_avg = [sum(x) / len(x) for x in zip(*b_inf_list)]

plt.plot(a_inf_avg, label='Group A')

# Crear gráfica para el número de nodos B infectados
plt.plot(b_inf_avg, label='Group B')

# Configurar la leyenda y los títulos de la gráfica
plt.legend()
plt.xlabel('Model step')
plt.ylabel('Number of infected nodes')

plt.savefig("b-b y a-b igual 0,5, a-a 0,1 con algortimo.png")

