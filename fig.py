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
                    # Eliminamos la arista que conecta al nodo focal con neighbor_node, si existe
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



sa=0.75
steps=350000
p=0.2
c=0.75
n=100



grafo_conexo = False

while not grafo_conexo:
    # Generar el grafo utilizando la distribución de Erdos-Renyi
    h = nx.erdos_renyi_graph(n, p)
    
    # Verificar si el grafo es conexo
    grafo_conexo = nx.is_connected(h)

 # Asignar grupos a los nodos
node_labels = ["A"] * (n//2) + ["B"] * (n//2)
random.shuffle(node_labels)
for i, node in enumerate(h.nodes()):
    group = node_labels[i]
    h.nodes[node]["group"] = group



pos = nx.spring_layout(h)
color_map = {'A': 'red', 'B': 'blue'}
colors = [color_map[h.nodes[node]['group']] for node in h.nodes()]

plt.figure(figsize=(10, 8))
nx.draw_networkx(h, pos, node_color=colors, edge_color='black', width=0.5, with_labels=False, node_size=50)
plt.axis('off')

# Guardamos el grafo
plt.savefig('grafo1.png')

# Calcula el grado medio
grado_medio = h.number_of_edges()

print("Grado medio del grafo:", grado_medio)

g, Taa = calculate_taa(n, h, c, sa, steps)

grado_medio = g.number_of_edges()

print("Grado medio del grafo:", grado_medio)

# Dibujamos el grafo
pos1 = nx.spring_layout(g)

colors1 = [color_map[g.nodes[node]['group']] for node in g.nodes()]

plt.figure(figsize=(10, 8))
nx.draw_networkx(g, pos1, node_color=colors1, edge_color='black', width=0.5, with_labels=False, node_size=50)
plt.axis('off')

# Guardamos el grafo
plt.savefig('grafo2.png')

