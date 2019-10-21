# -*- coding: utf-8 -*-
import math 
import operator
from pprint import pprint
import numpy as np # for Floyd-Warshall matrices
 

def create_graph(directed = True, weighted = False): # TP1
	g = { 'directed': directed, 'weighted': weighted, 'nb_nodes': 0, 'nb_edges': 0, 'weight_attribute': None , 'nodes': {}, 'nodes_index' : {}, 'edges': {}  } # le sous dictionnaire 'nodes_index' a été ajouté et est un duplicata du dictionnaire 'nodes' et ayant comme attribut l'index du sommet. Ceci afin de trier le dictionnaire.
	return g


#Cette fonction attribut un index à chaques sommet, elle permetra a posteriori de trier le dictionnaire de manière répétitive.
def trigraphe(g):
	i = 0
	for u in sorted(g['nodes']):
		g['nodes_index'][u] = g['nodes'][u]
	for u in sorted(g['nodes_index']):
		g['nodes_index'][u] = i
		i = i+1
	return g


def add_node(g, n, attributes = None): # TP1
	if n not in g['nodes']: # ensure node does not already exist
		if attributes is None: # create empty attributes if not provided
			attributes = {}
		g['nodes'][n] = attributes
		g['edges'][n] = {} # init outgoing edges
		g['nb_nodes'] += 1
	return g['nodes'][n] # return node attributes


def add_edge(g, n1, n2, attributes = None): # TP1
	# create nodes if they do not exist
	if n1 not in g['nodes']: add_node(g, n1) # ensure n1 exists
	if n2 not in g['nodes']: add_node(g, n2) # ensure n2 exists
	# add edge(s) only if they do not exist
	if n2 not in g['edges'][n1]:
		if attributes is None:
			attributes = {}
		g['edges'][n1][n2] = attributes
		if not g['directed']:
			g['edges'][n2][n1] = g['edges'][n1][n2] # share the same attributes as n1->n2
		g['nb_edges'] += 1
	return g['edges'][n1][n2] # return edge atributes

def load_SIF(filename, directed=True): # TP1
	# line syntax: nodeD <relationship type> nodeE nodeF nodeB
	g = create_graph(directed) # new empty graph
	with open(filename) as f: # OPEN FILE
		# PROCESS THE REMAINING LINES
		row = f.readline().rstrip() # read next line and remove ending whitespaces
		while row:
			vals = row.split('\t') # split line on tab
			for i in range(2, len(vals)):
				att = { 'type': vals[1] } # set edge type
				add_edge(g, vals[0], vals[i], att)
			row = f.readline().rstrip() # read next line
	return g # return graph

def load_TAB(filename, directed=True): 
	g = create_graph(directed)
	with open(filename) as f: 
		# GET COLUMNS NAMES
		tmp = f.readline().rstrip()
		attNames= tmp.split('\t')
		# REMOVES FIRST TWO COLUMNS WHICH CORRESPONDS TO THE LABELS OF THE CONNECTED VERTICES
		attNames.pop(0)
		attNames.pop(0)
		# PROCESS THE REMAINING LINES
		row = f.readline().rstrip()
		while row:
			vals = row.split('\t')
			v1 = vals.pop(0)
			v2 = vals.pop(0)
			att = {}
			for i in range(len(attNames)):
				att[ attNames[i] ] = vals[i]
			add_edge(g, v1, v2, att)
			row = f.readline().rstrip() # NEXT LINE
	return g

#Plus court chemin d'un sommet à tous les sommets.
def BFS(g, s): 
	g['BFS'] = {'color':{},'d':{},'pi': {},'Q':[]}#Cree un sous dictionnaire 'BFS' attaché au dictionnaire g
	for u in sorted(g['nodes_index']):#Initialisation des sommets
		g['BFS']['color'][u]= 'WHITE' 
		g['BFS']['d'][u]= float('inf')
		g['BFS']['pi'][u] = None
	g['BFS']['color'][s] = 'GRAY' #Le sommet dont on commence le traitement est grisé
	g['BFS']['d'][s] = 0 #On inititalise la distance parcourue
	g['BFS']['Q'].append(s) #On enfile le noeud dans le Q
	while g['BFS']['Q'] != []:
		u = g['BFS']['Q'].pop(0)# On défile la liste en FIFO ce qui en fait une file
		for v in sorted(g['edges'][u]):#Parcours du graphe de manière ordonnée 
			if g['BFS']['color'][v]=='WHITE':#si un sommet est non traité
				g['BFS']['color'][v]='GRAY'#On le traite
				g['BFS']['d'][v] = g['BFS']['d'][u] + 1#On incrémente la distance
				g['BFS']['pi'][v] = u#Le sommet d'ou l'on viens est ajouté à la liste des prédécesseurs
				g['BFS']['Q'].append(v)# on enfile le sommet v dans la file
		g['BFS']['color'][u] = 'BLACK'# On colore en noir le sommet dont le traitement est terminé
	return(g)


def DFS(g):
	g['DFS']={'color':{},'pred':{},'time': 0,'class':{},'d':{},'f':{}, 'Cyclique' : False, 'tritopo' : [] }
	for u in sorted(g['nodes_index']):#initialisation des sommets
		g['DFS']['color'][u]='WHITE'
		g['DFS']['pred'][u]= None
	g['DFS']['time']= 0
	for u in sorted(g['nodes_index'].keys(), reverse = False, key=lambda t: t[0]):#Tri du dictionnaire selon les vertex
		if g['DFS']['color'][u]=='WHITE':#Si le vertex n'est pas traité,
			DFS_VISIT(g,u)#On le visite
	if g['DFS']['Cyclique']==False: #Si le graphe est acyclique
		for u in sorted(g['DFS']['f'].items(), reverse = True, key=lambda t: t[1]):#On tri selon les valeurs de temps de fin de traitement de manière décroissante
			g['DFS']['tritopo'].append(u)#on ajoute au sous dictionnaire 'tritopo' qui s'avère etre une liste qui sera ordonnée sous la forme du tri topologique de notre graphe apres DFS	
	return(g)
	
def DFS_VISIT(g,u):
	g['DFS']['color'][u] = 'GRAY'
	g['DFS']['time'] = g['DFS']['time'] + 1
	g['DFS']['d'][u] = g['DFS']['time']
	for v in sorted(g['edges'][u].keys(), key=lambda t: t[0]):#Tri des voisins de u selon leur ordre d'indexation
		if g['DFS']['color'][v]=='WHITE':
			g['DFS']['pred'][v]= u#si on est arrivé au sommet v en passant par u,
			DFS_VISIT(g,v)#récursion afin de visiter l'ensemble du graphe
			g['DFS']['class'][(u,v)] = 'TREE EDGE'#Voir compte rendu
		elif g['DFS']['color'][v]=='GRAY':
			g['DFS']['class'][(u,v)] = 'BACK EDGE'#Voir compte rendu
			g['DFS']['Cyclique']=True#Ajout de l'information sur le caractère cyclique du graphe
		elif g['DFS']['d'][u] > g['DFS']['d'][v]:
			g['DFS']['class'][(u,v)] = 'CROSS EDGE'#Voir compte rendu
		elif g['DFS']['d'][v] > g['DFS']['d'][u]:
			g['DFS']['class'][(u,v)] = 'FORWARD EDGE'#Voir compte rendu
	g['DFS']['color'][u]='BLACK'
	g['DFS']['time'] = g['DFS']['time'] + 1
	g['DFS']['f'][u]=g['DFS']['time']
	return(g)


#g = load_TAB('testDFS.tab')
#g = DFS(g)
#pprint(g)	


def initialize_single_source(g, s):
	g['BellFord'] = {'d' : {} , 'pi' : {}}
	for v in g['nodes_index']:
		g['BellFord']['d'][v]= float('inf')
		g['BellFord']['pi'][v] = None
	g['BellFord']['d'][s]=0
	return(g)

def relax(g , u, v):
	if g['BellFord']['d'][v] > g['BellFord']['d'][u] + int(g['edges'][u][v]['weight']):
		g['BellFord']['d'][v] = g['BellFord']['d'][u] + int(g['edges'][u][v]['weight'])
		g['BellFord']['pi'][v] = u 

def bellman_ford(g , s):
	initialize_single_source(g,s)
	for i in range(len(g['nodes']) - 1):
		for u  in g['edges']:
			for v in g['edges'][u]:			
				relax(g,u , v)
	return(g)
	
def adjacency_matrix(g):
	"""
	Retourne la matrice d'adjacence associée au graphe avec 1 si aucune valeur de poids d'arc n'est définis
	
	"""
	matrix_size = int(g['nb_nodes'])
	mat = np.full( (matrix_size, matrix_size), np.inf)#Création d'une matrice carrée rempli d'infini
	for u in sorted(g['nodes_index']):
		for v in sorted(g['edges'][u]):#pour tous les voisins de u.
			if g['edges'][u][v]['weight'] == None:
				mat[g['nodes_index'][u],g['nodes_index'][v]]=1#le point correspondant aux coordonnées des vertex prend la valeur 1 
			else:
				mat[g['nodes_index'][u],g['nodes_index'][v]]=g['edges'][u][v]['weight']#ou il prend la valeur du poids de l'arc
	return(mat)

	
def successeur_matrix(g):
	"""
	Retourne la matrice de successeur de chaques sommets.
	"""
	matrix_size = g['nb_nodes']
	N = np.full( (matrix_size, matrix_size), None)
	adj = adjacency_matrix(g)
	for i in range(g['nb_nodes']):
		for j in range(g['nb_nodes']):
			if adj[i,j]!=np.inf:
				N[i,j]=int(i)
			if i == j:
				N[i,j]= None
	return(N)

#pour lindex des sommets :
#print('index des sommets')
#	for u in g['nodes']:
#		print(u, g['nodes'][u])
def floyd_warshall(g):
	g['Floyd'] = {'D' : adjacency_matrix(g), 'N' : successeur_matrix(g), 'path' : [], 'diameter' : 0 }
	for i in range(len(g['Floyd']['D'])):
		for j in range(len(g['Floyd']['D'])):
			g['Floyd']['D'][i,i]=int(0)
	for k in range(len(g['nodes'])):
		for i in range(len(g['nodes'])):
			for j in range(len(g['nodes'])):
				if g['Floyd']['D'][i,k] + g['Floyd']['D'][k,j] < g['Floyd']['D'][i,j]:
					g['Floyd']['D'][i,j] = g['Floyd']['D'][i,k]+g['Floyd']['D'][k,j]
					g['Floyd']['N'][i,j] = g['Floyd']['N'][k,j]
	lst = []
	for i in range(len(g['Floyd']['D'])):
		lst.append(max(g['Floyd']['D'][i]))
	g['Floyd']['diameter'] = int((max(lst)))
	return(g)


#a modifier
def shortest_path(g,i,j):
	if i != j:
		if g['Floyd']['D'][i,j] == np.inf:
			print("il n'y a pas de chemins")
		else : 
			g['Floyd']['path'].append(i)
			k = int(g['Floyd']['N'][i,j])
			km= -1
			print(k)
			g['Floyd']['path'].append(k)
			while math.isnan(g['Floyd']['N'][k,i])==False and k != km and k != j:			
				k = int(g['Floyd']['N'][k,i])
				g['Floyd']['path'].append(k)
				km = k 
			g['Floyd']['path'].append(j)
		return(g)
	else : 
		return(g)
	
if __name__ == "__main__":	
	
	print('________________________________\n#### create_graph ####')
	print('________________________________\n')
	g = create_graph()
	pprint(g)
 
	print('________________________________\n## Adding node A')
	print('________________________________\n')
	add_node(g,'A')
	pprint(g)
 
	print('________________________________\n## Adding edge A -> B')
	print('________________________________\n')
	add_edge(g, 'A', 'B')
	pprint(g)
 
	print('________________________________\n## #Loading SIF file dressing.sif')
	print('________________________________\n')
	G = load_SIF('dressing.sif')
	pprint(G)
	
	print('________________________________\n	Test BFS : ')
	print('________________________________\n')

	g1 = load_TAB('testBFS.tab')
	g1 = trigraphe(g1)
	g1 = BFS(g1, 's')
	pprint(g1)
	print('________________________________\n	Test DFS : ')
	print('________________________________\n')
	g2 = load_TAB('testDFS.tab')
	g2 = trigraphe(g2)
	g2 = DFS(g2)
	pprint(g2)
	print('________________________________\n	Test Bellman-Ford : ')
	print('________________________________\n')
	g3 = load_TAB('M1BBS_Graphe_Bellman-Ford.tab')
	g3 = trigraphe(g3)
	g3 = bellman_ford(g3, 'A')
	pprint(g3)

	print('________________________________\n	Test Floyd-Warshall : ')
	print('________________________________\n')
	g4 = load_TAB('M1BBS_Graphe_Floyd-Warshall.tab')
	g4=trigraphe(g4)
	g4 = floyd_warshall(g4)
	pprint(g4)








