#!/usr/bin/python3
# -*- coding: utf-8 -*-

import TPGraphlib as gr # Graph library from part 1 of the project
import math
from pprint import pprint


def load_OBO(filename):
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""

	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'): return
		# otherwise create node
		id = lines.pop(0)[4:].rstrip()
		term = gr.add_node(g,id)
		term['id'] = id
		term['type'] = 'GOTerm'
		for line in lines:
			# attributes (name, namespace, def)
			if line.startswith('name: '): term['name'] = line[6:]
			elif line.startswith('namespace: '): term['namespace'] = line[11:]
			elif line.startswith('def: '): term['def'] = line[5:]
			elif line.startswith('alt_id: '): g['alt_id'][ line[8:] ] = id # alternate ids
			# relationships
			elif line.startswith('is_a:'): # is_a
				parent = line[6:line.index('!')].rstrip()
				e = gr.add_edge(g,id, parent)
				e['type'] = 'is_a'
			elif line.startswith('relationship: part_of '): # part_of
				line = line[line.index('GO:'):]
				dest = line[:line.index(' ')]
				e = gr.add_edge(g,id, dest)
				e['type'] = 'part_of'
	#
	g=gr.create_graph(directed=True, weighted=False)
	g['alt_id'] = {} # alternate GO ids
	with open(filename) as f:
		line = f.readline().rstrip()
		# skip header to reach 1st Term
		while not line.startswith('[Term]'):
			line = f.readline().rstrip()
		buff = []
		line = f.readline()
		stop = False
		while line and not stop:
			# buffer lines until the next Term is found
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				# next Term found: create corresponding node and edges in parseTerm and empty buffer
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	return g

def load_GOA(go, filename):
	"""
	parse GOA file and add annotated gene products to previsouly loaded graph go

	Extract of a file to be parsed:
	!gaf-version: 2.1
	!GO-version: http://purl.obolibrary.org/obo/go/releases/2016-10-29/go.owl
	UniProtKB  A5A605  ykfM      GO:0006974  PMID:20128927   IMP              P  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20100901  EcoCyc
	UniProtKB  A5A605  ykfM      GO:0016020  GO_REF:0000037  IEA              C  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20161029  UniProt
	UniProtKB  P00448  sodA      GO:0004784  GO_REF:0000003  IEA  EC:1.15.1.1 F  Superoxide dismutase [Mn]       SODM_ECOLI|sodA|JW3879|b3908  protein taxon:83333  20161029  UniProt
	UniProtKB  P00393  ndh  NOT  GO:0005737  PMID:6784762    IDA              C  NADH dehydrogenase              DHNA_ECOLI|ndh|JW1095|b1109   protein taxon:83333  20100621  EcoliWiki
	    0        1       2   3       4             5          6        7      8             9                              10
	             id    name        go_id               evidence-codes                     desc                           aliases
	"""

	names = {}
	go['names'] = names # gene names or gene product names (column 3)
	with open(filename) as f:
		line = f.readline()
		while line:
			if not line.startswith('!'):
				cols = line.rstrip().split('\t')
				id = cols[1]
				go_id = cols[4]
				if go_id not in go['nodes']: # GOTerm not found search alternate ids
					if go_id in go['alt_id']: # success
						go_id = go['alt_id'][go_id] # replace term
					else: # warn user
						print('Warning: could not attach a gene product (%s) to a non existing GO Term (%s)' % (id, go_id))
				if go_id in go['nodes']:
					# create node for gene product if not already present
					if id not in go['nodes']:
						g = gr.add_node(go,id)
						g['id'] = id
						g['type'] = 'GeneProduct'
						names[cols[2]] = id
					# create or update gene product attributes
					gp = go['nodes'][id]
					gp['name'] = cols[2]
					gp['desc'] = cols[9]
					gp['aliases'] = cols[10]
					# attach gene product to GOTerm
					go_term = go['nodes'][go_id]
					e = gr.add_edge(go, id, go_id)
					e['type'] = 'annotation'
					if 'evidence-codes' not in e: e['evidence-codes'] = []
					e['evidence-codes'].append( cols[6] )
				else: # go_id or alt_id not found in GOTerms
					print('Error: could not attach a gene product (%s) to non existing GO Term (%s)' % (id, go_id))
			line = f.readline()

def load_GOAMD(go, filename):
	"""
	parse GOA file and add annotated gene products to previsouly loaded graph go

	Extract of a file to be parsed:
	!gaf-version: 2.1
	!GO-version: http://purl.obolibrary.org/obo/go/releases/2016-10-29/go.owl
	UniProtKB  A5A605  ykfM      GO:0006974  PMID:20128927   IMP              P  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20100901  EcoCyc
	UniProtKB  A5A605  ykfM      GO:0016020  GO_REF:0000037  IEA              C  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20161029  UniProt
	UniProtKB  P00448  sodA      GO:0004784  GO_REF:0000003  IEA  EC:1.15.1.1 F  Superoxide dismutase [Mn]       SODM_ECOLI|sodA|JW3879|b3908  protein taxon:83333  20161029  UniProt
	UniProtKB  P00393  ndh  NOT  GO:0005737  PMID:6784762    IDA              C  NADH dehydrogenase              DHNA_ECOLI|ndh|JW1095|b1109   protein taxon:83333  20100621  EcoliWiki
	    0        1       2   3       4             5          6        7      8             9                              10
	             id    name        go_id               evidence-codes                     desc                           aliases
	"""
	nbgp = 0
	erreurgp = 0
	names = {}
	go['names'] = names # gene names or gene product names (column 3)
	with open(filename) as f:
		line = f.readline()
		while line:
			if not line.startswith('!'):
				cols = line.rstrip().split('\t')
				id = cols[1]
				go_id = cols[4]
				if go_id not in go['nodes']: # GOTerm not found search alternate ids
					if go_id in go['alt_id']: # success
						go_id = go['alt_id'][go_id] # replace term
					else: # warn user
						print('Warning: could not attach a gene product (%s) to a non existing GO Term (%s)' % (id, go_id))
						erreurgp = erreurgp + 1
				if go_id in go['nodes']:
					# create node for gene product if not already present
					if id not in go['nodes']:
						nbgp = nbgp + 1
						g = gr.add_node(go,id)
						g['id'] = id
						g['type'] = 'GeneProduct'
						names[cols[2]] = id
					# create or update gene product attributes
					gp = go['nodes'][id]
					gp['name'] = cols[2]
					gp['desc'] = cols[9]
					gp['aliases'] = cols[10]
					# attach gene product to GOTerm
					go_term = go['nodes'][go_id]
					e = gr.add_edge(go, go_id, id)


					e['type'] = 'annotation'
					if 'evidence-codes' not in e: e['evidence-codes'] = []
					e['evidence-codes'].append( cols[6] )
				else: # go_id or alt_id not found in GOTerms
					print('Error: could not attach a gene product (%s) to non existing GO Term (%s)' % (id, go_id))
					erreurgp = erreurgp + 1
			line = f.readline()
	print(nbgp, erreurgp)

def load_OBOMD(filename):#Charge la GeneOntology avec les arc inversés (utile pour la fonction Max_Depth)
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""

	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'): return
		# otherwise create node
		id = lines.pop(0)[4:].rstrip()
		term = gr.add_node(g,id)
		term['id'] = id
		term['type'] = 'GOTerm'
		for line in lines:
			# attributes (name, namespace, def)
			if line.startswith('name: '): term['name'] = line[6:]
			elif line.startswith('namespace: '): term['namespace'] = line[11:]
			elif line.startswith('def: '): term['def'] = line[5:]
			elif line.startswith('alt_id: '): g['alt_id'][ line[8:] ] = id # alternate ids
			# relationships
			elif line.startswith('is_a:'): # is_a
				parent = line[6:line.index('!')].rstrip()
				e = gr.add_edge(g,parent, id)#ordre inversé
				e['type'] = 'is_a'
			elif line.startswith('relationship: part_of '): # part_of
				line = line[line.index('GO:'):]
				dest = line[:line.index(' ')]
				e = gr.add_edge(g,dest ,id)#ordre inversé
				e['type'] = 'part_of'
	#
	nbgo = 0
	g=gr.create_graph(directed=True, weighted=False)
	g['alt_id'] = {} # alternate GO ids
	with open(filename) as f:
		line = f.readline().rstrip()
		# skip header to reach 1st Term
		while not line.startswith('[Term]'):
			line = f.readline().rstrip()
		buff = []
		line = f.readline()
		stop = False
		while line and not stop:
			# buffer lines until the next Term is found
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				nbgo = nbgo + 1
				# next Term found: create corresponding node and edges in parseTerm and empty buffer
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	print(nbgo)
	return g

def allgp(go):

	lst = []
	for u in go['nodes']:
		if go['nodes'][u]['type']=='GeneProduct':
			lst.append(go['nodes'][u]['id'])
	return lst

def allgo(go):

	lst = []
	for u in go['nodes']:
		if go['nodes'][u]['type']=='GOTerm':
			lst.append(go['nodes'][u]['id'])
	return lst

def GOTerms(go, gp, goterm, all=True, evidence_code=None):
	"""
	return the GOTerms associated to the provided gene product (gp)

	go: Gene Ontology graph
	gp: gene product
	goterm: dictionnaire des résultats
	all: if True, all the GOTerms and their ancestors will be return, otherwise only the GOTerms directly associated to the gene product will be returned.
	evidence_code: ignored for the moment

	Returns a list of GOTerms identifiers, e.g. ['GO:0005215','GO:0005515','GO:0006810','GO:0006974','GO:0008643']
	"""
	#goterm = { 'nodes' : {} }

	if gp in go['nodes']:
		for u in go['edges'][gp]:#Pour tous les voisins du vertex choisi
			if go['nodes'][u]['type'] == 'GOTerm':
				if u not in goterm:#Evsite les doublons
					goterm[u]=go['nodes'][u]
			if all:
				for v in go['edges'][u]:
					if v not in goterm:
						goterm[v] = GOTerms(go , v, goterm, all = True)#Entre dans une boucle récursive jusqu'a épuisement des voisins
	else :
		print('Numero inexistant')
	return (goterm)


def recherche_max(g, gp):
	"""
	Recherche du chemin le plus mong d'un poin à tous les autres points du graphe, la fonction s'inspire du BFS mais choisira le chemin le plus long de tous les chemins accessibles via le sommet entré en argument. Explications supplémentaires sur le compte rendu.
	"""

	g['max'] = {'d':{},'Q':[], 'lst' : [] }
	for u in sorted(g['nodes_index']):
		g['max']['d'][u]= 0
	g['max']['d'][gp] = 0
	g['max']['Q'].append(gp)
	while g['max']['Q'] != []:
		i = g['max']['Q'].pop(0)#Mise en place d'une file.
		for v in sorted(g['edges'][i]):
			g['max']['d'][v] = max((g['max']['d'][i] + 1), g['max']['d'][v])#choix du max entre la valeur du sommet et celle atteignable en ajoutant 1 au sommet parent
			g['max']['Q'].append(v)#Engendre une récursivité jusqu'a épuisement des noeuds atteignable
			g['max']['lst'].append(g['max']['d'][v])
	maxd= max(g['max']['lst'])#Retour de la distance du chemin maximal
	return(maxd)

def puits(go):
	"""
	Détermination des puits du graphe, à savoir les racines des trois GeneOntology : Biological Precess, Molecular Fuction, Cellullar Component.
	"""

	go['MaxDepth'] = { 'lst' : {}}
	for u in go['nodes']:
		if not(go['edges'][u]):# Le condition pour être un pit est de ne pas avoir d'arc sortant (pour le graph chargé normalement)
			go['MaxDepth']['lst'][u] = go['nodes'][u]
	return go

def Max_Depth(go, goinv):
	"""
	Détermine les profondeurs des trois GeneOntology. Part des puits puis pour le graphe inversé, applique la recherche de la distance maximale.
	"""

	for u in go['MaxDepth']['lst']:
		if go['nodes'][u]['namespace']=='biological_process':
			print('Profondeur max biological process')
			print(recherche_max(goinv, u))
		if go['nodes'][u]['namespace']=='cellular_component':
			print('Profondeur max cellular component')
			print(recherche_max(goinv, u))
		if go['nodes'][u]['namespace']=='molecular_function':
			print('Profondeur max molecular function')
			print(recherche_max(goinv, u))
	return(go)

def GeneProducts(go, term, geneproduct, all=True, evidence_code=None):
	"""
	return the gene products anotated by the provided GOTerm

	go: Gene Ontology graph
	term: GOTerm id
	all: if True, all the gene products directly and undirectly annotated (linked to a descendant of GOTerm) will be return, otherwise only the gene products directly associated to the GOTerm will be returned.
	evidence_code: ignored for the moment

	Returns a list of gene products identifiers, e.g. ['P0AAG5', 'P0AFY6', 'P10907', 'P16676', 'P23886']
	"""

	if all:
		for u in go['edges']:
			if go['nodes'][u]['type'] == 'GeneProduct':
				goterm = {}
				got = GOTerms(go , u, goterm, all = True)#Recherche dans tous les arc reliant un GeneProduct à des GOTerm, si le GOTerm choisi est présent. Recherche exhaustive pour les descendants
				if term in got:
					geneproduct[u] = go['nodes'][u]#Renvoi le GeneProduct qui est lié au GOTerm
	else:
		for u in go['edges']:
			if go['nodes'][u]['type'] == 'GeneProduct':
				goterm = {}
				got = GOTerms(go , u, goterm, all = False)#Recherche dans tous les arc reliant un GeneProduct à des GOTerm, si le GOTerm choisi est présent. Recherche restreint aux descendants directs
				if term in got:
					geneproduct[u] = go['nodes'][u]
	return (geneproduct)


##### lib tests #####
if __name__ == "__main__":

	print('GeneOntology lib tests')
	print('________________________________\n	Chargement des graphes pour les tests : ')
	print('________________________________\n')
	go = load_OBO('go-basic.obo')
	load_GOA(go , 'H_influenzae_ATCC_51907.goa')
	go=gr.trigraphe(go)
	goinv=load_OBOMD('go-basic.obo')
	goinv = gr.trigraphe(goinv)
	print()
	choix=str(input('Voulez vous visualiser le graphe ? O/N'))
	if choix == 'O' or choix == 'o':
		pprint(go)
	print('________________________________\n	Test de la fonction GOTerm : ')
	print('________________________________\n')
	goterm = {}
	print('Recherche resterinte : \n')
	pprint(GOTerms(go, 'P45092', goterm, all = False))
	print('Recherche exhaustive : \n')
	pprint(GOTerms(go, 'P45092', goterm, all = True))
	print('________________________________\n	Test de la fonction GeneProduct : ')
	print('________________________________\n')
	geneproduct = {}
	print('Recherche resterinte : \n')
	pprint(GeneProducts(go , 'GO:0003333' , geneproduct, all = False))
	print('Recherche exhaustive : \n')
	pprint(GeneProducts(go , 'GO:0003333' , geneproduct, all = True))
	print('________________________________\n	Test de la fonction Max_Depth : ')
	print('________________________________\n')
	goinv=load_OBOMD('go-basic.obo')
	goinv = gr.trigraphe(goinv)
	go=puits(go)
	print('Sources de la Gene Ontology')
	pprint(go['MaxDepth'])
	Max_Depth(go, goinv)



