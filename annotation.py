#!/usr/bin/python
# -*- coding: utf-8 -*-

###Berlin Tristan
###Programme permettant d'obtenir les fichier sets à partir de la
###Gene Ontology et d'un fichier GOA lié

import GeneOntology as gr
import datetime

spc = input("Veuillez entrer le nom du fichier :")
num = input ("Quelle annotation? Direct[1] ou Implicit[2] :")

###Récupérer la liste des gènes à la manière de la fonction GeneProducts de l'année dernière
def recupGene(Got, GP, all = True):
    if all:
        for t in OBO["edges"][Got]:
            if OBO["nodes"][t]["type"] == "GeneProduct" and (OBO["nodes"][t]["id"]) not in GP:
                GP.append((OBO["nodes"][t]["id"]))
            elif OBO["nodes"][t]["type"] == "GOTerm":
                recupGene(OBO["nodes"][t]["id"], GP)
    else:
        for t in OBO["edges"][Got]:
            if OBO["nodes"][t]["type"] == "GeneProduct" and (OBO["nodes"][t]["id"]) not in GP:
                GP.append((OBO["nodes"][t]["id"]))
    return(GP)

if int(num) == 1 or str(num) == "Direct":
    ##Création du fichier sets
    nomFichier = str(spc) + "_direct.sets"
    fichier = open(nomFichier, "x")
    ##Chargement de la Gene Ontology et du fichier GOA
    OBO = gr.load_OBOMD("go.obo")
    gr.load_GOAMD(OBO,spc)
    exit
    fichier.write("# format: sets")
    fichier.write("\n# version: 1.0")
    fichier.write("\n# strain: "+ str(spc))
    fichier.write("\n# date: "+ str(datetime.datetime.now()))
    fichier.write("\n# comment: Gene Ontology terms direct")
    ###Pour les directs
    for i in OBO["nodes"]:
        if OBO["nodes"][i]["type"] == "GOTerm":
            listeGP = []
            Got = OBO["nodes"][i]["id"]
            Desc = OBO["nodes"][i]["name"]
            recupGene(Got, listeGP, all = False)
            if len(listeGP) != 0:
                listeGP = "\t".join(listeGP)
                fichier.write("\n" +str(Got) + "\t" +  str(Desc) + "\t" + str(listeGP))
    fichier.close()


###Pour les indirects
elif int(num) == 2 or str(num) == "Implicit":
    ##Création du fichier
    nomFichier = str(spc) + "_implicit.sets"
    fichier = open(nomFichier, "x")
    #Chargement de la Gene Ontology et du fichier GOA
    OBO = gr.load_OBOMD("go.obo")
    gr.load_GOAMD(OBO, spc)
    fichier.write("# format: sets")
    fichier.write("\n# version: 1.0")
    fichier.write("\n# strain: "+ str(spc))
    fichier.write("\n# date: "+ str(datetime.datetime.now()))
    fichier.write("\n# comment: Gene Ontology terms direct")
    for i in OBO["nodes"]:
        if OBO["nodes"][i]["type"] == "GOTerm":
           # GP = gr.GeneProducts(OBO, OBO["nodes"][i]["id"])
           Got = OBO["nodes"][i]["id"]
           Desc = OBO["nodes"][i]["name"]
           GP = []
           recupGene(Got, GP)
           if len(GP) != 0:
               GP = "\t".join(GP)
               fichier.write("\n" +str(Got) + "\t" +  str(Desc) + "\t" + str(GP))
    fichier.close()
