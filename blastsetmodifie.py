#!/usr/bin/env python
# Copyright 2018 BARRIOT Roland
# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import numpy as np
import random as rd
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency

# SCRIPT PARAMETERS
# e.g. ./blastset.py --sets EcolA.biocyc.sets --query 'ALAS ARGS ASNS ASPS CYSS GLTX GLYQ GLYS HISS ILES'
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')
parser.add_argument('-t', '--sets', required=True, help='Target sets (categories).')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
parser.add_argument('-c', '--adjust', required=False, action="store_true", help='Adjust for multiple testing (FDR).')
parser.add_argument('-m', '--measure', required=False, default='binomial', help='Dissimilarity index: binomial (default), hypergeometric, chi2 or coverage')
parser.add_argument('-l', '--limit', required=False, type=int, default=0, help='Maximum number of results to report.')
parser.add_argument('-r', '--random', required=False, type=int, default=10, help='Number of randomized genes. Doesn''t work if there is already a query. 10 by default' )
parser.add_argument('-R', '--round', required=False, type=int, default=1, help='How much do you want to run the program. Not necessary with a query')
param = parser.parse_args()

class ComparedSet(object):
    def __init__(self, id, name = '', common = 0, size = 0, pvalue = 1, elements = [], c = []):
        self.id = id
        self.name = name
        self.common = common
        self.size = size
        self.pvalue = pvalue
        self.elements = elements
        self.c = c

# LOAD REFERENCE SETS
def load_sets(filename):
    sets = {}
    ids = set()
    listgene = []
    with open( filename ) as f:
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            words = l.split('\t')
            if len(words) > 2 and not words[0].startswith('#'):
                id = words.pop(0)
                name = words.pop(0)
                words = set(words)
                sets[ id ] = { 'name': name, 'elements': words}
                for i in words:
                    listgene.append(i)
                ids |= words
    return [ sets, len( ids ), listgene ]
(sets, population_size, listgene) = load_sets(param.sets)
j = 0
while j < param.round:
    # LOAD QUERY
    text = param.query
    query = set()
    if isfile(text):
        with open(text) as f:
            content = f.read()
            lines = content.split('\n')
            for l in lines:
                if l!='':
                    query |= set(l.split())
    elif text == "random":
        i =0
        hasard = ""
        while i < param.random:
            ah = rd.choice(listgene)
            if ah not in query:
                hasard = hasard + " " +ah
                i = i+1
        query |= set(hasard.split())
    else: # parse string
        query |= set(text.split())

    # EVALUATE SETS
    results = []
    q = len(query)
    g = population_size
    for id in sets:
        elements = sets[ id ][ 'elements' ]
        t = len(elements)
        common_elements = elements.intersection( query )
        c = len(common_elements)
        if param.measure=='binomial': # binom.cdf(>=success, attempts, proba)
            # p_success = 384/2064, 152 attempts, 61 success
            #~ pval = binom.pmf(61, 152, 384.0/2064)
            pval = binom.cdf( q - c, q, 1 - float(t)/g)
        elif param.measure=='hypergeometric': # hypergeom.sf(common-1, population, target, query) = 1-p( X <= x-1 ) = p( X >= x )
            pval = hypergeom.sf(c-1, g, t, q)
        elif param.measure=='chi2': #chi2_contingency
            chi2, pval, a, e = chi2_contingency(np.array([[c,q-c,q],[t-c,g-q-t+c,g-q],[t,g-t,g]]))
        elif param.measure=='coverage': #Coverage
            if c == 0:
                pval = 1
            elif c==q and c==t:
                pval = 0
            else:
                pval = 1 - ((c / q) * (c / t))
        else:
            print('sorry,-m  %s not (yet) implemented' % ( param.measure ))
            exit(1)
        r = ComparedSet( id, sets[id]['name'], c, t, pval, elements, common_elements)
        results.append( r )

    # PRINT SIGNIFICANT RESULTS
    results.sort(key=lambda an_item: an_item.pvalue)
    i=1
    for r in results:
        # FDR
        if param.adjust and r.pvalue > param.alpha * i / len(results): break
        # limited output
        if param.limit > 0 and i>param.limit: break
        # alpha threshold
        elif r.pvalue > param.alpha : break
        # OUTPUT
        print("%s\t%s\t%s/%s\t%s\t%s" % ( r.id, r.pvalue, r.common, r.size, r.name, ', '.join(r.c)))
        i+=1
    j += 1
