import json
import re
from NetBiolGraph import NetBiolGraph
from itertools import chain
import argparse
import sys
import os
from time import strftime

def parse_args(args):
    help_text = \
    """
    === SalmoNet 2.0 Metabolic Layer construction script ===

    **Description**

    This script needs a locustag mapping file with a specific uniprot-like format, 
    see uni2stm.tab, and combines it with (.json) metabolic models, such as the ones
    from BIGG to output indirect Enzyme-Enzyme interactions.

    ** Parameters: **

    -i, --input-file <path>          : Path to .json file [mandatory]

    **Exit codes**
    Exit code1: The specified input file doesn't exist!


    """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
        help="<Path to .json file [mandatory]>",
        type=str,
        dest="input_file",
        action="store",
        required=True)
    parser.add_argument("-o", "--output-file",
        help="<path to an output mitab file> [mandatory]",
        type=str,
        dest="output_file",
        action="store",
        required=True)

    results = parser.parse_args(args)

    return results.input_file, results.output_file

def main():
    input_file, output_file = parse_args(sys.argv[1:])
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Parsing model')
    abspath_output_file = os.path.abspath(output_file)

    uniprot_mapping = {}
    #open mapping file
    with open("uni2stm.tab") as f:
        for row in f:
            row = row.strip()
            w = row.split("\t")
            olns = w[3].replace(" ", "")
            for oln in olns.split(";"):
                if oln != "":
                    if oln not in uniprot_mapping:
                        uniprot_mapping[oln] = {"uniprot": [], "name": "", "alternatives": [], "taxid": ""}
                    uniprot_mapping[oln]["uniprot"].append(w[0])
                    names = w[1].replace(" ", "")
                    for name in names.split(";"):
                        if name != "":
                            if uniprot_mapping[oln]["name"] == "":
                                uniprot_mapping[oln]["name"] = name
                            uniprot_mapping[oln]["alternatives"].append(name)
                    uniprot_mapping[oln]["taxid"] = w[4]
    uniprot_node_data = {}                
    
    for u in uniprot_mapping:
        if uniprot_mapping[u]["name"] == "":
            uniprot_mapping[u]["name"] = u
        for ac in uniprot_mapping[u]["uniprot"]:
            uniprot_node_data[ac] = {
                "label": uniprot_mapping[u]["name"],
                "type": "enzyme",
                "aliases": {
                    "uniprot_ac": [ac,],
                    "oln": [u,],
                    "alternative_names": uniprot_mapping[oln]["alternatives"],
                },
                "species": int(uniprot_mapping[oln]["taxid"])
            }

    '''
    Get a list of keys from dictionary which has the given value
    '''
    def getKeysByValue(dictOfElements, valueToFind):
        listOfKeys = []
        listOfItems = dictOfElements.items()
        for item  in listOfItems:
            if item[1] == valueToFind:
                listOfKeys.append(item[0])
        return  listOfKeys

    with open(input_file,'r') as access_json:
        salmonellaJson = json.load(access_json)

    enzym = {}
    metabolit = {}
    F_NET = {}
    F_graph = NetBiolGraph.NetBiolGraph()
    F_nodes = []



    for react in salmonellaJson['reactions']:

        #map yara's IDs to locus tag
        yaraIdmap = {}
        for i in salmonellaJson['genes']:
            yID = i['id']
            try:
                lt = i['notes']['original_bigg_ids'][0]
            except KeyError:
                next
            if yID not in yaraIdmap:
                try:
                    yaraIdmap[i['id']] = i['notes']['original_bigg_ids'][0]
                except KeyError:
                    next

        Mr = []
        Mp = []

        #get reactants
        Mr = getKeysByValue(react['metabolites'], -1)
        
        #get products
        Mp = getKeysByValue(react['metabolites'], -1)
        E1 = []
        
        #get and / or out of the lists
        if react['gene_reaction_rule'] not in E1:
            if E1 != "":
                E1.append(re.split(' or | and ',react['gene_reaction_rule']))
        
        E2 = []
        
        for yaraIDlist in E1:
            for yaraID in yaraIDlist:
                if yaraIdmap.get(yaraID) != '':
                    try:
                        E2.append(yaraIdmap[yaraID]) 
                    except KeyError:
                        next
                

        is_low = False
        is_upp = False
        #unlist the resulting nested list
        #E = list(chain.from_iterable(E2))    
        #print(E)
        E = E2
        
        is_low = False
        is_upp = False

        if react['lower_bound'] < 0.0:
            is_low = True
        if react['upper_bound'] > 0.0:
            is_upp = True
        
        for e in E:
            if e not in enzym:
                enzym[e] = {'in': [], 'out': []}
            if is_upp:
                for m in Mr:
                    if m not in enzym[e]['in']:
                        enzym[e]['in'].append(m)
                for m in Mp:
                    if m not in enzym[e]['out']:
                        enzym[e]['out'].append(m)
            if is_low:
                for m in Mr:
                    if m not in enzym[e]['out']:
                        enzym[e]['out'].append(m)
                for m in Mp:
                    if m not in enzym[e]['in']:
                        enzym[e]['in'].append(m)

        if is_upp:
            for m in Mr:
                if m not in metabolit:
                    metabolit[m] = {'in': [], 'out': []}
                for e in E:
                    if e not in metabolit[m]:
                        metabolit[m]['out'].append(e)
            for m in Mp:
                if m not in metabolit:
                    metabolit[m] = {'in': [], 'out': []}
                for e in E:
                    if e not in metabolit[m]:
                        metabolit[m]['in'].append(e)
        if is_low:
            for m in Mr:
                if m not in metabolit:
                    metabolit[m] = {'in': [], 'out': []}
                for e in E:
                    if e not in metabolit[m]:
                        metabolit[m]['in'].append(e)
            for m in Mp:
                if m not in metabolit:
                    metabolit[m] = {'in': [], 'out': []}
                for e in E:
                    if e not in metabolit[m]:
                        metabolit[m]['out'].append(e)
        if len(E) > 1:
            eu = []
            for e in E:
                if e in uniprot_mapping:
                    for u in uniprot_mapping[e]["uniprot"]:
                        eu.append(u)
            if len(eu) > 1:
                for e1 in eu:
                    for e2 in eu:
                        if e1 != e2:
                            edge = "%s;%s" % (e1, e2)
                            edge_r = "%s;%s" % (e2, e1)
                            mo = "co-catalysis"
                            if edge in F_NET:
                                if mo not in F_NET[edge]:
                                    F_NET[edge].append(mo)
                                    F_graph.edit_edge(e1, e2, common_metabolit=[mo,])
                            elif edge_r in F_NET:
                                if mo not in F_NET[edge_r]:
                                    F_NET[edge_r].append(mo)
                                    F_graph.edit_edge(e2, e1, common_metabolite=[mo,])
                            else:
                                F_NET[edge_r] = []
                                F_NET[edge_r].append(mo)
                                F_graph.add_edge(e1, e2)
                                F_graph.edit_edge(e1, e2, xrefs={"db": ["STM_v1.0",], "pubmed": []})
                                F_graph.edit_edge(e1, e2, common_metabolite=[mo,])
                                if e1 not in F_nodes:
                                    F_nodes.append(e1)
                                    n = F_graph.n_node(e1, **uniprot_node_data[e1])
                                if e2 not in F_nodes:
                                    F_nodes.append(e2)
                                    n = F_graph.n_node(e2, **uniprot_node_data[e2])

    #
    except_m = []
    NET = {}
    graph = NetBiolGraph.NetBiolGraph()
    nodes = []

    for e in enzym:
        for mo in enzym[e]['out']:
            if mo not in except_m:
                for ee in metabolit[mo]['out']:
                    if e != ee:
                        if e in uniprot_mapping and ee in uniprot_mapping:
                            for u in uniprot_mapping[e]["uniprot"]:
                                for uu in uniprot_mapping[ee]["uniprot"]:
                                    edge = "%s;%s" % (u, uu)
                                    edge_r = "%s;%s" % (uu, u)
                                    if edge in NET:
                                        if mo not in NET[edge]:
                                            NET[edge].append(mo)
                                            graph.edit_edge(u, uu, common_metabolit=[mo,])
                                    elif edge_r in NET:
                                        if mo not in NET[edge_r]:
                                            NET[edge_r].append(mo)
                                            graph.edit_edge(uu, u, common_metabolite=[mo,])
                                    else:
                                        NET[edge_r] = []
                                        NET[edge_r].append(mo)
                                        graph.add_edge(u, uu)
                                        graph.edit_edge(u, uu, xrefs={"db": ["STM_v1.0",], "pubmed": []})
                                        graph.edit_edge(u, uu, common_metabolite=[mo,])
                                        if u not in nodes:
                                            nodes.append(u)
                                            n = graph.n_node(u, **uniprot_node_data[u])
                                        if uu not in nodes:
                                            nodes.append(uu)
                                            n = graph.n_node(uu, **uniprot_node_data[uu])

    for e in enzym:
        for mo in enzym[e]['in']:
            if mo not in except_m:
                for ee in metabolit[mo]['in']:
                    if e != ee:
                        if e in uniprot_mapping and ee in uniprot_mapping:
                            for u in uniprot_mapping[e]["uniprot"]:
                                for uu in uniprot_mapping[ee]["uniprot"]:
                                    edge = "%s;%s" % (u, uu)
                                    edge_r = "%s;%s" % (uu, u)
                                    if edge in NET:
                                        if mo not in NET[edge]:
                                            NET[edge].append(mo)
                                            graph.edit_edge(u, uu, common_metabolit=[mo,])
                                    elif edge_r in NET:
                                        if mo not in NET[edge_r]:
                                            NET[edge_r].append(mo)
                                            graph.edit_edge(uu, u, common_metabolite=[mo,])
                                    else:
                                        NET[edge_r] = []
                                        NET[edge_r].append(mo)
                                        graph.add_edge(u, uu)
                                        graph.edit_edge(u, uu, xrefs={"db": ["STM_v1.0",], "pubmed": []})
                                        graph.edit_edge(u, uu, common_metabolite=[mo,])
                                        if u not in nodes:
                                            nodes.append(u)
                                            n = graph.n_node(u, uniprot_node_data[u])
                                        if uu not in nodes:
                                            nodes.append(uu)
                                            n = graph.n_node(uu, uniprot_node_data[uu])

    ##########################################################
    co_metabolite = {}
    for edge in NET:
        for m in NET[edge]:
            if m not in co_metabolite:
                co_metabolite[m] = 0
            co_metabolite[m] += 1
    for m in co_metabolite:
        if co_metabolite[m] > 10:
            except_m.append(m)
    ##########################################################

    for e in enzym:
        for mo in enzym[e]['out']:
            if mo not in except_m:
                for ee in metabolit[mo]['out']:
                    if e != ee:
                        if e in uniprot_mapping and ee in uniprot_mapping:
                            for u in uniprot_mapping[e]["uniprot"]:
                                for uu in uniprot_mapping[ee]["uniprot"]:
                                    edge = "%s;%s" % (u, uu)
                                    edge_r = "%s;%s" % (uu, u)
                                    if edge in F_NET:
                                        if mo not in F_NET[edge]:
                                            F_NET[edge].append(mo)
                                            F_graph.edit_edge(u, uu, common_metabolit=[mo,])
                                    elif edge_r in F_NET:
                                        if mo not in F_NET[edge_r]:
                                            F_NET[edge_r].append(mo)
                                            F_graph.edit_edge(uu, u, common_metabolite=[mo,])
                                    else:
                                        F_NET[edge_r] = []
                                        F_NET[edge_r].append(mo)
                                        F_graph.add_edge(u, uu)
                                        F_graph.edit_edge(u, uu, xrefs={"db": ["STM_v1.0",], "pubmed": []})
                                        F_graph.edit_edge(u, uu, common_metabolite=[mo,])
                                        if u not in F_nodes:
                                            F_nodes.append(u)
                                            n = F_graph.n_node(u, **uniprot_node_data[u])
                                        if uu not in F_nodes:
                                            F_nodes.append(uu)
                                            n = F_graph.n_node(uu, **uniprot_node_data[uu])

    for e in enzym:
        for mo in enzym[e]['in']:
            if mo not in except_m:
                for ee in metabolit[mo]['in']:
                    if e != ee:
                        if e in uniprot_mapping and ee in uniprot_mapping:
                            for u in uniprot_mapping[e]["uniprot"]:
                                for uu in uniprot_mapping[ee]["uniprot"]:
                                    edge = "%s;%s" % (u, uu)
                                    edge_r = "%s;%s" % (uu, u)
                                    if edge in F_NET:
                                        if mo not in F_NET[edge]:
                                            F_NET[edge].append(mo)
                                            F_graph.edit_edge(u, uu, common_metabolit=[mo,])
                                    elif edge_r in F_NET:
                                        if mo not in F_NET[edge_r]:
                                            F_NET[edge_r].append(mo)
                                            F_graph.edit_edge(uu, u, common_metabolite=[mo,])
                                    else:
                                        F_NET[edge_r] = []
                                        F_NET[edge_r].append(mo)
                                        F_graph.add_edge(u, uu)
                                        F_graph.edit_edge(u, uu, xrefs={"db": ["STM_v1.0",], "pubmed": []})
                                        F_graph.edit_edge(u, uu, common_metabolite=[mo,])
                                        if u not in F_nodes:
                                            F_nodes.append(u)
                                            n = F_graph.n_node(u, uniprot_node_data[u])
                                        if uu not in F_nodes:
                                            F_nodes.append(uu)
                                            n = F_graph.n_node(uu, uniprot_node_data[uu])

    with open("%s.csv" % output_file, "w") as f:
        print(f'MESSAGE [{strftime("%H:%M:%S")}]: Writing results to the output files: {abspath_output_file}')
        for edge in F_NET:
            f.write("%s;%s;%s;enzym-enzym interaction\n" % (edge, ",".join(F_NET[edge]), len(F_NET[edge])))

    #with open("Metabolite_%s.csv" % output_file, "w") as f:
    #    co_metabolite = {}
    #    for edge in F_NET:
    #        for m in F_NET[edge]:
    #            if m not in co_metabolite:
    #                co_metabolite[m] = 0
    #            co_metabolite[m] += 1
    #    for m in co_metabolite:
    #        f.write("%s;%s\n" % (m, co_metabolite[m]) )


if __name__ == '__main__':
    main()