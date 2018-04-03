# ==============================
# Athout: Eugene Seo
# Date: 10.21.2016
# Description: CS519 Homework 1
# ==============================
import pandas
hs = pandas.read_csv("hsmetnet.txt", sep="\t", names=["v_left", "v_right"])

from igraph import *
from igraph import summary
meta_graph = Graph.TupleList(hs.values.tolist(), directed=True)
summary(meta_graph)

from collections import defaultdict
metabolite_set = set()
reaction_set = set()
metabolite_degree = defaultdict(int)
metabolite_idx = []
for v in meta_graph.vs:
	if 'REACTION' in v['name']:
		reaction_set.add(v)
		continue
	else:
		metabolite_degree[v['name']] = v.degree()
		metabolite_idx.append(v.index)
		metabolite_set.add(v)

print "A. number of distinct metabolities:", len(metabolite_set)
print "A. number of distinct reactions:", len(reaction_set)
print "A. number of edges:", hs.shape[0]

import numpy as np
from operator import itemgetter
print sorted(metabolite_degree.items(), key=itemgetter(1), reverse=True)[0:6]

import matplotlib.pyplot as plt
degree, frac_metabolites = zip(*[(left, count) for left, _, count in meta_graph.degree_distribution(1, vertices=metabolite_idx).bins()])
plt.loglog(degree, frac_metabolites)
plt.xlabel("degree")
plt.ylabel("fraction of metabolites")
plt.show()

print power_law_fit(np.log(frac_metabolites)).alpha

shortest_paths = meta_graph.shortest_paths(source=metabolite_idx, target=metabolite_idx, mode=ALL)
inf_idx = np.where(np.isinf(shortest_paths[0])==True)
print shortest_paths[0][inf_idx[0][0]]

metabolite_cluster = meta_graph.clusters(mode='weak')
giant_cluster = metabolite_cluster.giant()
giant_memnbers = []
for v in giant_cluster.vs:
	if 'REACTION' in v['name']:
		continue
	else:
		giant_memnbers.append(v.index)

short_path_giant = np.array(giant_cluster.shortest_paths(source=giant_memnbers, target=giant_memnbers, mode=ALL))
short_path_giant = short_path_giant[np.isfinite(short_path_giant)]
print np.average(short_path_giant)

print np.amax(short_path_giant)

betweeness = np.array(meta_graph.betweenness(vertices=metabolite_idx, directed=True))
degree = np.array(meta_graph.degree(vertices=metabolite_idx))
inds_keep = np.where(betweeness > 0)
ax = plt.gca()
ax.scatter(degree[inds_keep], betweeness[inds_keep])
ax.loglog()
plt.xlabel("degree")
plt.ylabel("betweeness centrality")
plt.show()

metabolites_degree_two = np.where(degree == 2)[0]
max_arg = np.argmax(betweeness[metabolites_degree_two])
metabolites_idx_max_betweenness = metabolite_idx[metabolites_degree_two[max_arg]]
print meta_graph.vs[metabolites_idx_max_betweenness]["name"]