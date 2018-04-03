# =================================
# Athout: Eugene Seo
# Date: 12.02.2016
# Description: CS519 Final Project
# =================================
import pandas
import numpy as np
import matplotlib.pyplot as plt
from igraph import *
from scipy.spatial.distance import *
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import fcluster
from scipy.stats import entropy
import pdb
from collections import OrderedDict

def get_color(color):
    '''
    b---blue   c---cyan  g---green    k----black
    m---magenta r---red  w---white    y----yellow'''
    color_set = ['g','r','c','m', 'b','k','y']
    return color_set[color % 7]

def plot_inoutdegree(top_max, cluster_member_num,
								indegreeByCluster, outdegreeByCluster, network_name):
	plt.figure(1, figsize=(10,5))
	ax = plt.gca()
	color = 0
	j = 1
	for i in range(top_max):
		if cluster_member_num[i] > 0:
			ax.scatter(indegreeByCluster[i], outdegreeByCluster[i],
			color=get_color(color), label="Cluster %d" %(j))#, alpha=.5)
			color += 1
			j += 1
	plt.title("clusters of %s gene regulatory network" %(network_name))
	plt.xlabel("indegree")
	plt.ylabel("outdegree")
	plt.legend(loc='upper right')


def plot_dendrogram(Z, threshold, linkage, network_name):
	plt.figure(3, figsize=(10,5))
	plt.title("Dendrogram on %s gene regulatory network" %(network_name))

	dn = hierarchy.dendrogram(Z, color_threshold=threshold)#, p=12, truncate_mode='lastp' , labels=np.array(vertex_name))
	'''
	print "TYPE!!:", type(dn['color_list'])
	a= np.array(dn[
	'color_list'])
	print "a!!:", a
	print 'g', len(np.where(a=='g')[0])
	print 'r', len(np.where(a=='r')[0])
	print 'b', len(np.where(a=='b')[0])
	print 'c', len(np.where(a=='c')[0])
	print 'm', len(np.where(a=='m')[0])
	print 'y', len(np.where(a=='y')[0])
	print 'k', len(np.where(a=='k')[0])
	'''
	plt.axhline(y=threshold, c='k')


def plot_heatmap(cluster_matrix, cluster_num):
	fig, ax = plt.subplots()
	im = ax.imshow(cluster_matrix, interpolation='nearest', cmap=plt.cm.Blues)
	plt.colorbar(im)
   	ax.set_xticks(np.arange(cluster_num), minor=False)
	ax.set_xticklabels(np.arange(cluster_num)+1, minor=False)
	ax.set_yticks(np.arange(cluster_num), minor=False)
	ax.set_yticklabels(np.arange(cluster_num)+1, minor=False)
	ax.tick_params(labelbottom='off',labeltop='on', labelleft="on")
   	plt.tight_layout()

def compute_entropy(cluster_matrix, cluster_num):
	entropyByCluster = []
	norm_value = np.log2(cluster_num)
	for i in range(cluster_num):
		entropyByCluster.append(entropy(cluster_matrix[i], base=2)/norm_value)

	return np.sum(entropyByCluster)/cluster_num

import pylab
def plot_together(Z, g_matrix, cluster_num):
	# Compute and plot first dendrogram.
	plt.figure()
	fig = pylab.figure(figsize=(8,8))
	ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
	#Y = sch.linkage(D, method='centroid')
	Z1 = hierarchy.dendrogram(Z, orientation='left')
	ax1.set_xticks([])
	ax1.set_yticks([])

	g_matrix = g_matrix.T
	g_dist = pdist(g_matrix)
	Y = hierarchy.linkage(g_dist, 'complete') #complete, average, single

	# Compute and plot second dendrogram.
	ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
	#Y = sch.linkage(D, method='single')
	Z2 = hierarchy.dendrogram(Y)
	ax2.set_xticks([])
	ax2.set_yticks([])

	# Plot distance matrix.
	axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
	idx1 = Z1['leaves']
	idx2 = Z2['leaves']
	g_matrix = g_matrix.T
	g_matrix = g_matrix[idx1,:]
	g_matrix = g_matrix[:,idx2]
	im = axmatrix.imshow(g_matrix, aspect='auto', origin='lower', cmap=pylab.cm.Blues)
   	axmatrix.set_xticks([])
	axmatrix.set_yticks([])

	# Plot colorbar.
	axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
	pylab.colorbar(im, cax=axcolor)
	fig.show()
	fig.savefig('dendrogram1.png')

def Clustering(Z, g, g_matrix, threshold, network_name, show_boolean=False, linkage="complete", Ecoli=True):
	clusters = fcluster(Z, threshold, criterion='distance')
	cluster_num = np.max(clusters)
	print '\nclusters_num: %d\n' %cluster_num
	cluster_member_num = []
	cluster_member= []
	transformed_g_matrix = np.zeros((cluster_num, len(g_matrix)))
	cluster_matrix = np.zeros((cluster_num,cluster_num))

	if Ecoli == True:
		print "Transfer"
		indegree_num = g.outdegree()
		outdegree_num = g.indegree()
	else:
		indegree_num = g.indegree()
		outdegree_num = g.outdegree()

	indegreeByCluster = []
	outdegreeByCluster = []

	for i in range(1, cluster_num+1):
		idx = np.where(clusters==i)
		print "[", i, "]", idx[0], "len", len(idx[0])
		print g.vs[2], g.vs[82]
		transformed_g_matrix[i-1] = g_matrix[idx].sum(0)
		cluster_member.append(idx)
		cluster_member_num.append(len(idx[0]))
		indegreeByCluster.append([indegree_num[i] for i in idx[0]])
		outdegreeByCluster.append([outdegree_num[i] for i in idx[0]])

	transformed_g_matrix = transformed_g_matrix.T
	for i in range(1, cluster_num+1):
		idx = np.where(clusters==i)
		cluster_matrix[i-1] = transformed_g_matrix[idx].sum(0)

	print "\ncluster_matrix\n", cluster_matrix.T
	cluster_matrix = cluster_matrix / cluster_matrix.sum(axis=0)
	cluster_matrix = cluster_matrix.T
	print "\ncluster_matrix\n", cluster_matrix
	plot_together(Z, g_matrix, cluster_num)
	top_max_num = cluster_num
	#top_max = np.array(cluster_member_num).argsort()[-top_max_num:][::-1]
	top_max = cluster_num
	plot_inoutdegree(top_max, cluster_member_num, indegreeByCluster, outdegreeByCluster, network_name)

	if show_boolean:
		plot_dendrogram(Z, threshold, linkage, network_name)
		plot_heatmap(cluster_matrix, cluster_num)

	return round(compute_entropy(cluster_matrix, cluster_num),2), cluster_num

def Dendrogram(g_matrix, linkageType, g, network_name, threshold, show_boolean=False, Ecoli=True):
	g_dist = pdist(g_matrix)
	Z = hierarchy.linkage(g_dist, linkageType) #complete, average, single
	return Clustering(Z, g, g_matrix, threshold, network_name, show_boolean, linkageType, Ecoli)

def Analysis(testfile, linkageType, threshold, network_name, show_boolean=False, Ecoli=True):
	edge_list_neph = pandas.read_csv(testfile,
                                 sep="\t",
                                 names=["regulator","target"])

	g = Graph.TupleList(edge_list_neph[["target","regulator"]].values.tolist(),
					directed=True)

	g_matrix= np.matrix(g.get_adjacency().data)
	if Ecoli == True:
		print "Transfer"
		g_matrix = g_matrix.T

	return Dendrogram(g_matrix, linkageType, g, network_name, threshold, show_boolean, Ecoli)

def RunMain():
	# neph_gene_network (538), ecolitfnet (133)
	filepath = "" #"../data/"
	testfile = ["ecolitfnet", "neph_gene_network"]
	linkageType = ["single","complete","average"]
	network_name = ["E. coli", "Human"]
	#print Analysis(filepath+testfile[0]+".txt", linkageType[1], 4.5, network_name[0], True, True) #4.5 True 3.2 False
	print Analysis(filepath+testfile[1]+".txt", linkageType[1], 16.7, network_name[1], True, False) #14.1 True 16.7 False
	plt.show()

if __name__ == "__main__":
    RunMain()