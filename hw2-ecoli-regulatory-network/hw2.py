# ==============================
# Athout: Eugene Seo
# Date: 11.21.2016
# Description: CS519 Homework 2
# ==============================

# coding: utf-8

# Step 1. load in the E. coli operon-operon network (gene regulatory network) data in edge-list format; display the first six lines of data as a sanity check. Use pandas.read_csv and pandas.DataFrame.head.

# In[1]:


import pandas
edge_list_ecoli = pandas.read_csv("ecolitfnet.txt", sep="\t", names=["source", "target"])
edge_list_ecoli.head(n=6)


# Step 2. Make an igraph directed graph from the network; print a graph summary as a sanity check.

# In[2]:


from igraph import Graph
from igraph import summary
ecoli_graph = Graph.TupleList(edge_list_ecoli.values.tolist(), directed=True)
summary(ecoli_graph)


# Q1. Which one of connected 3-vertex motifs is most fequent in the E. coli regularoty network?

# In[3]:


import numpy as np
three_vertex_motifs_counts = ecoli_graph.motifs_randesu(size=3)
print np.nanargmax(three_vertex_motifs_counts)


# Q2. Which one of these motifs has a count of 47 in the regularoty network? (FFL)

# In[4]:


print three_vertex_motifs_counts.index(47)


# Step 3. Compute the count of this FFL motif in the real network and random networks, and obtain the mean and standard deviation of the FFL counts for the random networks. Finally, obtain the Z-score for the count on the real network vs. the random networks. (Write a function motif_analysis that takes as its inputs motif_size, random_num, rewire_prob, and query motif name and its idx, and that prints out the counts of the motif on the real network, the mean and standard deviation of the motif on random networks, and the Z-score for the count on the real network vs. the random networks.)

# In[5]:


import timeit
def motif_analysis(motif_size=3, random_num=100, rewire_prob=1000, query='FFLs', idx=7):
    start_time = timeit.default_timer()
    real_motif_counts = ecoli_graph.motifs_randesu(size=motif_size)
    matrix_real_motif_counts = np.array(real_motif_counts)
    print 'Counts of all {0}-vertex motifs in the graph:\n{1}'.format(motif_size, real_motif_counts)
    print '\nNumber of {0} on the real network: {1}'.format(query, real_motif_counts[idx])


    random_motif_counts = []
    copied_ecoli_graph = ecoli_graph.copy()
    for i in range(random_num):
        if (i % 1000 == 0):
            print "iteration:", i
        copied_ecoli_graph.rewire(n=rewire_prob)
        random_motif_counts.append(copied_ecoli_graph.motifs_randesu(size=motif_size))
        matrix_random_motif_counts = np.array(random_motif_counts)
        mean_random_motif = np.mean(matrix_random_motif_counts, 0)
        std_random_motif = np.std(matrix_random_motif_counts, 0)
        z_scores = (matrix_real_motif_counts-mean_random_motif)/std_random_motif

    print 'Mean of {0} on random networks: {1}'.format(query, mean_random_motif[idx])
    print 'Standard deviation of {0} on random networks: {1}'.format(query, std_random_motif[idx])
    print 'Z-score for {0} count on the real network vs. the random networks: {1}'.format(query, z_scores[idx])

    anlysis_elapsed = timeit.default_timer() - start_time
    print "\nanlysis_elapsed:", anlysis_elapsed

    return real_motif_counts[idx]/mean_random_motif[idx]


# Step 4. Analyze the FFL motif. (Three-vertex motif analysis)

# In[6]:


FFL_ratio = motif_analysis(3, 10000, 1000, 'FFL', 7)


# Q3. Is this Z-score statistically significant? This Z-score is not statistically significant as it only shows about 1 standard deviation less than the mean.
#
# Q4. What is the ratio of the FFL count for the real network to the average FFL count for the random networks? How does ratio compare to the same ratio for the data in Table 1 in Shen-Orr et al., Nature Genetis, 2002?

# In[7]:


print "Ratio of the FFL for the real network to random networks:", FFL_ratio
print "Ratio of the FFL in Shen-Orr et al., Nature Genetis, 2002:", 40/6.9


# Q5. Given the modest ratio of the FFL frequency in the real network vs. randomly shuffled network, should we entertain the possibility that the high frequency of FFLs in the real network could be a consequence of the degree distribution rather than evolution specifically “favoring” FFLs as a building block for gene regulatory networks? According to the difference of two comparing ratios, we might not be able to assert that the FFL motif is the favoring feature of the E. coli gene regulatory network. There could be a possibility to have a higher frequency of FFLs from other factors including degree distribution.

# Step 5. Figure out the isomorphism class of the network shown in the assignment PDF.

# In[8]:


from igraph import plot
tuple_list = []
tuple_list.append([1,3])
tuple_list.append([1,4])
tuple_list.append([2,3])
tuple_list.append([2,4])
v4_motif_graph = Graph.TupleList(tuple_list, directed=True)
plot(v4_motif_graph, layout="sugiyama", bbox=(100,100), vertex_color = "pink")


# In[9]:


print v4_motif_graph.isoclass()


# Step 6. Analyze this 4-vertex DOR motif. (Four-vertex motif analysis)

# In[10]:


DOR_ratio = motif_analysis(4, 10000, 1000, 'DOR', 19)


# Q6. Is this Z-score statistically significant? This Z-score is still not that statistically significant as it shows about 3 standard deviation greater than the mean.
#
# Q7. Compute the ratio of the DOR count in the real network to the average DOR count for the random networks. How does this ratio compare to the same ratio for the data in Table 1 from Shen-Orr et al.? Are they consistent?

# In[23]:


print "Ratio of the DOR for the real network to random networks:", DOR_ratio
print "Ratio of the DOR in Shen-Orr et al., Nature Genetis, 2002:", 203/(57.0)


# Q8. Does this suggest that Shen-Orr's actual network randomization procedure is possibly not consistent with their description in the Methods section of their paper, i.e., that it may have had some kind of error? Both FFL and DOR's Z-scores and ratios are not consistent with Shen-Orr's results, and it may happen due to wrong random network generation and motif enumeration methods, according to the paper by Konagurthu and Lesk.

# (Bonus Question) Write a program that generates a 20-vertex directed graph using the Barabasi-Albert model in igraph (with m = 2 edges per step), and then count the number of appearances of each type of connected three-vertex motif (note: there are five connected three-vertex digraph isomorphism classes that do not contain any double-arrows). Your program should print out the count of each of the five connected isomorphism classes.

# Step 1. Generates a 20-vertex directed graph using the Barabasi-Albert model in igraph (with m = 2 edges)

# In[15]:


v_num = 20
g = Graph.Barabasi(v_num, 2, directed=True)
summary(g)


# Step 2. Plot the 20-vertex directed graph

# In[16]:


i = 0
for v in g.vs:
	v["label"] = i
	i = i + 1
plot(g, layout = g.layout("kk"), vertex_color = "pink", bbox = (300, 300))


# Step 3. Store each classe's indegree and outdegree in one list to distinguish motifs

# In[17]:


from collections import defaultdict
class_dic = defaultdict(int)
class2  = [2,0,0,0,1,1] # [indegree, outdegree] of 3 vertices motifs
class4  = [1,1,0,0,1,1]
class6  = [1,1,0,0,0,2]
class7  = [2,1,0,0,1,2]
class11 = [1,1,1,1,1,1]


# For each vertex, determine the corresposing class by comparing indegree & outdegree numbers

# In[18]:


import itertools
import operator
for x in itertools.combinations(range(v_num), 3): # for all combinations
	pruned_vs = g.vs.select(x[0], x[1], x[2])
	sub_g = g.subgraph(pruned_vs) # build a subgraph with 3 vertices
	degree_num = sub_g.indegree()
	degree_num.extend(sub_g.outdegree()) # store its indegree and outdegree in order

	if sub_g.is_connected(mode='WEAK'): # ignore not-connected three-vertex motifs
		if all(v == 0 for v in map(operator.sub, degree_num, class2)): # compare its degrees with each motif's degrees
			class_dic['class2'] += 1
		elif all(v == 0 for v in map(operator.sub, degree_num, class4)):
			class_dic['class4'] += 1
		elif all(v == 0 for v in map(operator.sub, degree_num, class6)):
			class_dic['class6'] += 1
		elif all(v == 0 for v in map(operator.sub, degree_num, class7)):
			class_dic['class7'] += 1
		elif all(v == 0 for v in map(operator.sub, degree_num, class11)):
			class_dic['class11'] += 1


# Print out the counts of each of the five connected isomorphism classes from my function and compare to the results from motifs_randesu function

# In[19]:


print "\nResult from my function:\n", [(key, value) for (key, value) in sorted(class_dic.items())]
motif_num = g.motifs_randesu(size=3)
print "\nResult from motif_num function:\n", motif_num

