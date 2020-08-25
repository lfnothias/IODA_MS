import logging
import operator
import sys

import numpy as np
from scipy.stats import multivariate_normal

logger = logging.getLogger("path_finder.curve")


def CentroidSampleControl(centers, intensity_threshold, intensity_ratio):
    centroid_dic = {}
    num_center = 1
    center_intensity_rt_charge = {}
    for i in range(len(centers)):
        if (centers[i, 4] > intensity_threshold) and (
            centers[i, 4] / (centers[i, 3] + 1e-4) > intensity_ratio
        ):
            centroid_dic[(centers[i, 1], centers[i, 0])] = num_center
            center_intensity_rt_charge[num_center] = (
                centers[i, 4],
                centers[i, 1],
                centers[i, 2],
            )
            num_center += 1
    return centroid_dic, num_center - 1, center_intensity_rt_charge


def GMMCluster(data, centroid_dic, restriction, set_restriction):
    if set_restriction:
        labels = np.array([-1] * len(data))
        for j in centroid_dic.keys():
            ind = np.where(
                (data[:, 0] >= j[0] - restriction[0])
                & (data[:, 0] <= j[0] + restriction[0])
                & (data[:, 1] >= j[1] - restriction[1])
                & (data[:, 1] <= j[1] + restriction[1])
            )
            labels[ind[0]] = 1
    else:
        labels = np.array([1] * len(data))
    return labels


# --------------------------------------------------------------
# -----------------GMM clustering Module------------------------
# --------------------------------------------------------------


def GMM(data, centers, centroid_dic, n_iter, Var_init, Var_max):
    #     ----------------------init--------------------
    n_clusters = len(centroid_dic)
    n_data = len(data)
    labels = [-1] * n_data
    Mu = [
        i for i in centroid_dic.keys()
    ]  # Mu is set as a rt_mz tuple, given apex it should not be changed
    Var = [[Var_init, Var_init]] * n_clusters
    Pi = [1 / n_clusters] * n_clusters
    W = np.ones((n_data, n_clusters)) / n_clusters
    #     ----------------------iter--------------------
    for i in range(n_iter):
        W = UpdateW(data, Mu, Var, Pi)
        if i == n_iter - 1:
            break
        Pi = UpdatePi(W)
        Var = UpdateVar(data, Mu, W, Var_max)
    for i in range(W.shape[0]):
        if W[i, :].sum() != 0:
            labels[i] = np.argmax(W[i, :]) + 1
    return labels


def UpdateW(data, Mu, Var, Pi):
    n_data, n_clusters = len(data), len(Pi)
    pdfs = np.zeros((n_data, n_clusters))
    for i in range(n_clusters):
        if Var[i][0] != 0 and Var[i][1] != 0:
            pdfs[:, i] = Pi[i] * multivariate_normal.pdf(data, Mu[i], np.diag(Var[i]))
    W = pdfs / np.sum(pdfs, axis=1).reshape(-1, 1)
    return W


def UpdatePi(W):
    Pi = W.sum(axis=0) / W.sum()
    return Pi


def UpdateVar(data, Mu, W, Var_max):
    n_clusters = W.shape[1]
    Var = np.zeros((n_clusters, 2))
    for i in range(n_clusters):
        if W[:, i].sum() != 0:
            Var[i] = np.average((data - Mu[i]) ** 2, axis=0, weights=W[:, i])
            Var[i][0] = min(Var[i][0], Var_max)
            Var[i][1] = min(Var[i][1], Var_max)
    return Var


# --------------------------------------------------------------
# ------------------------Class DEFs----------------------------
# --------------------------------------------------------------
class Cluster:
    def __init__(self, node, num):
        self.nodes = [node]
        self.intensity = node.intensity
        self.num = num

    def addNode(self, node):
        if node.cluster != self.num:
            logger.error("Node into different cluster")
        self.nodes.append(node)
        self.intensity += node.intensity


class Node:
    def __init__(self, rawSig, cluster, num):
        self.rt = rawSig[0]
        self.mz = [rawSig[1], rawSig[1]]
        self.intensity = rawSig[2]
        self.cluster = cluster
        self.num = num

    def addRawSig(self, rawSig, cluster):
        if rawSig[0] != self.rt:
            logger.error("Different rt, should not merge into same node")
        if cluster != self.cluster:
            logger.error("Different cluster, should not merger into same node")

        if rawSig[1] < self.mz[0]:
            self.mz[0] = rawSig[1]
        elif rawSig[1] > self.mz[1]:
            self.mz[1] = rawSig[1]

        self.intensity += rawSig[2]


class Graph:
    def __init__(self, vertices, node_cluster):
        self.V = vertices
        self.graph = {}
        self.inEdge = {}
        self.outEdge = {}

    def addEdge(self, edge):
        if edge[0] not in self.graph.keys():
            self.graph[edge[0]] = [(edge[1], edge[2])]
        else:
            self.graph[edge[0]].append((edge[1], edge[2]))

        if edge[0] not in self.inEdge.keys() and edge[0] != 0:
            self.inEdge[edge[0]] = {}

        tmp = {edge[0]: True}
        if edge[1] not in self.inEdge.keys():
            self.inEdge[edge[1]] = tmp
        else:
            self.inEdge[edge[1]].update(tmp)

        tmp = {edge[1]: True}
        if edge[0] not in self.outEdge.keys():
            self.outEdge[edge[0]] = tmp
        else:
            self.outEdge[edge[0]].update(tmp)

    def topologicalSort(self, curr, stack):
        while len(curr) != 0:
            tmp_val = curr[0]
            curr = curr[1:]
            if tmp_val in self.outEdge.keys():
                for i in self.outEdge[tmp_val]:
                    del self.inEdge[i][tmp_val]
                    if len(self.inEdge[i]) == 0:
                        curr.append(i)
            stack.append(tmp_val)

    def shortestPath(self, s, t):
        stack = []
        curr = [s]
        self.topologicalSort(curr, stack)
        dist = [1] * (self.V)
        dist[s] = 0
        ancestor = [None] * (self.V)
        while stack:
            i = stack[0]
            stack = stack[1:]
            if i in self.graph.keys():
                for node, weight in self.graph[i]:
                    if dist[node] > dist[i] + weight:
                        dist[node] = dist[i] + weight
                        ancestor[node] = i

        return dist, ancestor


# --------------------------------------------------------------
# -------------------Other Helper Func--------------------------
# --------------------------------------------------------------
def NodeCreate(data, labels):
    nodes = []
    rt_label_node_dic = {}
    node_cluster = {}
    node_index = {}
    num = 1
    for i in range(len(labels)):
        if labels[i] != 0:
            if (data[i, 0], labels[i]) not in rt_label_node_dic.keys():
                rt_label_node_dic[(data[i, 0], labels[i])] = Node(
                    data[i, :], labels[i], num
                )
                node_cluster[num] = labels[i]
                num += 1
            else:
                rt_label_node_dic[(data[i, 0], labels[i])].addRawSig(
                    data[i, :], labels[i]
                )

    for i in rt_label_node_dic:
        nodes.append(rt_label_node_dic[i])
    return nodes, num, node_cluster


def ClusterCreate(nodes, num_centers, int_max):
    clusters = [None] * num_centers
    num = 1
    for i in nodes:
        if clusters[i.cluster - 1] is None:
            clusters[i.cluster - 1] = Cluster(i, i.cluster)
        else:
            clusters[i.cluster - 1].addNode(i)

    ret = []
    for i in clusters:
        if i is None or i.intensity < int_max:
            continue
        i.nodes.sort(key=operator.attrgetter("rt"))
        ret.append(i)
    return ret


def EdgeCreate(clusters, int_max, num_node, delay, min_scan, max_scan):
    edges = []
    # create weight 1 edges
    for i in clusters:
        for j in range(len(i.nodes)):
            sum_int = 0
            for k in range(j, len(i.nodes)):
                sum_int += i.nodes[k].intensity
                if (
                    sum_int > int_max
                    and i.nodes[k].rt > i.nodes[j].rt + min_scan
                    and i.nodes[k].rt < i.nodes[j].rt + max_scan
                ):
                    edges.append([i.nodes[j].num, i.nodes[k].num, -1])
                    break

    # create weight 0 edges
    for i in range(len(clusters)):
        for j in clusters[i].nodes:
            edges.append([0, j.num, 0])  # add start node
            edges.append([j.num + num_node + 1, num_node + 1, 0])  # add end node
            for k in range(i + 1, len(clusters)):
                for m in clusters[k].nodes:
                    if m.rt > j.rt + delay:
                        edges.append([j.num, m.num, 0])
                if clusters[k].nodes[0].rt > j.rt + delay:
                    break
    return edges


# Avoid resample same cluster
def AddPrimeNode(num_node, edge, node_cluster):
    num_edge = len(edge)
    modified_node = {}
    for i in range(num_edge):
        if (
            edge[i][0] != 0
            and edge[i][1] != num_node + 1
            and node_cluster[edge[i][0]] == node_cluster[edge[i][1]]
        ):
            modified_node[edge[i][1]] = True
            edge[i][1] += num_node + 1

    for i in range(num_edge):
        if edge[i][0] in modified_node.keys() and num_node >= edge[i][1]:
            edge[i][0] += num_node + 1

    return edge


# Extract Path from ancestors
def PathExtraction(dist, ancestors, num_node):
    idx = ancestors[num_node]
    path = []
    while idx is not None and idx != 0:
        if dist[idx] - dist[ancestors[idx]] == -1:
            path.append([ancestors[idx], idx])
        idx = ancestors[idx]
    path.reverse()
    return path


# Map the path to index
def IndexHis(path, num_node, nodes, center_intensity_rt_charge, node_cluster):
    index_his = []
    for i in path:
        if i[1] > num_node:
            i[1] -= num_node
        if node_cluster[i[0]] != node_cluster[i[1]]:
            logger.error("Not collecting same sample")
        intensity_rt_charge = center_intensity_rt_charge[node_cluster[i[0]]]
        index_his.append(
            [
                (nodes[i[0]-1].rt, nodes[i[1]-1].rt),
                (
                    min(nodes[i[0]-1].mz[0], nodes[i[1]-1].mz[0]),
                    max(nodes[i[0]-1].mz[1], nodes[i[1]-1].mz[1]),
                ),
                intensity_rt_charge[0],
                intensity_rt_charge[1],
                intensity_rt_charge[2],
            ]
        )
    return index_his


def ClusterRemove(path, num_node, node_cluster, clusters):
    ind = {}
    newClusters = []
    for i in range(len(path)):
        ind[node_cluster[path[i][0]]] = True

    for i in range(len(clusters)):
        if clusters[i].num not in ind.keys():
            newClusters.append(clusters[i])
    return newClusters


def WriteFile(file_name, indice_his, restriction, delay, isolation, min_time, max_time):
    restriction_rt = restriction[0]
    restriction_mz = restriction[1]
    text_file = open(file_name, "wt")
    for i in range(len(indice_his)):
        n = text_file.write("path" + str(i) + "\t")
        for j in range(len(indice_his[i])):
            mz_index = (indice_his[i][j][1][0] + indice_his[i][j][1][1]) / 2.0
            iso = max(indice_his[i][j][1][1] - mz_index, isolation)
            start = indice_his[i][j][0][0]
            end = indice_his[i][j][0][1]
            intensity = indice_his[i][j][2]
            rt = indice_his[i][j][3]
            charge = indice_his[i][j][4]
            dur = end - start
            if j != len(indice_his[i])-1 and end > indice_his[i][j+1][0][0] - delay:
                logger.error("Not enough time for delay. start: %.4f  end: %.4f", end, indice_his[i][j+1][0][0] - delay)
            if dur < min_time or dur > max_time:
                logger.error("Too long / short for scan period: dur %.4f", dur)
            n = text_file.write(
                "{:.4f}".format(mz_index)
                + " "
                + "{:.4f}".format(iso)
                + " "
                + "{:.4f}".format(dur)
                + " "
                + "{:.4f}".format(start)
                + " "
                + "{:.4f}".format(end)
                + " "
                + "{:.4f}".format(intensity)
                + " "
                + "{:.4f}".format(rt)
                + " "
                + "{:.4f}".format(charge)
                + "\t"
            )
        n = text_file.write("\n")
    text_file.close()


def PathGen(
    infile_raw,
    infile_feature,
    intensity_threshold,
    intensity_ratio,
    intensity_accu,
    restriction,
    num_path,
    delay,
    min_scan,
    max_scan,
):
    n_iter = 2
    Var_init = restriction[0]
    Var_max = restriction[0]

    try:
        data = np.genfromtxt(infile_raw, skip_header=12)
        data = data[~np.isnan(data)]
        data = data[np.nonzero(data)]
        data = data.reshape(-1, 3)
        centers = np.genfromtxt(infile_feature, delimiter=",", skip_header=1)
    except:
        logger.error("error in reading data from input file", exc_info=sys.exc_info())
        sys.exit()

    try:
        data = data[data[:, 1].argsort()]  # First sort doesn't need to be stable.
        data = data[data[:, 0].argsort(kind="mergesort")]
        centers = centers[centers[:, 0].argsort()]  # First sort doesn't need to be stable.
        centers = centers[centers[:, 1].argsort(kind="mergesort")]
    except:
        logger.error("error is sorting data", exc_info=sys.exc_info())
        sys.exit()

    logger.info("=============")
    logger.info("File Read")
    logger.info("=============")

    try:
        centroid_dic, num_center, center_intensity_rt_charge = CentroidSampleControl(
            centers, intensity_threshold, intensity_ratio
        )
        labels = GMMCluster(data, centroid_dic, restriction, True)
        data_clean = data[labels != -1]
        labels_clustered = GMM(
            data_clean[:, :2], centers, centroid_dic, n_iter, Var_init, Var_max
        )
    except:
        logger.error("error in clustering data", exc_info=sys.exc_info())
        sys.exit()

    logger.info("Begin Finding Path")
    logger.info("=============")
    indice_his = []
    try:
        nodes, num_node, node_cluster = NodeCreate(data_clean, labels_clustered)
        clusters = ClusterCreate(nodes, num_center, intensity_threshold)
    except:
        logger.error("error in creating nodes and clusters", exc_info=sys.exc_info())
        sys.exit()
    try:
        for i in range(num_path):
            edges = EdgeCreate(
                clusters, intensity_threshold, num_node, delay, min_scan, max_scan
            )
            edges = AddPrimeNode(num_node, edges, node_cluster)
            g = Graph(num_node * 2 + 2, node_cluster)
            for j in edges:
                g.addEdge(j)
            s = 0
            t = num_node + 1
            dist, ancestors = g.shortestPath(s, t)
            if dist[num_node + 1] >= 0:
                break
            logger.info(
                "[%d/%d]: features: %d, rest: %d"
                % (i + 1, num_path, -dist[num_node + 1], len(clusters))
            )
            path = PathExtraction(dist, ancestors, num_node + 1)
            index_his = IndexHis(
                path, num_node + 1, nodes, center_intensity_rt_charge, node_cluster
            )
            indice_his.append(index_his)
            clusters = ClusterRemove(path, num_node, node_cluster, clusters)
    except:
        logger.error("error in finding paths", exc_info=sys.exc_info())
        sys.exit()

    return indice_his
