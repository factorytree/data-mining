import numpy as np
import math
from utils.linked_list import Node, LinkedList

def cal_affinity_matrix(X,K,num_graph):
    affinity = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            X_ij = X[i][j].ravel()
            affinity[i][j] = np.matmul(np.matmul(X_ij.T, K[i][j]), X_ij)
    # 两种归一化方法，实现时采用了第二种
    # affinity = affinity / np.max(affinity)
    # affinity=(affinity-np.min(affinity))/(np.max(affinity)-np.min(affinity))
    return affinity

def cal_pairwise_consistency_matrix(X,num_graph,num_node):
    C = np.zeros((num_graph, num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            sum = 0
            for k in range(num_graph):
                sum += np.linalg.norm(X[i][j] - np.matmul(X[i][k], X[k][j]))
            C[i][j] = 1 - sum / (2* num_graph * num_node)
    return C

def single_affinity(X,K):
    X_ij=X.ravel()
    return np.matmul(np.matmul(X_ij.T, K), X_ij)

def cal_unary_consistency_matrix(X,num_graph,num_node):
    C=np.zeros(num_graph)
    for k in range(num_graph):
        sum = 0
        for i in range(num_graph-1):
            for j in range(i+1,num_graph):
                # F范式： 绝对值平方和再开方
                sum += np.linalg.norm(X[i][j] - np.matmul(X[i][k], X[k][j]))
        sum /= (num_node*num_graph*(num_graph-1))
        C[k]=1-sum
    return C


def mgm_floyd(X, K, num_graph, num_node):
    """
    :param K: affinity matrix, (num_graph, num_graph, num_node^2, num_node^2)
    :param num_graph: number of graph, int
    :param num_node: number of node, int
    :return: matching results, (num_graph, num_graph, num_node, num_node)
    """

    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    max=np.max(affinity_matrix)
    min=np.min(affinity_matrix)

    linked_k=LinkedList()
    for i in range(num_graph):
        linked_k.add(i)
    linked_nodes = LinkedList()
    while not linked_k.isEmpty():
        v=linked_k.pop()
        for i in range(v):
            linked_nodes.add(i)
        for i in range(v+1,num_graph):
            linked_nodes.append(i)
        while not linked_nodes.isEmpty():
            x = linked_nodes.pop()
            for y in range(num_graph):
                # calculate S_org
                J_xy_ori=single_affinity(X[x][y],K[x][y])
                J_xy=(J_xy_ori-min)/(max-min)
                S_org=J_xy

                # calculate S_opt
                X_xv_vy = np.matmul(X[x][v],X[v][y])
                J_xv_vy=(single_affinity(X_xv_vy,K[x][y])-min)/(max-min)
                S_opt=J_xv_vy

                # compare and update
                if S_org < S_opt:
                    X[x][y]=np.matmul(X[x][v],X[v][y])
                    X[y][x]=np.matmul(X[y][v],X[v][x])
                    if linked_nodes.search(y):
                        linked_nodes.remove(y)
                        linked_nodes.append(y)
                    if linked_k.search(y):
                        linked_k.remove(y)
                        linked_k.append(y)
                    if linked_k.search(x):
                        linked_k.remove(x)
                        linked_k.append(x)
                if J_xy_ori<min:
                    min=J_xy_ori
                elif J_xy_ori>max:
                    max=J_xy_ori

    # set lambda and repeat above process
    Lambda=0.45
    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    max = np.max(affinity_matrix)
    min = np.min(affinity_matrix)
    # use unary consistency to speed up
    consistency_matrix = cal_unary_consistency_matrix(X, num_graph, num_node)

    flag=False # use flag to check whether X is updated
    for i in range(num_graph):
        linked_k.add(i)
    while not linked_k.isEmpty():
        v=linked_k.pop()
        for i in range(v):
            linked_nodes.add(i)
        for i in range(v+1,num_graph):
            linked_nodes.add(i)
        if flag:
            consistency_matrix = cal_unary_consistency_matrix(X,num_graph,num_node)
            flag=False
        while not linked_nodes.isEmpty():
            x=linked_nodes.pop()
            for y in range(num_graph):
                J_xy_ori = single_affinity(X[x][y], K[x][y])
                J_xy = (J_xy_ori - min) / (max - min)
                Cp_xy = consistency_matrix[y]
                S_org = (1 - Lambda) * J_xy + Lambda * Cp_xy

                X_xv_vy = np.matmul(X[x][v], X[v][y])
                J_xv_vy=(single_affinity(X_xv_vy,K[x][y])-min)/(max-min)
                C_xv_vy = consistency_matrix[v]
                S_opt = (1 - Lambda) * J_xv_vy + Lambda * C_xv_vy

                if S_org < S_opt:
                    X[x][y] = np.matmul(X[x][v], X[v][y])
                    X[y][x] = np.matmul(X[y][v],X[v][x])
                    flag=True
                    if linked_nodes.search(y):
                        linked_nodes.remove(y)
                        linked_nodes.append(y)
                    if linked_k.search(y):
                        linked_k.remove(y)
                        linked_k.append(y)
                    if linked_k.search(x):
                        linked_k.remove(x)
                        linked_k.append(x)
                if J_xy_ori<min:
                    min=J_xy_ori
                elif J_xy_ori>max:
                    max=J_xy_ori
    return X
