import numpy as np
import math

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
                # sum += np.sum(np.fabs(X[i][j] - np.matmul(X[i][k], X[k][j])))
            C[i][j] = 1 - sum / (2* num_graph * num_node)
    return C

def cal_unary_consistency_matrix(X,num_graph,num_node):
    C=np.zeros(num_graph)
    for k in range(num_graph):
        sum = 0
        for i in range(num_graph-1):
            for j in range(i+1,num_graph):
                # F范式： 绝对值平方和再开方
                sum += np.linalg.norm(X[i][j] - np.matmul(X[i][k], X[k][j]))
                # 另一种方法： 直接求绝对值的和
                # sum += np.sum(np.fabs(X[i][j] - np.matmul(X[i][k], X[k][j])))
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
    Lambda=0
    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    #consistency_matrix = cal_unary_consistency_matrix(X,num_graph,num_node)
    consistency_matrix = cal_pairwise_consistency_matrix(X, num_graph, num_node)

    for v in range(num_graph):
        for x in range(num_graph):
            for y in range(x+1,num_graph):
                # calculate S_org
                J_xy=(affinity_matrix[x][y]-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
                Cp_xy=consistency_matrix[x][y]
                S_org=(1-Lambda)*J_xy+Lambda*Cp_xy

                # calculate S_opt
                X_xv_vy = np.matmul(X[x][v],X[v][y]).ravel()
                J_xv_vy=(np.matmul(np.matmul(X_xv_vy.T, K[x][y]), X_xv_vy)-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
                #C_xv_vy = consistency_matrix[v]
                C_xv_vy=math.sqrt(consistency_matrix[x][v]*consistency_matrix[v][y])
                S_opt=(1-Lambda)*J_xv_vy+Lambda*C_xv_vy

                # compare and update
                if S_org < S_opt:
                    X[x][y]=np.matmul(X[x][v],X[v][y])
                    X[y][x]=np.matmul(X[y][v],X[v][x])

    # set lambda and repeat above process
    Lambda=0.45
    consistency_matrix1 = cal_pairwise_consistency_matrix(X, num_graph, num_node)

    for v in range(num_graph):
        affinity_matrix = cal_affinity_matrix(X, K, num_graph)
        # use unary consistency to speed up
        consistency_matrix = cal_unary_consistency_matrix(X,num_graph,num_node)
        for x in range(num_graph):
            for y in range(x+1,num_graph):
                J_xy = (affinity_matrix[x][y] -np.min(affinity_matrix))/ (np.max(affinity_matrix)-np.min(affinity_matrix))
                Cp_xy = consistency_matrix1[x][y]
                S_org = (1 - Lambda) * J_xy + Lambda * Cp_xy

                X_xv_vy = np.matmul(X[x][v], X[v][y]).ravel()
                J_xv_vy = (np.matmul(np.matmul(X_xv_vy.T, K[x][y]), X_xv_vy)-np.min(affinity_matrix)) / (np.max(affinity_matrix)-np.min(affinity_matrix))
                C_xv_vy = consistency_matrix[v]
                #C_xv_vy = math.sqrt(consistency_matrix[x][v] * consistency_matrix[v][y])
                S_opt = (1 - Lambda) * J_xv_vy + Lambda * C_xv_vy

                if S_org < S_opt:
                    X[x][y] = np.matmul(X[x][v], X[v][y])
                    X[y][x] = np.matmul(X[y][v],X[v][x])

    return X