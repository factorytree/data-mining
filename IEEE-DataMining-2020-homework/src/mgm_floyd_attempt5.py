import numpy as np
import math
from sklearn.cluster import SpectralClustering
import warnings
warnings.filterwarnings("ignore")

def cal_affinity_matrix(X,K,num_graph):
    affinity = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            X_ij = X[i][j].ravel()
            affinity[i][j] = np.matmul(np.matmul(X_ij.T, K[i][j]), X_ij)
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
    Lambda=0
    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    max=np.max(affinity_matrix)
    min=np.min(affinity_matrix)
    affinity_matrix=(affinity_matrix-min)/(max-min)


    clu_number=2
    cluster = SpectralClustering(n_clusters=clu_number, affinity='precomputed')
    labels_ = cluster.fit_predict(affinity_matrix)
    # print(labels_)
    clusters = [[] for i in range(clu_number)]
    for i in range(num_graph - 1):
        clusters[labels_[i]].append(i)
    index=[0]
    tmp=0
    for i in range(len(clusters)):
        tmp+=len(clusters[i])
        index.append(tmp)
    graph_rearrange=[]
    for item in clusters:
        graph_rearrange.extend(item)
    # raise Exception('STOP!')

    for i in range(len(clusters)):
        for v in clusters[i]:
            begin = index[i]
            end = index[i + 1]
            rearrange = graph_rearrange[begin:end] + graph_rearrange[:begin] + graph_rearrange[end:]
            for x in rearrange:
                for y in rearrange:
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
                    if J_xy_ori<min:
                        min=J_xy_ori
                    elif J_xy_ori>max:
                        max=J_xy_ori

    # set lambda and repeat above process
    Lambda=0.45
    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    max = np.max(affinity_matrix)
    min = np.min(affinity_matrix)
    # consistency_matrix = cal_pairwise_consistency_matrix(X, num_graph, num_node)
    # use unary consistency to speed up
    consistency_matrix = cal_unary_consistency_matrix(X, num_graph, num_node)

    flag=False # use flag to check whether X is updated
    for i in range(len(clusters)):
        for v in clusters[i]:
            begin = index[i]
            end = index[i + 1]
            rearrange = graph_rearrange[begin:end] + graph_rearrange[:begin] + graph_rearrange[end:]
            if flag:
                # consistency_matrix = cal_pairwise_consistency_matrix(X, num_graph, num_node)
                # use unary consistency to speed up
                consistency_matrix = cal_unary_consistency_matrix(X,num_graph,num_node)
                flag=False
            for x in rearrange:
                for y in rearrange:
                    J_xy_ori = single_affinity(X[x][y], K[x][y])
                    J_xy = (J_xy_ori - min) / (max - min)
                    Cp_xy = consistency_matrix[y]
                    S_org = (1 - Lambda) * J_xy + Lambda * Cp_xy

                    X_xv_vy = np.matmul(X[x][v], X[v][y])
                    J_xv_vy=(single_affinity(X_xv_vy,K[x][y])-min)/(max-min)
                    # if use unary consistency
                    C_xv_vy = consistency_matrix[v]
                    # if use pairwise consistency
                    # C_xv_vy = math.sqrt(consistency_matrix[x][v] * consistency_matrix[v][y])
                    S_opt = (1 - Lambda) * J_xv_vy + Lambda * C_xv_vy

                    if S_org < S_opt:
                        X[x][y] = np.matmul(X[x][v], X[v][y])
                        X[y][x] = np.matmul(X[y][v],X[v][x])
                        flag=True
                    if J_xy_ori<min:
                        min=J_xy_ori
                    elif J_xy_ori>max:
                        max=J_xy_ori
    return X

# result
# 1. a.从Gx到Gy遍历 -> b.遍历Gx, Gy从Gx开始往后遍历 -> c.第一次采用a.,第二次采用b. -> d.第一次采用b.,第二次采用a.
#      class car    Motorbike Face   Winebottle Duck
#   (a)acc   0.8435 0.9020    0.9525 0.7776     0.7884
#   (a)time  2.8891 3.0251    3.0453 3.1512     3.1565
#   (b)acc   0.8087 0.8736    0.9523 0.7611     0.7569
#   (b)time  2.4123 2.6173    2.6173 2.7567     2.6516 (最快)
#   (c)acc   0.8226 0.9003    0.9530 0.7639     0.7756
#   (c)time  2.7339 2.7534    2.8449 2.8325     2.9590
#   (d)acc   0.8503 0.9074    0.9560 0.7863     0.7862
#   (d)time  2.7719 2.9673    3.0984 3.2205     3.1821 (最准确)
#attempt2  a.第一步都遍历，改变第二步 b.遍历Gx,Gy从Gx开始，改变第二步 e.第一步第二步都采取新方法
#一次错误的代码
#   (a)acc   0.8243 0.9001    0.9521 0.7813     0.7675
#   (a)time  0.5456 0.5165    0.5195 0.5076     0.5225
#   (b)acc   0.8295 0.9092    0.9522 0.7844     0.7702
#   (b)time  0.5766 0.5217    0.5270 0.5560     0.5676
#   (c)acc   0.7544 0.8749    0.9499 0.7388     0.7252
#   (c)time  0.3232 0.2980    0.3335 0.3527     0.3319
#修正代码后 a.第一步第二步都采取新方法
#   (a)acc   0.8509 0.9027    0.9550 0.7715     0.7897
#   (a)time  2.7420 3.0015    2.8579 2.9757     2.8539
#further
#   (a)acc   0.8605 0.9152    0.9559 0.7831     0.7863
#   (a)time  2.6616 2.7430    2.7477 2.9322     2.9592
#attempt3
#   (a)acc   0.8290 0.9065    0.9542 0.7593     0.7825
#   (a)time  3.2250 3.2635    3.1534 3.4088     3.3669
#attempt4
#   (a)acc   0.8623 0.9179    0.9574 0.7717     0.7725
#   (a)time  3.1819 3.2954    3.1382 3.3588     3.2585


