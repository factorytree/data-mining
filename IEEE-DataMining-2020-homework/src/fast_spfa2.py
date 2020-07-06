import queue
import random
import numpy as np

def affinity_score(X_ij, K_ij):
    vec_X = X_ij.reshape(-1, 1)
    vec_XT = vec_X.transpose()
    affinity = np.matmul(np.matmul(vec_XT, K_ij), vec_X)[0][0]
    return affinity

def get_affinity_score(X, K, num_graph):
    affi_m = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            X_ij = X[i][j]
            K_ij = K[i][j]
            affi_m[i][j] = affinity_score(X_ij,K_ij)
    return  affi_m

def pairwise_score(X, i,j):
    num_graphs, a,num_nodes,c= X.shape
    sum_ =sum([ np.linalg.norm(X[i][k] - np.matmul( X[i][k],X[k][j]))/2 for k in range(num_graphs)])
    return (1 - sum_/ (2*num_graphs*num_nodes))

def get_pairwise_socre(X,num_graph):
    pc_m = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            pc_m[i][j] = pairwise_score(X, i, j)
    return  pc_m



def fast_spfa(K, X, num_graph, num_node):
    """
    :param K: affinity matrix, (num_graph, num_graph, num_node^2, num_node^2)
    :param X: matching results, X[:-1, :-1] is the matching results obtained by last iteration of MGM-SPFA,
              X[num_graph,:] and X[:,num_graph] is obtained via two-graph matching solver(RRWM), We suppose the last
              graph is the new coming graph. (num_graph, num_graph, num_node, num_node)
    :param num_graph: number of graph, int
    :param num_node: number of node, int
    :return: X, matching results, match graph_m to {graph_1, ... , graph_m-1)
    """
    Lambda = 0.4
    c_min = 2


    affi_ma = get_affinity_score(X, K, num_graph)
    pc_ma= get_pairwise_socre(X,num_graph)
    # 归一化
    affi_min = np.min(affi_ma)
    affi_max = np.max(affi_ma)
    nor_func = lambda  x : (x - affi_min) / (affi_max - affi_min)

    # processing graph n
    G_n = num_graph - 1

    # calculate the clusters
    clu_number = max(1, num_graph//c_min)
    all_graph = list(range(num_graph - 1))
    random.shuffle(all_graph)
    clusters = [all_graph[i:i + c_min] for i in range(0, len(all_graph), c_min)]
    if len(clusters) > clu_number:
        lst_clus = clusters[-1]
        for i in range(0, len(lst_clus), clu_number):
            for j in range(clu_number):
                if (i + j) <= len(lst_clus) - 1:
                    clusters[j].append(lst_clus[i + j])
        clusters = clusters[0:-1]


    # for each clusters do spfa
    for now_clus in clusters:
        sub_allgraph = now_clus + [num_graph - 1]
        graph_queue = queue.Queue()
        for gra in now_clus:
            graph_queue.put(gra)
        # update X
        upate_times = 0
        while(not (graph_queue.empty())):
            upate_times +=1
            G_x= graph_queue.get()
            for G_y in now_clus :
                if G_x == G_y :
                    continue
                # use equation 7
                J_yn = nor_func(affinity_score(X[G_y][G_n], K[G_y][G_n]))
                pc_yy = pc_ma[G_y][G_y]
                pc_yn = pc_ma[G_y][G_n]
                C_yn = np.sqrt(pc_yy * pc_yn)
                s_org = (1 - Lambda) * J_yn + Lambda * C_yn
                # use equation 7 too
                X_yxxn = np.matmul(X[G_y][G_x],X[G_x][G_n])
                J_yxxn = nor_func(affinity_score(X_yxxn, K[G_y][G_n]))
                pc_yx = pc_ma[G_y][G_x]
                pc_xn = pc_ma[G_x][G_n]
                C_yxxn = np.sqrt(pc_yx * pc_xn)
                s_opt = (1 - Lambda) * J_yxxn + Lambda * C_yxxn
                if s_org < s_opt :
                    X[G_y][G_n] = X_yxxn
                    X[G_n][G_y] = X_yxxn.T
                    pc_ma[G_y][G_n] = pairwise_score(X, G_y, G_n)
                    pc_ma[G_n][G_y] = pairwise_score(X, G_n, G_y)
                    graph_queue.put(G_y)
                # every two times, update pc_ma
                # if upate_times % 2 == 0:
                #     pc_ma = get_pairwise_socre(X, num_graph)
                # limit updating
                if upate_times > num_graph * num_graph :
                    break

        # last update
        for G_x in now_clus:
            for G_y in now_clus:
                if G_x == G_y :
                    continue
                    # use equation 7
                J_xy = nor_func(affinity_score(X[G_x][G_y], K[G_x][G_y]))
                pc_xx = pc_ma[G_x][G_x]
                pc_xy = pc_ma[G_x][G_y]
                C_xy = np.sqrt(pc_xx * pc_xy)
                s_org = (1 - Lambda) * J_xy + Lambda * C_xy
                # use equation 7 too
                X_xnny = np.matmul(X[G_x][G_n], X[G_n][G_y])
                J_xnny = nor_func(affinity_score(X_xnny, K[G_x][G_y]))
                pc_xn = pc_ma[G_x][G_n]
                pc_ny = pc_ma[G_n][G_y]
                C_xnny = np.sqrt(pc_xn * pc_ny)
                s_opt = (1 - Lambda) * J_xnny + Lambda * C_xnny
                if s_org < s_opt:
                    X[G_x][G_y] = X_xnny

    # out of clusters do general update
    pc_ma = get_pairwise_socre(X, num_graph)
    for G_x in range(num_graph -1 ):
        for G_y in range(num_graph - 1):
            if G_x == G_y:
                continue
                # use equation 7
            J_xy = nor_func(affinity_score(X[G_x][G_y], K[G_x][G_y]))
            pc_xx = pc_ma[G_x][G_x]
            pc_xy = pc_ma[G_x][G_y]
            C_xy = np.sqrt(pc_xx * pc_xy)
            s_org = (1 - Lambda) * J_xy + Lambda * C_xy
            # use equation 7 too
            X_xnny = np.matmul(X[G_x][G_n], X[G_n][G_y])
            J_xnny = nor_func(affinity_score(X_xnny, K[G_x][G_y]))
            pc_xn = pc_ma[G_x][G_n]
            pc_ny = pc_ma[G_n][G_y]
            C_xnny = np.sqrt(pc_xn * pc_ny)
            s_opt = (1 - Lambda) * J_xnny + Lambda * C_xnny
            if s_org < s_opt:
                X[G_x][G_y] = X_xnny
    return X









