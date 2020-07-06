import numpy as np
import math

def cal_affinity_matrix(X,K,num_graph):
    affinity = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            X_ij = X[i][j].ravel()
            affinity[i][j] = np.matmul(np.matmul(X_ij.T, K[i][j]), X_ij)
    # affinity = affinity / np.max(affinity)
    #affinity=(affinity-np.min(affinity))/(np.max(affinity)-np.min(affinity))
    return affinity

def cal_pairwise_consistency_matrix(X,num_graph,num_node):
    C = np.zeros((num_graph, num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            sum = 0
            for k in range(num_graph):
                sum += np.linalg.norm(X[i][j] - np.matmul(X[i][k], X[k][j]))
                #sum += np.sum(np.fabs(X[i][j] - np.matmul(X[i][k], X[k][j])))
            C[i][j] = 1 - sum / (2* num_graph * num_node)
    return C

def cal_unary_consistency_matrix(X,num_graph,num_node):
    C=np.zeros(num_graph)
    for k in range(num_graph):
        sum = 0
        for i in range(num_graph-1):
            for j in range(i+1,num_graph):
                sum += np.linalg.norm(X[i][j] - np.matmul(X[i][k], X[k][j]))
                #sum += np.sum(np.fabs(X[i][j] - np.matmul(X[i][k], X[k][j])))
        sum /= (num_node*num_graph*(num_graph-1))
        C[k]=1-sum
    return C
def mgm_spfa(K, X, num_graph, num_node):
    """
    :param K: affinity matrix, (num_graph, num_graph, num_node^2, num_node^2)
    :param X: matching results, X[:-1, :-1] is the matching results obtained by last iteration of MGM-SPFA,
              X[num_graph,:] and X[:,num_graph] is obtained via two-graph matching solver(RRWM), We suppose the last
              graph is the new coming graph. (num_graph, num_graph, num_node, num_node)
    :param num_graph: number of graph, int
    :param num_node: number of node, int
    :return: X, matching results, match graph_m to {graph_1, ... , graph_m-1)
    """
    q=[]
    Lambda=0.3
    count=0
    N=num_graph-1
    for i in range(N):
        q.append(i)

    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    #consistency_matrix = cal_unary_consistency_matrix(X,num_graph,num_node)
    consistency_matrix = cal_pairwise_consistency_matrix(X, num_graph, num_node)
    #consistency_matrix2 = cal_unary_consistency_matrix(X,num_graph,num_node)
    while(len(q)!=0):
        x=q[0]
        del q[0]
        count += 1
        for y in range(N):
            if y!=x:

                J_yN=(affinity_matrix[y][N]-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
                Cp_yN=consistency_matrix[y][N]
                S_org=(1-Lambda)*J_yN+Lambda*Cp_yN

                # calculate S_opt
                X_yxX_xN = np.matmul(X[y][x],X[x][N]).ravel()
                J_yx_xN=(np.matmul(np.matmul(X_yxX_xN.T, K[y][N]), X_yxX_xN)-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
                #C_yx_xN = consistency_matrix2[x]
                C_yx_xN=math.sqrt(consistency_matrix[y][x]*consistency_matrix[x][N])
                S_opt=(1-Lambda)*J_yx_xN+Lambda*C_yx_xN

                #S_org=(1-Lambda)*get_affinity_score(X_yN,K_yN,max_affinity)+Lambda*np.sqrt(pairwise_consistency[Gy][Gy]*pairwise_consistency[Gy][num_graph-1])
                #S_opt=(1-Lambda)*get_affinity_score(X_yxX_xN,K_yN,max_affinity)+Lambda*np.sqrt(pairwise_consistency[Gy][Gx]*pairwise_consistency[Gx][num_graph-1])
                if S_org<S_opt:
                    X[y][N]=np.matmul(X[y][x],X[x][N])
                    X[N][y]=X[y][N].transpose()

                    q.append(y)
        if count>num_graph*num_graph:
                break
        if count%2==0:
                consistency_matrix=cal_pairwise_consistency_matrix(X,num_graph,num_node)

    consistency_matrix2 = cal_unary_consistency_matrix(X,num_graph,num_node)
    flag=True
    for x in range(num_graph-1):
        if flag==False:
            break
        for y in range(x+1, num_graph-1):
            count+=1
            N=num_graph-1
            J_xy=(affinity_matrix[x][y]-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
            Cp_xy=consistency_matrix[x][y]
            S_org=(1-Lambda)*J_xy+Lambda*Cp_xy

            # calculate S_opt
            X_xNX_Ny = np.matmul(X[x][N],X[N][y]).ravel()
            J_xN_Ny=(np.matmul(np.matmul(X_xNX_Ny.T, K[y][N]), X_xNX_Ny)-np.min(affinity_matrix))/(np.max(affinity_matrix)-np.min(affinity_matrix))
            #C_xv_vy = consistency_matrix[v]
            C_xN_Ny=consistency_matrix2[N]
            S_opt=(1-Lambda)*J_xN_Ny+Lambda*C_xN_Ny

            if S_org<S_opt:
                X[x][y]=np.matmul(X[x][N],X[N][y])
                X[y][x]=X[x][y].transpose()

    return X
