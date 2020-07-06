import numpy as np
import math
import time
def cal_affinity_matrix(X,K,num_graph):
    affinity = np.zeros((num_graph,num_graph))
    for i in range(num_graph):
        for j in range(num_graph):
            X_ij = X[i][j].ravel()
            affinity[i][j] = np.matmul(np.matmul(X_ij.T, K[i][j]), X_ij)
    # affinity = affinity / np.max(affinity)
    #affinity=(affinity-np.min(affinity))/(np.max(affinity)-np.min(affinity))
    return affinity
def single_affinity(X,K):
    X_ij=X.ravel()
    return np.matmul(np.matmul(X_ij.T, K), X_ij)
def pairwise_score(X, i,j):
    num_graphs, a,num_nodes,c= X.shape
    sum_ =sum([ np.linalg.norm(X[i][k] - np.matmul( X[i][k],X[k][j]))/2 for k in range(num_graphs)])
    return (1 - sum_/ (2*num_graphs*num_nodes))
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
    Lambda=0.3
    count=0
    N=num_graph-1
    q=[i for i in range(N)]
    affinity_matrix = cal_affinity_matrix(X, K, num_graph)
    max=np.max(affinity_matrix)
    min=np.min(affinity_matrix)

    consistency_matrix = cal_pairwise_consistency_matrix(X, num_graph, num_node)
    updatetime=0
    total=0
    flag1=0
    flag2=0
    #print("q :",q)
    while len(q)!=0:
        x=q[0]
        q.remove(x)
        count+=1
        #print("use ",x,"to update ")
        flag1=1-flag1
        if flag1==1:
            flag2=0
        #print("flag1",flag1)
        for y in range(N):
            if y!=x:

                # calculate S_org
                J_yN=(single_affinity(X[y][N],K[y][N])-min)/(max-min)
                Cp_yN=np.sqrt(consistency_matrix[y][y]*consistency_matrix[y][N])
                S_org=(1-Lambda)*J_yN+Lambda*Cp_yN

                # calculate S_opt
                X_yxX_xN=np.matmul(X[y][x],X[x][N])
                J_yx_xN=(single_affinity(X_yxX_xN,K[y][N])-min)/(max-min)
                C_yx_xN=np.sqrt(consistency_matrix[y][x]*consistency_matrix[x][N])
                S_opt=(1-Lambda)*J_yx_xN+Lambda*C_yx_xN

                if S_org<S_opt:
                    X[y][N]=X_yxX_xN
                    X[N][y]=X_yxX_xN.transpose()

                    flag2=1
                    #if y not in q:
                        #print("y not in q")
                    q.append(y)

        if count>num_graph*num_graph:
            q=[]
            break
        #time1=time.time()
        if count%2==0 and count!=0 and flag2==1:
            updatetime+=1
            #print("update")
            consistency_matrix=cal_pairwise_consistency_matrix(X,num_graph,num_node)
        #time2=time.time()
        #total+=time2-time1
    #print(updatetime,total)


    consistency_matrix=cal_pairwise_consistency_matrix(X,num_graph,num_node)
    for x in range(num_graph-2):
        for y in range(x+1,num_graph-1):

            N=num_graph-1
            J_xy=(single_affinity(X[x][y],K[x][y])-min)/(max-min)
            Cp_xy=np.sqrt(consistency_matrix[x][x]*consistency_matrix[x][y])
            S_org=(1-Lambda)*J_xy+Lambda*Cp_xy

            # calculate S_opt
            X_xNX_Ny = np.matmul(X[x][N],X[N][y])
            J_xN_Ny=(single_affinity(X_xNX_Ny,K[x][y])-min)/(max-min)

            C_xN_Ny=np.sqrt(consistency_matrix[x][N]*consistency_matrix[N][y])
            S_opt=(1-Lambda)*J_xN_Ny+Lambda*C_xN_Ny

            if S_org<S_opt:
                X[x][y]=X_xNX_Ny
                X[y][x]=X_xNX_Ny.transpose()
    return X


