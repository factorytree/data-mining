clusters=[
    [2,3,4,5],
    [1,6],
    [7,8,9]
]
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
        begin=index[i]
        end=index[i+1]
        rearrange=graph_rearrange[begin:end]+graph_rearrange[:begin]+graph_rearrange[end:]
        print(rearrange)