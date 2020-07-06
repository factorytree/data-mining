import matplotlib.pyplot as plt
import numpy as np

file=open('mgm_floyd_result','r')
graph_number=[4,8,12,16,20,24,28,32]
accuracy={'Car':[],'Motorbike':[],'Face':[],'Winebottle':[],'Duck':[]}
time={'Car':[],'Motorbike':[],'Face':[],'Winebottle':[],'Duck':[]}
line=file.readline()
for line in file.readlines():
    line=line.split()
    accuracy[line[1]].append(float(line[2]))
    time[line[1]].append(float(line[3]))
file.close()


def acc_graph():
    plt.plot(graph_number, accuracy['Car'],'-o', color='skyblue', label='Car')
    plt.plot(graph_number, accuracy['Motorbike'],'-o', color='blue', label='Motorbike')
    plt.plot(graph_number, accuracy['Face'],'-o', color='pink', label='Face')
    plt.plot(graph_number, accuracy['Winebottle'],'-o', color='lightgreen', label='Winebottle')
    plt.plot(graph_number, accuracy['Duck'],'-o', color='red', label='Duck')
    plt.legend()

    plt.xlabel('Graph Number')
    plt.ylabel('Accuracy')
    plt.show()

def time_graph():
    plt.plot(graph_number, time['Car'],'-o', color='skyblue', label='Car')
    plt.plot(graph_number, time['Motorbike'],'-o', color='blue', label='Motorbike')
    plt.plot(graph_number, time['Face'],'-o', color='pink', label='Face')
    plt.plot(graph_number, time['Winebottle'],'-o', color='lightgreen', label='Winebottle')
    plt.plot(graph_number, time['Duck'],'-o', color='red', label='Duck')
    plt.legend()

    plt.xlabel('Graph Number')
    plt.ylabel('Time')
    plt.show()

def double_yaxis():
    kind=['Car','Motorbike','Face','Winebottle','Duck']
    baseline_acc=[0.6246, 0.8251, 0.9308, 0.6830, 0.5969]
    baseline_time=[4.0835, 3.9270, 3.9204, 4.0389, 3.9087]
    acc=[0.8087,0.8736,0.9523,0.7611,0.7569]
    time=[2.5691,2.6323,2.5678,2.9049,2.7001]
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(kind, baseline_acc,'-o', color='skyblue',label='Baseline Accuracy')
    ax1.plot(kind,acc, '-o',color='blue',label='My Accuracy')
    ax2.plot(kind, baseline_time, '-o',color='pink',label='Baseline Time')
    ax2.plot(kind,time,'-o',color='red',label='My Time')

    ax1.set_xlabel("Class Name")
    ax1.set_ylabel("Accuracy")
    ax2.set_ylabel("Time")

    ax1.legend()
    ax2.legend()
    plt.show()

time_graph()