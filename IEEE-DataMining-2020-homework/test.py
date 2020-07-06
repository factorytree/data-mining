import matplotlib.pyplot as plt



Class=['Car','Motorbike','Face','Winebottle','Duck']
acc={'a':[0.8435,0.9020,0.9525,0.7776,0.7884],
     'b':[0.8087,0.8736,0.9523,0.7611,0.7569],
     'c':[0.8226,0.9003,0.9530,0.7639,0.7756],
     'd':[0.8503,0.9074,0.9560,0.7863,0.7862],
     'idea1':[0.8509,0.9027,0.9550,0.7715,0.7897],
     'idea1_further':[0.8605,0.9152,0.9559,0.7831,0.7863],
     'attempt3':[0.8290,0.9065,0.9542,0.7593,0.7825],
     'attempt4':[0.8623,0.9179,0.9574,0.7717,0.7725],
     'attempt5_last':[0.8280,0.8864,0.9437,0.7770,0.7598],
     'attempt5_first':[0.8086 ,0.8893 ,0.9438 ,0.7705 ,0.7603]}
time={'a':[2.8891,3.0251,3.0453,3.1512,3.1565],
      'b':[2.4123,2.6173,2.6173,2.7567,2.6516],
      'c':[2.7339,2.7534,2.8449,2.8325,2.9590],
      'd':[2.7719,2.9673,3.0984,3.2205,3.1821],
      'idea1':[2.7420,3.0015,2.8579,2.9757,2.8539],
      'idea1_further':[2.6616,2.7430,2.7477,2.9322,2.9592],
      'attempt3':[3.2250,3.2635,3.1534,3.4088,3.3669],
      'attempt4':[3.1819,3.2954,3.1382,3.3588,3.2585],
      'attempt5_last':[2.8704,2.8551,2.9279,2.9547,3.3872],
      'attempt5_first':[3.1195,3.0367,3.0393,3.1420,3.1060]}

def acc_graph():
    plt.plot(Class, acc['a'],'-o', color='skyblue', label='a')
    plt.plot(Class, acc['b'],'-o', color='blue', label='b')
    plt.plot(Class, acc['c'],'-o', color='pink', label='c')
    plt.plot(Class, acc['d'],'-o', color='red', label='d')
    plt.plot(Class, acc['idea1_further'], '-o', color='green', label='attempt 2 further')
    plt.plot(Class, acc['attempt5_last'], '-o', color='lightgreen', label='attempt 5 last')

    plt.legend()

    plt.xlabel('Class')
    plt.ylabel('Accuracy')
    plt.savefig('../improvement_or_changes/acc_attempt5_last')
    plt.show()

def time_graph():
    plt.plot(Class, time['a'],'-o', color='skyblue', label='a')
    plt.plot(Class, time['b'],'-o', color='blue', label='b')
    plt.plot(Class, time['c'],'-o', color='pink', label='c')
    plt.plot(Class, time['d'],'-o', color='red', label='d')
    plt.plot(Class, time['idea1_further'], '-o', color='green', label='attempt 2 further')
    plt.plot(Class, time['attempt5_last'], '-o', color='lightgreen', label='attempt 5 last')

    plt.legend()

    plt.xlabel('Class')
    plt.ylabel('Time')
    plt.savefig('../improvement_or_changes/time_attempt5_last')
    plt.show()

#acc_graph()
time_graph()