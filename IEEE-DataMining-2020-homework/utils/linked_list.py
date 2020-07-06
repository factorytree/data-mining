# 先定一个node的类
class Node():  # value + next
    def __init__(self, value=None, next=None):
        self._value = value
        self._next = next
        self.times = 0

    def getValue(self):
        return self._value

    def getNext(self):
        return self._next

    def setValue(self, new_value):
        self._value = new_value

    def setNext(self, new_next):
        self._next = new_next

    def addtimes(self):
        self.times+=1

    def gettimes(self):
        return self.times


# 实现Linked List及其各类操作方法
class LinkedList():
    def __init__(self):  # 初始化链表为空表
        self._head = Node()
        self._tail = None
        self._length = 0

    # 检测是否为空
    def isEmpty(self):
        return self._head.getValue() == None

        # add在链表前端添加元素:O(1)

    def add(self, value):
        newnode = Node(value, None)  # create一个node（为了插进一个链表）
        newnode.setNext(self._head)
        self._head = newnode
        if self._length==0:
            self._tail=self._head
        self._length+=1

    def pop(self):
        if self._length==1:
            self._tail=None
        value=self._head.getValue()
        self._head=self._head.getNext()
        self._length-=1
        return value

    # append在链表尾部添加元素
    def append(self, value):
        #newnode = Node(value)
        if self.isEmpty():
            self.add(value)  # 若为空表，将添加的元素设为第一个元素
        else:
            newnode = Node(value)
            newnode.setNext(self._tail.getNext())
            self._tail.setNext(newnode)
            self._tail=newnode
            self._length+=1
        # else:
        #     current = self._head
        #     self._length += 1
        #     while current.getNext() != None:
        #         current = current.getNext()  # 遍历链表
        #     current.setNext(newnode)  # 此时current为链表最后的元素


    # search检索元素是否在链表中
    def search(self, value):
        current = self._head
        foundvalue = False
        while current != None and not foundvalue:
            if current.getValue() == value:
                foundvalue = True
            else:
                current = current.getNext()
        return foundvalue

    # index索引元素在链表中的位置
    def index(self, value):
        current = self._head
        count = 0
        found = None
        while current != None and not found:
            count += 1
            if current.getValue() == value:
                found = True
            else:
                current = current.getNext()
        if found:
            return count
        else:
            raise ValueError('%s is not in linkedlist' % value)

    # remove删除链表中的某项元素
    def remove(self, value):
        current = self._head
        if self._length==1:
            self._tail=None
        self._length-=1
        pre = None
        while current != None:
            if current.getValue() == value:
                if not pre:
                    self._head = current.getNext()
                else:
                    pre.setNext(current.getNext())
                    if current==self._tail:
                        self._tail=pre
                break
            else:
                pre = current
                current = current.getNext()

    # insert链表中插入元素
    def insert(self, pos, value):
        if pos <= 1:
            self.add(value)
        elif pos > self.size():
            self.append(value)
        else:
            temp = Node(value)
            count = 1
            pre = None
            current = self._head
            while count < pos:
                count += 1
                pre = current
                current = current.getNext()
            pre.setNext(temp)
            temp.setNext(current)

    def addtimes(self,value):
        current=self._head
        pre=None
        pre_value=None
        while current!=None:
            if current.getValue() == value:
                current.addtimes()
                if not pre_value:
                    if not pre:
                        return
                    pre.setNext(current.getNext())
                    current.setNext(self._head)
                    self._head=current
                elif current.gettimes()>pre.gettimes():
                    pre.setNext(current.getNext())
                    current.setNext(pre_value.getNext())
                    pre_value.setNext(current)
                return
            else:
                pre=current
                current = current.getNext()
                if pre.gettimes()>current.gettimes():
                    pre_value=pre

    def print(self):
        current=self._head
        while current!=None:
            print(current,current.getValue(),current.gettimes())
        print('\n')
