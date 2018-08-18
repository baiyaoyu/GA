import math
import numpy as np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

P=800
def f(x):
	return 20.0/2047*x-10.0
#计算适应度#
def cal(x):
	s1=0
	s2=0
	for i in range(1,6):
		s1=s1+i*math.sin((1+i)*x[0]+i)
		s2=s2+i*math.sin((1+i)*x[1]+i)
	s=s1*s2
	return s
#解码将基因型解码为表现型#
def Decode(bits):
	x1=int(bits[0])
	x2=int(bits[11])
	for x in range(1,11):
		x1=int(bits[x])+(x1<<1)
		x2=int(bits[x+11])+(x2<<1)
	t1=f(x1)
	t2=f(x2)
	return [t1,t2]
#编码将表现型编码为基因型#
def Encode(t):
	x1=int((t[0]+10.0)/20.0*2047)
	x2=int((t[1]+10.0)/20.0*2047)
	bits1=[]
	bits2=[]
	for x in range(0,11):
		if x1%2==0:
			bits1.append('0')
		else:
			bits1.append('1')
		x1=x1//2
		if x2%2==0:
			bits2.append('0')
		else:
			bits2.append('1')
		x2=x2//2
	bits1.reverse()
	bits2.reverse()
	bits=bits1+bits2
	return bits
#随机生成100个个体#
def ran():
	array=np.random.rand(P,2)*20-10
	return array

def getindex(x,percentage):
	for i in percentage:
		if i>x:
			t =percentage.index(i)
			return t
#生成下一代种群个体
def getNextGen(array):
	degree=[]
	percentage=[]
	sum=0
	s=0
	for x in array:
		t=cal(x)+225
		degree.append(t)
		sum=sum+t
	for x in degree:
		s=s+x
		percentage.append(s/sum)
	#print(percentage)
	r=np.random.rand(P)
	next_array=[]
	for x in r:
		index=getindex(x,percentage)
		next_array.append(array[index].copy())
	for i in range(0,P):
		array[i]=next_array[i]
#模拟交叉互换
def getSwap(array):
	percentage=np.random.rand(P//2)                  
	for i in percentage:
		if i < 0.03:
			pos=np.where(percentage==i)[0][0]*2
			bits_A=Encode(array[pos])
			bits_B=Encode(array[pos+1])
			pos_bit=np.random.randint(0,22)
			value=copy.deepcopy(bits_A[pos_bit])
			bits_A[pos_bit]=copy.deepcopy(bits_B[pos_bit])
			bits_B[pos_bit]=copy.deepcopy(value)
			array[0]=Decode(bits_A)
			array[0]=Decode(bits_B)
#模拟基因突变
def getChange(array):
	percentage=np.random.rand(P)
	for i in percentage:
		if i < 0.01:
			pos=np.where(percentage==i)[0][0]
			bits_C=Encode(array[pos])
			pos_bit=np.random.randint(0,22)
			if bits_C[pos_bit]=='1':
				bits_C[pos_bit]='0'
			else:
				bits_C[pos_bit]='1'
			array[pos]=Decode(bits_C)

def getMax(array):
	max=0
	t=np.array(2)
	for i in array:
		temp=cal(i)
		if(max<temp):
			max=temp
			t=copy.deepcopy(i)
	return max,t

def drawf(t,max):
	fig=plt.figure()
	ax=Axes3D(fig)
	X=np.arange(-10,10,0.1)
	Y=np.arange(-10,10,0.1)
	X,Y=np.meshgrid(X,Y)
	S1=np.sin((2)*X+1)+2*np.sin((3)*X+2)+3*np.sin((4)*X+3)+4*np.sin((5)*X+4)+5*np.sin((6)*X+5)
	S2=np.sin((2)*Y+1)+2*np.sin((3)*Y+2)+3*np.sin((4)*Y+3)+4*np.sin((5)*Y+4)+5*np.sin((6)*Y+5)
	Z=S1*S2
	ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=plt.get_cmap('rainbow'))
	ax.contourf(X,Y,Z,zdir='z',offset=-2,cmap='rainbow')
	ax.set_zlim(-200,200)
	ax.scatter(t[0], t[1], max, c = 'g', marker = '.')
	plt.show()

def main():
	array=ran()
	for g in range(0,400):
		getSwap(array)
		getChange(array)
		getNextGen(array)
		getMax(array)
		#print(array)
		max,t=getMax(array)
		print("%d代: %f" %(g,max))
	print(t)
	#for x in array:
		#print(Encode(x))
	drawf(t,max)

if __name__ == '__main__':
	main()