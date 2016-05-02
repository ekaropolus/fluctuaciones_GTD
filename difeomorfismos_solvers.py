from numpy import *
from pylab import *
from sympy import Matrix
from sympy import *
from sympy.abc import beta
from itertools import *
import time
import gc


def platexini00(file_name):
    filename = file_name
    FILE = open(filename,'w')                   
    FILE.writelines('\\begin{document}')
    FILE.close()

def platexend00(file_name):
    filename = file_name
    FILE = open(filename,'a')                   
    FILE.writelines('\\end{document}')
    FILE.close()
    
def platex00(expresion,file_name,mode):
    filename = file_name
    FILE = open(filename,'a')                   
    FILE.writelines(latex(expresion,mode=mode))
    FILE.close()
	
def send_D_to_latex(equations,determinants,integrability,ZX,file_name):
	to_latex = []
	clean_D = []
	
	platexini00(file_name)
	for x in equations:
		tabix = 0
		clean_eq = []
		for z in x:
			equations=Eq(ZX[tabix],z)
			clean_eq.append(equations)
			tabix = tabix + 1
		clean_D.append(clean_eq)
	tabix = 0
	for x in clean_D:
		platex00('Difeomorfismo '+str(tabix)+':',file_name,'plain')
		platex00(Matrix(x),file_name,'equation')
		platex00(determinants[tabix],file_name,'equation')
		platex00(Matrix(integrability[tabix]),file_name,'equation')
		tabix = tabix + 1
	
	platexend00(file_name)

	

def set_coordinates(n):
	Z = []
	B = []
	
	
	n_index = range(1,n+1)
	
	B.append(Symbol('F'))
	for x in n_index:
	   E = 'X_' + str(x)
	   B.append(Symbol(E))
	   
	for x in n_index:
		I = 'Y_' + str(x)
		B.append(Symbol(I))
		
	Z.append(Symbol('\Phi'))
	for x in n_index:
	   X = Symbol('E_' + str(x))
	   Z.append(X)
	   
	for x in n_index:
		Y = Symbol('I_' + str(x))
		Z.append(Y)
						
	return Z,B

def set_constants_a(a,letter,value):
	A = []
	c_list = []
	
	for x in range(0,a):
		coefficient = Symbol(str(letter)+'_'+str(x+1))
		if value != 's':
			coefficient = value
		A.append(coefficient)
		c_list.append(coefficient)
			
	return A,c_list
	
def set_constants_axb(a,b,letter,value):
	A = []
	c_list = []
	
	for x in range(0,a):
		B = []
		for y in range(0,b):
			coefficient = Symbol(str(letter)+'_'+str(x+1)+str(y+1))
			if value != 's':
				coefficient = value
			B.append(coefficient)
			c_list.append(coefficient)
		A.append(B)
		
	return A,c_list
	
def set_constants_axbxc00(a,b,c,letter,value):
	A = []
	c_list = []
	
	for x in range(0,a):
		B = zeros(b,c)
		for y in range(0,b):
			for z in range(0,c):
				if 0==0:
					B[y,z] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1))
				if value != 's':
					B[y,z] = value
				c_list.append(B[y,z])
			
		A.append(B)
		
	return A,c_list


	A = []
	c_list = []
	
	for x in range(0,a):
		B = zeros(b,c)
		for y in range(0,b):
			for z in range(0,c):
				if y!=z:
					B[y,z] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1))
				if value != 's':
					B[y,z] = value
				c_list.append(B[y,z])
			
		A.append(B)
		
	return A,c_list
	
def set_constants_axbxc(a,b,c,letter,value,oper):
	A = []
	c_list = []
	
	for x in range(0,a):
		B = zeros(b,c)
		for y in range(0,b):
			for z in range(0,c):
				if oper(y,z):
					B[y,z] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1))
				if value != 's':
					B[y,z] = value
				c_list.append(B[y,z])
			
		A.append(B)
		
	return A,c_list
	
def set_constants_axbxc00(a,b,c,letter,value):
	A = []
	c_list = []
	
	for x in range(0,a):
		B = zeros(b,c)
		for y in range(0,b):
			for z in range(0,c):
				B[y,z] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1))
				if value != 's':
					B[y,z] = value
				c_list.append(B[y,z])
			
		A.append(B)
		
	return A,c_list

def set_constants_axbxcxd(a,b,c,d,letter,value):
	A = []
	c_list = []
	C = []
	for x in range(0,a):
		A = []
		for y in range(0,b):
			B = zeros(c,d)
			for z in range(0,c):
				for w in range(0,d):
					B[z,w] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1)+str(w+1))
					if value != 's':
						B[z,w] = value
					c_list.append(B[z,w])
			A.append(B)
		C.append(A)
	return C
		
def set_constants_axbxcxdxe(a,b,c,d,e,letter,value):
	
	c_list = []
	A = []
	for x in range(0,a):
		B = []
		for y in range(0,b):
			C = []
			for z in range(0,c):
				D = zeros(c,d)
				for w in range(0,d):
					for r in range(0,e):
						D[w,r] = Symbol(str(letter)+'_'+str(x+1)+str(y+1)+str(z+1)+str(w+1)+str(r+1))
						if value != 's':
							D[w,r] = value
						c_list.append(D[w,r])
				C.append(D)
			B.append(C)
		A.append(B)
		
	return A,c_list
	
def linear_solution(A,B,C,D,G,H,n):
	EQ = []
	PM = zeros(2*n+1)
	
	e_index = range(1,n+1)
	i_index = range(n+1,2*n+1)
	
	oper = 0
	for x in range(0,2*n+1):
		for y in range(x,2*n+1):
			if x == 0 and y == 0:
				oper = wise_mul(A,D,[],n)
				EQ.append(oper)
			if x == 0 and y in e_index:
				b = y - 1
				oper = wise_mul(A,G,[b],n)
				EQ.append(oper)
			if x == 0 and y in i_index:
				b = y - 1 - n
				oper = wise_mul(A,H,[b],n)
				EQ.append(oper)
			if x in e_index and y in e_index:
				a = x - 1
				b = y - 1
				oper = wise_mul(B,G,[a,b],n)
				EQ.append(oper)
			if x in e_index and y in i_index:
				a = x - 1
				b = y - 1 - n
				
				oper = wise_mul(B,H,[a,b],n)#-wise_mul(C,G,[a,b],n)
				if a != b: EQ.append(oper)
			if x in i_index and y in i_index:
				a = x - 1 - n
				b = y - 1 - n
				oper = wise_mul(C,H,[a,b],n)
				EQ.append(oper)
			PM[x,y] = oper
		
	return EQ, PM
	
def matrix_gtdielineal(B,C,G,H,ZE):

	EQ = []
	
	n = size(ZE)
	m = int((n - 1)/2)
	m_index =  range (0,m)
	
	h = zeros(m)
	f = zeros(m)
	
	for a in m_index:
		for b in m_index:
			h[a,b] = wise_mul(B,H,[a,b],m)+wise_mul(C,G,[a,b],m)
			#if a!=b:
				#EQ.append(h[a,b])
			if a==b:
				f[a,b] = wise_mul(B,H,[a,b],m)-wise_mul(C,G,[a,b],m)
	

	return f,h,EQ
	
def matrix_gtdiemixed(B,C,H,G,M,N,ZE):

	EQ = []
	eq = []
	
	m = size(ZE)
	n = int((m - 1)/2)
	n_index = range (0,n)
	
	h = zeros(n)
	f = zeros(n)
	
	
	for a in n_index:
		for b in n_index:
			if a == b:
				o1 = wise_mul(B,H,[a,b],n)
				o2 = wise_mul(C,G,[b,a],n)
				h[a,b] = h[a,b] + o1 + o2
				f[a,b] = f[a,b] + o1 - o2
				EQ.append(o1+o2)
				for z in range(0,n):
					o1 = wise_mul02(B,N,[a,z,b],n)
					o2 = wise_mul02(G,M,[a,z,b],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+1]
					EQ.append(o1+o2)
					o1 = wise_mul02(H,M,[b,a,z],n)
					o2 = wise_mul02(C,N,[b,a,z],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+n+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+n+1]
					EQ.append(o1+o2)
					for w in range(0,n):
						o1 = wise_mul02(M,N,[a,w,z,b],n)
						o2 = wise_mul02(M,N,[z,b,a,w],n)
						h[a,b] = h[a,b] + (o1+o2)*ZE[z+1]*ZE[w+n+1]
						f[a,b] = f[a,b] + (o1-o2)*ZE[z+1]*ZE[w+n+1]
						if z == w: 
							EQ.append(o1+o2-1) 
						else: 
							EQ.append(o1+o2)
		eq.append(EQ)
		EQ = []
	EQ = eq
	return f,h,EQ
		
def wise_mul(X,Y,index,n):
	oper = 0
	
	if len(index) == 0:
		for x in range(0,n):
			oper = oper + X[x]*Y[x]
	elif len(index) == 1:
		for x in range(0,n):
			oper = oper + X[x]*Y[x][index[0]]
	elif len(index) == 2:
		for x in range(0,n):
			oper = oper + X[x][index[0]]*Y[x][index[1]]
		

	return oper
	
def wise_mul00(X,Y,index,n):
	oper = 0
	
	if len(index) == 0:
		for x in range(0,n):
			oper = oper + X[x]*Y[x]
	elif len(index) == 1:
		for x in range(0,n):
			oper = oper + X[x]*Y[index[0]][x]
	elif len(index) == 2:
		for x in range(0,n):
			oper = oper + X[x][index[0]]*Y[x][index[1]]
		

	return oper
	
def wise_mulxx(X,Y,Z,index,n,power,option):
	oper = 0
	PX = []
	PY = []
	if option == 0:
		for x in range(0,n):
			for y in range(0,n):
				for z in range(0,n):
					oper = oper + pow(X[x],power[0])*pow(Y[y],2-power[1])*Y[z]*Z[index[0]][x][y,z]
	if option == 1:
		for x in range(0,n):
			oper = oper + X[x]*Z[index[0]][x]
	if option == 2:
		for x in range(0,n):
			oper = oper + Y[x]*Z[index[0]][x]
	if option == 3:
		for x in range(0,n):
			for y in range(0,n):
				oper = oper + pow(X[x],2-power[0])*pow(Y[y],power[1])*Z[index[0]][x,y]
	if option == 4:
		for x in range(0,n):
			for y in range(0,n):
				for z in range(0,n):
					for w in range(0,n):
						PE = pow(X[x],power[0])*pow(X[y],power[1])*pow(Y[z],(1-power[2]))*pow(Y[w],(1-power[3]))
						if PE not in PX:
							oper = oper + PE*Z[index[0]][x][y][z,w]
							PX.append(PE)
	if option == 5:
		for x in range(0,n):
			for y in range(0,n):
				for z in range(0,n):
					for w in range(0,n):
						PE = pow(X[x],(1-power[0]))*pow(X[y],(1-power[1]))*pow(Y[z],power[2])*pow(Y[w],power[3])
						if PE not in PX:
							oper = oper + PE*Z[index[0]][x][y][z,w]
							PX.append(PE)
	return oper
	

		
		

	return oper
	
def wise_mul01(E,I,Z,index,n):
	oper = 0
	
	if len(index) == 0:
		pass
	elif len(index) == 2:
		pass
	elif len(index) == 1:
		for x in range(0,n):
			for y in range(0,n):
				oper = oper + E[x]*I[y]*Z[index[0]][x,y]
	elif len(index)==3:
		for x in range(0,n):
			for y in range(0,n):
				for z in range(0,n):
					oper = oper + E[x]*I[y]*I[z]*Z[index[0]][x][y,z]
		

	return oper
	
def wise_mul03(E,I,Z,index,n):
	oper = 0
	
	if len(index) == 0:
		pass
	elif len(index) == 2:
		pass
	elif len(index) == 3:
		for x in range(0,n):
			for y in range(0,n):
				oper = oper + E[x]*I[y]*Z[index[0]][x,y]
	elif len(index) == 1:
		for x in range(0,n):
			for y in range(0,n):
				for z in range(0,n):
					oper = oper + E[x]*I[y]*I[z]*Z[index[0]][x][y,z]
		

	return oper
	
def wise_mul05(E,I,Z,index,n):
	oper = 0
	
	if len(index) == 0:
		pass
	elif len(index) == 2:
		pass
	elif len(index) == 1:
		for x in range(0,n):
			for y in range(0,n):
				oper = oper + E[y]*I[x]*Z[index[0]][x,y]
		

	return oper
	
def wise_mul02(X,Y,index,n):
	oper = 0
	
	if len(index) == 0:
		pass
	elif len(index) == 2:
		for x in range(0,n):
			oper = oper + X[x]*Y[x][index[0],index[1]]
	elif len(index) == 3:
		for x in range(0,n):
			oper = oper + X[x][index[0]]*Y[x][index[1],index[2]]
	elif len(index) == 4:
		for x in range(0,n):
			oper = oper + X[x][index[0],index[1]]*Y[x][index[2],index[3]]
	elif len(index) == 1:
		for x in range(0,n):
			for y in range(0,n):
				oper = oper + E[y]*I[x]*Z[index[0]][x,y]

		

	return oper

def Jacobian00(Z0,ZE,n):
	G = zeros(2*n+1)
	for x in range(0,2*n+1):
		for y in range(0,2*n+1):
			if x == 0:
				G[x,y] = Z0[x][y]
			else:
				G[x,y] = diff(Z0[x],ZE[y])
			
	d = G.det()
	return G,d
	
def mixed_solution(A,B,C,D,G,H,M,N,n):
	EQ = []
	eq = []
	PM = zeros(2*n+1)
	
	e_index = range(1,n+1)
	i_index = range(n+1,2*n+1)
	
	oper = 0
	for x in range(0,2*n+1):
		for y in range(x,2*n+1):
			EQ = []
			if x == 0 and y == 0:
				oper = wise_mul(A,D,[],n)
				EQ.append(oper)
			if x == 0 and y in e_index:
				b = y - 1
				oper = wise_mul(A,G,[b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(A,N,[b,z],n)
					EQ.append(oper)
			if x == 0 and y in i_index:
				b = y - 1 - n
				oper = wise_mul(A,H,[b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(A,N,[z,b],n)
					EQ.append(oper)
			if x in e_index and y in e_index:
				a = x - 1
				b = y - 1
				oper = wise_mul(B,G,[a,b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(B,N,[a,b,z],n)
					EQ.append(oper)
					oper = wise_mul02(G,M,[b,a,z],n)
					EQ.append(oper)
					for w in range(0,n):
						oper = wise_mul02(M,N,[a,z,b,w],n)
						EQ.append(oper)
			if x in e_index and y in i_index:
				a = x - 1
				b = y - 1 - n
				oper = wise_mul(B,H,[a,b],n)#-wise_mul(C,G,[b,a],n)
				if a != b: EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(B,N,[a,z,b],n)#-wise_mul02(G,M,[a,z,b],n)
					if a != b: EQ.append(oper)
					oper = wise_mul02(H,M,[b,a,z],n)#-wise_mul02(C,N,[b,a,z],n)
					if a != b: EQ.append(oper)
					for w in range(0,n):
						oper = wise_mul02(M,N,[a,w,z,b],n)#-wise_mul02(M,N,[z,b,a,w],n)
						if a != b: 
							#EQ.append((x,y))
							EQ.append(oper)
			if x in i_index and y in i_index:
				a = x - 1 - n
				b = y - 1 - n
				oper = wise_mul(C,H,[a,b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(C,N,[a,z,b],n)
					EQ.append(oper)
					oper = wise_mul02(H,M,[b,z,a],n)
					EQ.append(oper)
					for w in range(0,n):
						oper = wise_mul02(M,N,[w,a,z,b],n)
						EQ.append(oper)
			#PM[x,y] = EQ
			eq = eq + EQ
		
	return eq, PM
	
def mixed_difeomorphism(A,B,C,D,G,H,M,N,n,ZX,ZE):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	for x in range(0,n):
		a = A[x]*ZE[0]+wise_mul00(ZE[1:n+1],B,[x],n)+wise_mul00(ZE[n+1:2*n+1],C,[x],n)+wise_mul01(ZE[1:n+1],ZE[n+1:2*n+1],M,[x],n)
		difeo.append(a)
		
	for x in range(0,n):
		a = D[x]*ZE[0]+wise_mul00(ZE[1:n+1],G,[x],n)+wise_mul00(ZE[n+1:2*n+1],H,[x],n)+wise_mul01(ZE[1:n+1],ZE[n+1:2*n+1],N,[x],n)
		difeo.append(a)
		
	return difeo
	
def GI_difeomorphism(A,B,C,D,G,H,M,N,J,K,n,ZX,ZE):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	for x in range(0,n):
		a = A[x]*ZE[0]+wise_mul00(ZE[1:n+1],B,[x],n)+wise_mul00(ZE[n+1:2*n+1],C,[x],n)+wise_mul01(ZE[1:n+1],ZE[n+1:2*n+1],M,[x],n)+wise_mul01(ZE[1:n+1],ZE[1:n+1],J,[x],n)
		difeo.append(a)
		
	for x in range(0,n):
		a = D[x]*ZE[0]+wise_mul00(ZE[1:n+1],G,[x],n)+wise_mul00(ZE[n+1:2*n+1],H,[x],n)+wise_mul05(ZE[1:n+1],ZE[n+1:2*n+1],N,[x],n)+wise_mul01(ZE[n+1:2*n+1],ZE[n+1:2*n+1],K,[x],n)
		difeo.append(a)
		
	return difeo
	
def pow_difeomorphism(cof,n,ZX,ZE,power):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	for x in range(0,n):
		a = 0
		for y in range(0,len(cof[0])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[0][y],[x],n,power[0][y],4)
		for y in range(0,len(cof[2])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[2][y],[x],n,power[0][y],1)
		for y in range(0,len(cof[4])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[4][y],[x],n,power[0][y],2)
		difeo.append(a)
		
	for x in range(0,n):
		a = 0
		for y in range(0,len(cof[1])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[1][y],[x],n,power[1][y],5)
		for y in range(0,len(cof[3])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[3][y],[x],n,power[0][y],2)
		for y in range(0,len(cof[5])):
			a = a + wise_mulxx(ZE[1:n+1],ZE[n+1:2*n+1],cof[5][y],[x],n,power[0][y],1)
		difeo.append(a)
		
	return difeo
	
def linear_difeomorphism(A,B,C,D,G,H,n,ZX,ZE):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	for x in range(0,n):
		a = A[x]*ZE[0]+wise_mul00(ZE[1:n+1],B,[x],n)+wise_mul00(ZE[n+1:2*n+1],C,[x],n)
		difeo.append(a)
		
	for x in range(0,n):
		a = D[x]*ZE[0]+wise_mul00(ZE[1:n+1],G,[x],n)+wise_mul00(ZE[n+1:2*n+1],H,[x],n)
		difeo.append(a)
		
	return difeo
	
def get_replacement00(sol,case):
	replacement = []
	
	for x in sol[1]:
		tabix = 0
		for y in x:
			if y == case and case == case:
				replacement.append((sol[0][tabix],y))
			elif y != 0 and case == 'n0':
				replacement.append((sol[0][tabix],y))
			tabix = tabix+1
	return replacement
	

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))
	
def one_set_solutions(p_set,EQ):
	all_solutions = set()

	for x in p_set:
		sol = solve(EQ,x,set= True)
		if sol != []:
			replacementno = get_replacement00(sol,'n0')
			replacement0 = get_replacement00(sol,0)
			all_solutions.add(tuple(replacementno+replacement0))
	l = list(all_solutions)
	return l
	
def one_set_solutions00(EQ):
	all_solutions = set()
	sol = solve(EQ,set= True)
	if sol != []:
		replacementno = get_replacement00(sol,'n0')
		replacement0 = get_replacement00(sol,0)
		all_solutions.add(tuple(replacementno+replacement0))
	l = list(all_solutions)
	return l

def subs_solutions(LD,subs):
	LR = []
	for x in LD:
		LR.append(x.subs(subs))
	return LR
	
def subs_solutions00(M,subs,n):
	N = zeros(n)
	for x in range(0,n):
		for y in range (0,n):
			N[x,y] = M[x,y].subs(subs)
		
	return N
	
def assert_jacobian(LR,ZE):
	return 1 if Matrix(LR).jacobian(ZE).det() == 0 else 0
	
def assert_matrix(M,n):
	Zero = zeros(n)
	if M == Zero:
		return 1
	else:
		return 0

def seek_mixed(LD,h_LD,f_LD,EQ_LD,D_LD,s,n,ZE,toler,max_sol):
	s_SD = []
	LR = LD[0]
	h = h_LD[0]
	f = f_LD[0]
	EQ = EQ_LD[0]
	LD = []
	h_LD=[]
	f_LD=[]
	EQ_LD = []
	LD.append(LR)
	h_LD.append(h)
	f_LD.append(f)
	EQ_LD.append(EQ)
	deti = D_LD[0]
	tabix = 0
	sabix = 0
	zero_list = [0]*(2*n+1)
	for z in s:
		if tabix > toler:
			break
		tabix = tabix + 1
			
		EQ = EQ_LD[0]
		LR = LD[0]
		h = h_LD[0]
		f = f_LD[0]
		det_subs = deti
		for x in z:
			k = det_subs.subs(x)
			det_subs = k
		det_subs = simplify(det_subs)
		if det_subs != 0:
			for y in z:
				EQ = subs_solutions(EQ,y)
				#h = subs_solutions00(h,x,n)
				#f = subs_solutions00(f,x,n)
				LR = subs_solutions(LR,y)
			LD.append(LR)
			h_LD.append(h)
			f_LD.append(f)
			EQ_LD.append(EQ)
			D_LD.append(det_subs)
			s_SD.append(z)
			if sabix > max_sol:
				break
			sabix = sabix + 1
			
	return LD, h_LD, f_LD,EQ_LD,D_LD,s_SD
	
def mixed_solution00(A,B,C,D,G,H,M,N,n):
	EQ = []
	eq = []
	PM = zeros(2*n+1)
	
	e_index = range(1,n+1)
	i_index = range(n+1,2*n+1)
	
	oper = 0
	for x in range(0,2*n+1):
		for y in range(x,2*n+1):
			EQ = []
			if x == 0 and y == 0:
				oper = wise_mul(A,D,[],n)
				EQ.append(oper)
			if x == 0 and y in e_index:
				b = y - 1
				oper = wise_mul(A,G,[b],n)
				EQ.append(oper)
			if x == 0 and y in i_index:
				b = y - 1 - n
				oper = wise_mul(A,H,[b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(A,N,[z,b],n)
					oper1 = wise_mul02(A,N,[b,z],n)
					EQ.append(oper+oper1)
			if x in e_index and y in e_index:
				a = x - 1
				b = y - 1
				oper = wise_mul(B,G,[a,b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(G,M,[b,a,z],n)
					oper1 = wise_mul02(G,M,[b,z,a],n)
					EQ.append(oper+oper1)
			if x in e_index and y in i_index:
				a = x - 1
				b = y - 1 - n
				oper = wise_mul(B,H,[a,b],n)
				if a != b: EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(B,N,[a,z,b],n)
					oper1 = wise_mul02(B,N,[a,b,z],n)
					if a != b: EQ.append(oper+oper1)
					oper = wise_mul02(H,M,[b,a,z],n)
					oper1 = wise_mul02(H,M,[b,z,a],n)
					if a != b: EQ.append(oper+oper1)
					for w in range(0,n):
						oper = wise_mul02(M,N,[a,z,b,w],n)
						oper1 = wise_mul02(M,N,[a,z,w,b],n)
						oper2 = wise_mul02(M,N,[z,a,b,w],n)
						oper3 = wise_mul02(M,N,[z,a,w,b],n)
						if a != b:
							EQ.append(oper+oper2+oper3+oper1)
			if x in i_index and y in i_index:
				a = x - 1 - n
				b = y - 1 - n
				oper = wise_mul(C,H,[a,b],n)
				EQ.append(oper)
				for z in range(0,n):
					oper = wise_mul02(C,N,[a,z,b],n)
					EQ.append(oper)

			#PM[x,y] = EQ
			eq = eq + EQ
		
	return eq, PM
	
def matrix_gtdiemixed00(B,C,H,G,M,N,ZE):

	EQ = []
	eq = []
	m = size(ZE)
	n = int((m - 1)/2)
	n_index = range (0,n)
	
	h = zeros(n)
	f = zeros(n)
	
	
	for a in n_index:
		for b in n_index:
			if a == b:
				o1 = wise_mul(B,H,[a,b],n)
				o2 = wise_mul(C,G,[b,a],n)
				h[a,b] = h[a,b] + o1 + o2
				f[a,b] = f[a,b] + o1 - o2
				EQ.append(o1+o2)
				for z in range(0,n):
					o1 = wise_mul02(B,N,[a,z,b],n)
					o2 = wise_mul02(B,N,[a,b,z],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+n+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+n+1]
					EQ.append(o1+o2)
					o1 = wise_mul02(H,M,[b,a,z],n)
					o2 = wise_mul02(H,M,[b,z,a],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+1]
					EQ.append(o1+o2)
					for w in range(0,n):
						o1 = wise_mul02(M,N,[a,z,b,w],n)
						o2 = wise_mul02(M,N,[a,z,w,b],n)
						o3 = wise_mul02(M,N,[z,a,b,w],n)
						o4 = wise_mul02(M,N,[z,a,w,b],n)
						h[a,b] = h[a,b] + (o1+o2)*ZE[z+n+1]*ZE[w+1]
						f[a,b] = f[a,b] + (o1-o2)*ZE[z+n+1]*ZE[w+1]
						if z == w: 
							EQ.append(ZE[z+n+1]*ZE[w+1])
							EQ.append(o1+o2+o3+o4-1)
						else:
							EQ.append(ZE[z+n+1]*ZE[w+1])
							EQ.append(o1+o2+o3+o4)
		eq.append(EQ)
		EQ = []
	EQ = eq	
	return f,h,EQ
	
def matrix_gtdieGIII(B,C,H,G,M,N,ZE):

	EQ = []
	eq = []
	m = size(ZE)
	n = int((m - 1)/2)
	n_index = range (0,n)
	
	h = zeros(n)
	f = zeros(n)
	
	
	for a in n_index:
		for b in n_index:
			if a == b:
				o1 = wise_mul(B,H,[a,b],n)
				o2 = wise_mul(C,G,[b,a],n)
				h[a,b] = h[a,b] + o1 + o2
				f[a,b] = f[a,b] + o1 - o2
				EQ.append(o1+o2)
				for z in range(0,n):
					o1 = wise_mul02(B,N,[a,z,b],n)
					o2 = wise_mul02(B,N,[a,b,z],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+n+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+n+1]
					EQ.append(o1+o2)
					o1 = wise_mul02(H,M,[b,a,z],n)
					o2 = wise_mul02(H,M,[b,z,a],n)
					h[a,b] = h[a,b] + (o1+o2)*ZE[z+1]
					f[a,b] = f[a,b] + (o1-o2)*ZE[z+1]
					EQ.append(o1+o2)
					for w in range(0,n):
						o1 = wise_mul02(M,N,[a,z,b,w],n)
						o2 = wise_mul02(M,N,[a,z,w,b],n)
						o3 = wise_mul02(M,N,[z,a,b,w],n)
						o4 = wise_mul02(M,N,[z,a,w,b],n)
						h[a,b] = h[a,b] + (o1+o2)*ZE[z+n+1]*ZE[w+1]
						h[a,b] = h[a,b]
						f[a,b] = f[a,b] + (o1-o2)*ZE[z+n+1]*ZE[w+1]
						if z == w: 
							EQ.append(ZE[z+n+1]*ZE[w+1])
							EQ.append(o1+o2+o3+o4-1)
						else:
							EQ.append(ZE[z+n+1]*ZE[w+1])
							EQ.append(o1+o2+o3+o4)
		eq.append(EQ)
		EQ = []
	EQ = eq	
	return f,h,EQ
	
def mixed_difeomorphism00(A,B,C,D,G,H,M,N,n,ZX,ZE):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	for x in range(0,n):
		a = A[x]*ZE[0]+wise_mul00(ZE[1:n+1],B,[x],n)+wise_mul00(ZE[n+1:2*n+1],C,[x],n)+wise_mul01(ZE[1:n+1],ZE[1:n+1],M,[x],n)
		difeo.append(a)
		
	for x in range(0,n):
		a = D[x]*ZE[0]+wise_mul00(ZE[1:n+1],G,[x],n)+wise_mul00(ZE[n+1:2*n+1],H,[x],n)+wise_mul01(ZE[n+1:2*n+1],ZE[n+1:2*n+1],N,[x],n)
		difeo.append(a)
		
	return difeo
	
def GIII_difeomorphism(A,B,C,D,G,H,M,N,n,ZX,ZE):
	difeo = []
	
	a= 0
	
	a = Function('F')(*ZE)
	difeo.append(a)
	
	k = Symbol('k')
	b = k+1
	ZE1 = []
	
	tabix = 0
	for x in ZE:
		if tabix == 0:
			ZE1.append(x)
		else:
			ZE1.append(pow(x,b)/pow(2*b,1/2))
		tabix = tabix + 1
	
	for x in range(0,n):
		a = A[x]*ZE1[0]+wise_mul00(ZE1[1:n+1],B,[x],n)+wise_mul00(ZE1[n+1:2*n+1],C,[x],n)+wise_mul01(ZE1[1:n+1],ZE1[1:n+1],M,[x],n)
		difeo.append(a)
		
	for x in range(0,n):
		a = D[x]*ZE1[0]+wise_mul00(ZE1[1:n+1],G,[x],n)+wise_mul00(ZE1[n+1:2*n+1],H,[x],n)+wise_mul01(ZE1[n+1:2*n+1],ZE1[n+1:2*n+1],N,[x],n)
		difeo.append(a)
		
	return difeo
	
def seek_mixed00(LD,h_LD,f_LD,EQ_LD,D_LD,s,n,ZE,toler,max_sol,postpro):
	s_LD = []
	LR = LD[0]
	h = h_LD[0]
	f = f_LD[0]
	EQ = EQ_LD[0]
	LD = []
	h_LD=[]
	f_LD=[]
	EQ_LD = []
	LD.append(LR)
	h_LD.append(h)
	f_LD.append(f)
	EQ_LD.append(EQ)
	deti = D_LD[0]
	tabix = 0
	sabix = 0
	
	print ('Comienza substitucion')
	for z in s:
		merged_subs = merge_dicts(*z)
		if tabix >= toler:
			break
		tabix = tabix + 1
		EQ = EQ_LD[0]
		LR = LD[0]
		#h = h_LD[0]
		#f = f_LD[0]
		det_subs = deti
		for x in z:
			det_subs = det_subs.subs(x)
		
		
		if det_subs == 0:
			continue
			
		for x in z:
			LR = subs_solutions(LR,x)
			
			#EQ = subs_solutions(EQ,x)
		if LR not in LD: #and det_subs != 0:
			LD.append(LR)
			#	h_LD.append(h)
			#	f_LD.append(f)
			EQ_LD.append(EQ)
			s_LD.append(z)#(merged_subs)
			D_LD.append(det_subs)
			if sabix >= max_sol:
				break
			sabix = sabix + 1
		else:
			tabix = tabix - 1
	
	print ('Finaliza substitucion')		
	if postpro[0] == 1:
		print ('Comienza limpieza')
		dmy_LD = []
		dmy_D = []
		tabix = 0
		for x in LD:
			deta = D_LD[tabix]
			if 0 not in x and deta !=0:
				dmy_LD.append(x)
				dmy_D.append(D_LD[tabix])
			tabix = tabix + 1
		LD = dmy_LD
		D_LD = dmy_D
		print ('Finaliza substitucion')
			
	return LD, h_LD, f_LD,EQ_LD,D_LD,s_LD
	
def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def consistency_subs(AE):
    
    EQ = list(AE)
    tabix = 0
    for x in AE:
        rabix = 0
        for y in AE:
            if rabix != tabix: 
                f = y.subs([(x,0)])
                AE[rabix] = f
            rabix = rabix + 1
        tabix = tabix +1
        
    return EQ
	
def EQ_sub(EQ,sub_EQ):
	C_EQ = []
	for x in EQ:
		step = x
		for y in sub_EQ:
			step = step.subs([(y,0)])
		C_EQ.append(step)
		
	return C_EQ
	
def send_DH_to_latex(equations,determinants,integrability,equations_matrix,ZX,ini,file_name):
    to_latex = []
    clean_D = []
    
    platexini00(file_name)
    for x in equations:
        tabix = 0
        clean_eq = []
        for z in x:
            equations=Eq(ZX[tabix],z)
            clean_eq.append(equations)
            tabix = tabix + 1
        clean_D.append(clean_eq)
    tabix = 0
    for x in clean_D:
        platex00('Difeomorfismo '+str(ini + tabix)+':',file_name,'plain')
        platex00(Matrix(x),file_name,'equation')
        platex00(determinants[tabix],file_name,'equation')
        platex00(Matrix(integrability[tabix]),file_name,'equation') 
        platex00(Matrix(equations_matrix[tabix]),file_name,'equation')
        tabix = tabix + 1

    platexend00(file_name)
	
def solve_part(equation,sim,exclude):
    ex = []
    for x in exclude:
        sx = list(x.free_symbols)
        for y in sx:
            if y not in ex:
                ex.append(y)
    all_sol = {}
    ambar = list(equation.values())
    for z in ambar:
                fs = list(z.free_symbols)
                for x in ex:
                    if x in fs: fs.remove(x)
                for w in fs:
                    sol = solve(z,w,dict=True)
                    if len(sol) == 0: break
						
                    if list(sol[0].keys())[0] not in all_sol.keys():
                        #if list(sol[0].values())[0] != 0:
                        all_sol.update(sol[0])
                        break
    soli = {}
    tabix = 0
    for x in all_sol:
        soli[x] = sim(all_sol[x]) 
        tabix = tabix + 1
        #if tabix >= 1:
        #    break
    return soli
	
def solve_pair(equation1,equation2,sim):
    
    equation = {}
    for x in equation1:
        if x in list(equation2.keys()):
            equation[x] = equation1[x] - equation2[x]
                

    all_sol = {}
    ambar = list(equation.values())
    for z in ambar:
        fs = list(z.free_symbols)
        for w in fs:
            sol = solve(z,w,dict=True)
            if list(sol[0].keys())[0] not in all_sol.keys():
                all_sol.update(sol[0])
                #break
    soli = {}
    tabix = 0
    for x in all_sol:
        soli[x] = sim(all_sol[x]) 
        tabix = tabix + 1
        #if tabix >= 1:
        #    break
    return soli
	

	
def clean_system00(EQS,n):
    EQC = EQS[:,:]
    for x in range(0,2*n+1):
        for y in range(x,2*n+1):
            if EQC[x,y] != {} and EQC[x,y] != 0:
                for w in EQC[x,y]:
                    s = factor(EQC[x,y][w])
                    fac = 1
                    for r in s.args:
                        if type(r) == Integer:
                            fac = r
                            break
                    EQC[x,y][w] = expand(EQC[x,y][w]/fac)
                    
                for w in EQC[x,y]:
                    deli = {}
                    deli[EQC[x,y][w]] = 0
                    for r in EQC[x,y]:
                        if r != w:
                            EQC[x,y][r] = EQC[x,y][r].subs(deli)
                                
                inv_map = {v: k for k, v in EQC[x,y].items()}
                inv_map = dict.fromkeys(inv_map, 0)
                for a in range(0,2*n+1):
                    for b in range(a,2*n+1):
                        if a == x and b == y: continue                               
                        for z in EQC[a,b]:
                            EQC[a,b][z] = EQC[a,b][z].subs(inv_map)
    return EQC
	
def seek_cohe(check,m,n,price,sim):
	lars = set()
	for x in range(-4,4):
		fito = {}
		fito[n] = x
		ew = check
		for w in ew:
			y = powsimp(w.subs(fito))
			if y == price:
				z = sim(ew[w])
				lars.add((x,m,(w,z)))

	for x in range(-4,4):
		fito = {}
		fito[m] = x
		ew = check
		for w in ew:
			y = powsimp(w.subs(fito))
			if y == price:
				z = sim(ew[w])
				lars.add((n,x,(w,z)))
	
	fito = {}
	fito[n] = 0
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((0,0,(x,z)))
	
	fito = {}
	fito[n] = 1
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((1,0,(x,z)))
				
	fito = {}
	fito[n] = 2
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((2,0,(x,z)))
				
	fito = {}
	fito[n] = 3
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((3,0,(x,z)))
				
	fito = {}
	fito[n] = 4
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((4,0,(x,z)))
				
	fito = {}
	fito[n] = -1
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((-1,0,(x,z)))
				
	fito = {}
	fito[n] = -2
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((-2,0,(x,z)))
				
	fito = {}
	fito[n] = -3
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((-3,0,(x,z)))
				
	fito = {}
	fito[n] = -4
	fito[m] = 0
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((-4,0,(x,z)))
	
	fito = {}
	fito[n] = -1
	fito[m] = 1
	ew = check
	for x in ew:
			y = powsimp(x.subs(fito))
			if y == price:
				z = sim(ew[x])
				lars.add((-1,1,(x,z)))

	
	return lars

