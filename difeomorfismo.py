from numpy import *
from pylab import *
#from mpmath import *
from sympy import Matrix
from sympy import *
from sympy.abc import beta
from itertools import *
import time
import gc


def platexini(file_name):
    filename = file_name
    FILE = open(filename,'w')                   
    FILE.writelines('\\begin{document}')
    FILE.close()

def platexend(file_name):
    filename = file_name
    FILE = open(filename,'a')                   
    FILE.writelines('\\end{document}')
    FILE.close()
    
def platex(expresion,file_name):
    filename = file_name
    FILE = open(filename,'a')                   
    FILE.writelines(latex(expresion,mode='equation'))
    FILE.close()

def n_grados(n):
	A = Symbol('\Lambda')
	k = Symbol('k')
	Z = []
	B = []
	
	B.append(Symbol('\Phi '))
	n_index = range(1,n+1)
	
	
	for x in n_index:
	   E = 'E_' + str(x)
	   B.append(Symbol(E))
	   
	for x in n_index:
		I = 'I_' + str(x)
		B.append(Symbol(I))
		
	Z.append(B[0]+Function('F')(*B[1:2*n+1]))
	for x in n_index:
	   X = Function('X_' + str(x))(*B)
	   Z.append(X)
	   
	for x in n_index:
		Y = Function('Y_' + str(x))(*B)
		Z.append(Y)
						
	return Z,B
	

def n_grados001(n):
	A = Symbol('\Lambda')
	k = Symbol('k')
	Z = []
	B = []
	
	B.append(Symbol('F'))
	n_index = range(1,n+1)
	
	
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
	
def S_component(a,b,Z,B):
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	S = 0
	
	for d in e_index:
		S = S + (diff(Z[d],B[a])*diff(Z[d+m],B[b]))+(diff(Z[d],B[b])*diff(Z[d+m],B[a]))
		
	return S

def S_component001(a,b,Z,B):
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	S = 0
	
	for d in e_index:
		S = S + (diff(Z[d],B[a])*diff(Z[d+m],B[b]))+(diff(Z[d],B[b])*diff(Z[d+m],B[a]))
		
	S = S/2
		
		
	return S
	
def S_component002(a,b,Z,B):
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	S = 0
	
	for d in e_index:
		S = S + expand(diff(Z[d],B[a])*diff(Z[d+m],B[b]))
		
	return S

def S_component003(a,b,Z,B):
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	S = 0
	
	for d in e_index:
		S = S + diff(Z[d],B[a])*diff(Z[d+m],B[b])+ diff(Z[d],B[b])*diff(Z[d+m],B[a])
		
	return S
	
def H_component001(a,Z,B,f,sim):
    
	if f == 0:
		f = Symbol('f')
		
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	i_index = range (m+1,2*m+1)
	
	H2 = diff(Z[0],B[a])
	
	H1 = 0
	for d in e_index:
		H1 = H1 + f[d-1,d-1]*(diff(Z[d],B[a])*Z[d+m])
	
	if a == 0:
		H1 = H1 + 1
	
	if a in e_index:
		H1 = H1 - B[a+m]
		
	H1 = H1/f[0,0]
		

	H = (H2,H1,sim(integrate(H1,B[a])))
	
	
	return H

def H_component003(a,Z,B,f,sim):
    	
	n = size(Z)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
		
	H1 = 0
	for d in e_index:
		H1 = H1 + (diff(Z[d],B[a])*Z[d+m])
	
	if a == 0:
		H1 = H1+1
	
	if a in e_index:
		H1 = H1 - (f[a-1,a-1])*B[a+m]
		

	H = sim(H1)
	
	return H

	
def H_component002(a,f,Z,B):
    
			
	n = size(Z)
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	i_index = range (m+1,2*m+1)
	
	H2 = f[0,0]*diff(Z[0],B[a])
	
	H1 = 0
	for d in e_index:
		H1 = H1 + f[d-1,d-1]*(diff(Z[d],B[a])*Z[d+m])
		
	H = H2-H1
		
	return H
	
def H_component(a,b,f,Z,B):
	H = H_component002(a,f,Z,B)*H_component002(b,f,Z,B)
	return H
	
def difeomorphism_equations(Z,B,sim,k):
	
	if k == 0:
		k = Symbol('k')#(Z[0],Z[1],Z[2])
	
	n = size(Z)
	n_index = range(0,n)
	H = zeros(n)
	S = zeros(n)
	simbolo = 0
	for a in n_index:
		for b in n_index:
			simbolo = H_component(a,b,Z,B)
			H[a,b] = k*sim(simbolo)
			simbolo = S_component(a,b,Z,B)
			S[a,b] = sim(simbolo)
	
	return H, S	
def difeomorphism_equations004(Z,B,sim,k,N):
	
	if k == 0:
		k = Symbol('k')#(Z[0],Z[1],Z[2])
	
	start_time = time.time()
	n = size(B)
	n_index = range(0,n)
	H = zeros(n)
	S = zeros(n)
	simbolo = 0
	for a in n_index:
		for b in n_index:
			simbolo = H_component(a,b,Z,B)
			H[a,b] = sim(simbolo)
			simbolo = S_component001(a,b,Z,B,N)
			S[a,b] = sim(simbolo)
	
	print("--- %s seconds ---" % (time.time() - start_time))
	return H, S
	
def compatibility_equation(Z,B,sim,N):
	
	
	#start_time = time.time()
	
	R =[]
	
	n = size(B)
	m = int((n-1)/2)
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	for e in e_index:
		i = e + m
		SFX = S_component001(0,e,Z,B,N)
		SXX = S_component001(e,e,Z,B,N)
		SYY = S_component001(i,i,Z,B,N)
		SXY = S_component001(e,i,Z,B,N)
		SFY = S_component001(0,i,Z,B,N)
		SFF = S_component001(0,0,Z,B,N)
		X = B[e]
		#r = ((X+SFX)*SYY)+(((1/2)-SXY)*SFY)
		r = ((pow(B[e],2)-SXX)*pow(SFY,2))-(1-SFF)*((1/2)-SXY)
		if sim != 0:
			r = sim(r)
		R.append(r)
	
	
	#print("--- %s seconds ---" % (time.time() - start_time))
	return R
	
def difeomorphism_equations005(Z,B,sim,f,h,a_index,b_index):
	#start_time = time.time()
	
	n = size(B)
	H = zeros(n)
	S = zeros(n)
	simbolo = 0
	for a in range(a_index,b_index):
		for b in range(a,b_index):
			simbolo = H_component(a,b,f,Z,B)
			H[a,b] = sim(simbolo)
			simbolo = S_component002(a,b,Z,B)+S_component002(b,a,Z,B)
			#if a+int((n-1)/2) == b and a !=0:
			#	simbolo = simbolo - h[a-1,a-1]
			S[a,b] = sim(simbolo)
			H[b,a] = H[a,b]
			S[b,a] = S[a,b]
			
	
	#print("--- %s seconds ---" % (time.time() - start_time))
	return H, S
	
def difeomorphism_equations001(Z,B,sim,N):
	n = size(Z)
	if k == 0:
		k = Function('k')(*B[1:2*n+1])
	n_index = range(0,n)
	H = zeros(n)
	S = zeros(n)
	simbolo = 0
	for a in n_index:
		for b in n_index:
			simbolo = H_component(a,b,Z,B)
			H[a,b] = k*sim(simbolo)
			simbolo = S_component001(a,b,Z,B)
			S[a,b] = sim(simbolo)
	
	return H, S	
	
def difeomorphism_equations002(Z,B,sim,f,h,a_index,b_index):
	
	n = size(Z)
	if f == 0:
		f = Symbol('f')
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	i_index = range (1+m,2*m+1)
	
	H = []
	S = zeros(n)
	simbolo = 0
	for a in range(a_index,b_index):
		H.append(H_component001(a,Z,B,f,sim))
		for b in range(a,b_index):
			simbolo = S_component001(a,b,Z,B)
			#if a+int((n-1)/2) == b and a !=0:
			#	simbolo = simbolo - h[a-1,a-1]
			
				
			
			S[a,b] = sim(simbolo)
			#S[b,a] = sim(simbolo)
	
	return H, S

def poly_args(x,dc,dZE,sim):
	dEQ = {}
	if type(x) == Add:
		for y in x.args:
			sZ = y.subs(dc)
			l = sZ.subs(dZE)
			simplify(l)
			sZ = powsimp(sZ/l)
			if sZ not in dEQ.keys():
				dEQ[sZ] = 0
			sc = y.subs(dZE)
			dEQ[sZ] = dEQ[sZ] + sc
	elif type(x) == Mul:
		sZ = x.subs(dc)
		l = sZ.subs(dZE)
		sZ = powsimp(sZ/l)
		if sZ not in dEQ.keys():
			dEQ[sZ] = 0
		sc = x.subs(dZE)
		dEQ[sZ] = dEQ[sZ] + sc
		
	return dEQ
	

def clean_factor(equation):
	s = factor(equation)
	fac = 1
	for r in s.args:
		if type(r) == Integer:
			fac = r
			break
	equation = equation/fac
	return equation
	
def get_zero_eq(equation):
	deli = {}
	deli[equation] = 0
	return deli

def get_zero_co(equation,dZE,dc,sim):
	mali = {}
	mali = poly_args(equation,dZE,dc,sim)
	mali = dict.fromkeys(mali,0)
	return mali
def set_zero_by_eq(EQA,deli,mali,n):
	EQD = EQA[:,:]
	for a in range(0,2*n+1):
		for b in range(0,2*n+1):
			if EQD[a,b] != 0 and EQD != {}:
				for z in EQD[a,b]:
					if EQD[a,b][z] != 0:
						clean_counter = len(mali)
						resto = expand(EQD[a,b][z])
						tabix = 0
						r1 = resto
						for x in mali:
							deni = {}
							deni[x] = mali[x]
							r2 = r1
							r1 = r1.subs(deni)
							if r1 != r2:
								tabix = tabix + 1
						if tabix == clean_counter:
							resto = r1
						fijo = factor(EQD[a,b][z] - resto)
						fijo = expand(clean_factor(fijo)).subs(deli)
						EQD[a,b][z] = fijo + resto
	return EQD
	
def clean_system02(EQS,n,dZE,dc,sim,exclude):
	elias = []
	complice = []
	EQC = EQS[:,:]
	for x in range(0,2*n+1):
		for y in range(0,2*n+1):
			if EQC[x,y] != {} and EQC[x,y] != 0:
				for w in EQC[x,y]:
					if (x,y) not in exclude:
						if EQC[x,y][w] != 0:
							equation = expand(clean_factor(EQC[x,y][w]))
							deli = get_zero_eq(equation)
							mali = get_zero_co(equation,dZE,dc,sim)
							EQB = EQC[:,:]
							EQB[x,y][w] = 0
							EQC = set_zero_by_eq(EQB,deli,mali,n)
							elias.append(EQC)
							complice.append((x,y,w,equation))
							EQC[x,y][w] = equation
						
					
	return EQC, elias, complice
	
def poly_difeomorphism(Z,B,dc,dZE,h,sim):
	
	n = size(Z)
	m = int((n-1)/2)
		
	S = zeros(n)
	EQ = zeros(n)
	simbolo = 0
	
	uni = [1]
	for a in range(0,n):
		for b in range(0,n):
			
			simbolo = S_component002(a,b,Z,B)
			if a in range(1,m+1) and b in range(m+1,n):
				u = a-1
				v = b-1-m
				simbolo = simbolo + S_component002(b,a,Z,B)
				if u == v:
				   simbolo = simbolo - h[a-1,a-1]
			if a in range(m+1,n) and b in range(1,m+1):
				u = a-1-m
				v = b-1
				simbolo =  2*(S_component002(b,a,Z,B)-simbolo)
				simbolo = expand(simbolo)
			
			EQ[a,b] = poly_args(simbolo,dc,dZE,sim)
			S[a,b] = simbolo
			
			
			#S[b,a] = sim(simbolo)
	
	return S,EQ



	
def contactomorphism_equations(Z,B,sim,f):
	
	n = size(Z)
	if f == 0:
		f = Symbol('f')
	
	
	H = []
	for a in range(0,n):
		H.append(H_component001(a,Z,B,f,sim)[1])
			
	return H

def difeomorphism_equations003(Z,B,sim,k):
	n = size(Z)
	if k == 0:
		k = Function('k')(*B[1:2*n+1])
	n_index = range(0,n)
	m = int((n - 1)/2)
	e_index = range (1,m+1)
	i_index = range (1+m,2*m+1)
	#print e_index, i_index
	H = zeros(n)
	S = zeros(n)
	simbolo = 0
	for a in n_index:
		for b in n_index:
			if a == 0 or b == 0:
				simbolo = H_component(a,b,Z,B)
				H[a,b] = k*sim(simbolo)
			else:
				if a in e_index and b in e_index:
					simbolo = (B[b+m]+S_component(0,b,Z,B))*(B[a+m]+S_component(0,a,Z,B))/H_component(0,0,Z,B)
					H[a,b] = sim(simbolo)/k
				if a in i_index and b in i_index:
					simbolo = S_component(0,b,Z,B)*S_component(0,a,Z,B)/H_component(0,0,Z,B)
					H[a,b] = sim(simbolo)/k
				if a in e_index and b in i_index:
					simbolo = (S_component(0,b,Z,B))*(B[a+m]+S_component(0,a,Z,B))/H_component(0,0,Z,B)
					H[a,b] = sim(simbolo)/k
				if a in i_index and b in e_index:
					simbolo = (B[b+m]+S_component(0,b,Z,B))*(S_component(0,a,Z,B))/H_component(0,0,Z,B)
					H[a,b] = sim(simbolo)/k
			
			simbolo = S_component(a,b,Z,B)
			S[a,b] = sim(simbolo)
	
	return H, S	

	
def matrix_gtdie(B):
	n = size(B)
	m = int((n - 1)/2)
	n_index = range (1,m+1)
	N = []
	N.append(Function('h0')(*B[1:2*m+1]))
	for d in n_index:
		N.append(Function('h'+str(d))(*B[1:2*m+1]))
		
	return N
	
def matrix_gtdieI(B,Z):
	n = size(B)
	m = (n - 1)/2
	n_index = range (1,m+1)
	N = []
	N.append(1)
	
	sum = 0
	for x in range(1,m+1):
		sum = sum + Z[x]*Z[x+m]
	
	for d in n_index:
		N.append(sum)
		
	return N
	

	
def sod_creator(Z,B,seed):
	a = product(seed,repeat=size(B))
	b =islice(a, 1, None)
	s = permutations(b,size(Z))
	
	return s
	
	
def difeomorphism_001(n):
	
	b = beta
	Z = []
	
	P = Symbol('\Phi ')
	Z.append(P)
	n_index = range(1,n+1)
	
	
	for x in n_index:
	   E = 'E_' + str(x)
	   E = Symbol(E)
	   E = Pow(E,b)/b
	   Z.append(E)
	   
	for x in n_index:
		I = 'I_' + str(x)
		I = Symbol(I)
		I = Pow(I,b)/b
		Z.append(I)
	
	
						
	return Z

	
def difeomorphism_002(n,B):
	
	b = beta
	Z = []
	
	P = Symbol('\Phi ')
	Z.append(P)
	n_index = range(1,n+1)
	
	L = 0
	for x in n_index:
	   E = 'E_' + str(x)
	   E = Symbol(E)
	   #E = Pow(E,b)/b
	   I = 'I_' + str(x)
	   I = Symbol(I)
	   #I = Pow(I,b)/b
	   L = L + E*I
	
	s = size(B)
	m = int((s - 1)/2)
	e_index = range (1,m+1)
	
	base = 1	
	for e in e_index:
		i = e+m
		if e == base:
			Z.append(Pow(L,2)*B[i])
		else:
			Z.append(1/B[e])
	
	base = base + 1
	for e in e_index:
		i = e+m
		if e == base:
			Z.append(Pow(L,2)*B[e])
		else:
			Z.append(1/B[i])
		
						
	return Z


def difeomorphism_003(n,B):

	
	b = beta
	Z = []
	
	P = Symbol('\Phi ')
	
	n_index = range(1,n+1)
	
	
	
	
		
	s = size(B)
	m = (s - 1)/2
	e_index = range (1,m+1)
	
	L = Function('L')(B[m+1])#(*B[m+1:2*m+1])
	P = P
	Z.append(P)
	
	for e in e_index:
		i = e+m
		variables = B[1:m+1]
		variables.remove(B[e])
		Z.append(Function('f'+str(e))(*variables))
		#Z.append(Pow(L,2)*B[e])
		
	for e in e_index:
		i = e+m
		#Z.append(B[i])
		variables = B[1+m:2*m+1]
		variables.remove(B[i])
		ariables = B[1:2*m+1]
		ariables.remove(B[i])
		Z.append(Function('g'+str(e))(*ariables)+Function('f'+str(i))(*variables))
		
						
	return Z
	
def difeomorphism_100(n):

	A = Symbol('\Lambda')
	
	Z = []
	B = []
	
	B.append(Symbol('\Phi '))
	n_index = range(1,n+1)
	
	
	for x in n_index:
	   E = 'E_' + str(x)
	   B.append(Symbol(E))
	   
	for x in n_index:
		I = 'I_' + str(x)
		B.append(Symbol(I))
		
	k = Function('k')(*B[0:2*n+1])
	k = k*k
		
	Z.append(Function('F')(*B[0:2*n+1]))
	for x in n_index:
	   X = Function('X_' + str(x))(*B[0:2*n+1])
	   Z.append(X)
	   
	for x in n_index:
		Y = Function('Y_' + str(x))(*B[0:2*n+1])
		Z.append(Y)
						
	return Z,B
	
	
def difeomorphism_200(Z,B):
	Z0 = []
	
	n = size(B)
	
	m = int((n-1)/2)
	Z0.append(Function(str(Z[0]))(*B[0:m+1]))
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	for e in e_index:
		variables = B[1:m+1]
		variables.remove(B[e])
		valor = Function(str(Z[e]))(*variables)
		Z0.append(valor)
		
	for i in i_index:
		variables = B[0:2*m+1]
		variables.remove(B[i])
		valor = Function(str(Z[i]))(*variables)
		Z0.append(valor)
		
	
	return Z0
	
def difeomorphism_300(Z,B):

	Z0 = []
	
	n = size(B)
	
	m = int((n-1)/2)
	Z0.append(Function(str(Z[0]))(*B[0:2*m+1]))
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	for e in e_index:
		variables = B[0:2*m+1]
		valor = Function(str(Z[e]))(*variables)
		Z0.append(valor)
		
	for i in i_index:
		variables = B[0:2*m+1]
		valor = Function(str(Z[i]))(*variables)
		Z0.append(valor)
		
	
	return Z0
	
def difeomorphism_sod(Z,B,s):

	Z0 = []
	
	n = size(B)
	m = int((n-1)/2)
	variables = B[0:2*m+1]
	
	
	
	d_index = range(0,n)
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	
	
	eriables = []
	iriables = []
	priables = []
	
	priables = list(variables)
	for d in d_index:
		if s[0][d] != 1:
			priables.remove(B[d])
	valor = Function(str(Z[0]))(*priables)
	Z0.append(valor)
	
	for e in e_index:
		eriables = list(variables)
		for d in d_index:
			if s[e][d] != 1:
				eriables.remove(B[d])
		valor = Function(str(Z[e]))(*eriables)
		Z0.append(valor)
		
	for i in i_index:
		iriables = list(variables)
		for d in d_index:
			if s[i][d] != 1:
				iriables.remove(B[d])
		valor = Function(str(Z[i]))(*iriables)
		Z0.append(valor)
		
	return Z0
	
def difeomorphism_soa(Z,B,s):

	Z0 = []
	c = []
	
	n = size(B)
	m = int((n-1)/2)
	variables = B[0:2*m+1]
	
	
	
	d_index = range(0,n)
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	
	
	eriables = []
	iriables = []
	priables = []
	
	priables = list(variables)
	as_sum = 0
	as_pro = 1
	for d in d_index:
		if s[0][d] == 0:
			priables.remove(B[d])
		elif s[0][d] == 2:
			priables.remove(B[d])
			as_sum = as_sum + Function(str(Z[0])+str(d))(B[d])
		elif s[0][d] == 4:
			priables.remove(B[d])
			as_sum = as_sum + Symbol(str(Z[0])+str(d))*B[d]
			c.append(Symbol(str(Z[0])+str(d)))
		elif s[0][d] == 3:
			priables.remove(B[d])
			as_pro = as_pro*Function(str(Z[0])+str(d))(B[d])
	if priables == []:
		as_pro = 0 if as_pro == 1  else as_pro
		valor = as_sum + as_pro
	else:
		if as_sum == 0 and as_pro == 1:
			valor = as_sum + as_pro*Function(str(Z[0]))(*priables)
		else:
			valor = as_sum + as_pro*Function(str(Z[0])+'_r')(*priables)
	Z0.append(valor)
	
	for e in e_index:
		eriables = list(variables)
		as_sum = 0
		as_pro = 1
		for d in d_index:
			if s[e][d] == 0:
				eriables.remove(B[d])
			elif s[e][d] == 2:
				eriables.remove(B[d])
				as_sum = as_sum + Function(str(Z[e])+str(d))(B[d])
			elif s[e][d] == 3:
				eriables.remove(B[d])
				as_pro = as_pro*Function(str(Z[e])+str(d))(B[d])
			elif s[e][d] == 4:
				eriables.remove(B[d])
				as_sum = as_sum + Symbol(str(Z[e])+str(d))*B[d]
				c.append(Symbol(str(Z[e])+str(d)))
		if eriables == []:
			as_pro = 0 if as_pro == 1  else as_pro
			valor = as_sum + as_pro
		else:
			if as_sum == 0 and as_pro == 1:
				valor = as_sum + as_pro*Function(str(Z[e]))(*eriables)
			else:
				valor = as_sum + as_pro*Function(str(Z[e])+'_r')(*eriables)
		Z0.append(valor)
		
	for i in i_index:
		iriables = list(variables)
		as_sum = 0
		as_pro = 1
		for d in d_index:
			if s[i][d] == 0:
				iriables.remove(B[d])
			elif s[i][d] == 2:
				iriables.remove(B[d])
				as_sum = as_sum + Function(str(Z[i])+str(d))(B[d])
			elif s[i][d] == 3:
				iriables.remove(B[d])
				as_pro = as_pro*Function(str(Z[i])+str(d))(B[d])
			elif s[i][d] == 4:
				iriables.remove(B[d])
				as_sum = as_sum + Symbol(str(Z[i])+str(d))*B[d]
				c.append(Symbol(str(Z[i])+str(d)))
		if iriables == []:
			as_pro = 0 if as_pro == 1  else as_pro
			valor = as_sum + as_pro
		else:
			if as_sum == 0 and as_pro == 1:
				valor = as_sum + as_pro*Function(str(Z[i]))(*iriables)
			else:
				valor = as_sum + as_pro*Function(str(Z[i])+'_r')(*iriables)
		Z0.append(valor)
		
	return Z0,c
	
def difeomorphism_sop(Z,B,s):

	Z0 = []
	
	n = size(B)
	m = int((n-1)/2)
	variables = B[0:2*m+1]
	
	
	
	d_index = range(0,n)
	e_index = range(1,m+1)
	i_index = range(m+1,2*m+1)
	
	
	
	eriables = []
	iriables = []
	priables = []
	c = []
	
		
	
	valor = get_poly(s[0],n,B,c,B[0])
	Z0.append(valor)
	
	
	for e in e_index:
		valor = get_poly(s[e],n,B,c,B[e])
		Z0.append(valor)
	
	
	for i in i_index:
		valor = get_poly(s[i],n,B,c,B[i])
		Z0.append(valor)
		
	return Z0,c
	
def get_poly(s,n,B,c,v):
	d_index = range(0,n)
	terms = []
	valor = 0
	for d in d_index:
		if s[d] != N:
			for grade in range(1,s[d]+1):
					terms.append(pow(B[d],grade))
	terms.append(1)

	poly_size = size(terms)
	
	coefficients = zeros(poly_size)
	
	
	for x in range(0,poly_size):
		for y in range(x,poly_size):
			coefficients[x,y] = Symbol(str(v)+str(x)+str(y))
			c.append(coefficients[x,y])
	
	enco = 0
	anco = 0
	for d in d_index:
		if s[d] != N or s[d] != 0:
			anco = enco
			enco = enco + s[d]			
			coefficients[anco:enco,anco:enco] = zeros(s[d])
		
	
	for i in range(0,poly_size):
		for j in range (0, poly_size):
			if i < j:
				valor = valor+coefficients[i,j]*terms[i]*terms[j]
	
	return valor

def assert_ii(H,n):
	if H[n+1:2*n+1,n+1:2*n+1] == zeros(n): 
		return 0
	else: 
		return 1

def assert_00(H,n):
	if H[0,0] == 0: 
		return 1
	else: 
		return 0
		
def assert_0i(H,n):
	if H[0,n+1:2*n+1] == zeros(1,n): 
		return 0
	else: 
		return 1
		
def assert_0e(H,n):
	if H[0,1:n+1] == zeros(1,n): 
		return 1
	else: 
		return 0

def assert_ee(H,n):
	if H[1:n+1,1:n+1] == zeros(n): 
		return 1
	else: 
		return 0

def assert_ee00(H,n):
	E = H[1:n+1,1:n+1]
	O = H[0,1:n+1]
	for i in range(0,n):
		for j in range(i,n):
			l = E[i,j] - O[i]*O[j]
			if l != 0:
				return 1
	return 0

def assert_ie(H,n):
	M = H[n+1:2*n+1,1:n+1]
	v = 1
	for i in range(0,n):
		if M[i,i] == 0:
			return 1
	
	return 0

def assert_co(H,n):
	if H[0] == 0:
		return 0
	return 1

def assert_all(H,n):
	return 0
def assert_zero(H,n):
	if H != 0:
		return 0
	return 1
	
def get_difeomorphisms(n,s,which_difeo,ini_tabix,end_tabix):
	
	
	
	
	correct_difeomorphism = []
	R = []	
	s2 = []	
	Z = []
	B = []
	Z,B = n_grados001(n)
	N = matrix_gtdie(Z)
	#print ('trabajamos con: ' + str(size(s)))
	tabix = ini_tabix
	
	start_time = time.time()
	while tabix <= end_tabix:
		try:
			s1 = s[tabix]
		except IndexError:
			tabix = end_tabix + 1
			break
		else:
			am_i_good = 0
			Z0 = which_difeo(Z,B,s1)
			H = zeros(n)
			S = zeros(n)
			H,S = difeomorphism_equations004(Z0,B,factor,0,N)
			D = H + S
			
			am_i_good = assert_00(D,n) + assert_0e(D,n) + assert_ee00(D,n) + assert_ie(D,n) + assert_ii(D,n) + assert_0i(D,n)
			if am_i_good == 0:
				correct_difeomorphism.append(Z0)
				R.append(D)
				s2.append(s1)
				#print (tabix, "--- %s seconds ---" % (time.time() - start_time))
				#s1 = []	
		tabix = tabix + 1
		gc.collect()
	#print("--- %s seconds ---" % (time.time() - start_time))
			#if tabix % 10000 == 0:
			#	print tabix
			#	
			#	s1 = [] 
	print("--- %s seconds ---" % (time.time() - start_time))	
	return correct_difeomorphism, R, s2
	
def get_difeomorphisms00(n,s,which_difeo,ini_tabix,end_tabix,a_index,b_index,asserts):
	import time

	correct_difeomorphism = []
	R = []	
	s2 = []	
	Z = []
	B = []
	Z,B = n_grados001(n)
	N = matrix_gtdie(Z)
	#print ('trabajamos con: ' + str(size(s)))
	tabix = 0
	
	start_time = time.time()
	for s1 in s:
		if tabix < ini_tabix:
			tabix = tabix + 1
			continue
			
		am_i_good = 0
		Z0 = which_difeo(Z,B,s1)
		H = []
		S = zeros(n)
		H,S = difeomorphism_equations002(Z0,B,factor,0,N,a_index,b_index)
		D = H + S
		
		for asserto in range(0,size(asserts)):
			am_i_good = am_i_good + asserts[asserto](D,n)
			
		if am_i_good == 0:
			correct_difeomorphism.append(Z0)
			R.append(D)
			s2.append(s1)
		
		if tabix >= end_tabix:
			break
		
		tabix = tabix + 1
	print("--- %s seconds ---" % (time.time() - start_time))		
	return correct_difeomorphism, R, s2
	
def get_difeomorphisms01(n,s,which_difeo,ini_tabix,end_tabix,a_index,b_index,asserts):
	import time

	correct_difeomorphism = []
	R = []	
	s2 = []	
	Z = []
	B = []
	B,Z = n_grados001(n)
	N = matrix_gtdie(Z)
	#print ('trabajamos con: ' + str(size(s)))
	tabix = 0
	
	start_time = time.time()
	for s1 in s:
		if tabix < ini_tabix:
			tabix = tabix + 1
			continue
			
		am_i_good = 0
		Z0 = which_difeo(Z,B,s1)
		#D = compatibility_equation(Z0,B,simplify,N)
		h = Matrix(Z0)
		D = h.jacobian(B).det()
		
		
		for asserto in range(0,size(asserts)):
			am_i_good = am_i_good + asserts[asserto](D,n)
			
		if am_i_good == 0:
			correct_difeomorphism.append(Z0)
			R.append(D)
			s2.append(s1)
		
		if tabix >= end_tabix:
			break
		
		tabix = tabix + 1
	print("--- %s seconds ---" % (time.time() - start_time))		
	return correct_difeomorphism, R, s2
	

def set_replacement(n,Z):
	replacement = [(Z[i],1) for i in range(0,2*n+1)]
	return replacement
	
def set_EQ(n,S,replacement,Z):
	EQ=[]
	e_index = range(1,n+1)
	i_index = range(n+1,2*n+1)

	for x in range(0,2*n+1):
		for y in range(x,2*n+1):
			if x in e_index and y in i_index:
				pass
			elif x in i_index and y in e_index:
				pass
			else:
				EQ.append(Eq(expand(S[x,y].subs(replacement)),0))
	return EQ
	
def set_c(n,S,replacement):
	c=set()
	e_index = range(1,n+1)
	i_index = range(n+1,2*n+1)
	for x in range(0,2*n+1):
		for y in range(x,2*n+1):
			if x in e_index and y in i_index:
				pass
			elif x in i_index and y in e_index:
				pass
			else:
				z = expand(S[x,y]).subs(replacement)
				w = z.atoms()
				for i in w:
					if type(i) == Symbol:
						c.add(i)
	return c
	
def get_replacement(sol):
	replacement = []
	tabix = 0
	for x in sol[1]:
		for y in x:
			replacement.append((sol[0][tabix],y))
			tabix = tabix+1
	return replacement
	
def solve_s(n,Z,S):
	replacement = set_replacement(n,Z)
	EQ=set_EQ(n,S,replacement,Z)
	c=set_c(n,S,replacement)
	sol = solve(EQ,c, set= True, manual = True)
	replacement = get_replacement(sol)
	
	return EQ,c,replacement

	
def check_me(ew,EQS,n,dZE,dc):
    keye = list(ew.keys())
    matchb = []

    for a in range(0,len(keye)):
        i = keye[a]
        x = ew[i]
        if x != 0:
            argse = poly_args(expand(x),dZE,dc,factor)
            key_argse = list(argse.keys())
            matcha = []
            for x1 in range(0,2*n+1):
                for y1 in range(x1,2*n+1):
                    aw = EQS[x1,y1]
                    if aw == 0: continue
                    keya = list(aw.keys())
                    for b in range(0,len(keya)):
                        j = keya[b]
                        y = aw[j]
                        if y != 0:
                            argsa = poly_args(expand(y),dZE,dc,factor)
                            key_argsa = list(argsa.keys())
                            if key_argsa == key_argse:
                                matcha.append([argsa,x1,y1,i,j])
            matchb.append(matcha)
    return matchb