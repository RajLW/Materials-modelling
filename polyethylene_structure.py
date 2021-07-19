import numpy as np
from math import cos, sin, pi
import os

tol = 3

lc = 1.5260
lh = 1.0900
dqc = -0.120
dqh = 0.060
mol = 1
atomc = 1
atomh = 2

bondcc = 1
bondch = 2

angleccc = 1
anglehch = 2

dihedralcccc = 1
dihedralhcch = 2

thetac = (109.5)*(pi/180)

def read_data():
	mol = 0
	while mol <1:
		flag = 0
		while flag!=1:
			monolen = int(input("\nPlease enter number of monomer units you want"))
			if monolen <=0:
				print("\nInvalid value, please enter again")
			if monolen >150:
				print("\nVery large value, please enter again (current limit is <=150 monomer units")
			else:
				flag = 1
		flag1 = 0
		while flag1 != 1:
			flag = int(input("\nNow to enter co-ordinates of simulation box that you want press 1, else to use default values (80x60x50) press 2"))
			if flag ==2 or flag==1:
				flag1 = 1
			else:
				print("\nInvalid values, enter again")
				
		xlo = 0
		xhi = 83
		ylo = 0
		yhi = 63
		zlo = 0
		zhi = 55
		
		flag1 = 0
		if flag == 1:
			while flag1 != 1:
				#xlo = float(input("\nXlo (in A):"))
				xhi = float(input("\nXhi (in A):"))
				
				#ylo = float(input("\nYlo (in A):"))
				yhi = float(input("\nYhi (in A):"))
				
				#zlo = float(input("\nZlo (in A):"))
				zhi = float(input("\nZhi (in A):"))
			
				if abs(xhi-xlo) < 180 and abs(yhi-ylo) < 60 and abs(zhi-zlo) < 50:
					flag1 = 1
				else:
					print("\nVery large simulation box, current limit is (80x60x50), enter values again or you can change values in the code")
		#x0 = x0 + ((monolen)*lc*sin(thetac/2)) + (tol)
		molx = int((abs(xhi-xlo))/(((monolen)*lc*sin(thetac/2))+tol))
		moly = int((abs(yhi-ylo))/((2*(lh*sin(pi/3)))+tol))
		molz = int((abs(zhi-zlo))/((1*((lc*cos(thetac/2))+(2*lh*cos(pi/3))))+ tol))
		mol = molx*moly*molz
		if mol <1:
			print("\nEither simulation box is too small or no. of monomer units are too high, please enter the values again")
	return mol, molx, moly, molz, xlo, xhi, ylo, yhi, zlo, zhi, monolen
	

def atomsf(f, monolen, mol, xc, yc, zc, k1):
	k= k1
	l=3
	m=3
	xcn = xc
	zcn = zc
	nofatoms = 0
	yh1 = (yc + (lh*sin(pi/3)))
	yh2 = (yc - (lh*sin(pi/3)))
	
	for n in range(1,monolen+1):
		xcn = xcn + (lc*sin(thetac/2))
		zcn = zcn + (((-1)**(n))*(lc*cos(thetac/2)))
		ycn = yc
		
		xh1n = xcn
		yh1n = yh1
		zh1n = zcn +(1*((-1)**(n))*(lh*cos(pi/3)))
		
		xh2n = xcn
		yh2n = yh2
		zh2n = zh1n
		
		tempc = np.array([xcn, ycn, zcn])
		temph1 = np.array([xh1n, yh1n, zh1n])
		temph2 = np.array([xh2n, yh2n, zh2n])
		
		
		tu = "\n" + str(k) +" 	" + str(mol)+" 	" + str(atomc)+" 	" + str(dqc)+" 	" + str(tempc[0])+" 	" + str(tempc[1])+" 	" + str(tempc[2])
		f.write(tu)
		tu = "\n"+ str(k+1) +" 	" + str(mol)+" 	" + str(atomh)+" 	" + str(dqh)+"	 " + str(temph1[0])+" 	" + str(temph1[1])+" 	" + str(temph1[2])
		f.write(tu)
		tu = "\n"+ str(k+2) +" 	" + str(mol)+" 	" + str(atomh)+" 	" + str(dqh)+" 	" + str(temph2[0])+" 	" + str(temph2[1])+" 	" + str(temph2[2])
		nofatoms = k+2
		f.write(tu)
		k=k+3
		
		
	return nofatoms, k
	



def bondsf(f, monolen, mol, xc, yc, zc, k1, l1):
	xcn = xc
	k= k1
	l= l1
	m=3
	zcn = zc
	nofbonds = 0
	yh1 = (yc + (lh*sin(pi/3)))
	yh2 = (yc - (lh*sin(pi/3)))
	
	for n in range(1,monolen+1):

		xcn = xcn + (lc*sin(thetac/2))
		zcn = zcn + ((-1**(n))*(lc*cos(thetac/2)))
		ycn = yc
		
		
		xh1n = xcn
		yh1n = yh1
		zh1n = zcn +(1*(-1**(n))*(lh*cos(pi/3)))
		
		xh2n = xcn
		yh2n = yh2
		zh2n = zh1n
		
		tempc = np.array([xcn, ycn, zcn])
		temph1 = np.array([xh1n, yh1n, zh1n])
		temph2 = np.array([xh2n, yh2n, zh2n])
		
		 
		tu =  "\n"+ str(l) +" 	" + str(bondcc)+" 	" + str(k)+" 	" + str(k-3)
		f.write(tu)
		tu = "\n"+ str(l+1) +" 	" + str(bondch)+" 	" + str(k)+" 	" + str(k+1)
		f.write(tu)
		tu = "\n"+ str(l+2) +" 	" + str(bondch)+" 	" + str(k)+" 	" + str(k+2)
		f.write(tu)
		
		nofbonds = l+2
		k=k+3
		l=l+3
		
	return nofbonds, k, l
		



def anglesf(f, monolen, mol, xc, yc, zc, k1, m1):
	k= k1
	l=3
	m= m1
	xcn = xc
	zcn = zc
	nofangles = 0
	yh1 = (yc + (lh*sin(pi/3)))
	yh2 = (yc - (lh*sin(pi/3)))
	
	for n in range(1,monolen):

		xcn = xcn + (lc*sin(thetac/2))
		zcn = zc + ((-1**(n))*(lc*cos(thetac/2)))
		ycn = yc
		
		
		xh1n = xcn
		yh1n = yh1
		zh1n = zcn +(1*(-1**(n))*(lh*cos(pi/3)))
		
		xh2n = xcn
		yh2n = yh2
		zh2n = zh1n
		
		tempc = np.array([xcn, ycn, zcn])
		temph1 = np.array([xh1n, yh1n, zh1n])
		temph2 = np.array([xh2n, yh2n, zh2n])
		
		tu = "\n" + str(m) +" 	" + str(angleccc)+" 	" + str(k-6)+" 	" + str(k-3)+" 	" + str(k)
		f.write(tu)
		tu = "\n"+ str(m+1) +" 	" + str(anglehch)+" 	" + str(k+1)+" 	" + str(k)+" 	" + str(k+2)
		f.write(tu)
		nofangles = m+1
		m=m+2
			
		k=k+3
		l=l+3
		
		
	return nofangles, k, m
		
def dihedralsf(f, monolen, mol, xc, yc, zc, k1, z1):
	k= k1
	l=3
	m=3
	z= z1
	xcn = xc
	zcn = zc
	nofangles = 0
	yh1 = (yc + (lh*sin(pi/3)))
	yh2 = (yc - (lh*sin(pi/3)))
	
	for n in range(1,monolen-1):
		k=k+3
		xcn = xcn + (lc*sin(thetac/2))
		zcn = zc + ((-1**(n))*(lc*cos(thetac/2)))
		ycn = yc
		
		
		xh1n = xcn
		yh1n = yh1
		zh1n = zcn +(1*(-1**(n))*(lh*cos(pi/3)))
		
		xh2n = xcn
		yh2n = yh2
		zh2n = zh1n
		
		tempc = np.array([xcn, ycn, zcn])
		temph1 = np.array([xh1n, yh1n, zh1n])
		temph2 = np.array([xh2n, yh2n, zh2n])
		
						#[(z+4), mol, (k-9), (k-6), (k-3), k]])
		tu = "\n" + str(z) +" 	" + str(dihedralhcch)+" 	" + str(k-2)+" 	" + str(k-3)+" 	" + str(k)+" 	" + str(k+1)
		f.write(tu)
		tu = "\n" + str(z+1) +" 	" + str(dihedralhcch)+" 	" + str(k-2)+" 	" + str(k-3)+" 	" + str(k)+" 	" + str(k+2)
		f.write(tu)
		tu = "\n" + str(z+2) +" 	" + str(dihedralhcch)+" 	" + str(k-1)+" 	" + str(k-3)+" 	" + str(k)+" 	" + str(k+1)
		f.write(tu)
		tu = "\n" + str(z+3) +" 	" + str(dihedralhcch)+" 	" + str(k-1)+" 	" + str(k-3)+" 	" + str(k)+" 	" + str(k+2)
		f.write(tu)
		tu = "\n" + str(z+4) +" 	" + str(dihedralcccc)+" 	" + str(k-9)+" 	" + str(k-6)+" 	" + str(k-3)+" 	" + str(k)
		f.write(tu)
		nofdihedrals = z+4
		z = z+5
		
		l=l+3
		
	
	return nofdihedrals, k, z


def write_tempdata_atoms(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1):
	
	xc = x0
	yc = y0
	zc = z0

	xh1 = xc
	yh1 = (y0 + (lh*sin(pi/3)))
	zh1 = (zc + (lh*cos(pi/3)))

	xh2 = xc
	yh2 = (y0 - (lh*sin(pi/3)))
	zh2 = zh1

	c1 = np.array([[xc, yc, zc]])
	h1 = np.array([[xh1, yh1, zh1]])
	h2 = np.array([[xh2, yh2, zh2]])
	if mol==1:
		f.write("\n"+" Atoms\n")
	tu = "\n" + str(k1) +" 	" + str(mol)+" 	" + str(atomc)+" 	" + str(dqc)+" 	" + str(c1[0][0])+" 	" + str(c1[0][1])+" 	" + str(c1[0][2])
	f.write(tu)
	tu = "\n"+ str(k1+1) +" 	" + str(mol)+" 	" + str(atomh)+" 	" + str(dqh)+"	 " + str(h1[0][0])+" 	" + str(h1[0][1])+" 	" + str(h1[0][2])
	f.write(tu)
	tu = "\n"+ str(k1+2) +" 	" + str(mol)+" 	" + str(atomh)+" 	" + str(dqh)+" 	" + str(h2[0][0])+" 	" + str(h2[0][1])+" 	" + str(h2[0][2])
	f.write(tu)
	k1 = k1+3
	
	nofatoms, k = atomsf(f, monolen, mol, xc, yc, zc, k1)
	
	return nofatoms,k

def write_tempdata_bonds(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1):	
	xc = x0
	yc = y0
	zc = z0

	xh1 = xc
	yh1 = (y0 + (lh*sin(pi/3)))
	zh1 = (zc + (lh*cos(pi/3)))

	xh2 = xc
	yh2 = (y0 - (lh*sin(pi/3)))
	zh2 = zh1

	c1 = np.array([[xc, yc, zc]])
	h1 = np.array([[xh1, yh1, zh1]])
	h2 = np.array([[xh2, yh2, zh2]])
	if mol==1:
		f.write("\n"+" Bonds\n")
		
	tu = "\n" + str(l1) +" 	" + str(bondch)+" 	" + str(k1+1)+" 	" + str(k1)
	f.write(tu)
	tu = "\n"+ str(l1+1) +" 	" + str(bondch)+" 	" + str(k1+2)+" 	" + str(k1)
	f.write(tu)
	l1 = l1+2
	k1 = k1+3
	nofbonds, k, l = bondsf(f, monolen, mol, xc, yc, zc, k1, l1)
	
	return nofbonds, k, l

def write_tempdata_angles(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1):	
	xc = x0
	yc = y0
	zc = z0

	xh1 = xc
	yh1 = (y0 + (lh*sin(pi/3)))
	zh1 = (zc + (lh*cos(pi/3)))

	xh2 = xc
	yh2 = (y0 - (lh*sin(pi/3)))
	zh2 = zh1

	c1 = np.array([[xc, yc, zc]])
	h1 = np.array([[xh1, yh1, zh1]])
	h2 = np.array([[xh2, yh2, zh2]])
	if mol==1:
		f.write("\n"+" Angles\n")
	
	tu = "\n" + str(m1) +" 	" + str(anglehch)+" 	" + str(k1+1)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	k1 = k1+3
	
	tu = "\n" + str(m1+1) +" 	" + str(anglehch)+" 	" + str(k1+1)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	k1 = k1+3
	m1 = m1+2
	
	nofangles, k, m = anglesf(f, monolen, mol, xc, yc, zc, k1, m1)
	
	return nofangles, k, m

def write_tempdata_dihedrals(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1):	
	xc = x0
	yc = y0
	zc = z0
	mol1 = mol
	xh1 = xc
	yh1 = (y0 + (lh*sin(pi/3)))
	zh1 = (zc + (lh*cos(pi/3)))

	xh2 = xc
	yh2 = (y0 - (lh*sin(pi/3)))
	zh2 = zh1

	c1 = np.array([[xc, yc, zc]])
	h1 = np.array([[xh1, yh1, zh1]])
	h2 = np.array([[xh2, yh2, zh2]])
	if mol==1:
		f.write("\n"+" Dihedrals\n")
	if mol1 > 1:
		k1 = k1+3		
			
	tu = "\n" + str(z1) +" 	" + str(dihedralhcch)+" 	" + str(k1-2)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+1)
	f.write(tu)
	tu = "\n" + str(z1+1) +" 	" + str(dihedralhcch)+" 	" + str(k1-2)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	tu = "\n" + str(z1+2) +" 	" + str(dihedralhcch)+" 	" + str(k1-1)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+1)
	f.write(tu)
	tu = "\n" + str(z1+3) +" 	" + str(dihedralhcch)+" 	" + str(k1-1)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	z1 = z1+4
	k1 = k1+3
	tu = "\n" + str(z1) +" 	" + str(dihedralhcch)+" 	" + str(k1-2)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+1)
	f.write(tu)
	tu = "\n" + str(z1+1) +" 	" + str(dihedralhcch)+" 	" + str(k1-2)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	tu = "\n" + str(z1+2) +" 	" + str(dihedralhcch)+" 	" + str(k1-1)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+1)
	f.write(tu)
	tu = "\n" + str(z1+3) +" 	" + str(dihedralhcch)+" 	" + str(k1-1)+" 	" + str(k1-3)+" 	" + str(k1)+" 	" + str(k1+2)
	f.write(tu)
	z1 = z1+4

	nofdihedrals, k, z = dihedralsf(f, monolen, mol, xc, yc, zc, k1, z1)
	
	return nofdihedrals, k, z

def write_data(g, nofatoms, nofbonds, nofangles, nofdihedrals, xlo, xhi, ylo, yhi, zlo, zhi):
	
	f1 = open("Polyethylene.data", "a")
	f2 = open("OPLS.txt", "r")
	f1.write("#Polyethylene\n")
	f1.write(str(nofatoms)+" atoms\n")
	f1.write(str(nofbonds)+" bonds\n")
	f1.write(str(nofangles)+" angles\n")
	f1.write(str(nofdihedrals)+" dihedrals\n")
	f1.write(str(0)+" impropers\n"+"\n")
	
	f1.write(str(2)+" atom types\n")
	f1.write(str(2)+" bond types\n")
	f1.write(str(2)+" angle types\n")
	f1.write(str(2)+" dihedral types\n")
	f1.write(str(0)+" improper types\n")
	
	f1.write(str(xlo)+"	"+str(xhi)+"	xlo xhi\n")
	f1.write(str(ylo)+"	"+str(yhi)+"	ylo yhi\n")
	f1.write(str(zlo)+"	"+str(zhi)+"	zlo zhi\n")
	for line in f2:
		f1.write(line)
	for line in g:
		f1.write(line)
	f1.close()
	f2.close()

def write_simubox(f, decision):
	mol1, molx1, moly1, molz1, xlo1, xhi1, ylo1, yhi1, zlo1, zhi1, monolen = read_data()
	mol = 1
	flag = 1
	k1 = 1
	l1 = 1
	m1 = 1
	z1 = 1
	x0 = xlo1 + 3
	y0 = ylo1 + 3
	z0 = zlo1 + 5
	if decision == 1:
		molx1 = 1
		moly1 = 1
		molz1 = 1
	for moly in range(1, moly1+1):
		x0 = xlo1 + 3
		z0 = zlo1 + 5
		for molz in range(1, molz1+1):
			x0 = xlo1 + 3
			for molx in range(1, molx1+1):
				nofatoms, k = write_tempdata_atoms(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1)
				
				mol = mol +1
				x0 = x0 + ((monolen)*lc*sin(thetac/2)) + (tol)
				k1 = k
				
		
			z0 = z0 + (2*((lc*cos(thetac/2)))) + tol
		y0 = y0 + (2*(lh*sin(pi/3))) + tol
	mol = 1
	flag = 1
	k1 = 1
	l1 = 1
	m1 = 1
	z1 = 1
	x0 = xlo1 + 3
	y0 = ylo1 + 3
	z0 = zlo1 + 5
	for moly in range(1, moly1+1):
		x0 = xlo1 + 3
		z0 = zlo1 + 5
		for molz in range(1, molz1+1):
			x0 = xlo1 + 3
			for molx in range(1, molx1+1):
				
				nofbonds, k, l = write_tempdata_bonds(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1)
				
				mol = mol +1
				x0 = x0 + ((monolen)*lc*sin(thetac/2)) + (tol)
				k1 = k
				l1 = l
				
		
			z0 = z0 + (2*((lc*cos(thetac/2)))) + tol
		y0 = y0 + (2*(lh*sin(pi/3))) + tol
	mol = 1
	flag = 1
	k1 = 1
	l1 = 1
	m1 = 1
	z1 = 1
	x0 = xlo1 + 3
	y0 = ylo1 + 3
	z0 = zlo1 + 5
	for moly in range(1, moly1+1):
		x0 = xlo1 + 3
		z0 = zlo1 + 5
		for molz in range(1, molz1+1):
			x0 = xlo1 + 3
			for molx in range(1, molx1+1):
				
				nofangles, k, m = write_tempdata_angles(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1)
				
				mol = mol +1
				x0 = x0 + ((monolen)*lc*sin(thetac/2)) + (tol)
				k1 = k
				m1 = m
		
			z0 = z0 + (2*((lc*cos(thetac/2)))) + tol
		y0 = y0 + (2*(lh*sin(pi/3))) + tol
	mol = 1
	flag = 1
	k1 = 4
	l1 = 1
	m1 = 1
	z1 = 1
	x0 = xlo1 + 3
	y0 = ylo1 + 3
	z0 = zlo1 + 5
	for moly in range(1, moly1+1):
		x0 = xlo1 + 3
		z0 = zlo1 + 5
		for molz in range(1, molz1+1):
			x0 = xlo1 + 3
			for molx in range(1, molx1+1):
				
				nofdihedrals, k, z = write_tempdata_dihedrals(f, mol, monolen, x0, y0, z0, k1, l1, m1, z1)
				mol = mol +1
				x0 = x0 + ((monolen)*lc*sin(thetac/2)) + (tol)
				k = k + 3
				k1 = k
				z1 = z
		
			z0 = z0 + (2*((lc*cos(thetac/2)))) + tol
		y0 = y0 + (2*(lh*sin(pi/3))) + tol
	return nofatoms, nofbonds, nofangles, nofdihedrals, xlo1, xhi1, ylo1, yhi1, zlo1, zhi1

f = open("name.data", "a")
condition = True
while condition:
 decision = int(input(" To generte single polymer chain press 1, to generate multiple polymer chains press 2"))
 if decision <1 or decision >2:
 	print(" Enter valid values")
 else:
 	condition = False
nofatoms, nofbonds, nofangles, nofdihedrals, xlo1, xhi1, ylo1, yhi1, zlo1, zhi1 = write_simubox(f, decision)	
f.close()	

g = open("name.data", "r")
write_data(g, nofatoms, nofbonds, nofangles, nofdihedrals, xlo1, xhi1, ylo1, yhi1, zlo1, zhi1)
g.close()	
os.remove("name.data")	
	
	
	
	
	
	
	
	
	
	
