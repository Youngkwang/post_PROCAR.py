#!/Users/yj359/opt/anaconda3/bin/python

#First Created : 2015/01/23 by Youngkwang
#Last Modified : 2022/10/05 by Youngkwang

import YK_vasp as vasp
import math
import os

wave_suffix = input("please type suffix for igor wave name: ")

PROCAR = open("PROCAR", "r")
OUTCAR = open("OUTCAR" , "r")

PROCAR_lines = PROCAR.readlines()
OUTCAR_lines = OUTCAR.readlines()

num_IONS = int(PROCAR_lines[1].split()[11])
num_KPOINTS = int(PROCAR_lines[1].split()[3])
num_BANDS = int(PROCAR_lines[1].split()[7])
ref_E = -0.06365 #float(vasp.get_keyword_from("E-fermi", "OUTCAR").split()[2])

starting_KPOINT = 17 ##This option is for choosing kpoint range (Defualt: 1)
ending_KPOINT = num_KPOINTS ##This option is for choosing kpoint range (Defualt: num_KPOINTS)
starting_BAND = 1 ##This option is for choosing band range (Defualt: 1)
ending_BAND = num_BANDS ##This option is for choosing band range (Defialt: num_BANDS)
starting_ION = 1 ##This option is for choosing ion range (Defualt: 1)
ending_ION = num_IONS ##This option is for choosing ion range (Defialt: num_IONS)

check_SOC = vasp.get_keyword_from("LSORBIT", "OUTCAR").split()[2]

pi = 3.14159265

kvec_loc = int(OUTCAR_lines.index(vasp.get_keyword_from("2pi/SCALE", "OUTCAR")))
kvec_loc_bef = kvec_loc + 1

temp = float(OUTCAR_lines[kvec_loc_bef].split()[0])
value = temp * 2 * pi
kvec_X_bef = value

temp = float(OUTCAR_lines[kvec_loc_bef].split()[1])
value = temp * 2 * pi
kvec_Y_bef = value

temp = float(OUTCAR_lines[kvec_loc_bef].split()[2])
value = temp * 2 * pi
kvec_Z_bef = value

kvec_K_bef = 0

print ("K-vector is being generated")
kvec_mat = [float(0) for x in range(0, num_KPOINTS)]

for i in range(starting_KPOINT, num_KPOINTS + 1):
	j = kvec_loc + i

	temp = float(OUTCAR_lines[j].split()[0])
	value = temp * 2 * pi
	kvec_X = value

	temp = float(OUTCAR_lines[j].split()[1])
	value = temp * 2 * pi
	kvec_Y = value

	temp = float(OUTCAR_lines[j].split()[2])
	value = temp * 2 * pi
	kvec_Z = value

	kvec_K = kvec_K_bef + math.sqrt((kvec_X - kvec_X_bef) ** 2 + (kvec_Y - kvec_Y_bef) ** 2 + (kvec_Z - kvec_Z_bef) ** 2)
	kvec_mat[i - 1] = kvec_K

	kvec_X_bef = kvec_X
	kvec_Y_bef = kvec_Y
	kvec_Z_bef = kvec_Z
	kvec_K_bef = kvec_K

eigen_val_mat = [[float(0) for x in range(0, num_KPOINTS)] for y in range(0, num_BANDS)]
elec_occ_mat = [[float(0) for x in range(0, num_KPOINTS)] for y in range(0, num_BANDS)]
pro_band_mat = [[[[float(0) for x in range(0, num_KPOINTS)] for y in range(0, num_BANDS)] for z in range(0, num_IONS + 1)] for t in range(0, 4)]

for m in range(starting_BAND, ending_BAND + 1):
	print ("Band number %03d is being generated" % m)

	for o in range(starting_ION, ending_ION + 2):

		for k in range(starting_KPOINT, num_KPOINTS + 1):
			if check_SOC == "F":
				h = 4 + (k - 1) * (num_BANDS * (num_IONS + 5) + 3) + (2 + (m - 1) * (num_IONS + 5))
			elif check_SOC == "T":
				h = 4 + (k - 1) * (num_BANDS * 4 * (num_IONS + 2) + 3) + (2 + (m - 1) * 4 * (num_IONS + 2))


			eigen_VAL = float(PROCAR_lines[h - 1].split()[4])
			eigen_vbm_VAL = eigen_VAL - ref_E
			elec_occ = float(PROCAR_lines[h - 1].split()[7])

			eigen_val_mat[m - 1][k - 1] = eigen_vbm_VAL
			elec_occ_mat[m - 1][k - 1] = elec_occ

			if check_SOC == "F":
				l =  4 + (k - 1) * (num_BANDS * (num_IONS + 5) + 3) + (4 + (m - 1) * (num_IONS + 5)) + o
			elif check_SOC == "T":
				l =  4 + (k - 1) * (num_BANDS * 4 * (num_IONS + 2) + 3) + (4 + (m - 1) * 4 * (num_IONS + 2)) + o

			pro_s = float(PROCAR_lines[l - 1].split()[1])
			pro_p = float(PROCAR_lines[l - 1].split()[2])
			pro_d = float(PROCAR_lines[l - 1].split()[3])
			pro_t = float(PROCAR_lines[l - 1].split()[4])


			pro_band_mat[0][o - 1][m - 1][k - 1] = pro_s
			pro_band_mat[1][o - 1][m - 1][k - 1] = pro_p
			pro_band_mat[2][o - 1][m - 1][k - 1] = pro_d
			pro_band_mat[3][o - 1][m - 1][k - 1] = pro_t

band_for_write = open("band_structure.itx", "w")

band_for_write.write("IGOR\n")
band_for_write.write("WAVES/D	")
band_for_write.write("Kvector_%s	" % wave_suffix)

for y in range(starting_BAND - 1, ending_BAND):
	band_for_write.write("E_band_%03d_%s	" % ((y + 1), wave_suffix))
		
band_for_write.write("\n")
band_for_write.write("BEGIN\n")

for x in range(starting_KPOINT - 1, ending_KPOINT):
	band_for_write.write("%.10f	" % kvec_mat[x])

	for y in range(starting_BAND - 1, ending_BAND):
		band_for_write.write("%.8f	" % eigen_val_mat[y][x])

	band_for_write.write("\n")

band_for_write.write("END\n")

band_for_write.write("X Display	")
band_for_write.write("as	")
band_for_write.write('"band_structure_%s"' % wave_suffix)
band_for_write.write(";DelayUpdate\n")

count_band = ending_BAND - starting_BAND + 1
count_start = starting_BAND - 1
count_end = starting_BAND + 4

while count_band > 5:
	band_for_write.write("X AppendToGraph	")

	for y in range(count_start, count_end):
		band_for_write.write("E_band_%03d_%s	" % ((y + 1), wave_suffix))

	band_for_write.write("vs	")
	band_for_write.write("Kvector_%s	" % wave_suffix)
	band_for_write.write(";DelayUpdate\n")

	count_band -= 5
	count_start += 5
	count_end += 5

if count_start != ending_BAND:
	band_for_write.write("X AppendToGraph	")

	for y in range(count_start, ending_BAND):
		band_for_write.write("E_band_%03d_%s	" % ((y + 1), wave_suffix))

	band_for_write.write("vs	")
	band_for_write.write("Kvector_%s	" % wave_suffix)
	band_for_write.write(";DelayUpdate\n")

for y in range(starting_BAND - 1, ending_BAND):
		band_for_write.write("X ModifyGraph rgb(E_band_%03d_%s)=(65535,43690,0),lsize(E_band_%03d_%s)=2;DelayUpdate\n" % ((y + 1), wave_suffix, (y + 1), wave_suffix))

band_for_write.write("X ModifyGraph width=226.772,height=340.157;DelayUpdate\n")
band_for_write.write("X ModifyGraph tick=2;DelayUpdate\n")
band_for_write.write("X ModifyGraph mirror=1;DelayUpdate\n")
band_for_write.write("X ModifyGraph nticks(left)=5;DelayUpdate\n")
band_for_write.write("X ModifyGraph fSize=24;DelayUpdate\n")
band_for_write.write("X ModifyGraph lblMargin=15;DelayUpdate\n")
band_for_write.write("X ModifyGraph standoff=0;DelayUpdate\n")
band_for_write.write("X ModifyGraph axisOnTop=1;DelayUpdate\n")
band_for_write.write("X ModifyGraph axThick=2;DelayUpdate\n")
band_for_write.write("X ModifyGraph zero(left)=8,zeroThick(left)=2;DelayUpdate\n")
band_for_write.write("X ModifyGraph grid(bottom)=1,gridHair(bottom)=0,gridStyle(bottom)=5,gridRGB(bottom)=(0,0,0);DelayUpdate\n")
band_for_write.write('X Label left "Energy (eV)";DelayUpdate\n')
band_for_write.write('X Label bottom "Wavevector";DelayUpdate\n')
band_for_write.write("X SetAxis left -3,3\n")

oband_for_write = open("oband_structure.itx", "w")

oband_for_write.write("IGOR\n")
oband_for_write.write("WAVES/D	")
oband_for_write.write("Kvector_occ_%s	" % wave_suffix)

for y in range(starting_BAND - 1, ending_BAND):
	oband_for_write.write("E_band_%03d_occ_%s	" % ((y + 1), wave_suffix))

for y in range(starting_BAND - 1, ending_BAND):
	oband_for_write.write("O_band_%03d_occ_%s	" % ((y + 1), wave_suffix))
		
oband_for_write.write("\n")
oband_for_write.write("BEGIN\n")

for x in range(starting_KPOINT - 1, ending_KPOINT):
	oband_for_write.write("%.10f	" % kvec_mat[x])

	for y in range(starting_BAND - 1, ending_BAND):
		oband_for_write.write("%.8f	" % eigen_val_mat[y][x])

	for y in range(starting_BAND - 1, ending_BAND):
		oband_for_write.write("%.8f	" % elec_occ_mat[y][x])

	oband_for_write.write("\n")

oband_for_write.write("END\n")

oband_for_write.write("X Display	")
oband_for_write.write("as	")
oband_for_write.write('"oband_structure_occ_%s"' % wave_suffix)
oband_for_write.write(";DelayUpdate\n")

count_band = ending_BAND - starting_BAND + 1
count_start = starting_BAND - 1
count_end = starting_BAND + 4

while count_band > 5:
	oband_for_write.write("X AppendToGraph	")

	for y in range(count_start, count_end):
		oband_for_write.write("E_band_%03d_occ_%s	" % ((y + 1), wave_suffix))

	oband_for_write.write("vs	")
	oband_for_write.write("Kvector_occ_%s	" % wave_suffix)
	oband_for_write.write(";DelayUpdate\n")

	count_band -= 5
	count_start += 5
	count_end += 5

if count_start != ending_BAND:
	oband_for_write.write("X AppendToGraph	")

	for y in range(count_start, ending_BAND):
		oband_for_write.write("E_band_%03d_occ_%s	" % ((y + 1), wave_suffix))

	oband_for_write.write("vs	")
	oband_for_write.write("Kvector_occ_%s	" % wave_suffix)
	oband_for_write.write(";DelayUpdate\n")

for y in range(starting_BAND - 1, ending_BAND):
	oband_for_write.write("X ModifyGraph mode(E_band_%03d_occ_%s)=0,lsize=2,zColor(E_band_%03d_occ_%s)={O_band_%03d_occ_%s,0,1,BlueGreenOrange,1};DelayUpdate\n" % ((y + 1), wave_suffix, (y + 1), wave_suffix, (y + 1), wave_suffix))
	#oband_for_write.write("X ModifyGraph mode(E_band_%03d_occ_%s)=3,marker(E_band_%03d_occ_%s)=19,msize(E_band_%03d_occ_%s)=4.5,zColor(E_band_%03d_occ_%s)={O_band_%03d_occ_%s,0,1,BlueGreenOrange,1};DelayUpdate\n" % ((y + 1), wave_suffix, (y + 1), wave_suffix, (y + 1), wave_suffix, (y + 1), wave_suffix, (y + 1), wave_suffix))

oband_for_write.write("X ModifyGraph width=226.772,height=340.157;DelayUpdate\n")
oband_for_write.write("X ModifyGraph tick=2;DelayUpdate\n")
oband_for_write.write("X ModifyGraph mirror=1;DelayUpdate\n")
oband_for_write.write("X ModifyGraph nticks(left)=5;DelayUpdate\n")
oband_for_write.write("X ModifyGraph fSize=24;DelayUpdate\n")
oband_for_write.write("X ModifyGraph lblMargin=15;DelayUpdate\n")
oband_for_write.write("X ModifyGraph standoff=0;DelayUpdate\n")
oband_for_write.write("X ModifyGraph axisOnTop=1;DelayUpdate\n")
oband_for_write.write("X ModifyGraph axThick=2;DelayUpdate\n")
oband_for_write.write("X ModifyGraph zero(left)=8,zeroThick(left)=2;DelayUpdate\n")
oband_for_write.write("X ModifyGraph grid(bottom)=1,gridHair(bottom)=0,gridStyle(bottom)=5,gridRGB(bottom)=(0,0,0);DelayUpdate\n")
oband_for_write.write('X Label left "Energy (eV)";DelayUpdate\n')
oband_for_write.write('X Label bottom "Wavevector";DelayUpdate\n')
oband_for_write.write("X SetAxis left -3,3\n")


orbital_type = ["1_s", "2_p", "3_d", "4_t"]

for t in range(0, len(orbital_type)):
	for z in range(starting_ION - 1, ending_ION):
	#for z in range(starting_ION - 1, ending_ION + 1): ## turn on when you need sum of pbands for all elements (starting_ION, ending_ION must be defualt)
		pband_for_write = open("%02dion_pband_%s.itx" % ((z + 1), orbital_type[t]), "w")

		pband_for_write.write("IGOR\n")
		pband_for_write.write("WAVES/D	")
		pband_for_write.write("Kvector_%02dion_%s_%s	" % ((z + 1), orbital_type[t], wave_suffix))

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("E_band_%03d_%02dion_%s_%s	" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix))

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("P_band_%03d_%02dion_%s_%s	" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix))

		pband_for_write.write("\n")
		pband_for_write.write("BEGIN\n")

		for x in range(starting_KPOINT - 1, ending_KPOINT):
			pband_for_write.write("%.10f	" % kvec_mat[x])

			for y in range(starting_BAND - 1, ending_BAND):
				pband_for_write.write("%.8f	" % eigen_val_mat[y][x])

			for y in range(starting_BAND - 1, ending_BAND):
				pband_for_write.write("%.8f	" % pro_band_mat[t][z][y][x])

			pband_for_write.write("\n")

		pband_for_write.write("END\n")

		pband_for_write.write("X Display	")
		pband_for_write.write("as	")
		pband_for_write.write('"pband_%02dion_%s_%s"' % ((z + 1), orbital_type[t], wave_suffix))
		pband_for_write.write(";DelayUpdate\n")

		count_band = ending_BAND - starting_BAND + 1
		count_start = starting_BAND - 1
		count_end = starting_BAND + 4

		while count_band > 5:
			pband_for_write.write("X AppendToGraph	")

			for y in range(count_start, count_end):
				pband_for_write.write("E_band_%03d_%02dion_%s_%s	" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix))

			pband_for_write.write("vs	")
			pband_for_write.write("Kvector_%02dion_%s_%s	" % ((z + 1), orbital_type[t], wave_suffix))
			pband_for_write.write(";DelayUpdate\n")

			count_band -= 5
			count_start += 5
			count_end += 5

		if count_start != ending_BAND:
			pband_for_write.write("X AppendToGraph	")

			for y in range(count_start, ending_BAND):
				pband_for_write.write("E_band_%03d_%02dion_%s_%s	" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix))

			pband_for_write.write("vs	")
			pband_for_write.write("Kvector_%02dion_%s_%s	" % ((z + 1), orbital_type[t], wave_suffix))
			pband_for_write.write(";DelayUpdate\n")

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("X ModifyGraph mode(E_band_%03d_%02dion_%s_%s)=0,lsize=2,zColor(E_band_%03d_%02dion_%s_%s)={P_band_%03d_%02dion_%s_%s,0,0.6,Red,1};DelayUpdate\n" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix))
			#pband_for_write.write("X ModifyGraph mode(E_band_%03d_%02dion_%s_%s)=3,marker(E_band_%03d_%02dion_%s_%s)=19,zmrkSize(E_band_%03d_%02dion_%s_%s)={P_band_%03d_%02dion_%s_%s,0,1,2,6},zColor(E_band_%03d_%02dion_%s_%s)={P_band_%03d_%02dion_%s_%s,0,1,Red,1};DelayUpdate;DelayUpdate\n" % ((y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix, (y + 1) ,(z + 1), orbital_type[t], wave_suffix))

		pband_for_write.write("X ModifyGraph width=226.772,height=340.157;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph tick=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph mirror=1;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph nticks(left)=5;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph fSize=24;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph lblMargin=15;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph standoff=0;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph axisOnTop=1;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph axThick=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph zero(left)=8,zeroThick(left)=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph grid(bottom)=1,gridHair(bottom)=0,gridStyle(bottom)=5,gridRGB(bottom)=(0,0,0);DelayUpdate\n")
		pband_for_write.write('X Label left "Energy (eV)";DelayUpdate\n')
		pband_for_write.write('X Label bottom "Wavevector";DelayUpdate\n')
		pband_for_write.write("X SetAxis left -3,3\n")

element_type = vasp.get_elements_type()
element_num = vasp.get_elements_num()
element_info = vasp.get_elements_info()

print ("num_IONS: %d" % num_IONS)
print ("num_KPOINTS: %d" % num_KPOINTS)
print ("num_BANDS: %d" % num_BANDS)
print ("ref_E: %.4f eV" % ref_E)

print (element_info)

pro_band_elements_mat = [[[[float(0) for x in range(0, num_KPOINTS)] for y in range(0, num_BANDS)] for z in range(0, len(element_type))] for t in range(0, 4)]

for t in range(0, len(orbital_type)):
	element_num_start = 0
	for z in range(0, len(element_type)):
		element_num_end = element_num_start + int(element_num[z])
		for e in range(element_num_start, element_num_end):
			#print "t: %d, z: %d, e: %d" % (t, z, e) ##for checking whether iterations are going well or not
			for y in range(starting_BAND - 1, ending_BAND):
				for x in range(starting_KPOINT - 1, ending_KPOINT):
					pro_band_elements_mat[t][z][y][x] += pro_band_mat[t][e][y][x]
		element_num_start = element_num_end

for t in range(0, len(orbital_type)):
	for z in range(0, len(element_type)):
		pband_for_write = open("%s_pband_%s.itx" % (element_type[z], orbital_type[t]), "w")

		pband_for_write.write("IGOR\n")
		pband_for_write.write("WAVES/D	")
		pband_for_write.write("Kvector_%s_%s_%s	" % (element_type[z], orbital_type[t], wave_suffix))

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("E_band_%03d_%s_%s_%s	" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix))

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("P_band_%03d_%s_%s_%s	" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix))

		pband_for_write.write("\n")
		pband_for_write.write("BEGIN\n")

		for x in range(starting_KPOINT - 1, ending_KPOINT):
			pband_for_write.write("%.10f	" % kvec_mat[x])

			for y in range(starting_BAND - 1, ending_BAND):
				pband_for_write.write("%.8f	" % eigen_val_mat[y][x])

			for y in range(starting_BAND - 1, ending_BAND):
				pband_for_write.write("%.8f	" % pro_band_elements_mat[t][z][y][x])

			pband_for_write.write("\n")

		pband_for_write.write("END\n")

		pband_for_write.write("X Display	")
		pband_for_write.write("as	")
		pband_for_write.write('"pband_%s_%s_%s"' % (element_type[z], orbital_type[t], wave_suffix))
		pband_for_write.write(";DelayUpdate\n")

		count_band = ending_BAND - starting_BAND + 1
		count_start = starting_BAND - 1
		count_end = starting_BAND + 4

		while count_band > 5:
			pband_for_write.write("X AppendToGraph	")

			for y in range(count_start, count_end):
				pband_for_write.write("E_band_%03d_%s_%s_%s	" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix))

			pband_for_write.write("vs	")
			pband_for_write.write("Kvector_%s_%s_%s	" % (element_type[z], orbital_type[t], wave_suffix))
			pband_for_write.write(";DelayUpdate\n")

			count_band -= 5
			count_start += 5
			count_end += 5

		if count_start != ending_BAND:
			pband_for_write.write("X AppendToGraph	")

			for y in range(count_start, ending_BAND):
				pband_for_write.write("E_band_%03d_%s_%s_%s	" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix))

			pband_for_write.write("vs	")
			pband_for_write.write("Kvector_%s_%s_%s	" % (element_type[z], orbital_type[t], wave_suffix))
			pband_for_write.write(";DelayUpdate\n")

		for y in range(starting_BAND - 1, ending_BAND):
			pband_for_write.write("X ModifyGraph mode(E_band_%03d_%s_%s_%s)=0,lsize=2,zColor(E_band_%03d_%s_%s_%s)={P_band_%03d_%s_%s_%s,0,0.6,Red,1};DelayUpdate\n" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix))
			#pband_for_write.write("X ModifyGraph mode(E_band_%03d_%s_%s_%s)=3,marker(E_band_%03d_%s_%s_%s)=19,zmrkSize(E_band_%03d_%s_%s_%s)={P_band_%03d_%s_%s_%s,0,1,2,6},zColor(E_band_%03d_%s_%s_%s)={P_band_%03d_%s_%s_%s,0,1,Red,1};DelayUpdate;DelayUpdate\n" % ((y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix, (y + 1) ,element_type[z], orbital_type[t], wave_suffix))

		pband_for_write.write("X ModifyGraph width=226.772,height=340.157;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph tick=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph mirror=1;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph nticks(left)=5;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph fSize=24;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph lblMargin=15;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph standoff=0;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph axisOnTop=1;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph axThick=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph zero(left)=8,zeroThick(left)=2;DelayUpdate\n")
		pband_for_write.write("X ModifyGraph grid(bottom)=1,gridHair(bottom)=0,gridStyle(bottom)=5,gridRGB(bottom)=(0,0,0);DelayUpdate\n")
		pband_for_write.write('X Label left "Energy (eV)";DelayUpdate\n')
		pband_for_write.write('X Label bottom "Wavevector";DelayUpdate\n')
		pband_for_write.write("X SetAxis left -3,3\n")
