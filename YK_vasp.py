#!/usr/bin/python

#First Created : 2015/01/23 by Youngkwang
#Last Modified : 2017/12/22 by Youngkwang

import math
import numpy as np
from operator import itemgetter


#### General

def get_keyword_from(keyword, target):
	with open(target, "r") as target:
		for line in target:
			if keyword in line:
				return line

def get_keywords_from(keyword, target):
	with open(target, "r") as target:
		lines = []
		for line in target:
			if keyword in line:
				lines.append(line)
		return lines

def get_line_from(line, target):
	with open(target, "r") as target:
		target_lines = target.readlines()
		return target_lines[line - 1]	


#### CONTCAR related

def get_elements_type(target = "CONTCAR"):
	with open(target, "r") as CONTCAR:
		CONTCAR_lines = CONTCAR.readlines() 
		elements_type = CONTCAR_lines[5].split()

		return elements_type

def get_elements_num(target = "CONTCAR"):
	with open(target, "r") as CONTCAR:
		CONTCAR_lines = CONTCAR.readlines() 
		elements_num = CONTCAR_lines[6].split()

		return elements_num

def get_elements_info(target = "CONTCAR"):
	elements_info = {}
	for i in range(0, len(get_elements_type())):
		elements_info[get_elements_type()[i]] = get_elements_num()[i]

	return elements_info

def get_cell_vectors(target = "CONTCAR"):
	with open(target, "r") as CONTCAR:
		CONTCAR_lines = CONTCAR.readlines()

		scale_factor = float(CONTCAR_lines[1].split()[0])

		cell_vector_a_temp = CONTCAR_lines[2].split()
		cell_vector_b_temp = CONTCAR_lines[3].split()
		cell_vector_c_temp = CONTCAR_lines[4].split()

		cell_vector_a = [float(0) for x in range(0, len(cell_vector_a_temp))]
		cell_vector_b = [float(0) for x in range(0, len(cell_vector_b_temp))]
		cell_vector_c = [float(0) for x in range(0, len(cell_vector_c_temp))]

		for i in range(0, len(cell_vector_a)):
			cell_vector_a[i] = float(cell_vector_a_temp[i]) * scale_factor
		for i in range(0, len(cell_vector_b)):
			cell_vector_b[i] = float(cell_vector_b_temp[i]) * scale_factor
		for i in range(0, len(cell_vector_c)):
			cell_vector_c[i] = float(cell_vector_c_temp[i]) * scale_factor

		cell_vectors = {'a' : cell_vector_a, 'b' : cell_vector_b, 'c' :cell_vector_c}

		return cell_vectors

def get_lattice_constants(target = "CONTCAR"):
	with open(target, "r") as CONTCAR:
		CONTCAR_lines = CONTCAR.readlines()

		scale_factor = float(CONTCAR_lines[1].split()[0])

		cell_vector_a_temp = CONTCAR_lines[2].split()
		cell_vector_b_temp = CONTCAR_lines[3].split()
		cell_vector_c_temp = CONTCAR_lines[4].split()

		cell_vector_a = [float(0) for x in range(0, len(cell_vector_a_temp))]
		cell_vector_b = [float(0) for x in range(0, len(cell_vector_b_temp))]
		cell_vector_c = [float(0) for x in range(0, len(cell_vector_c_temp))]

		for i in range(0, len(cell_vector_a)):
			cell_vector_a[i] = float(cell_vector_a_temp[i]) * scale_factor
		for i in range(0, len(cell_vector_b)):
			cell_vector_b[i] = float(cell_vector_b_temp[i]) * scale_factor
		for i in range(0, len(cell_vector_c)):
			cell_vector_c[i] = float(cell_vector_c_temp[i]) * scale_factor

		lattice_constant_a = math.sqrt(cell_vector_a[0] ** 2 + cell_vector_a[1] ** 2 + cell_vector_a[2] ** 2)
		lattice_constant_b = math.sqrt(cell_vector_b[0] ** 2 + cell_vector_b[1] ** 2 + cell_vector_b[2] ** 2)
		lattice_constant_c = math.sqrt(cell_vector_c[0] ** 2 + cell_vector_c[1] ** 2 + cell_vector_c[2] ** 2)

		lattice_constants = {'a0' : lattice_constant_a, 'b0' : lattice_constant_b, 'c0' : lattice_constant_c}

		return lattice_constants

def get_atoms_position(target = "CONTCAR", position_type = "defualt"):
	with open(target, "r") as CONTCAR:
		CONTCAR_lines = CONTCAR.readlines()

		elements_type = get_elements_type()
		elements_num = get_elements_num()
		elements_info = get_elements_info()
		cell_vectors = get_cell_vectors()

		if list(CONTCAR_lines[7].split()[0])[0] == "S" or list(CONTCAR_lines[7].split()[0])[0] == "s":
			first_atom_line = 9
		else:
			first_atom_line = 8

		atoms_position={}

		for i in range(0, len(elements_type)):
			atoms_position[elements_type[i]] = [[float(0) for x in range(0, 3)] for y in range(0, int(elements_num[i]))]

		element_position_line = first_atom_line

		if position_type == "defualt":
			for i in range(0, len(elements_type)):
				for j in range(element_position_line, element_position_line + int(elements_num[i])):
					atoms_position[elements_type[i]][j - element_position_line] = [float(CONTCAR_lines[j].split()[k]) for k in range (0, 3)]
				element_position_line = element_position_line + int(elements_num[i])

		elif position_type == "c":
			if list(CONTCAR_lines[first_atom_line - 1].split()[0])[0] == "D" or list(CONTCAR_lines[first_atom_line - 1].split()[0])[0] == "d":
				for i in range(0, len(elements_type)):
					for j in range(element_position_line, element_position_line + int(elements_num[i])):
						atoms_position[elements_type[i]][j - element_position_line] = np.asarray(cell_vectors["a"]) * float(CONTCAR_lines[j].split()[0]) + np.asarray(cell_vectors["b"]) * float(CONTCAR_lines[j].split()[1]) + np.asarray(cell_vectors["c"]) * float(CONTCAR_lines[j].split()[2])
						atoms_position[elements_type[i]][j - element_position_line] = atoms_position[elements_type[i]][j - element_position_line].tolist()
					element_position_line = element_position_line + int(elements_num[i])

		for i in atoms_position:
			atoms_position[i] = sorted (atoms_position[i], key = itemgetter(2,0,1))



		return atoms_position

#print get_atoms_position(position_type = "c")


#### OUTCAR related

def get_VBM_band_num(target = "OUTCAR"):
	with open(target, "r") as OUTCAR:
		VBM_band_num = int(float(get_keyword_from("NELEC", target).split()[2]) / 2)
		return VBM_band_num
	
def get_CBM_band_num(target = "OUTCAR"):
	CBM_band_num = int(get_VBM_band_num(target) + 1)
	return CBM_band_num


#def get_VBM_value(target = "OUTCAR"):


#### Igor related

def turn_G_into_Gamma(axis_let_mat):
	for i in range(0, len(axis_let_mat)):
		if axis_let_mat[i] == 'G':
			axis_let_mat[i] = " \F'Symbol'" + axis_let_mat[i] + " "
		else:
			axis_let_mat[i] = " " + axis_let_mat[i] + "\F'Symbol' "

	return axis_let_mat











