from sys import path
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description='Convert Obj to .off and generate Graspit xml files')
parser.add_argument('--inputPath', type=str, help='Path to the maxRange_Scale.txt file')

args = parser.parse_args()
maxRange_scale_path = args.inputPath
path_elements = maxRange_scale_path.split('/')
path_elements[-1] = path_elements[-2] + '.obj'
inputPath = '/'.join(path_elements)

vertices = []
faces = ''
with open(inputPath, 'r') as fin:
    for line in fin:
        if line[0] == 'v' or line[0:2]=='vt':
            try:
                vertices.append([*map(lambda x: -1 * float(x[1:]) if x[0]=="-" else float(x), [i for i in line.replace('\n', '').split(" ")[1:] if i != ''])])
            except IndexError:
                print(line)
                print([i for i in line.replace('\n', '').split(" ")[1:] if i != ''])
                exit(1)
        if line[0] == "f":
            faces += line
try:
    with open(maxRange_scale_path, 'r') as f:
        content = f.readlines()
        scale = float(content[1])
        if path_elements[-2] == 'ShadowHand10' or path_elements[-2]== 'BarrettHand10':
            scale=1.0
        maxRange = np.array([*map(lambda x: -1 * float(x[1:]) if x[0]=="-" else float(x), content[0].split(" ")[:-1])])
        # print(scale)
        # print(maxRange)
except FileNotFoundError:
    print(f"PLEASE RUN \'{path_elements[-2]+'.sh'}\' FIRST TO GET THE \'maxRange_Scale.txt\' FILE!")
    exit(1)
v = np.array(vertices)
# print(np.max(v[:,0]), np.min(v[:,0]))
# print(np.max(v[:,1]), np.min(v[:,1]))
# print(np.max(v[:,2]), np.min(v[:,2]))
v *= scale

centroid_before = np.array([np.max(v[:,0]), np.max(v[:,1]), np.max(v[:,2])])
centroid_after = maxRange
movement = centroid_after - centroid_before

v += movement
# print(f'movement = {movement}')
# print(np.max(v[:,0]), np.min(v[:,0]))
# print(np.max(v[:,1]), np.min(v[:,1]))
# print(np.max(v[:,2]), np.min(v[:,2]))  

path_elements[-1] = path_elements[-2] + '_small_tmp.obj'
if path_elements[-2] == 'ShadowHand2':
    # print(np.array([np.max(v[:,0]), np.max(v[:,1]), np.max(v[:,2])]))
    v -= np.array([np.max(v[:,0]), np.max(v[:,1]), np.max(v[:,2])])
# print(f'movement = {movement}')
# print(np.max(v[:,0]), np.min(v[:,0]))
# print(np.max(v[:,1]), np.min(v[:,1]))
# print(np.max(v[:,2]), np.min(v[:,2]))    

obj_output_path = '/'.join(path_elements)
with open(obj_output_path, 'w') as fout:
    for i in range(v.shape[0]):
        line = f'v {v[i, 0]} {v[i, 1]} {v[i, 2]}\n'
        fout.write(line)
    fout.write(faces)
v *= 1000

path_elements[-1] = path_elements[-2] + '_small.obj'
obj_output_path = '/'.join(path_elements)
with open(obj_output_path, 'w') as fout:
    for i in range(v.shape[0]):
        line = f'v {v[i, 0]} {v[i, 1]} {v[i, 2]}\n'
        fout.write(line)
    fout.write(faces)
print("SCALED .OBJ FILE GENERATED!")
def obj2off(objpath, offpath):
    '''
    Convert obj file to off file
         :param objpath: path to the .obj file
         :param offpath: the save address of the path of the .off file
         :return: none
    '''
    line = ""

    vset = []
    fset = []
    with open(objpath,'r') as f:
        lines = f.readlines()
    p = re.compile(r'/+')
    space = re.compile(r' +')

    for line in lines:
                 #Get a line in the obj file as a string
        tailMark = " "
        line = line+tailMark
        if line[0]!='v' and line[0]!='f' :
            continue

        parameters = space.split(line.strip())
        if parameters[0] == "v": #if it is a vertex
            Point = []
            Point.append(eval( parameters[1]) )
            Point.append(eval( parameters[2]) )
            Point.append(eval( parameters[3]) )
            vset.append(Point)

        elif parameters[0] == "f": #if it is a face, the index of the vertices is stored
            vIndexSets = [] # collection of temporary storage points
            for i in range(1,len(parameters) ):
                x = parameters[i]
                ans = p.split(x)[0]
                index = eval(ans)
                index -= 1 #Because the vertex index starts at 1 in the obj file, and the vertex we store starts at 0, so we want to subtract 1
                vIndexSets.append(index)

                fset.append(vIndexSets)

    with open(offpath, 'w') as out:
        out = open(offpath, 'w')
        out.write("OFF\n")
        out.write(str(vset.__len__()) + " " + str(fset.__len__()) + " 0\n")
        for j in range(len(vset)):
            out.write(str(vset[j][0]) + " " + str(vset[j][1]) + " " + str(vset[j][2]) + "\n")

        for i in range(len(fset)):
            s = str(len( fset[i] ))
            for j in range( len( fset[i] ) ):
                s = s+ " "+ str(fset[i][j])
            s += "\n"
            out.write(s)

        print("{} converts to {} success!".format( p.split(objpath)[-1], p.split(offpath)[-1] ))

path_elements[-1] = path_elements[-2] + '_small.off'
off_path ='/'.join(path_elements[:-3]) + '/GraspItDataset/' + path_elements[-2] + '_small.off'
obj2off(obj_output_path, off_path)
xml_path = '/'.join(path_elements[:-3]) + '/GraspItDataset/' + path_elements[-2] + '_small.xml'
path_elements[-1] = 'initialParameters.txt'
init_param_path1 = '/'.join(path_elements)
init_param_path2 = '/'.join(path_elements[:-3]) + '/GraspItDataset/' + path_elements[-2] + '_small_initialParameters.txt'

with open(xml_path, 'w') as xml:
    xml.write(f'<root>\n\t<material>plastic</material>/n/t<mass>100</mass>\n\t<cog>0 0 0</cog>\n\t<geometryFile type="off">{path_elements[-2]+"_small.off"}</geometryFile>\n</root>')
with open(init_param_path1, 'r') as params:
    content = params.readlines()
    with open(init_param_path2, 'w') as out_params:
        out_params.write(content[0])

print("GRASPIT REQUIRED FILES GENERATED!")