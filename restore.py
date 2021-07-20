from transforms3d.quaternions import rotate_vector
import xml.etree.ElementTree as ET
import argparse 
import transforms3d
import os

scale = 0.001
parser = argparse.ArgumentParser(description='Grasp World Files.')
parser.add_argument('--path', type=str, help='Path to the grasp file', default='/home/jiangtang/IRC/EBM_Hand/grasp_generation/output/ShadowHand2_small.xml/grasp2.xml')
parser.add_argument('--hand_type', type = str, help='Type of the hand', default='ShadowHand')
args = parser.parse_args()
path =  args.path
hand = args.hand_type
tree = ET.parse(path)
root = tree.getroot()
print(root.tag)

transform = []
for trans in root.iter('fullTransform'):
    transform.append(trans.text)
    print(trans.text)
robot_position = transform[1]
robot_trans = [i * scale for i in [*map(lambda x: float(x[1:]) if x[0]=="+" else -1 * float(x[1:]), robot_position[robot_position.find("[")+1: robot_position.find("]")].split(" "))]]
print(robot_trans)
if 'BarrettHand' in hand:
    error_in_z_axis = 0.0787185
    robot_trans[2] -= error_in_z_axis
else:
    error_in_z_axis = 0.198529 - 0.19494474812631868
    error_in_x_axis = 0.111251 - 0.11216547500617721
    error_in_y_axis = 0.011497 -  0.014274999999999998
    robot_trans[0] -= error_in_x_axis
    robot_trans[1] -= error_in_y_axis
    robot_trans[2] -= error_in_z_axis
robot_rots_quat = [*map(lambda x: float(x[1:]) if x[0]=="+" else -1 * float(x[1:]), robot_position[robot_position.find("(")+1: robot_position.find(")")].split(" "))]
print(robot_rots_quat)
robot_rots_ax = transforms3d.euler.quat2euler(robot_rots_quat, 'sxzy')
print(robot_rots_ax)
dofs = [*map(lambda x: float(x[1:]) if x[0]=="+" else -1 * float(x[1:]), [*filter(None, [*root.iter("dofValues")][0].text.split(" "))])]
print(dofs)
if 'BarrettHand' in hand: 
    A = 0.32142929 
    dofs  = [dofs[0], dofs[1], dofs[1] * A, dofs[0],  dofs[4], dofs[4] * A, dofs[7], dofs[7] * A]
else:
    A = 0.8
    dofs = [dofs[10], dofs[11], dofs[12], dofs[12] * A, 
            dofs[7], dofs[8], dofs[9], dofs[9] * A, 
            dofs[4], dofs[5], dofs[6], dofs[6] * A,
            dofs[0], dofs[1], dofs[2], dofs[3], dofs[3] * A,
            dofs[13], dofs[14], dofs[15], dofs[16], dofs[17]]
    dofs = [i * -1 for i in dofs]
elements = path.split('/')
obj_name = elements[-2].replace('_small.xml', '')
print(obj_name)
paramters = []
if obj_name == 'ShadowHand2':
    robot_trans[0] += 6.49825
    robot_trans[1] += 9.24166 
    robot_trans[2] += 0.872342

paramters.extend(robot_trans)
paramters.extend(robot_rots_ax)
paramters.extend(dofs)
result_path = path[:-3] + 'txt'
print(result_path)
paramters = " ".join([str(i) for i in paramters])
with open(result_path, 'w') as f:
    f.write(paramters)


dat_root = '/home/jiangtang/IRC/Bilevel-Differential-Grasp-Optimization/GraspDataset/' + obj_name
for file in os.listdir(dat_root):
    if file.endswith(".dat"):
        dat_path = os.path.join(dat_root, file)
        dat_scale = file.split('_')[-1].replace('.dat', '')
        print(dat_path)
if 'BarrettHand' in hand:
    cmd_line = '/home/jiangtang/IRC/Bilevel-Differential-Grasp-Optimization/data/BarrettHand/bh280.urdf 200 '+dat_path+' '+obj_name+' '+dat_scale+' 2 1 ./ ' + result_path
else:
    cmd_line = '/home/jiangtang/IRC/Bilevel-Differential-Grasp-Optimization/data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 '+dat_path+' '+obj_name+' '+dat_scale+' 2 1 ./ ' + result_path

print(cmd_line)