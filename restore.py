import xml.etree.ElementTree as ET
import argparse 
import transforms3d

scale = 0.001
A = 0.32142929
parser = argparse.ArgumentParser(description='Grasp World Files.')
parser.add_argument('--path', type=str, help='Path to the grasp file', default='/home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand2_small.xml/BarrettHand2_small.xml/grasp0.xml')

args = parser.parse_args()
path =  args.path

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
robot_rots_quat = [*map(lambda x: float(x[1:]) if x[0]=="+" else -1 * float(x[1:]), robot_position[robot_position.find("(")+1: robot_position.find(")")].split(" "))]
print(robot_rots_quat)
robot_rots_ax, theta = transforms3d.quaternions.quat2axangle(robot_rots_quat)
print(theta, robot_rots_ax)
dofs = [*map(lambda x: float(x[1:]) if x[0]=="+" else -1 * float(x[1:]), [*filter(None, [*root.iter("dofValues")][0].text.split(" "))])]
print(dofs)
dofs  = [dofs[0], dofs[1], dofs[1] * A, dofs[0],  dofs[4], dofs[4] * A, dofs[7], dofs[7] * A]
paramters = []
paramters.extend(robot_trans)
paramters.extend(robot_rots_ax)
paramters.extend(dofs)
result_path = path[:-3] + 'txt'
print(result_path)
paramters = " ".join([str(i) for i in paramters])
with open(result_path, 'w') as f:
    f.write(paramters)