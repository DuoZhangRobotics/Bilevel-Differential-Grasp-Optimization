import os,random,transforms3d,math,re,sys,subprocess
from shutil import copyfile
import numpy as np
import argparse

#from itertools import imap

def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)


def anyTrue(predicate, sequence):
    #return True in imap(predicate, sequence)
    for s in sequence:
        if predicate(s):
            return True
    return False    


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
       os.makedirs(path)


def filterFiles(folder, exts, list):
    for fileName in os.listdir(folder):
        if os.path.isdir(folder + '/' + fileName):
            filterFiles(folder + '/' + fileName, exts, list)
        elif anyTrue(fileName.endswith, exts):
            list.append(folder + '/' + fileName)


def find_running():
    username = subprocess.check_output(['whoami']).decode().splitlines()[0].replace('\"', '').replace(' ', '')
    running = []
    lines = subprocess.check_output(['squeue', '-u', username, '--format=\"%.100j\"']).decode().splitlines()
    for l in lines:
        running.append(l.replace('\"', '').replace(' ', ''))
    return running


def gen_script_for_each_xml(testRun, root_path, robot_xml_path, xml_name, built_binary_path, gen_scripts_path, result_path):
    robot_ori_quat = [1, 0, 0, 0]
    # robot_ori_tran = [1000, 0, 0]
    scale=0.001
    # scale = 1
    elements = re.split('/', xml_name)
    print(elements[len(elements)-1])
    element_name =  elements[len(elements)-1].split('.')[0]
    if 'BarrettHand' in element_name:
        error_in_z_axis = 0.0787185
    else:
        error_in_z_axis = 0.198529 - 0.19494474812631868
        error_in_x_axis = 0.111251 - 0.11216547500617721
        error_in_y_axis = 0.011497 -  0.014274999999999998
    Dataset_name = elements[-2]
    print("element name", element_name)
    print("dataset name", Dataset_name)
    initial_parameters = []
    with open(Dataset_name+'/'+element_name+'_initialParameters.txt') as f:
        initial_parameters.extend([*map(float, f.readline().split())])
    initial_parameters = [-0.00649268, -0.00219252, -0.141296, -0.0784382, 0.0587627, 0.813007, 0.357541, 1.03552, 0.332847, 0.357541, 1.05205, 0.33816, 0.991532, 0.318707 ]
    # initial_parameters = [0 for i in range(len(initial_parameters))]
    # ip = '0.0377626 0.0242435 0.0358731 1.29584 -0.376385 -1.10845 0.352237 -1.56747 -0.887228 -0.495569 0.314999 -1.4011 -0.816813 -0.35346 -0.0760167 -1.24467 -1.41401 -0.192404 -0.00649681 0.436332 -1.34866 -1.42993 -0.126254 -0.452586 -0.539669 -0.2618 -0.5237 -1.5708'
    
    # initial_parameters = [*map(float, ip.split(' '))]
    print(initial_parameters)
    if 'BarrettHand' in element_name:
        robot_ori_tran = initial_parameters[:3]
        robot_ori_tran[2] += error_in_z_axis
        robot_ori_tran = [i/scale for i in robot_ori_tran]
    else:
        robot_ori_tran = initial_parameters[:3]
        # robot_ori_tran = [0,0,0]
        robot_ori_tran[0] += error_in_x_axis
        robot_ori_tran[1] += error_in_y_axis
        robot_ori_tran[2] += error_in_z_axis

        robot_ori_tran = [i/scale for i in robot_ori_tran]
    # robot_ori_tran = [0.3, 0.1 , -0.2]
    sub_path = gen_scripts_path + '/' + elements[len(elements)-1]
    mkdir(root_path + '/' + sub_path)
    print(root_path, sub_path, elements)
    if os.path.exists(sub_path + '/run.sh'):
        print('Detected Script: ' + root_path + '/' + sub_path + '/run.sh')
    else:
        print('Generating Script: ' + root_path + '/' + sub_path + '/run.sh')

        f = open(root_path + '/' + sub_path + '/run.sh', 'w')

        script = '#!/bin/bash\n'
        # script += '#SBATCH --job-name=' + elements[len(elements)-1] + '\n'
        # script += '#SBATCH --cpus-per-task=2\n'
        # script += '#SBATCH --time=10-00:00:00\n'
        # script += '#SBATCH --output=%x.out\n'
        # script += '#SBATCH --error=%x.err\n'
        # script += '#SBATCH --ntasks=1\n'
        # script += '#SBATCH --mem=1G\n'



        if os.path.exists(root_path + '/' + sub_path + '/' + elements[len(elements)-1]):
            print('Detected Filefold: ' + root_path + '/' + sub_path + '/' + elements[len(elements)-1])
        else:
            mkdir(root_path + '/' + sub_path + '/' + elements[len(elements)-1])

            # robot_zrot = math.pi / (j + 1)  # radians
            # robot_yrot = math.pi / (k + 1)
            # robot_xrot = math.pi / (m + 1)
            if 'BarrettHand' in element_name:
                robot_zrot = math.pi + initial_parameters[5]
            else:
                robot_zrot = initial_parameters[5]
            robot_yrot = initial_parameters[4]
            robot_xrot = initial_parameters[3]
            if 'BarrettHand' in element_name:
                robot_dof = [initial_parameters[6], initial_parameters[7], initial_parameters[10], initial_parameters[12]]
            else:
                robot_dof = [initial_parameters[18], initial_parameters[19], initial_parameters[20], initial_parameters[21],
                             initial_parameters[14], initial_parameters[15], initial_parameters[16],
                             initial_parameters[10], initial_parameters[11], initial_parameters[12],
                             initial_parameters[6], initial_parameters[7], initial_parameters[8],
                             initial_parameters[23], initial_parameters[24], initial_parameters[25], initial_parameters[26], initial_parameters[27]]
                robot_dof = [-1 * i for i in robot_dof]
                # robot_xrot += math.pi
                # robot_yrot -= math.pi
            print("Rot = ", initial_parameters[3:6])
                    # print(robot_xrot, robot_yrot, robot_zrot)
                    # print('--------------------------------')
            R1 = transforms3d.taitbryan.euler2mat(0,0, robot_xrot)
            R2 = transforms3d.taitbryan.euler2mat(robot_zrot,0,0 )
            R3 = transforms3d.taitbryan.euler2mat(0, robot_yrot,0)
            robot_R = R3 @ R2 @ R1
            # robot_R = transforms3d.taitbryan.euler2mat(robot_zrot, robot_yrot, robot_xrot)  # rotations
            # robot_R = transforms3d.taitbryan.euler2mat(0, 0, 0)  # rotations
            print("robot_R = ", robot_R)
            if 'ShadowHand' in element_name:
                x_rot = math.pi/2
                y_rot = 0
                z_rot = math.pi / 2
                R1 = transforms3d.taitbryan.euler2mat(z_rot, 0 ,0)
                R2 = transforms3d.taitbryan.euler2mat(0, 0, x_rot)
                R = R2@R1
                robot_R = robot_R @ R
                # robot_xrot += x_rot
                # robot_yrot += y_rot

            print("robot_R = \n", robot_R)
            # robot_coor_T_R = np.matmul(robot_R, robot_ori_tran)
            robot_coor_T_R = robot_ori_tran

            robot_quat2mat = transforms3d.quaternions.quat2mat(robot_ori_quat)
            print("robot_quat2mat = \n", robot_quat2mat)
            robot_quat2mat_R = np.matmul(robot_R, robot_quat2mat)
            robot_mat2quat = transforms3d.quaternions.mat2quat(robot_quat2mat_R)
            print("robot_mat2quat = ", robot_mat2quat)
            # script += root_path + '/' + built_binary_path + ' '
            script +=  built_binary_path + ' '
            script += '--bodyFile=' + xml_name + ' '
            script += '--bodyRot=1,0,0,0 '
            script += '--bodyTrans=0,0,0 '
            if "BarrettHand" in element_name:
                script += '--robotFile=' + root_path + '/' + robot_xml_path + ' '
            else:
                script += '--robotFile=' + root_path + '/' + 'graspitmodified_lm/graspit/models/robots/ShadowHand/ShadowHandSimple.xml' + ' '
            script += '--robotRot=' + str(robot_mat2quat[0]) + ',' + str(robot_mat2quat[1]) + ',' + str(robot_mat2quat[2]) + ',' + str(robot_mat2quat[3]) + ' '
            script += '--robotTrans=' + str(robot_coor_T_R[0]) + ',' + str(robot_coor_T_R[1]) + ',' + str(robot_coor_T_R[2]) + ' '
            # script += '--robotDOF=0,0,0,0 '
            script += '--robotDOF='
            for value in robot_dof:
                script += str(value)+', '
            script += '--resultFile=' + root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '/\n'
            script += 'cp -r ' + root_path + '/' + sub_path + '/' + elements[len(elements)-1] + ' '
            script += root_path + '/' + result_path + '\n\n'
            print(script)
        f.write(script)
        f.close()

    # running = find_running()

    # if os.path.exists(root_path + '/' + result_path + elements[len(elements)-1] + '_' + str(4) + '_' + str(4) + '_' + str(4)):
    #     print('Dataset Generated!')
    # elif elements[len(elements)-1] in running:
    #     print('Running Script ' + elements[len(elements)-1] + '!')
    # elif not testRun:
    #     subprocess_cmd('cd ' + root_path + '/' + sub_path + '; sbatch run.sh;')
    # else:
    #     print('Test Run!')

    # print(sub_path)


def gen_all_scripts(testRun, root_path, object_xml_path, robot_xml_path, built_binary_path, gen_scripts_path, final_result_path):

    if not os.path.exists(root_path + '/' + gen_scripts_path):
        os.mkdir(root_path + '/' + gen_scripts_path)

    if not os.path.exists(root_path + '/' + final_result_path):
        os.mkdir(root_path + '/' + final_result_path)

    xml_file_list = []
    filterFiles(root_path + '/' + object_xml_path, '.xml', xml_file_list)
    # print(xml_file_list)

    # for i in range(len(xml_file_list)):
    for i in range(len(xml_file_list)):
        gen_script_for_each_xml(testRun, root_path, robot_xml_path, xml_file_list[i], built_binary_path,
                                gen_scripts_path, final_result_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=os.path.basename(__file__))
    parser.add_argument('--testRun', default=False,action='store_true')
    parser.add_argument('--root_path',default=os.path.dirname(os.path.realpath(__file__)))
    parser.add_argument('--object_xml_path',default='GraspItDataset')
    parser.add_argument('--robot_xml_path',default='graspitmodified_lm/graspit/models/robots/BarrettBH8_280/BarrettBH8_280.xml')
    parser.add_argument('--built_binary_path',default='/home/jiangtang/graspitmodified-build/graspit/graspit_cmdline')
    parser.add_argument('--gen_scripts_path',default='batch')
    parser.add_argument('--final_result_path',default='output')
    args = parser.parse_args()
    gen_all_scripts(args.testRun, 
                    args.root_path,  
                    args.object_xml_path,
                    args.robot_xml_path,
                    args.built_binary_path, 
                    args.gen_scripts_path, 
                    args.final_result_path)

