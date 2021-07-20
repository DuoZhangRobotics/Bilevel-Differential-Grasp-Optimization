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
    robot_ori_tran = [1000, 0, 0]

    elements = re.split('/', xml_name)
    # print(elements[len(elements)-1])
    sub_path = gen_scripts_path + '/' + elements[len(elements)-1]
    mkdir(root_path + '/' + sub_path)

    if os.path.exists(sub_path + '/run.sh'):
        print('Detected Script: ' + root_path + '/' + sub_path + '/run.sh')
    else:
        print('Generating Script: ' + root_path + '/' + sub_path + '/run.sh')

        f = open(root_path + '/' + sub_path + '/run.sh', 'w')

        script = '#!/bin/bash\n'
        script += '#SBATCH --job-name=' + elements[len(elements)-1] + '\n'
        script += '#SBATCH --cpus-per-task=2\n'
        script += '#SBATCH --time=10-00:00:00\n'
        script += '#SBATCH --output=%x.out\n'
        script += '#SBATCH --error=%x.err\n'
        script += '#SBATCH --ntasks=1\n'
        script += '#SBATCH --mem=1G\n'

        for j in range(5):
            for k in range(5):
                for m in range(5):

                    if os.path.exists(root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '_' + str(j) + '_' + str(k) + '_' + str(m)):
                        print('Detected Filefold: ' + root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '_' + str(j) + '_' + str(k) + '_' + str(m))
                    else:
                        mkdir(root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '_' + str(j) + '_' + str(k) + '_' + str(m))

                    robot_zrot = math.pi / (j + 1)  # radians
                    robot_yrot = math.pi / (k + 1)
                    robot_xrot = math.pi / (m + 1)
                    # print(robot_xrot, robot_yrot, robot_zrot)
                    # print('--------------------------------')
                    robot_R = transforms3d.taitbryan.euler2mat(robot_zrot, robot_yrot, robot_xrot)  # rotations
                    robot_coor_T_R = np.matmul(robot_R, robot_ori_tran)
                    robot_quat2mat = transforms3d.quaternions.quat2mat(robot_ori_quat)
                    robot_quat2mat_R = np.matmul(robot_R, robot_quat2mat)
                    robot_mat2quat = transforms3d.quaternions.mat2quat(robot_quat2mat_R)

                    script += root_path + '/' + built_binary_path + ' '
                    script += '--bodyFile=' + xml_name + ' '
                    script += '--bodyRot=1,0,0,0 '
                    script += '--bodyTrans=0,0,0 '
                    script += '--robotFile=' + root_path + '/' + robot_xml_path + ' '
                    script += '--robotRot=' + str(robot_mat2quat[0]) + ',' + str(robot_mat2quat[1]) + ',' + str(robot_mat2quat[2]) + ',' + str(robot_mat2quat[3]) + ' '
                    script += '--robotTrans=' + str(robot_coor_T_R[0]) + ',' + str(robot_coor_T_R[1]) + ',' + str(robot_coor_T_R[2]) + ' '
                    script += '--robotDOF=0,0,0,0 '
                    script += '--resultFile=' + root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '_' + str(j) + '_' + str(k) + '_' + str(m) + '/\n'
                    script += 'cp -r ' + root_path + '/' + sub_path + '/' + elements[len(elements)-1] + '_' + str(j) + '_' + str(k) + '_' + str(m) + ' '
                    script += root_path + '/' + result_path + '\n\n'

        f.write(script)
        f.close()

    running = find_running()

    if os.path.exists(root_path + '/' + result_path + elements[len(elements)-1] + '_' + str(4) + '_' + str(4) + '_' + str(4)):
        print('Dataset Generated!')
    elif elements[len(elements)-1] in running:
        print('Running Script ' + elements[len(elements)-1] + '!')
    elif not testRun:
        subprocess_cmd('cd ' + root_path + '/' + sub_path + '; sbatch run.sh;')
    else:
        print('Test Run!')

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
    parser.add_argument('--object_xml_path',default='good_shapes')
    parser.add_argument('--robot_xml_path',default='graspitmodified_lm/graspit/models/robots/BarrettBH8_280/BarrettBH8_280.xml')
    parser.add_argument('--built_binary_path',default='graspitmodified-build/graspit/graspit_cmdline')
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

