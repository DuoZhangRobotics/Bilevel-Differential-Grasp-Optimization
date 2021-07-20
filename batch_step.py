import os
import random
import transforms3d
import math
import re
import sys
import subprocess
from shutil import copyfile
import numpy as np
from itertools import imap


def subprocess_cmd(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)


def anyTrue(predicate, sequence):
    return True in imap(predicate, sequence)


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)


def filterFiles(folder, exts, list):
    for fileName in os.listdir(folder):
        if os.path.isdir(folder + '/' + fileName):
            filterFiles(folder + '/' + fileName, exts, list)
        elif anyTrue(fileName.endswith, exts):
            # print fileName
            list.append(folder + '/' + fileName)


running = []

main_path = './batch'
if not os.path.exists(main_path):
    os.mkdir(main_path)

obj_path = '/homes/mliu819/grasp_off_xml'
xml_file_list = []
filterFiles(obj_path, '.xml', xml_file_list)

print(xml_file_list)
index = 0
for i in range(len(xml_file_list)):
    xml_name = xml_file_list[i]

    robot_ori_quat = [1, 0, 0, 0]
    robot_ori_tran = [1000, 0, 0]

    for j in range(5):
        for k in range(5):
            for m in range(5):
                robot_zrot = math.pi/(j+1)  # radians
                robot_yrot = math.pi/(k+1)
                robot_xrot = math.pi/(m+1)
                # print(robot_xrot, robot_yrot, robot_zrot)
                # print('--------------------------------')
                robot_R = transforms3d.taitbryan.euler2mat(robot_zrot, robot_yrot, robot_xrot)  # rotations
                robot_coor_T_R = np.matmul(robot_R, robot_ori_tran)
                robot_quat2mat = transforms3d.quaternions.quat2mat(robot_ori_quat)
                robot_quat2mat_R = np.matmul(robot_R, robot_quat2mat)
                robot_mat2quat = transforms3d.quaternions.mat2quat(robot_quat2mat_R)
                elements = re.split('/', xml_name)
                print(elements[4])
                sub_path = main_path + '/' + elements[4] + '_' + str(j) + '_' + str(k) + '_' + str(m)

                # script
                target_path = '/homes/mliu819/output_data' + '/' + elements[4] + '_' + str(j) + '_' + str(k) + '_' + str(m) 
                folder = os.path.exists(target_path)
                print(i)
                if not folder:
                   index = index +1 
                   if index < 5000:
                       subprocess_cmd('cd ' + sub_path + '; sbatch run.sh; cd ../..;')
                       print(sub_path)
                   print('index: ',index)
                   print(sub_path)
                print(i)
