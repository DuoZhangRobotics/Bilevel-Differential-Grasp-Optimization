import os
import random
import transforms3d
import math
import re
import sys
import subprocess
from shutil import copyfile
import shutil
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

obj_path = '/homes/mliu819/output_data'
xml_file_list = os.listdir(obj_path)

#print(xml_file_list)
index = 0
for i in range(len(xml_file_list)):
    xml_name = xml_file_list[i]
    print(xml_name)
    folder = os.listdir(obj_path + '/' + xml_name)
    if len(folder) < 6:
       shutil.rmtree(obj_path + '/' + xml_name)
       index = index +1
       print(index)       
