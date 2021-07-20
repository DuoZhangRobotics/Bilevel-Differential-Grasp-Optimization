#!/bin/bash
#SBATCH --job-name=BarrettHand1_small.xml
#SBATCH --cpus-per-task=2
#SBATCH --time=10-00:00:00
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --ntasks=1
#SBATCH --mem=1G
/home/jiangtang/IRC/EBM_Hand/grasp_generation/~/graspitmodified-build/graspit/graspit_cmdline --bodyFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/GraspDataset/BarrettHand1_small.xml --bodyRot=1,0,0,0 --bodyTrans=0,0,0 --robotFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/graspitmodified_lm/graspit/models/robots/BarrettBH8_280/BarrettBH8_280.xml --robotRot=0.3938676265373032,-0.04245897418232365,-0.024391529091924546,-0.9178619621630555 --robotTrans=-6.49268,-2.19252,-62.57750000000001 --robotDOF=0.357541, 1.03552, 1.05205, 0.991532, --resultFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand1_small.xml/BarrettHand1_small.xml/
cp -r /home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand1_small.xml/BarrettHand1_small.xml /home/jiangtang/IRC/EBM_Hand/grasp_generation/output

