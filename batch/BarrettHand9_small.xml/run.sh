#!/bin/bash
#SBATCH --job-name=BarrettHand9_small.xml
#SBATCH --cpus-per-task=2
#SBATCH --time=10-00:00:00
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --ntasks=1
#SBATCH --mem=1G
/home/jiangtang/IRC/EBM_Hand/grasp_generation/~/graspitmodified-build/graspit/graspit_cmdline --bodyFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/GraspDataset/BarrettHand9_small.xml --bodyRot=1,0,0,0 --bodyTrans=0,0,0 --robotFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/graspitmodified_lm/graspit/models/robots/BarrettBH8_280/BarrettBH8_280.xml --robotRot=0.7071080798594737,0.0,0.0,-0.7071054825112364 --robotTrans=0.0,-200.0,78.71849999999999 --robotDOF=0.5, 0.0, 0.0, 0.0, --resultFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand9_small.xml/BarrettHand9_small.xml/
cp -r /home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand9_small.xml/BarrettHand9_small.xml /home/jiangtang/IRC/EBM_Hand/grasp_generation/output

