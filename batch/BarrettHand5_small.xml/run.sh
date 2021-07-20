#!/bin/bash
#SBATCH --job-name=BarrettHand5_small.xml
#SBATCH --cpus-per-task=2
#SBATCH --time=10-00:00:00
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --ntasks=1
#SBATCH --mem=1G
/home/jiangtang/IRC/EBM_Hand/grasp_generation/~/graspitmodified-build/graspit/graspit_cmdline --bodyFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/GraspDataset/BarrettHand5_small.xml --bodyRot=1,0,0,0 --bodyTrans=0,0,0 --robotFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/graspitmodified_lm/graspit/models/robots/BarrettBH8_280/BarrettBH8_280.xml --robotRot=0.010777637540842464,-0.3300423243118613,-0.6160658312270594,-0.7151355104358714 --robotTrans=41.357,-279.028,-32.8785 --robotDOF=0.763057, 0.00230305, 0.537668, 0.604613, --resultFile=/home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand5_small.xml/BarrettHand5_small.xml/
cp -r /home/jiangtang/IRC/EBM_Hand/grasp_generation/batch/BarrettHand5_small.xml/BarrettHand5_small.xml /home/jiangtang/IRC/EBM_Hand/grasp_generation/output

