#! /bin/bash
# check if directory ./limist exists. if so, delete it
if [ -d "./limits" ]
then
    rm -r ./limits
fi

# check if directory ./__pycache exists. if so, delete it
if [ -d "./__pycache__" ]
then
    rm -r ./__pycache__
fi
# check if dirctory ./.ipynb_checkpoints exists. if so, delete it
if [ -d "./.ipynb_checkpoints" ]
then 
    rm -r ./.ipynb_checkpoints
fi

git add --all
git commit -m "$1"
git push -u origin "$2"



