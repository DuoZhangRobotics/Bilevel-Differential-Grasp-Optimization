#! /bin/bash
rm -rdf *.png
rm -r ./limits
rm -r ./__pycache__
rm -r ./.ipynb_checkpoints
rm ./.project
rm ./.pydevproject

git add --all
git commit -m "$1"
git push -u origin "$2"



