rm -r ./__pycache__
rm -r ./limits
git add --all
git commit -m "$1"
git push -u origin master
