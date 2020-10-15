param($1, $2)
rm -r ./.idea
git add --all
git commit -m $1
git push -u origin $2