for name in BarrettHandd*
do
    newname=BarrettHand"$(echo "$name" | cut -c13-)"
    mv "$name" "$newname"
done
