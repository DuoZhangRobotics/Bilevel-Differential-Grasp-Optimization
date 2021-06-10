[ -z $iteration ] && iteration=100
[ -z $path ] && path="../../.././build3"
[ -z $initial ] && initial="./initialParameters.txt"

echo "density = "$1
echo "iteration = "$iteration
echo "path = "$path
echo "initial parameter path = "$initial

if [ -f "BarrettHand1_$1_0.300000.dat" ]; then
    echo "BarrettHand1_$1_0.300000.dat exists"
else 
    echo "Generating object"
    $path/mainPointCloudObject BarrettHand1.obj $1 0.3
    echo "Generating gripper"
    $path/mainGripper ../.././data/BarrettHand/bh280.urdf $1 BarrettHand1_$1_0.300000.dat
fi

if [ -f "noFGT$1.dat" ]; then
    echo "noFGT$1.dat exists"
else 
    echo "Running our method noFGT"
    $path/mainGraspPlan ../.././data/BarrettHand/bh280.urdf $1 BarrettHand1_$1_0.300000.dat BarrettHand1 0.3 0 $iteration profile $initial >> noFGT$1.dat
    echo "Running our method FGT"
    $path/mainGraspPlan ../.././data/BarrettHand/bh280.urdf $1 BarrettHand1_$1_0.300000.dat BarrettHand1 0.3 1 $iteration profile $initial >> FGT$1.dat
fi
