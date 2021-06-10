[ -z $iteration ] && iteration=4000
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

if [ -f "FGT$1.dat" ]; then
    echo "FGT$1.dat exists"
else 
    echo "Running our method FGT"
    $path/mainGraspPlan ../.././data/BarrettHand/bh280.urdf $1 BarrettHand1_$1_0.300000.dat BarrettHand1 0.3 1 $iteration ./ $initial 1e-9 >> FGT$1.dat
    mv afterOptimize_Q_INF_CONSTRAINT_FGT_BarrettHand_BarrettHand1_0.3 afterOptimize_Q_INF_CONSTRAINT_FGT_BarrettHand_BarrettHand1_0.3_$1
fi

echo "Computing QInf"
$path/mainGraspPlan ../.././data/BarrettHand/bh280.urdf $1 BarrettHand1_$1_0.300000.dat BarrettHand1 0.3 2 1 ./ afterOptimize_Q_INF_CONSTRAINT_FGT_BarrettHand_BarrettHand1_0.3_$1/parameters.txt >> QInf$1.dat
