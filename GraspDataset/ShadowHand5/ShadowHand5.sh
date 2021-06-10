Help()
{
   # Display Help
   echo "Grasp Planning for our method, Q1 metric and Object closeness."
   echo
   echo "Syntax: scriptTemplate [-r|i|h|p]"
   echo "Options:"
   echo "-r, --rounds			Set max iterations"
   echo "-h, --help				Print this Help."
   echo "-i, --initial-path		Path to the initial paramters of the hand"
   echo "-p, --path				Path to the build folder"
   echo
}

while :
do
   case $1 in
      	-h | --help) # display Help
         	Help
         	exit;;
		-r | --rounds) iteration="$2";
			#echo $iteration;
			shift 2;
			;;
		-i | --initial) initial="$2";
			shift 2;
			;;
		-p | --path) path="$2";
			shift 2;
			;;
		--) # End of all options
          	shift
          	break
          	;;
      	-*)
		     echo "Error: Unknown option: $1" >&2
		     ## or call function display_help
		     exit 1 
		     ;;
      	*)  # No more options
          	break
          	;;
   esac
done

[ -z $iteration ] && iteration=1000
[ -z $path ] && path="../../.././build3"
[ -z $initial ] && initial="./initialParameters.txt"
echo "iteration = "$iteration
echo "path = "$path
echo "initial parameter path = "$initial


echo "***************************************Generating object***************************************"
$path/mainPointCloudObject ShadowHand5.obj 200 0.2
echo "***************************************Generating gripper***************************************"
$path/mainGripper ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat
echo "***************************************Running our method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 1 $iteration ./ $initial
echo "***************************************Running  Q1 method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 0 $iteration ./ $initial
echo "***************************************Running Object closeness method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 3 $iteration ./ $initial
echo "***************************************Testing our method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 2 1 ./ ./afterOptimize_Q_INF_CONSTRAINT_FGT_BarrettHand_ShadowHand5_0.2/parameters.txt
echo "***************************************Testing Q1 method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 2 1 ./ ./afterOptimize_Q_1_BarrettHand_ShadowHand5_0.2/parameters.txt
echo "***************************************Testing OC method***************************************"
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand5_200_0.200000.dat ShadowHand5 0.2 2 1 ./ ./afterOptimize_No_Metric_OC_BarrettHand_ShadowHand5_0.2/parameters.txt


