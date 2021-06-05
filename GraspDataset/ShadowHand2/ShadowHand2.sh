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
[ -z $initial ] && initial="/home/jiangtang/IRC/build3/beforeOptimize_Shadowhand_21_0.2/initialParameters.txt"
echo "iteration = "$iteration
echo "path = "$path
echo "initial parameter path = "$initial

$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand2_200_0.200000.dat ShadowHand2 0.2 1 $iteration $path/ $initial 
$path/mainGraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand2_200_0.200000.dat ShadowHand2 0.2 0 $iteration $path/ $initial 
$path/mainraspPlan ../.././data/ShadowHand/shadowhand_noarm_noknuckle.urdf 200 ShadowHand2_200_0.200000.dat ShadowHand2 0.2 3 $iteration $path/ $initial 
