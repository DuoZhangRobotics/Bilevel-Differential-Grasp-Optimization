# Help()
# {
#    # Display Help
#    echo "Generate scaled .off file, .xml file, parameter file for GraspIt."
#    echo
#    echo "Syntax: scriptTemplate [-o|h]"
#    echo "Options:"
#    echo "-o --outputPath        Output Path for GraspIt"   
#    echo
# }

# while :
# do
#    case $1 in
#       	-h | --help) # display Help
#          	Help
#          	exit;;
#         -o | --outputPath) outputPath="$2";
#             shift 2;
#             ;;
# 		--) # End of all options
#           	shift
#           	break
#           	;;
#       	-*)
# 		     echo "Error: Unknown option: $1" >&2
# 		     ## or call function display_help
# 		     exit 1 
# 		     ;;
#       	*)  # No more options
#           	break
#           	;;
#    esac
# done

python3 ../../graspitGeneration.py --inputPath $PWD/maxRange_Scale.txt
