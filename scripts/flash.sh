#!/bin/bash

############################################################
# Help                                                     #
############################################################
function help
{
    # Display Help
    echo "This script flashes boards from Silabs RTL kits."
    echo
    echo "Usage:"
    echo "flash.sh [-h] [-f FIRMWARE]"
    echo
    echo "Options:"
    echo "-h    Print this help."
    echo "-f    Firmware path."
    echo
}

############################################################
# Variables                                                #
############################################################

commander=$HOME/SimplicityCommander/commander/commander
root_dir=$PWD

############################################################
# Process the input options.                               #
############################################################

while getopts "hf:" option; do
    case $option in
        h) # Display help
            help
            exit;;
        f) # Firmware path
            firmware="$(cd "$(dirname "$OPTARG")"; pwd)/$(basename "$OPTARG")"
            cd $root_dir;;
        \?) # Invalid option
            echo "Error: Invalid option"
            echo "Use -h for help"
            exit 1;;
    esac
done

############################################################
# Main program                                             #
############################################################

if [[ ! $firmware ]]; then
    echo -e "Missing firmware argument (-f)."
    echo "Use -h for help"
    exit 1
fi

$commander flash $firmware
if [[ $? -ne 0 ]]; then
    echo -e "Failed to flash firmware"
    exit 1
fi
