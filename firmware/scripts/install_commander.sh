#!/bin/bash

commander_path=$HOME/SimplicityCommander

echo -e "\n*** Downloading Simplicity Commander ***\n"
mkdir $commander_path
cd $commander_path
wget https://www.silabs.com/documents/public/software/SimplicityCommander-Linux.zip

echo -e "\n*** Installing Simplicity Commander ***\n"
unzip *.zip
mv $commander_path/SimplicityCommander-Linux/* $commander_path
rmdir $commander_path/SimplicityCommander-Linux
tar -xvf *.tar.bz

echo -e "\n*** Removing .tar.bz and .zip ***\n"
rm *.tar.bz
rm *.zip*
