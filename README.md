# Fluorescence Intermittency Trajectory Analyzer

This software package performs standard analysis techniques for intensity trajectories which exhibit the fluorescence intermittency phenomenon. It performs threshold, intensity distribution, and power spectrum analyses. This is a quick analysis method which may need further polishing for presentation.

![alt tag](https://raw.githubusercontent.com/aruth2/FITA/master/FITA.png)

Prerequisites:

    GTK+3
    gnuplot
    FFTW3


Installation Instructions:

Linux:
compile with this command

    gcc FITA.c -o FITA `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread -lfftw3 

To Run:

    ./FITA


Mac:
Use this script to install the prerequisites:
    
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install gtk+3
    brew install gnuplot
    BASEDIR="$( dirname "$0" )"
    cd "$BASEDIR"

Then compile

    gcc FITA.c -o FITA `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread -lfftw3

And Run
    
    ./FITA

Windows:
There are some issues with using GTK+ in Windows, but I have compiled GTK+ code in Windows before.
If you would like to run this program on a Windows machine please contact me.
