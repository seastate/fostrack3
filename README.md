# fostrack

This is a fork of a lightweight utility, fostrack, by @occeltic for assembling trajectories from files with coordinates extracted from videos of moving objects. It is updated to work in Python3, and to output a few additional movement statistics.

A simple particle-tracking program for use with processed .fos-part particle files

.fos-part particle files are comma-delimited text files generated from the fosica video processing software available from Wallingford Imaging Systems. They contain information on the size, shape, and location of particles detected in each frame of a video.

fostrack is a Python tracking library which can easily be used through a command line interface. With it you can perform simple path stitching in a conservative way, visualize the resulting paths as animations, save animations as .mp4s, and save your tracking results for quick access at a later time (.fos-trk extension).

To begin using it, navigate to the directory where you've stored fostrack.py and pull up the help menu by typing:
python fostrack.py -h

An example .fos-part file is included in your download so you can try out some basic commands. 

There is also an example .mp4 which was produced from this .fos-part file using the command:

python fostrack.py process example.fos-part -r 10 -s 3 -l 10 --fps-in 5.3 --fps-out 15 --width 1296 --height 972 -m ./ -t ./
