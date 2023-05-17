# fostrack3

A simple particle-tracking program for use with processed .fos-part particle files

This is a fork of a lightweight utility, fostrack, by @occeltic for assembling trajectories from files with coordinates extracted from videos of moving objects. It is updated to work in Python3, and to output additional movement statistics.

.fos-part particle files are comma-delimited text files generated from the fosica video processing software available from Wallingford Imaging Systems. They contain information on the size, shape, and location of particles detected in each frame of a video.

fostrack is a Python tracking library which can easily be used through a command line interface. With it you can perform simple path stitching in a conservative way, visualize the resulting paths as animations, save animations as .mp4s, and save your tracking results for quick access at a later time (.fos-trk extension).

To begin using it, navigate to the directory where you've stored fostrack.py and pull up the help menu by typing:
python fostrack.py -h

An example .fos-part file is included in your download so you can try out some basic commands. 

There is also an example .mp4 which was produced from this .fos-part file using the command:

python3 fostrack.py process example.fos-part -r 10 -s 3 -l 10 --fps-in 5.3 --fps-out 15 --width 1296 --height 972 -m ./ -t ./ [--interactive False] [--stats 10] [ | cat > BLAH.txt]

where optional arguments are:

### --interactive False

Do not stop and wait for the user to confirm parameters before every analysis. The default behavior is to wait for confirmation.

### --stats 10
Calculate path statistics. The 10 refers to the frame interval for statistics. That is, with --stats 10, velocity will be measured by taking positions at the 0th, 10th, 20th etc. frames of a path, and dividing the displacement by 10 for each of those intervals.

There is a tradeoff: longer intervals are less subject to noise, but have fewer samples (and may miss short timescale variations).

"directions" are the direction vectors associated with those velocities, speeds are the magnitudes of velocity vectors. "dir_corrs" are the dot products between successive direction observations. That is, the correlation between directions in successive observations within a path.

It also outputs the average number of particles and paths, and the number of frames at which particles and paths were sampled to compute those averages.


### | cat > BLAH.txt

Output statistics results to the csv-formatted file BLAH.txt. The statistics method currently has a very basic and pretty clunky output structure. It prints the output to the screen, and leaves you to cut and past the output or write a parser (e.g.in R). This flag facilitates sequestering this output into a file.

The "|" means the output of the first command will be used as input for the second command, which in this case writes output to BLAH.txt.

### Batch processing

The code supports batch processing by submitting a set of names in the command line. For example, to process all the fos-part files in the directory MyDir, use

python3 fostrack.py $(ls MyDir/.fos-part) -r 10 -s 3 -l 10 --fps-in 5.3 --fps-out 15 --width 1296 --height 972 -m ./ -t ./ [--interactive False] [--stats 10] [ | cat > BLAH.txt]

## Output path file
The code produces a csv-format file with the suffix "fos-trk" that summarizes the paths found by the tracking analysis. This file can be reloaded without rerunning the analysis.

In a fos-trk file, the first section recapitulates the frame-by-frame list of particles metrics. The second section (headed by the keywords "all paths = ") contains a sequence of particle numbers, with numbers in each line representing successive observations of a particle within a given path.
