"""
fostrack.py is a simple particle tracking program for use with .fos-part video files

Owen Coyle
ocoyle@uw.edu

July 13th, 2015
"""

import matplotlib
matplotlib.use('TKAgg')
from matplotlib.patches import Ellipse
from math import pi, sqrt
from pylab import array
from matplotlib import pyplot as plt
from matplotlib import animation, gridspec
import argparse
import time
import os
from numpy import interp

class Pos2D(object):
    """
    A class for representing 2-D coordinates
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return "(%f, %f)" % (self.x, self.y)

    def dist(self, pos):
        """
        Return the Euclidean distance between self and pos
        """
        return ((self.x - pos.x)**2 + (self.y-pos.y)**2)**0.5

    def equals(self, pos):
        """Check is two pos are equivalent"""
        if (self.x == pos.x) and (self.y == pos.y):
            return True
        else:
            return False


class Particle(object):
    """
    A class for representing particles from a fos-vid file
    """
    def __init__(self, frameNo, time, cameraNo, x, y, area, partNo, boundWidth, boundHeight, minDim, maxDim,
                 minDimAngle, maxDimAngle, max2min, length, lengthAngle, width, widthAngle, length2width):
        self.frameNo = frameNo
        self.time = time
        self.cameraNo = cameraNo
        self.pos = Pos2D(x, y)
        self.area = area
        self.partNo = partNo
        self.boundWidth = boundWidth
        self.boundHeight = boundHeight
        self.minDim = minDim
        self.maxDim = maxDim
        self.minDimAngle = minDimAngle
        self.maxDimAngle = maxDimAngle
        self.max2min = max2min
        self.length = length
        self.lengthAngle = lengthAngle
        self.width = width
        self.widthAngle = widthAngle
        self.length2width = length2width

    def plot(self, ax, color):
        """
        Return an ellipse object for plotting, add to ax, with color
        """
        e = Ellipse(xy = array([self.pos.x, self.pos.y]), width = self.width, height = self.length, angle = self.lengthAngle/pi*180+90)
        #e = Ellipse(xy = array([self.pos.x, self.pos.y]), width = self.width, height = self.length, angle = self.lengthAngle/pi*180)
        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_facecolor(color)
        return e

    def repr(self):
        """
        Return a string representation of the particle
        """
        return "%i,%f,%i,%f,%f,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f" % (self.frameNo, self.time, self.cameraNo, self.pos.x, self.pos.y, self.area, self.partNo, self.boundWidth, self.boundHeight, self.minDim, self.maxDim, self.minDimAngle, self.maxDimAngle, self.max2min, self.length, self.lengthAngle, self.width, self.widthAngle, self.length2width)

    def equals(self, particle):
        """Check if this particle and particle are equivalent"""
        if self.frameNo != particle.frameNo:
            return False
        if self.time != particle.time:
            return False
        if self.cameraNo != particle.cameraNo:
            return False
        if not self.pos.equals(particle.pos):
            return False
        if self.area != particle.area:
            return False
        if self.partNo != particle.partNo:
            return False
        if self.boundWidth != particle.boundWidth:
            return False
        if self.boundHeight != particle.boundHeight:
            return False
        if self.minDim != particle.minDim:
            return False
        if self.maxDim != particle.maxDim:
            return False
        if self.minDimAngle != particle.minDimAngle:
            return False
        if self.maxDimAngle != particle.maxDimAngle:
            return False
        if self.max2min != particle.max2min:
            return False
        if self.length != particle.length:
            return False
        if self.lengthAngle != particle.lengthAngle:
            return False
        if self.width != particle.width:
            return False
        if self.widthAngle != particle.widthAngle:
            return False
        if self.length2width != particle.length2width:
            return False
        return True

class Path(object):
    """
    A class representing a series of particles
    """
    def __init__(self):
        self.length = 0
        self.particleList = []
        self.frameNoList = []
    
    def addParticle(self, particle):
        """
        add a particle to this path
        """
        i = 0
        while i < self.length:
            if self.particleList[i].frameNo > particle.frameNo:
                break
            i += 1
        self.particleList.insert(i, particle)
        self.frameNoList.insert(i, particle.frameNo)
        self.length += 1

    def removeParticle(self, particle):
        """
        remove particle
        """
        try:
            i = self.particleList.index(particle)
            self.particleList.pop(i)
            self.frameNoList.pop(i)
            self.length -= 1
        except ValueError:
            print("Something went wrong: tried to remove a particle from a path that didn't contain it")


    def latestParticle(self):
        """
        return the most recently added particle in this list
        """
        return self.particleList[self.length - 1]

    def searchFrame(self, frameNo):
        """
        search if the particle list contains a particle with a matching frameNo. Return the first match
        or None if there are no matches.
        """
        try:
            i = self.frameNoList.index(frameNo)
            return self.particleList[i]
        except ValueError:
            return None


    def searchTime(self, time):
        """
        search if the particle list contains a partice with a matching time. Return the first match or
        None if there are no matches.
        """
        for particle in self.particleList:
            if particle.time == time:
                return particle
        return None

    def spanFrame(self, frameNo):
        """
        search if the particle list spans frameNo (i.e. it either contains frameNo directly, or has
        frames before and after. Returns True if the Path spans the frameNo, False if not.
        """
        seenBefore = False
        for particle in self.particleList:
            if particle.frameNo < frameNo:
                seenBefore = True
            elif particle.frameNo == frameNo:
                return True
            elif seenBefore:
                return True
        return False

    def spanTime(self, time):
        """
        search if the particle list spans time. Returns True if it does, False if not
        """
        seenBefore = False
        for particle in self.particleList:
            if particle.time < time:
                seenBefore = True
            elif particle.time == time:
                return True
            elif seenBefore:
                return True
        return False
            
    def checkRadiusAt(self, pos, radius, frameNo):
        """
        check if there is a particle within radius of pos for frame frameNo in the particleList.
        Return True if there is, False if not
        """
        framePart = self.searchFrame(frameNo)
        if framePart == None:
            return False
        elif framePart.pos.dist(pos) > radius:
            return False
        else:
            return True

    def plot(self, ax, color, frameStart, frameEnd):
        """plot a representation of this path on axis ax in color, only include frames from frameStart to frameEnd inclusive"""
        xs = []
        ys = []
        for p in self.particleList:
            if frameStart <= p.frameNo <= frameEnd:
                xs.append(p.pos.x)
                ys.append(p.pos.y)
        return ax.plot(xs, ys, color = color)[0]

    def equals(self, path):
        """Check if this path is equivalent to path"""
        if self.length != path.length:
            return False
        if len(self.particleList) != len(path.particleList):
            return False
        for particle1, particle2 in zip(self.particleList, path.particleList):
            if not particle1.equals(particle2):
                return False
        return True

    def avgParticleLength(self):
        """Return the average particle length in this path"""
        if self.length == 0:
            return None
        else:
            sumLength = 0.0
            for particle in self.particleList:
                sumLength += particle.length
            return sumLength/self.length

    def interpPos(self,interval):
        """return the i,j pixel positions of this path at spacing interval
           begining at the first particle in the path. If the path contains
           the required frame, its observed position is used. If the required
           frame is spanned by the path but a missing frame, linear interpolation 
           is used."""
        # Ingest particles in path, accounting for skipped frames
        present_frames=[p.frameNo for p in self.particleList]
        present_xs=[p.pos.x for p in self.particleList]
        present_ys=[p.pos.y for p in self.particleList]
        frameSamples=list(range(present_frames[0],present_frames[-1],interval))
        xSamples=interp(frameSamples,present_frames,present_xs)
        ySamples=interp(frameSamples,present_frames,present_ys)
        return frameSamples,xSamples,ySamples

class Frame(object):
    """
    A class for holding all of the particles from a video frame
    """
    def __init__(self, frameNo):
        self.frameNo = frameNo
        self.particleList = []
        self.pathList = []
    
    def addParticle(self, particle):
        """
        Add a particle to the list of particles in this frame
        """
        self.particleList.append(particle)
        self.pathList.append(None)

    def assignParticle(self, particle, path):
        """
        Assign the given particle (already in particleList) to path (doesn't add it to path, 
        just changes the reference in the pathList)
        """
        for i, part in enumerate(self.particleList):
            if part is particle:
                self.pathList[i] = path
                break

    def unassignedParticles(self):
        """
        Return a list of currently unassigned particles
        """
        unassignedList = []
        for particle, path in zip(self.particleList, self.pathList):
            if path == None:
                unassignedList.append(particle)
        return unassignedList

    def equals(self, frame):
        """Check if self and frame are equivalent"""
        if self.frameNo != frame.frameNo:
            return False
        if len(self.particleList) != len(frame.particleList):
            return False
        if len(self.pathList) != len(frame.pathList):
            return False
        for particle1, particle2 in zip(self.particleList, frame.particleList):
            if not particle1.equals(particle2):
                return False
        for path1, path2 in zip(self.pathList, frame.pathList):
            if path1 == None:
                if path2 != None:
                    return False
            elif path2 == None:
                return False
            elif not path1.equals(path2):
                return False
        return True



class Tracker(object):
    """
    A class for assembling particles into paths
    """
    def __init__(self):
        self.videofile = ""
        self.frameList = []
        self.frameNos = []
        self.pathList = []
    
    def equals(self, tracker):
        """Check if self and tracker are equivalent"""
        if self.videofile != tracker.videofile:
            return False
        if len(self.frameList) != len(tracker.frameList):
            return False
        if len(self.frameNos) != len(tracker.frameNos):
            return False
        if len(self.pathList) != len(tracker.pathList):
            return False
        for frame1, frame2 in zip(self.frameList, tracker.frameList):
            if not frame1.equals(frame2):
                return False
        for frameNo1, frameNo2 in zip(self.frameNos, tracker.frameNos):
            if frameNo1 != frameNo2:
                return False
        for path1, path2 in zip(self.pathList, tracker.pathList):
            if path1 == None:
                if path2 != None:
                    return False
            elif path2 == None:
                return False
            elif not path1.equals(path2):
                return False
        return True


    def loadVideo(self, videofile):
        """
        import particles from the video file
        """
        try:
            f = open(videofile, "r")
        except IOError as e:
            print("I/O error(%s): %s %s" % (e.errno, e.strerror, videofile))
            return

        self.videofile = videofile
        
        #skip the first line
        f.readline()
        #now read in each line and make a particle
        currentFrame = None
        for line in f:
            print(line)
            words = line.split()
            print(words)
            #parse the different arguments
            frameNo = int(words[0])
            time = float(words[1])
            cameraNo = int(words[2])
            x = float(words[3])
            y = float(words[4])
            area = float(words[5])
            #area = int(words[5])
            particleNo = int(words[6])
            boundWidth = int(words[7])
            boundHeight = int(words[8])
            minDim = float(words[9])
            maxDim = float(words[10])
            minDimAngle = float(words[11])
            maxDimAngle = float(words[12])
            max2min = float(words[13])
            length = float(words[14])
            lengthAngle = float(words[15])
            width = float(words[16])
            widthAngle = float(words[17])
            length2width = float(words[18])
            #skip gap lines
            if particleNo == 0:
                continue
            #now store them appropriately
            Part = Particle(frameNo, time, cameraNo, x, y, area, particleNo, boundWidth, boundHeight, minDim, maxDim, minDimAngle, maxDimAngle, max2min, length, lengthAngle, width, widthAngle, length2width)
            if currentFrame == None:
                currentFrame = Frame(frameNo)
                self.frameList.append(currentFrame)
                self.frameNos.append(frameNo)
            elif currentFrame.frameNo < frameNo:
                currentFrame = Frame(frameNo)
                self.frameList.append(currentFrame)
                self.frameNos.append(frameNo)
            currentFrame.addParticle(Part)
        f.close()

    def simpleStitchForward(self, radius, minPathLength = 3, skipFrames = 1):
        """
        A simple path stitching routine. Stitches particles together by looking for particles within
        a search radius of the previous position. Only stitches particles together if there are no conflicting
        matches within the search radius. skipframes defines the number of frames, in absolute terms, that may
        be skipped while still making a match, i.e. skipframes = 0 would only allow matches between successive frames.
        minPathLength defines the minimum number of links a path must have to be considered a viable path and retained.
        In the event of a skip frames conflict, the most recent match will be used. This algorithm starts at the first
        frame, and proceeds to the end in one pass.
        """
        for frame in self.frameList:
            frameNo = frame.frameNo
            frame = self.frameList[self.frameNos.index(frameNo)]
            #make a list to store over assigned paths
            conflictPaths = []
            for particle in frame.particleList:
                myPos = particle.pos
                #look for potential matches in the current paths
                potentialMatchFrames = []
                frameOffsets = range(1, 2+skipFrames)
                for frameOffset in frameOffsets:
                    searchFrameNo = frameNo - frameOffset
                    potentialMatches = []
                    for path in self.pathList:
                        if path.checkRadiusAt(myPos, radius, searchFrameNo):
                            potentialMatches.append(path)
                    potentialMatchFrames.append(potentialMatches)
                #evaluate potential matches for one that is acceptable
                foundMatch = False
                for potentialMatches in potentialMatchFrames:
                    if len(potentialMatches) == 1:
                        path = potentialMatches[0]
                        #check if this path already has a particle for this frameNo, if so, slate for removal
                        if path.searchFrame(frameNo):
                            #only add this path if someone else hasn't already
                            if path not in conflictPaths:
                                conflictPaths.append(path)
                        #if this path doesn't currently have a conflict, then assign me to it
                        else:
                            path.addParticle(particle)
                            frame.assignParticle(particle, path)
                            foundMatch = True
                            break
                #if you didn't find a match then make a new path
                if foundMatch is False:
                    path = Path()
                    path.addParticle(particle)
                    frame.assignParticle(particle, path)
                    self.pathList.append(path)
            #now deal with the conflict paths
            for badPath in conflictPaths:
                particle = badPath.searchFrame(frameNo)
                if particle == None:
                    print("Something went wrong, couldn't find a particle that should be in a path")
                    break
                newPath = Path()
                newPath.addParticle(particle)
                badPath.removeParticle(particle)
                self.pathList.append(newPath)
                frame.assignParticle(particle, newPath)
            #for performance, clean out old paths that aren't long enough
            frameNoCutoff = frameNo - skipFrames - 1
            for path in reversed(self.pathList):
                #is the path old enough?
                if path.latestParticle().frameNo < frameNoCutoff:
                    #is the path too short
                    if path.length < minPathLength:
                        #reset all particles in the path to be unassigned
                        for particle in path.particleList:
                            frame = self.frameList[self.frameNos.index(particle.frameNo)]
                            frame.assignParticle(particle, None)
                        #delete the path
                        self.pathList.remove(path)
        #remove all paths remaining less than the minimum path length
        for path in reversed(self.pathList):
            if path.length < minPathLength:
                for particle in path.particleList:
                    frame = self.frameList[self.frameNos.index(particle.frameNo)]
                    frame.assignParticle(particle, None)
                self.pathList.remove(path)
                    
        #run a check for conflicts, print(a warning if so
        if self.checkForConflicts() is True:
            print("Uh oh, after stitching some particles and paths seem to conflict")

    def checkForConflicts(self):
        """
        Check that all particles in paths are assigned to the proper path in my frame library.
        Return False if there are no conflicts, True if there are
        """
        for path in self.pathList:
            for particle in path.particleList:
                frameNo = particle.frameNo
                frame = self.frameList[self.frameNos.index(frameNo)]
                assignedPath = frame.pathList[frame.particleList.index(particle)]
                if path is not assignedPath:
                    return True
        return False

    def getFramePaths(self, frameNo):
        """
        Return a list of all paths spanning frameNo
        """
        paths = []
        for path in self.pathList:
            if path.spanFrame(frameNo):
                paths.append(path)
        return paths

    def countParts(self, frameNos, minParticleLength = 0.0, maxParticleLength = float("inf")):
        """return a list of the number of particles in each frameNo in frameNos with particle length between min and max"""
        numParts = []
        for frameNo in frameNos:
            if frameNo in self.frameNos:
                frame = self.frameList[self.frameNos.index(frameNo)]
                numParts.append(len(frame.particleList))
            else:
                numParts.append(0)
        return numParts                

    def countPaths(self, frameNos, minParticleLength = 0.0, maxParticleLength = float("inf"),
                   minPathLength = 1, maxPathLength = float("inf")):
        """return a list of the number of paths spanning each frameNo in frameNos with particle length 
           between min and max, and path length between min and max"""
        numPaths = []
        for frameNo in frameNos:
            paths = self.getFramePaths(frameNo)
            #check if each path matches
            numMatches = 0
            for path in paths:
                if minPathLength <= path.length <= maxPathLength:
                    if minParticleLength <= path.avgParticleLength() <= maxParticleLength:
                        numMatches += 1
            numPaths.append(numMatches)
        return numPaths

    def summaryStats(self, frameNoInterval, minParticleLength = 0.0, maxParticleLength = float("inf"),
                     minPathLength = 1, maxPathLength = float("inf")):
        """calculate summary statistics across the current set of paths, sampling at frameNoInterval,
        with particle length between min and max, and path length between min and max"""
        # Get average number of particles and paths across all frames
        frame_sample_list=list(range(self.frameNos[0],self.frameNos[-1],frameNoInterval))
        no_samples=len(frame_sample_list)
        avg_no_parts=array(self.countParts(frame_sample_list)).mean()
        avg_no_paths=array(self.countPaths(frame_sample_list)).mean()
        # Get movement stats, per path, weighted by path length
        vels=[]    # list of observed velocities
        speeds=[]  # list of observed speeds
        dirs=[]    # list of observed direction vectors
        dir_corrs=[]    # list of observed directonal persistences
        for path in self.pathList:
            # test if path meets criteria
            if path.length<minPathLength or path.length> maxPathLength or \
               path.avgParticleLength()<minParticleLength or path.avgParticleLength()>maxParticleLength:
                continue # test fails, skip path
            fs,xs,ys=path.interpPos(frameNoInterval)
            drs=[]
            for i in range(len(fs)-1):
                # Caclulate velocity and speed, add to lists
                vel=[(xs[i+1]-xs[i])/frameNoInterval,(ys[i+1]-ys[i])/frameNoInterval]
                vels.append(vel)
                spd=sqrt(vel[0]**2+vel[1]**2)
                #spd=sqrt(vel[0]**2+vel[1]**2)/frameNoInterval
                speeds.append(spd)
                # Calculate directions. If stationary, direction is None.
                # We compare directions only within a path, so use drs and add to corr_dirs afterwards
                if spd>0.:
                    dr=[vel[0]/spd,vel[1]/spd]
                else:
                    dr=None
                drs.append(dr)
            for i in range(len(drs)-1): 
                if drs[i]!=None and drs[i+1]!=None:
                    dir_corr=drs[i][0]*drs[i+1][0]+drs[i][1]*drs[i+1][1]
                    dir_corrs.append(dir_corr)
            dirs.extend(drs)
        print('Path statistics for file: ',self.filename)
        print('Number samples: ',no_samples)
        print('Avg. number particles: ',avg_no_parts)
        print('Avg. number paths: ',avg_no_paths)
        print('vels: ',vels)
        print('speeds: ',speeds)
        print('directions: ',dirs)
        print('direction correlations: ',dir_corrs)
        print('\n\n\n')

    def skippedFrames(self, startFrame = None, endFrame = None):
        """return a list of the skipped frames from startFrame and endFrame"""
        if startFrame == None:
            startFrame = self.frameNos[0]
        if endFrame == None:
            endFrame = self.frameNos[-1]
        missingFrameNos = []
        for frameNo in range(startFrame, endFrame + 1):
            if frameNo not in self.frameNos:
                missingFrameNos.append(frameNo)
        return missingFrameNos
        
    def save(self, filename):
        """
        save a representation of this tracker object to filename
        """
        f = open(filename, "w")
	    #save the filename
        f.write("videofile = %s" % self.videofile)
        #list all of the particles in the video
        f.write("\nall particles =")
        particles = []
        for frame in self.frameList:
            for particle in frame.particleList:
                f.write("\n%s" % particle.repr())
                particles.append(particle)
        #now write all of the paths, indicating which particle by it's index in the particle list
        f.write("\nall paths =")
        paths = []
        for path in self.pathList:
            paths.append(path)
            first = True
            for particle in path.particleList:
                i = particles.index(particle)
                if first:
                    f.write("\n%i" % i)
                    first = False
                else:
                    f.write(",%i" % i)
        f.close()

    def loadFromFile(self, filename):
        """
        load a tracker object from a previously saved tracker object
        """
        self.filename=filename
        f = open(filename, "r")
        self.videofile = f.readline().strip("videofile = ").strip()
        #skip the next line
        f.readline()
        particles = []
        currentFrame = None
        for line in f:
            if "all paths =" in line:
                break
            #it's a particle
            words = line.split(",")
            #parse the different arguments
            frameNo = int(words[0])
            time = float(words[1])
            cameraNo = int(words[2])
            x = float(words[3])
            y = float(words[4])
            area = int(words[5])
            particleNo = int(words[6])
            boundWidth = int(words[7])
            boundHeight = int(words[8])
            minDim = float(words[9])
            maxDim = float(words[10])
            minDimAngle = float(words[11])
            maxDimAngle = float(words[12])
            max2min = float(words[13])
            length = float(words[14])
            lengthAngle = float(words[15])
            width = float(words[16])
            widthAngle = float(words[17])
            length2width = float(words[18])
            particle = Particle(frameNo, time, cameraNo, x, y, area, particleNo, boundWidth, boundHeight, minDim, maxDim, minDimAngle, maxDimAngle, max2min, length, lengthAngle, width, widthAngle, length2width)
            particles.append(particle)
            #store in a frame
            if currentFrame == None:
                frame = Frame(particle.frameNo)
                frame.addParticle(particle)
                self.frameList.append(frame)
                self.frameNos.append(frameNo)
                currentFrame = frame
            elif currentFrame.frameNo < particle.frameNo:
                frame = Frame(particle.frameNo)
                frame.addParticle(particle)
                self.frameList.append(frame)
                self.frameNos.append(frameNo)
                currentFrame = frame
            else:
                currentFrame.addParticle(particle)
        #now read in the paths
        for line in f:
            path = Path()
            self.pathList.append(path)
            words = line.split(",")
            for word in words:
                partNo = int(word)
                particle = particles[partNo]
                frameNo = particle.frameNo
                frame = self.frameList[self.frameNos.index(frameNo)]
                path.addParticle(particle)
                frame.assignParticle(particle, path)
        f.close()
        # If frame insterval was passed as an argument, calculate statistics
        if self.stats>0:
            self.summaryStats(self.stats)

class Animator(object):
    """An object for animating tracker objects"""
    def __init__(self):
        """initiate animation object"""
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.ellipses = ()
        self.ellipsesDelete = ()
        self.tracks = ()
        self.tracksDelete = ()
        self.frameNoText = None
        self.timelineTracker = None
        self.frameNos = []
        self.skippedFrameNos = []
        self.frames = []
        self.pathsByFrame = []
        self.percTracked = []
        self.title = ''
        self.speedFactor = None

    def init(self):
        """intiate the animation"""
        #delete old ellipses
        self.ellipsesDelete = self.ellipses
        self.ellipses = ()
        for eDel in self.ellipsesDelete:
            self.ax1.artists.remove(eDel)
        #delete the old tracks
        self.tracksDelete = self.tracks
        self.tracks = ()
        for tDel in self.tracksDelete:
            self.ax1.lines.remove(tDel)
        #set static items
        self.ax1.title.set_text(self.title)
        self.ax2.set_ylabel("Tracks(#)", color = "red")
        self.ax2.set_xlabel("Frame No.")
        self.ax3.set_ylabel("% Tracked", color = "blue")
        speedText = self.ax1.text(0.72, 0.05, "Speed x %0.1f" % self.speedFactor, transform = self.ax1.transAxes)

        numPaths = []
        for paths in self.pathsByFrame:
            numPaths.append(len(paths))
        timeline = self.ax2.plot(self.frameNos, numPaths, color = "red")[0]
        timeline2 = self.ax3.plot(self.frameNos, self.percTracked, color = "blue")[0]
        skippedFrameMarkers = ()
        for frameNo in self.skippedFrameNos:
            skippedFrameMarkers += (self.ax2.text(frameNo, self.ax2.get_ylim()[0], "*"),)
        #set dynamic items
        self.frameNoText.set_text('')
        self.timelineTracker.set_data([], [])
        #return a tuple of all items that have been updated
        return self.ellipses + self.ellipsesDelete + (self.frameNoText,) + self.tracks + self.tracksDelete + (timeline,) + (timeline2,) + (self.timelineTracker,)

    def animate(self, i):
        """perform an animation step"""
        #delete old ellipses
        self.ellipsesDelete = self.ellipses
        self.ellipses = ()
        for eDel in self.ellipsesDelete:
            #self.ax1.artists.remove(eDel)
            eDel.remove()
        #delete the old tracks
        self.tracksDelete = self.tracks
        self.tracks = ()
        for tDel in self.tracksDelete:
            #self.ax1.lines.remove(tDel)
            tDel.remove()
        #draw the new particles and paths
        frame = self.frames[i]
        frameNo = self.frameNos[i]
        for particle in frame.particleList:
            self.ellipses += (particle.plot(self.ax1, "blue"),)
        for path in self.pathsByFrame[i]:
            self.tracks += (path.plot(self.ax1, "red", self.frameNos[0], frameNo),)
        #update the frame number text
        self.frameNoText.set_text("Frame No. %i" % frameNo)
        #update the timeline tracker
        numPaths = len(self.pathsByFrame[i])
        self.timelineTracker.set_data(frameNo, numPaths)
        #return a tuple of all items that have been updated
        return self.ellipses + self.ellipsesDelete + self.tracks + self.tracksDelete + (self.frameNoText,) + (self.timelineTracker,)

    def setupAnimation(self, tracker, title, fpsIn, fpsOut, startFrame, endFrame, width, height, xlim, ylim, invertX, invertY):
        """helper function, sets up an animation"""
        self.fig = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios = [3, 1])
        self.ax1 = plt.subplot(gs[0], aspect = "equal")
        self.ax2 = plt.subplot(gs[1])
        self.ax3 = self.ax2.twinx()
        self.ellipses = ()
        self.tracks = ()
        if startFrame == None:
            startFrame = tracker.frameNos[0]
        if endFrame == None:
            endFrame = tracker.frameNos[-1]
        if xlim == None:
            if invertX == 0:
                xlim = (0, width)
            else:
                xlim = (width, 0)
        if ylim == None:
            if invertY == 0:
                ylim = (height, 0)
            else:
                ylim = (0, height)
        self.ax1.set_xlim(xlim)
        self.ax1.set_ylim(ylim)
        self.ax2.set_xlim(startFrame, endFrame)
        self.ax3.set_ylim(50, 100)
        self.frameNoText = self.ax1.text(0.72, 0.95, '', transform = self.ax1.transAxes)
        self.timelineTracker = self.ax2.plot([], [], 'bo')[0]
        self.title = title
        self.speedFactor = round((fpsOut/fpsIn)*2/2)
        self.frameNos = []
        self.frames = []
        for i, frameNo in enumerate(tracker.frameNos):
            if startFrame <= frameNo <= endFrame:
                self.frameNos.append(frameNo)
                self.frames.append(tracker.frameList[i])
        self.skippedFrameNos = tracker.skippedFrames()
        self.pathsByFrame = []
        for frameNo in self.frameNos:
            self.pathsByFrame.append(tracker.getFramePaths(frameNo))
        self.percTracked = []
        numPaths = tracker.countPaths(tracker.frameNos)
        numParts = tracker.countParts(tracker.frameNos)
        for numPath, numPart in zip(numPaths, numParts):
            self.percTracked.append(100.0*numPath/numPart)
        #setup the animation
        ani = animation.FuncAnimation(self.fig, self.animate, frames = len(self.frameNos), interval = 1000.0/fpsOut, blit = True, init_func = self.init, repeat = False)
        return ani

    def playVideo(self, tracker, title, fpsIn, fpsOut, startFrame = None, endFrame = None, width = 2592, height = 1944, xlim = None, ylim = None, invertX = 0, invertY = 0):
        ani = self.setupAnimation(tracker, title, fpsIn, fpsOut, startFrame, endFrame, width, height, xlim, ylim, invertX, invertY)
        plt.show()
        plt.close(self.fig)

    def saveVideo(self, tracker, title, fpsIn, fpsOut, filename, startFrame = None, endFrame = None, width = 2592, height = 1944, xlim = None, ylim = None, bitrate = 1024, invertX = 0, invertY = 0):
        ani = self.setupAnimation(tracker, title, fpsIn, fpsOut, startFrame, endFrame, width, height, xlim, ylim, invertX, invertY)
        ani.save(filename, fps = fpsOut, bitrate = bitrate)
        plt.close(self.fig)

def inputYN(str):
    """Ask the user for a yes or no answer by printing str.
        Keep asking until you get a y or n. Return a boolean True for y, False for n
    """
    yes = ["y", "yes", "Y", "Yes", "YES"]
    no = ["n", "no", "N", "No", "NO"]
    userInput = input(str)
    #userInput = raw_input(str)
    if userInput in yes:
        keepAsking = False
        toReturn = True
    elif userInput in no:
        keepAsking = False
        toReturn = False
    else:
        keepAsking = True
    while keepAsking:
        print("Please enter y or n")
        userInput = input(str)
        #userInput = raw_input(str)
        if userInput in yes:
            keepAsking = False
            toReturn = True
        elif userInput in no:
            keepAsking = False
            toReturn = False
        else:
            keepAsking = True
    return toReturn

def processFosPart(parsedArgs, parser):
    """Helper method for the main method. Handles the command for processing .fos-part files"""
    #check if the input is formatted properly
    if parsedArgs.radius == None:
        print("Please supply a search radius for particle stitching using the '-r' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.skip_frames == None:
        print("Please supply a # of permissable skipped frames for particle stitching using the '-s' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.min_track_length == None:
        print("Please supply a minimum track length for particle stitching using the '-l' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.fps_in == None:
        print("Please supply the speed at which input video was shot (FPS) using the '--fps-in' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.fps_out == None:
        print("Please supply the speed at which output video should be saved and displayed using the '--fps-out' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    radius = parsedArgs.radius
    skipFrames = parsedArgs.skip_frames
    minTrackLength = parsedArgs.min_track_length
    filelist = parsedArgs.filelist
    movieDir = parsedArgs.movie_dir
    trackDir = parsedArgs.track_dir
    fpsIn = parsedArgs.fps_in
    fpsOut = parsedArgs.fps_out
    interactive = parsedArgs.interactive
    stats = parsedArgs.stats
    bitRate = parsedArgs.bit_rate
    width = parsedArgs.width
    height = parsedArgs.height
    startFrame = parsedArgs.start_frame
    endFrame = parsedArgs.end_frame
    invertX = parsedArgs.invert_x
    invertY = parsedArgs.invert_y

    #print(a helpful message
    print("------------------------------------------------------------------")
    print("Welcome to fostrack.py")
    print("Please review my settings before proceeding...")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("--------------------------")
    print("Files to be analyzed:")
    for path in filelist:
        file = os.path.split(path)[1]
        print(file)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("Analysis Parameters:")
    print("radius = %0.1f" % radius)
    print("skip frames = %i" % skipFrames)
    print("min. track length = %i" % minTrackLength)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("I/O parameters:")
    print("dir. to save .fos-trk files: %s" % trackDir)
    print("dir. to save .mp4 movie files: %s" % movieDir)
    print("video dimensions: %i x %i pixels" % (width, height))
    if invertX == 0:
        print("invert x axis: False")
    else:
        print("invert x axis True")
    if invertY == 0:
        print("invert y axis: False")
    else:
        print("invert y axis: True")
    print("input video speed: %0.1f FPS" % fpsIn)
    print("output video speed: %0.1f FPS" % fpsOut)
    print("output video bit rate: %i kilobits/s" % bitRate)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    if trackDir == None:
        print("### WARNING ###")
        print("You have not supplied a directory for me to save .fos-trk files to. This means I won't save results and you'll have to re-stitch paths later")
    if movieDir == None:
        print("### WARNING ###")
        print("You have not supplied a directory for me to save .mp4 files to. This means I won't save copies of particle movies")
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    proceed = True
    if interactive == True:
        proceed = inputYN("Would you like to proceed using these settings? (y/n):")
    if proceed == False:
        print("Ok, exiting...")
        return
    display = False
    if interactive == True:
        display = inputYN("Would you like to display preview videos as we go? (y/n, faster if n):")
    for path in filelist:
        file = os.path.split(path)[1]
        prefix = file.strip(".fos-part")
        tracker = Tracker()
        print("Loading %s..." % file)
        tracker.loadVideo(path)
        print("Stitching %s..." % file)
        tracker.simpleStitchForward(radius, minPathLength = minTrackLength, skipFrames = skipFrames)
        animator = Animator()
        if display == True:
            animator.playVideo(tracker, prefix, fpsIn, fpsOut, width = width, height = height, startFrame = startFrame, endFrame = endFrame, invertX = invertX, invertY = invertY)
        if trackDir != None:
            savePath = "%s%s.fos-trk" % (trackDir, prefix)
            print("Saving tracks to %s..." % savePath)
            tracker.save(savePath)
        if movieDir != None:
            savePath = "%s%s.mp4" % (movieDir, prefix)
            print("Saving movie to %s..." % savePath)
            animator.saveVideo(tracker, prefix, fpsIn, fpsOut, savePath, width = width, height = height, startFrame = startFrame, endFrame = endFrame, bitrate = bitRate, invertX = invertX, invertY = invertY)
        # If frame interval was passed as an argument, calculate statistics
        tracker.filename=file
        if stats>0:
            tracker.summaryStats(stats)
        
    print("Analysis complete, exiting")
    print("------------------------------------------------------------------")

def playFosTrk(parsedArgs, parser):
    """Helper method, handles the case where we just want to play previously analyzed .fos-trk files without any changes"""
    if parsedArgs.fps_in == None:
        print("Please supply the speed at which input video was shot (FPS) using the '--fps-in' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.fps_out == None:
        print("Please supply the speed at which output video should be saved and displayed using the '--fps-out' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    filelist = parsedArgs.filelist
    movieDir = parsedArgs.movie_dir
    fpsIn = parsedArgs.fps_in
    fpsOut = parsedArgs.fps_out
    interactive = parsedArgs.interactive
    stats = parsedArgs.stats
    bitRate = parsedArgs.bit_rate
    width = parsedArgs.width
    height = parsedArgs.height
    startFrame = parsedArgs.start_frame
    endFrame = parsedArgs.end_frame
    invertX = parsedArgs.invert_x
    invertY = parsedArgs.invert_y

    #print a helpful message
    print("------------------------------------------------------------------")
    print("Welcome to fostrack.py")
    print("Please review my settings before proceeding...")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("--------------------------")
    print("Files to be displayed:")
    for path in filelist:
        file = os.path.split(path)[1]
        print(file)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("I/O parameters:")
    print("dir. to save .mp4 movie files: %s" % movieDir)
    print("video dimensions: %i x %i pixels" % (width, height))
    if invertX == 0:
        print("invert x axis: False")
    else:
        print("invert x axis True")
    if invertY == 0:
        print("invert y axis: False")
    else:
        print("invert y axis: True")
    print("input video speed: %0.1f FPS" % fpsIn)
    print("output video speed: %0.1f FPS" % fpsOut)
    print("output video bit rate: %i kilobits/s" % bitRate)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    if movieDir == None:
        print("### WARNING ###")
        print("You have not supplied a directory for me to save .mp4 files to. This means I won't save copies of particle movies")
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    proceed = True
    if interactive == True:
        proceed = inputYN("Would you like to proceed using these settings? (y/n):")
    if proceed == False:
        print("Ok, exiting...")
        return
    display = False
    if interactive == True:
        display = inputYN("Would you like to display preview videos as we go? (y/n, faster if n):")
    for path in filelist:
        file = os.path.split(path)[1]
        prefix = file.strip(".fos-trk")
        tracker = Tracker()
        print("Loading %s..." % file)
        tracker.loadFromFile(path)
        animator = Animator()
        if display == True:
            animator.playVideo(tracker, prefix, fpsIn, fpsOut, width = width, height = height, startFrame = startFrame, endFrame = endFrame, invertX = invertX, invertY = invertY)
        if movieDir != None:
            savePath = "%s%s.mp4" % (movieDir, prefix)
            print("Saving movie to %s..." % savePath)
            animator.saveVideo(tracker, prefix, fpsIn, fpsOut, savePath, width = width, height = height, startFrame = startFrame, endFrame = endFrame, bitrate = bitRate, invertX = invertX, invertY = invertY)
    # If frame interval was passed as an argument, calculate statistics
    if stats>0:
          tracker.filename=path
          tracker.summaryStats(tracker.stats)
    print("Analysis complete, exiting")
    print("------------------------------------------------------------------")

def processFosTrk(parsedArgs, parser):
    """Process a previously saved .fos-trk file"""
    """Helper method for the main method. Handles the command for processing .fos-part files"""
    #check if the input is formatted properly
    if parsedArgs.radius == None:
        print("Please supply a search radius for particle stitching using the '-r' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.skip_frames == None:
        print("Please supply a # of permissable skipped frames for particle stitching using the '-s' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.min_track_length == None:
        print("Please supply a minimum track length for particle stitching using the '-l' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.fps_in == None:
        print("Please supply the speed at which input video was shot (FPS) using the '--fps-in' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    if parsedArgs.fps_out == None:
        print("Please supply the speed at which output video should be saved and displayed using the '--fps-out' flag. Type 'python %s -h' for more information" % parser.prog)
        return
    radius = parsedArgs.radius
    skipFrames = parsedArgs.skip_frames
    minTrackLength = parsedArgs.min_track_length
    filelist = parsedArgs.filelist
    movieDir = parsedArgs.movie_dir
    trackDir = parsedArgs.track_dir
    fpsIn = parsedArgs.fps_in
    fpsOut = parsedArgs.fps_out
    interactive = parsedArgs.interactive
    self.stats = parsedArgs.stats
    bitRate = parsedArgs.bit_rate
    width = parsedArgs.width
    height = parsedArgs.height
    startFrame = parsedArgs.start_frame
    endFrame = parsedArgs.end_frame
    invertX = parsedArgs.invert_x
    invertY = parsedArgs.invert_y

    #print(a helpful message
    print("------------------------------------------------------------------")
    print("Welcome to fostrack.py")
    print("Please review my settings before proceeding...")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("--------------------------")
    print("Files to be re-analyzed:")
    for path in filelist:
        file = os.path.split(path)[1]
        print(file)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("Analysis Parameters:")
    print("radius = %0.1f" % radius)
    print("skip frames = %i" % skipFrames)
    print("min. track length = %i" % minTrackLength)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    print("I/O parameters:")
    print("dir. to save .fos-trk files: %s" % trackDir)
    print("dir. to save .mp4 movie files: %s" % movieDir)
    print("video dimensions: %i x %i pixels" % (width, height))
    if invertX == 0:
        print("invert x axis: False")
    else:
        print("invert x axis True")
    if invertY == 0:
        print("invert y axis: False")
    else:
        print("invert y axis: True")
    print("input video speed: %0.1f FPS" % fpsIn)
    print("output video speed: %0.1f FPS" % fpsOut)
    print("output video bit rate: %i kilobits/s" % bitRate)
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    if trackDir == None:
        print("### WARNING ###")
        print("You have not supplied a directory for me to save .fos-trk files to. This means I won't save results and you'll have to re-stitch paths later")
    if movieDir == None:
        print("### WARNING ###")
        print("You have not supplied a directory for me to save .mp4 files to. This means I won't save copies of particle movies")
    print("--------------------------")
    if interactive == True:
        time.sleep(1.5) #introduce a small delay so the user can read text
    proceed = True
    if interactive == True:
        proceed = inputYN("Would you like to proceed using these settings? (y/n):")
    if proceed == False:
        print("Ok, exiting...")
        return
    display = False
    if interactive == True:
        display = inputYN("Would you like to display preview videos as we go? (y/n, faster if n):")
    for path in filelist:
        file = os.path.split(path)[1]
        prefix = file.strip(".fos-trk")
        tracker = Tracker()
        print("Loading %s..." % file)
        tracker.loadFromFile(path)
        print("Stitching %s..." % file)
        tracker.simpleStitchForward(radius, minPathLength = minTrackLength, skipFrames = skipFrames)
        animator = Animator()
        if display == True:
            animator.playVideo(tracker, prefix, fpsIn, fpsOut, width = width, height = height, startFrame = startFrame, endFrame = endFrame, invertX = invertX, invertY = invertY)
        if trackDir != None:
            savePath = "%s%s.fos-trk" % (trackDir, prefix)
            print("Saving tracks to %s..." % savePath)
            tracker.save(savePath)
        if movieDir != None:
            savePath = "%s%s.mp4" % (movieDir, prefix)
            print("Saving movie to %s..." % savePath)
            animator.saveVideo(tracker, prefix, fpsIn, fpsOut, savePath, width = width, height = height, startFrame = startFrame, endFrame = endFrame, bitrate = bitRate, invertX = invertX, invertY = invertY)
        # If frame interval was passed as an argument, calculate statistics
        if self.stats>0:
              self.filename=path
              self.summaryStats(self.stats)
    print("Analysis complete, exiting")
    print("------------------------------------------------------------------")


def main():
    """The main method to run fostrack on a file or batch of files"""
    #setup the argparser
    parser = argparse.ArgumentParser(description = "Process .fos-part files into tracks")
    parser.add_argument("command", choices = ["load", "process", "play"], help = "load a .fos-trk file for further tracking, process a .fos-part file, or play a previously saved .fos-trk file") 
    parser.add_argument("filelist", help="A .fos-part file, or list of .fos-part files to be analyzed", nargs = "+")
    parser.add_argument("-t", "--track-dir", "--track-directory", help = "The directory to write .fos-trk files", nargs = "?", default = None)
    parser.add_argument("-m", "--movie-dir", "--movie-directory", help = "The directory to write .mp4 movies", nargs = "?", default = None)
    parser.add_argument("-r", "--radius", "--search-radius", help = "The search radius (in pixels) used for simple path stitching", type = float)
    parser.add_argument("-l", "--min-track-length", help = "The minimum track length for a track to be considered valid", type = int)
    parser.add_argument("-s", "--skip-frames", help = "The max. number of frames a track may span w/o an observation", type = int)
    parser.add_argument("--fps-in", help = "The fps of input video", type = float)
    parser.add_argument("--fps-out", help = "The fps of output movies", type = float)
    parser.add_argument("--interactive", help = "Use the program in interactive mode", nargs = "?", default = True)
    #parser.add_argument("--stats", help = "Output path statistics", nargs = "?", default = True)
    parser.add_argument("--stats", help = "Output path statistics", nargs = "?", default = 0, type = int)
    parser.add_argument("--bit-rate", help = "Bit rate for output .mp4 files", nargs = "?", default = 1024, type = int)
    parser.add_argument("--width", help = "Video width (pixels)", nargs = "?", default = 2592, type = int)
    parser.add_argument("--height", help = "Video height (pixels)", nargs = "?", default = 1944, type = int)
    parser.add_argument("--start-frame", help = "Video start frame", nargs = "?", default = None, type = int)
    parser.add_argument("--end-frame", help = "Video end frame", nargs = "?", default = None, type = int)
    parser.add_argument("--invert-x", help = "Invert x axis? 0 = False, 1 = True", nargs = "?", default = 0)
    parser.add_argument("--invert-y", help = "Invert y axis? 0 = False, 1 = True", nargs = "?", default = 0, type = int)
    parsedArgs = parser.parse_args()

    if parsedArgs.command == "play":
        return playFosTrk(parsedArgs, parser)
    elif parsedArgs.command == "process":
        return processFosPart(parsedArgs, parser)
    elif parsedArgs.command == "load":
        return processFosTrk(parsedArgs, parser)



if __name__ == "__main__":
    main()
    """
    filepath = "./Tracks1/Halocline_09_29_30s_0015param.fos-trk"
    tracker = Tracker()
    tracker.loadFromFile(filepath)
    print(tracker.skippedFrames())
    """
        
