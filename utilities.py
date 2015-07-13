#!/usr/bin/env python
#ROLLING HOUGH TRANSFORM FIBER EXTRACTION
#Susan Clark, Lowell Schudel

#-----------------------------------------------------------------------------------------
#Imports
#-----------------------------------------------------------------------------------------
from __future__ import division #Must be first line of code in the file
import sys
import time
import string
import copy
import itertools
import operator
import collections
import networkx as nx
import datetime
import scipy.stats
import numpy as np

class ProgressBar:
    # Create progress meter that looks like: 
    # message + ' ' + '[' + '#'*p + ' '*(length-p) + ']' + time_message
    stream = sys.stdout
    TEXTWIDTH = 79
    length = int(0.55 * TEXTWIDTH)

    def append(self, child):
        def child_in_children(children):
            if len(children) == 0:
                return False
            elif child in children:
                return True
            else:
                return any(map(child_in_children, map(operator.attrgetter('_children'), children)))

        if not isinstance(child, ProgressBar):
            raise ValueError('You can only append some ProgressBar to existing ProgressBar!')
        elif child._parent is not None:
            raise ValueError('A ProgressBar object can only have one parent!')
        elif self == child:
            raise ValueError('ProgressBar cannot be appended to itself!')
        elif child_in_children(self._children):
            raise ValueError('ProgressBar cannot be appended to its own children!')
        else:
            self._children.append(child)
            child._parent = self

            self._stop_time = None #TODO Flags for recalculation
            self._problem_size += child._problem_size + 1
            self._dp = 1.0/self._problem_size

            if self._incrementing == child._incrementing:
                self._step += child._step
            else:
                self._step += child._problem_size - child._step
            if self._incrementing:
                self._prog = self._step*self._dp
            else:
                self._prog = 1.0 - self._step*self._dp
            self._last_update = time.time()            

    def update(self, progress=None):
        #return progress #This line will short-circuit the rest of the function and preserve the desired return value
        #if self._parent is not None:
        #self._parent.update()
        if len(self._children) > 0:
            return self._children[0].update(progress)
        #return progress

        if progress is None:
            if self._incrementing:
                self._step += 1
            else:
                self._step -= 1
            self._prog += self._dp
        elif isinstance(progress, int) and (self._problem_size >= progress >= 0):
            self._step = progress
            self._prog = self._step*self._dp
        elif isinstance(progress, float) and (1.0 >= progress >= 0.0):
            self._prog = progress
            self._step = int(self._problem_size*self._prog)
        else:
            raise ValueError('Invalid input to ProgressBar.update(): '+str(progress))
            #return progress #TODO

        now = time.time()
        if self._start_time is None:
            self._start_time = now
            return progress 

        elif (self._stop_time is None) or (self._prog == 0.0) or (self._prog == 1.0) or ((now - self._last_update) > 1.0):        
            self._stop_time = self._start_time + (now - self._start_time)/self._prog
        
            sec_remaining = int(self._stop_time - now)
            if sec_remaining >= 3600:
                time_message = ' < ' + str(sec_remaining//3600  +1) + 'hrs'
            elif sec_remaining >= 60:
                time_message = ' < ' + str(sec_remaining//60  +1) + 'min'
            else:
                time_message = ' < ' + str(sec_remaining +1) + 'sec'
            messlen = ProgressBar.TEXTWIDTH-(ProgressBar.length+3)-len(time_message)
            p = int(ProgressBar.length*self._prog)
            if 0 <= p < ProgressBar.length:
                ProgressBar.stream.write('\r{2} [{0}{1}]{3}'.format('#'*p, ' '*(ProgressBar.length-p), string.ljust(self._message, messlen)[:messlen], time_message))
                ProgressBar.stream.flush()
            else:
                final_offset = ProgressBar.TEXTWIDTH-len(time_message)
                final_message = string.ljust('Finished:', final_offset)[:final_offset]
                ProgressBar.stream.write('\r{0}{1}'.format(final_message, time_message))
                ProgressBar.stream.flush()
                print ''
        return progress
                
    def __init__(self, problem_size=0, message='ProgressBar:', incrementing=True):
        assert isinstance(message, str)
        assert isinstance(incrementing, bool)
        assert isinstance(problem_size, int)
        assert problem_size >= 0
        
        self._stop_time = None
        self._message = message
        self._incrementing = incrementing
        self._problem_size = problem_size
        self._prog = 0.0
        if self._incrementing:
            self._step = 0
        else:
            self._step = self._problem_size
        self._children = list()
        self._parent = None

        if problem_size == 0:
            self._start_time = time.time()
            self._dp = 0.0
        else:
            self._start_time = None #time.time()
            self._dp = 1.0/self._problem_size
        self._last_update = time.time()
