#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
'''
   __init__.py 
'''
import os, sys, time, re, glob, itertools
import subprocess, multiprocessing, shlex
import gzip
from Bio.Seq import Seq
from Bio import SeqIO
import pysam
from timeFunc import timeit
from argparse import ArgumentParser
from collections import defaultdict


__all__ = ["os","sys","time","re","glob","itertools","subprocess","multiprocessing","shlex","gzip","Seq","SeqIO","pysam","ArgumentParser","defaultdict","timeFunc","findVar"]
