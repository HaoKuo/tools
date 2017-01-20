#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''python timeit.py
'''
from __future__ import division
__author__ = 'Hao Guo'

import time
import functools
 
def timeit(func):
    @functools.wraps(func)
    def wrapper():
        start = time.clock()
        func()
        end =time.clock()
        print('Time used: {}s'.format(end - start))
    return wrapper

@timeit
def foo():
    print 'in foo()'

def main():
    foo()

if __name__ == '__main__':

    main()


