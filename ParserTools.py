#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 14:43:23 2018

@author: alex
"""

import math


def is_float(string):
    """Check if a string can be converted to float."""
    try:
        float(string)
        return True
    except ValueError:
        return False

def is_square(integer):
    root = math.sqrt(integer)
    if int(root + 0.5) ** 2 == integer: 
        return True
    else:
        return False
