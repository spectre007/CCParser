#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 14:43:23 2018

@author: alex
"""



def isfloat(string):
    """Check if a string can be converted to float."""
    try:
        float(string)
        return True
    except ValueError:
        return False