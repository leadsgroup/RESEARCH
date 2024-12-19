#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 15:12:01 2024

@author: nn18
"""

def RK4(x0, dx, y0, SlopeFn):
    m1 = SlopeFn(x0, y0)
    m2 = SlopeFn(x0+(dx/2), y0+(m1*dx/2))
    m3 = SlopeFn(x0+(dx/2), y0+(m2*dx/2))
    m4 = SlopeFn(x0+dx, y0+(m3*dx))
    
    change = (dx/6)*(m1 + 2*m2 + 2*m3 + m4)
    return y0 + change