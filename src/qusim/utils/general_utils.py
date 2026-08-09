# -*- coding: utf-8 -*-
"""
@author: Zihao Wang
Original source: https://github.com/Zihao96/qtrlb/blob/main/qtrlb/utils/general_utils.py
"""
import time

def make_it_list(thing, default: list = None):
    """
    A crucial, life-saving function.
    """
    if isinstance(thing, list):
        return thing
    elif thing is None:
        return [] if default is None else default
    else:
        return [thing]
    

def timeit(func):
    """
    A decorator to time a function.
    """
    def wrapper(*args, **kwargs):
        t_i = time.time()
        result = func(*args, **kwargs)
        print(f'{func.__name__}: {round(time.time()-t_i)}s')
        return result
    return wrapper