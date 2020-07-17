# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 09:44:35 2020

@author: guinetj
"""
import os

def test_path():
    """
    return the data folder
    :return: data path.
    :rtype: str
    """
    data_folder = os.path.join(os.environ["TESTPATH"])
    return data_folder
