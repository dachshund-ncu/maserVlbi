#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# iso: 2020-10-18T19:30:00
#      2020-10-19T03:30:00
import sys
sys.path.append('data')
from maserVlbi import maserVlbi
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    eee = maserVlbi(sys.argv[1])
    colors = eee.spots.getJetColors()
    print(f"Project code: {eee.project_code}")
    plt.scatter(eee.spots.dRA, eee.spots.dDEC, s=np.log(eee.spots.flux * 1000)**2.0 * 2, c = colors, marker='o', edgecolor='black')
    plt.gca().invert_xaxis()
    plt.show()