#!/usr/bin/env python3
"""Summarize conserved-energy drift in this bundle's traj.properties output.

Columns: step, time/fs, conserved/eV, potential/eV, kinetic_md/eV.
Reports diagnostics, not a universal acceptance threshold.
"""
import argparse
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('properties')
    p.add_argument('--waters',type=int,required=True)
    p.add_argument('--discard-fs',type=float,default=0)
    a=p.parse_args()
    if a.waters<1 or a.discard_fs<0:p.error('Invalid water count or discard time')
    data=np.loadtxt(a.properties,ndmin=2)
    data=data[data[:,1]>=a.discard_fs]
    if len(data)<3 or data.shape[1]!=5 or not np.isfinite(data).all():
        p.error('Need at least three finite rows in the generated five-column format')
    t=(data[:,1]-data[0,1])/1000
    if not np.all(np.diff(t)>0):p.error('Times must strictly increase; analyze separate segments separately')
    e=data[:,2]-data[0,2]
    slope,intercept=np.polyfit(t,e,1)
    print(f'Duration analyzed: {t[-1]:.6g} ps; frames: {len(t)}')
    print(f'Linear drift: {slope*1000/a.waters:+.6g} meV/water/ps')
    print(f'Energy peak-to-peak: {np.ptp(e)*1000/a.waters:.6g} meV/water')
    print(f'Energy residual RMS after linear fit: {np.std(e-(slope*t+intercept))*1000/a.waters:.6g} meV/water')
    print('Compare matched NVE segments at dt and dt/2; a short trace cannot certify production drift.')


if __name__=='__main__':main()
