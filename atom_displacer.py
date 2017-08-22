import datetime
import os
import json
import random as rnd

import numpy as np

import copy

from definitions import *

# parameters:
# info: json with a resume of all the simulation
# delta: displeacement that we want (in atomic unit)
# points: step to generate in each direction
# directions: for now only ('x', 'y', 'z' are supported)
# outdir: directory where to generate all the new configuration
info = json.load(open('./v0.2.1.json'))
delta = np.float32(0.01 * 1.889725989)
points = [-1, 0, 1]
outdir = './displeacement_out'

# key of the simulation we are using as reference
# or None (random choice by algorithm)
key = '40a8c754f14f9990828f265f45b1349c5c119fc90169cb00c6a9d471'

# useful simulations parameters
max_seconds = 3600
pseudo_dir = '~/QE/pseudo'

# debug flag
debug = False

# hardcoded for now
atom_description = {
    'C': ['1', '4.00', '12.', 'C.rpbe-n-kjpaw_psl.0.1.UPF'],
    'H': ['2', '1.00', '1.', 'H.rpbe-kjpaw_psl.0.1.UPF'],
    'N': ['4', '5.00', '14.', 'N.rpbe-n-kjpaw_psl.0.1.UPF'],
    'O': ['3', '6.00', '16.', 'O.rpbe-n-kjpaw_psl.0.1.UPF']}

# ------- end of parameters -------

if not key:
    key = rnd.choice(list(info.keys()))

print('simulation info:')
print(key)
for n, t, p, f in info[key]['atom']:
    print(t + str(p))

T = np.array(info[key]['cell_side'], dtype=np.float32)
T1 = np.linalg.inv(T)

dx = np.dot(np.array([delta, 0, 0], dtype=np.float32), T1)
dy = np.dot(np.array([0, delta, 0], dtype=np.float32), T1)
dz = np.dot(np.array([0, 0, delta], dtype=np.float32), T1)
# vector norm
print('T matrix')
print(T)
print('T-1 marix')
print(T1)
if debug:
    print('T-1 * T')
    print(T1 * T)

if not os.path.exists(outdir):
    os.makedirs(outdir)

if not os.path.exists(os.path.join(outdir, 'json')):
    os.makedirs(os.path.join(outdir, 'json'))

counter = 0
zero_flag = True
for n, t, p0, f in info[key]['atom']:
    # n = number
    # t = type
    # p0 = position
    # f = forces

    subdir = str(n) + '_' + t
    if not os.path.exists(os.path.join(outdir, subdir)):
        os.makedirs(os.path.join(outdir, subdir))

    p0 = np.array(p0, dtype=np.float32, copy=True)

    L = []
    for dp in points:
        if dp > 0:
            lx = ['{}x'.format(str(dp)), p0 + dp * dx]
            ly = ['{}y'.format(str(dp)), p0 + dp * dy]
            lz = ['{}z'.format(str(dp)), p0 + dp * dz]
            if debug:
                print('--')
                print(np.linalg.norm(np.dot(p0 - lx[1], T)) / 1.889725989)
                print(dp)
            L.extend([lx, ly, lz])
        elif dp < 0:
            lx = ['-{}x'.format(str(-dp)), p0 + dp * dx]
            ly = ['-{}y'.format(str(-dp)), p0 + dp * dy]
            lz = ['-{}z'.format(str(-dp)), p0 + dp * dz]
            L.extend([lx, ly, lz])
        elif dp == 0:
            if zero_flag:
                L.append(['p0', p0])
                zero_flag = False

    for subsubdir, pos in L:
        opath = os.path.join(outdir, subdir, subsubdir)

        if not os.path.exists(opath):
            os.makedirs(opath)

        target = copy.deepcopy(info[key])
        perdolapazienza = copy.deepcopy(info[key]['atom'])
        perdolapazienza[int(n) - 1][2] = ['{:10.9f}'.format(x) for
                                          x in pos.tolist()]
        target['atom'] = perdolapazienza
        target['first'] = key
        target['last'] = key
        target['next'] = key
        target['previous'] = key

        with open(os.path.join(outdir, 'json', '%s_%s_%s.simulation' %
                               (str(n), t, subsubdir)), 'w') as w:
            json.dump(target, w)

        # buid up simulation
        sim_in = {}
        sim_job = {}

        aspec = '\n'.join(
                ['{} {} {}'.format(k, str(v[2]), str(v[3]))
                 for k, v in atom_description.items()])

        apos = '\n'.join([x[1] + '   ' + ' '.join(x[2])
                          for x in target['atom']])

        sim_in = OUT_FILE.format(
            calculation='scf',
            prefix='force_test',
            outdir='./tmp',
            pseudo_dir=pseudo_dir,
            crystal_structure=0,
            e_cutoff=80,
            e_cutroh=560,
            ATOMIC_SPECIES=aspec,
            ATOMIC_POSITIONS=apos,
            nat=len(target['atom']),
            K_points='   4   4    3  0 0 0',
            cell='bohr',
            x=' '.join(info[key]['cell_side'][0]),
            y=' '.join(info[key]['cell_side'][1]),
            z=' '.join(info[key]['cell_side'][2]),
            max_seconds=max_seconds)

        sim_job = JOB_FILE.format(
            processName=subdir + '_' + subsubdir,
            wallTime=str(datetime.timedelta(seconds=max_seconds)),
            fFile='qe.in',
            out='qe.out',)

        with open(os.path.join(opath, 'qe.in'), 'w') as w:
            w.write(sim_in)

        with open(os.path.join(opath, 'qe.job'), 'w') as w:
            w.write(sim_job)
        if counter % 100 == 0:
            finalScript = open('FinalScript_{}.sh'.format(int(counter / 100) +
                                                          1), 'w')

        finalScript.write('cd {}/{}/{}\n'.format(outdir, subdir, subsubdir))
        finalScript.write('pwd >> ../../../process.sub\n')
        finalScript.write('qsub {} >> ../../../process.sub\n'.format('qe.job'))
        finalScript.write('cd ../../../\n')
        counter += 1

finalScript.close()
