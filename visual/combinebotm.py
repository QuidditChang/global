#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Symmetric counterpart of combinesurf.py for CitcomS bottom (CMB) data.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

'''
Cut-paste and combine Citcom botm data

usage: combinebotm.py modelname timestep nodex nodey nodez n_surf_proc nprocx nprocy nprocz
'''

class CombineBotm(object):


    def __init__(self, grid):
        self.saved = [''] * (grid['nox'] * grid['noy'])
        return


    def readData(self, crdfilename, botmfilename, grid, cap):
        fp1 = file(crdfilename, 'r')
        fp2 = file(botmfilename, 'r')

        nprocx = int(cap['nprocx'])
        nprocy = int(cap['nprocy'])
        nprocz = int(cap['nprocz'])

        nox = grid['nox']
        noy = grid['noy']
        noz = grid['noz']

        mynox = 1 + (nox-1)/nprocx
        mynoy = 1 + (noy-1)/nprocy
        mynoz = 1 + (noz-1)/nprocz
        snodes = mynox * mynoy

        # skip header
        fp1.readline()
        fp2.readline()

        botm = fp2.readlines()
        if not (len(botm) == snodes):
            print '"%s" file size incorrect' % botmfilename

        crd = fp1.readlines()
        if not (len(crd) == snodes*mynoz):
            print '"%s" file size incorrect' % crdfilename

        # paste data — botm uses first z-node (index i*mynoz), surf uses last
        data = range(snodes)
        for i in range(len(data)):
            x, y, z = crd[i*mynoz].split()
            data[i] = '%s %s %s' % (x, y, botm[i])

        return data


    def join(self, data, me, grid, cap):
        nprocx = cap['nprocx']
        nprocy = cap['nprocy']
        nprocz = cap['nprocz']

        mylocz = me % nprocz
        mylocx = ((me - mylocz) / nprocz) % nprocx
        mylocy = (((me - mylocz) / nprocz - mylocx) / nprocx) % nprocy

        nox = int(grid['nox'])
        noy = int(grid['noy'])
        noz = int(grid['noz'])

        mynox = 1 + (nox-1)/nprocx
        mynoy = 1 + (noy-1)/nprocy

        if not len(data) == mynox * mynoy:
            raise ValueError, "data size"

        mynxs = (mynox - 1) * mylocx
        mynys = (mynoy - 1) * mylocy

        n = 0
        for i in range(mynys, mynys+mynoy):
            for j in range(mynxs, mynxs + mynox):
                m = j + i * nox
                self.saved[m] = data[n]
                n += 1

        return


    def write(self, filename, grid, data):
        fp = file(filename, 'w')
        fp.write('%d x %d\n' % (grid['nox'], grid['noy']))
        fp.writelines(data)
        return



if __name__ == '__main__':

    import sys

    if not len(sys.argv) == 10:
        print __doc__
        sys.exit(1)

    prefix = sys.argv[1]
    step = int(sys.argv[2])

    grid = {}
    grid['nox'] = int(sys.argv[3])
    grid['noy'] = int(sys.argv[4])
    grid['noz'] = int(sys.argv[5])

    nprocxy = int(sys.argv[6])
    cap = {}
    cap['nprocx'] = int(sys.argv[7])
    cap['nprocy'] = int(sys.argv[8])
    cap['nprocz'] = int(sys.argv[9])

    nproc_per_cap = cap['nprocx'] * cap['nprocy'] * cap['nprocz']
    # botm: only bottom-layer procs (rank % nprocz == 0)
    for i in range(nprocxy):
        cb = CombineBotm(grid)
        for n in range(i * nproc_per_cap, (i+1) * nproc_per_cap, cap['nprocz']):
            crdfilename  = '%s.coord.%d' % (prefix, n)
            botmfilename = '%s.botm.%d.%d' % (prefix, n, step)
            print 'reading', botmfilename
            data = cb.readData(crdfilename, botmfilename, grid, cap)
            cb.join(data, n, grid, cap)

        filename = '%s.botm%02d.%d' % (prefix, i, step)
        print 'writing', filename
        cb.write(filename, grid, cb.saved)
