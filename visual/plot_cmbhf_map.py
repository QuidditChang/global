#!/usr/bin/env python
'''
Plot CitcomS CMB heat flux (bhflux) global map using GMT4.

usage: plot_cmbhf_map.py timestep age storedir

Reads global.botm{00..11}.timestep from the current directory.
Each file: header line "nox x noy", then lines of "colat_rad lon_rad bhflux".
Outputs: cmbhf_{age}Ma.ps and cmbhf_{age}Ma.pdf
'''

import sys
import os
import math

def read_botm_caps(timestep, ncap=12):
    '''Read all cap botm files, return list of (lon_deg, lat_deg, bhflux) tuples.'''
    points = []
    for i in range(ncap):
        fname = 'global.botm%02d.%d' % (i, timestep)
        if not os.path.exists(fname):
            print 'Warning: %s not found, skipping' % fname
            continue
        fp = open(fname, 'r')
        fp.readline()  # skip header "nox x noy"
        for line in fp:
            parts = line.split()
            if len(parts) < 3:
                continue
            colat_rad = float(parts[0])
            lon_rad   = float(parts[1])
            bhflux    = float(parts[2])
            lat_deg = 90.0 - math.degrees(colat_rad)
            lon_deg = math.degrees(lon_rad)
            if lon_deg < 0.0:
                lon_deg += 360.0
            points.append((lon_deg, lat_deg, bhflux))
        fp.close()
    return points


def write_xyz(points, filename):
    fp = open(filename, 'w')
    for lon, lat, val in points:
        fp.write('%f %f %f\n' % (lon, lat, val))
    fp.close()


def plot_cmbhf(timestep, age, storedir,
               region='0/360/-90/90',
               proj='N180/10.0i',
               grd_min=None, grd_max=None,
               res='0.5/0.5'):
    '''Grid and plot CMB heat flux with GMT4.'''

    xyzfile  = 'cmbhf.xyz'
    grdfile  = 'cmbhf.grd'
    cptfile  = 'cmbhf.cpt'
    psfile   = 'cmbhf_%dMa.ps' % age
    pdffile  = 'cmbhf_%dMa.pdf' % age

    # --- read and write combined xyz ---
    points = read_botm_caps(timestep)
    if not points:
        print 'Error: no botm data found for timestep %d' % timestep
        sys.exit(1)
    write_xyz(points, xyzfile)
    print 'Read %d points from 12 caps' % len(points)

    # --- auto range if not specified ---
    if grd_min is None or grd_max is None:
        vals = [p[2] for p in points]
        vmin, vmax = min(vals), max(vals)
        half = max(abs(vmin), abs(vmax))
        grd_min = -half
        grd_max =  half

    # --- GMT4 settings ---
    cmd = ('gmtset PAPER_MEDIA A2 MEASURE_UNIT inch '
           'ANNOT_FONT_SIZE_PRIMARY 14p ANNOT_FONT_SIZE_SECONDARY 14p '
           'ANNOT_OFFSET_PRIMARY 10p HEADER_FONT_SIZE 18p MAP_TITLE_OFFSET 10p')
    os.system(cmd)

    # --- grid: blockmean then surface ---
    cmd = 'blockmean %s -I%s -R%s -V > mean.xyz' % (xyzfile, res, region)
    os.system(cmd)

    cmd = ('surface mean.xyz -G%s -I%s -R%s -T0.25 -Ll%f -Lu%f -V'
           % (grdfile, res, region, grd_min, grd_max))
    os.system(cmd)

    # --- color palette: blue-white-red diverging ---
    cmd = 'makecpt -Cpolar -T%f/%f/%f -D > %s' % (grd_min, grd_max,
                                                    (grd_max - grd_min) / 100.0,
                                                    cptfile)
    os.system(cmd)

    # --- base map: grdimage ---
    title = 'CMB Heat Flux at %d Ma' % age
    cmd = ('grdimage %s -J%s -C%s -R%s -B30g30/15g15:."CMB Heat Flux %d Ma": '
           '-X2.5 -Y16.0 -V -K -P > %s'
           % (grdfile, proj, cptfile, region, age, psfile))
    os.system(cmd)

    # --- coastlines ---
    cmd = 'pscoast -J -R -B -Dc -A5000 -W0.5p,black -P -O -K >> %s' % psfile
    os.system(cmd)

    # --- plate boundaries from storedir if available ---
    plates_dir = os.path.join(storedir, 'GPlates')
    if os.path.isdir(plates_dir):
        trench_sL = os.path.join(plates_dir, 'topology_subduction_boundaries_sL_%.2fMa.xy' % age)
        trench_sR = os.path.join(plates_dir, 'topology_subduction_boundaries_sR_%.2fMa.xy' % age)
        ridges    = os.path.join(plates_dir, 'topology_ridge_transform_boundaries_%.2fMa.xy' % age)

        if os.path.exists(trench_sL):
            cmd = 'psxy %s -J -R -B -m -W1p,blue -Sf0.2/0.07+l+t -Gblue -P -O -K >> %s' % (trench_sL, psfile)
            os.system(cmd)
        if os.path.exists(trench_sR):
            cmd = 'psxy %s -J -R -B -m -W1p,blue -Sf0.2/0.07+r+t -Gblue -P -O -K >> %s' % (trench_sR, psfile)
            os.system(cmd)
        if os.path.exists(ridges):
            cmd = 'psxy %s -J -R -B -m -W1p,red -P -O -K >> %s' % (ridges, psfile)
            os.system(cmd)

    # --- color scale bar ---
    cmd = ('psscale -C%s -Ba0.1/0.05:"CMB Heat Flux (nondim)": '
           '-D5.0/-0.5/4.0/0.3h -V -O >> %s' % (cptfile, psfile))
    os.system(cmd)

    # --- convert to PDF ---
    cmd = 'ps2pdf %s %s' % (psfile, pdffile)
    os.system(cmd)
    print 'Output: %s  %s' % (psfile, pdffile)

    # --- cleanup ---
    os.system('rm -f %s mean.xyz %s %s' % (xyzfile, grdfile, cptfile))


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print __doc__
        sys.exit(1)

    timestep = int(sys.argv[1])
    age      = int(sys.argv[2])
    storedir = sys.argv[3]

    plot_cmbhf(timestep, age, storedir)
