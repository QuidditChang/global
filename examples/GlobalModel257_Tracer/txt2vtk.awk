#!/bin/awk
{
	#id[NR-1]=$1
	#type[NR-1]=$2
	#object[NR-1]=$3
	#pos[NR-1]=$4 " " $5 " " $6
	pos[NR-1]=6371000.0 * $3 * sin($1) * cos($2) " " 6371000.0 * $3 * sin($1) * sin($2)" " 6371000.0 * $3 * cos($1)
	#pos[NR-1]=$1 " " $2 " " $3
	#pos[NR-1]=-log((1+cos($1))/(1-cos($1)))/2.0 " " $2 " " $3
	#pos[NR-1]=90 - $1 / 3.14159 * 180 " " $2 / 3.14159 * 180 " " $3 * 6371.0
	colat[NR-1]=$1/3.1415926*180.0
	lon[NR-1]=$2/3.1415926*180.0
	depth[NR-1]=(1-$3)*6371.0
	#vel[NR-1]=$7 " " $8 " " $9
	#vel[NR-1]=$4 * 0.000495 " " $5 * 0.000495 " " $6 * 0.000495
	vel[NR-1]=sin($1)*cos($2)*$6+cos($1)*cos($2)*$4-sin($1)*sin($2)*$5 " " sin($1)*sin($2)*$6+$3*cos($1)*sin($2)*$4+sin($1)*cos($2)*$5 " " cos($1)*$6-sin($1)*$4
	#mass[NR-1]=$10
	#veloz[NR-1]=$6
	temp[NR-1]=$7
	#density[NR-1]=$11
	#pressure[NR-1]=$12
	viscosity[NR-1]=$8
	composition[NR-1]=$9
	numParts=NR
}

END {
	print "# vtk DataFile Version 2.0"
	print "Converted from TextWriter format"
	print "ASCII\nDATASET UNSTRUCTURED_GRID"

	printf "POINTS %u float\n", NR
	for (i=0; i < numParts; ++i)
		print pos[i]
	print "\n"

        nx=129
	#nx=257
        ny=129
	#ny=257
        #nz=65
	nz=65
	printf "CELLS %u %u\n", (nx-1)*(ny-1)*(nz-1), 9*(nx-1)*(ny-1)*(nz-1)
	for (i=0; i<nx-1; i++)
	    for (j=0; j<ny-1; j++)
		for (k=0; k<nz-1; k++)
		    print 8, j*nx*nz+i*nz+k, j*nx*nz+i*nz+k+1, j*nx*nz+(i+1)*nz+k, j*nx*nz+(i+1)*nz+k+1, (j+1)*nx*nz+i*nz+k, (j+1)*nx*nz+i*nz+k+1, (j+1)*nx*nz+(i+1)*nz+k, (j+1)*nx*nz+(i+1)*nz+k+1
	print "\n"

	#printf "CELLS %u %u\n", numParts, 2*numParts
	#for (i=0; i < numParts; ++i)
	#	printf "1 %u\n", i
	#print "\n"

	printf "CELL_TYPES %u\n", (nx-1)*(ny-1)*(nz-1)
	for (i=0; i < (nx-1)*(ny-1)*(nz-1); ++i)
		print "11"
	print "\n"

	printf "POINT_DATA %u\n", numParts

	print "VECTORS Velocity float"
	for (i=0; i < numParts; ++i)
		print vel[i]
	print "\n"

	print "SCALARS Viscosity float"
	print "LOOKUP_TABLE default"
	for (i=0; i < numParts; ++i)
		print viscosity[i]
	print "\n"

	print "SCALARS Colat float"
	print "LOOKUP_TABLE default"
	for (i=0; i < numParts; ++i)
		print colat[i]
	print "\n"

        print "SCALARS Lon float"
        print "LOOKUP_TABLE default"
        for (i=0; i < numParts; ++i)
                print lon[i]
        print "\n"

        print "SCALARS Depth float"
        print "LOOKUP_TABLE default"
        for (i=0; i < numParts; ++i)
                print depth[i]
        print "\n"

	print "SCALARS Temperature float"
	print "LOOKUP_TABLE default"
	for (i=0; i < numParts; ++i)
		print temp[i]
	print "\n"

	print "SCALARS Composition float"
        print "LOOKUP_TABLE default"
        for (i=0; i < numParts; ++i)
                print composition[i]
        print "\n"

        #print "SCALARS Velocity float"
        #print "LOOKUP_TABLE default"
        #for (i=0; i < numParts; ++i)
        #        print veloz[i]
        #print "\n"

	#print "SCALARS Type int"
	#print "LOOKUP_TABLE default"
	#for (i=0; i < numParts; ++i)
	#	print type[i]
	#print "\n"

	#print "SCALARS Object int"
	#print "LOOKUP_TABLE default"
	#for (i=0; i < numParts; ++i)
	#	print object[i]
	#print "\n"

	#print "SCALARS ParticleId int"
	#print "LOOKUP_TABLE default"
	#for (i=0; i < numParts; ++i)
	#	print id[i]
	#print "\n"
}
