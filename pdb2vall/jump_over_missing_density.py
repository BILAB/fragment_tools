#!/usr/bin/env python

from sys import argv, exit, stderr
from os.path import exists

def jump_over_missing_density( fasta_fn, pdb_fn, request ):
    missing_density_rsn_list = open( fasta_fn, "r").readline().split("missing_density_rsn:")[1:][0].strip().split()
    pdb_xyz_file     = []
    pdb_cb_file      = []
    pdb_cen_file     = []
    pdb_torsion_file = []

    masslst          = dict()
    masslst["C"] = 12.0107
    masslst["O"] = 15.9994
    masslst["S"] = 32.066
    masslst["N"] = 14.00674
    masslst["H"] = 1.00794
    

    ## RENUMBER IDEALIZED PDB BACK TO PDB_SEQRES
    pdb_fn = pdb_fn[:5]
    if exists( pdb_fn + '_0001_0001.pdb' ):

        file = open( pdb_fn + '_0001_0001.pdb', "r")
        line = file.readline()

        newnum_xyz     = 1
        newnum_cb      = 1
        newnum_cen     = 1
        newnum_torsion = 1
        
        cen_mass = 0.0
        cen_xyz = [0.0,0.0,0.0]

        while line:
            # obtain CA xyz coordinates
            # centroid calculation expects all sidechain atoms come after CA in the residue 
            if line[12:16] == " CA ":
                # push centroid coords
                if len(pdb_xyz_file) >= 1:
                    if pdb_xyz_file[-1][17:20] == "GLY":
                        cen_line_edit = "%s CEN%s" % (pdb_xyz_file[-1][:12], pdb_xyz_file[-1][16:])
                        pdb_cen_file.append( cen_line_edit )
                    elif pdb_xyz_file[-1][17:20] == "ALA":
                        cen_line_edit = "%s CEN%s" % (pdb_cb_file[-1][:12], pdb_cb_file[-1][16:])
                        pdb_cen_file.append( cen_line_edit )
                    elif cen_mass != 0.0:
                        coords = "%8.3f%8.3f%8.3f" % (cen_xyz[0]/cen_mass, cen_xyz[1]/cen_mass, cen_xyz[2]/cen_mass)
                        cen_line_edit = "%s CEN%s%s%s" % (pdb_xyz_file[-1][:12], pdb_xyz_file[-1][16:30], coords, pdb_xyz_file[-1][54:] )
                        pdb_cen_file.append( cen_line_edit )
                cen_mass = cen_xyz[0] = cen_xyz[1] = cen_xyz[2] = 0.0

                if str(newnum_xyz) in missing_density_rsn_list:
                    newnum_xyz = newnum_xyz + 1
                    newnum_cb = newnum_cb + 1
                    newnum_cen = newnum_cen + 1
                    continue

                newnum_cen += 1

                xcord = line[30:38]
                ycord = line[38:46]
                zcord = line[46:54]
                xyz_line_edit = "%s%4d %s" %( line[:22], newnum_xyz, line[27:] )
                pdb_xyz_file.append( xyz_line_edit )
                newnum_xyz = newnum_xyz + 1

            # Obtain CB xyz coordinates. if GLY, Obtain CA xyz instead.
            if ((line[17:20] != "GLY") and (line[12:16] == " CB ")):

                cbxcord = line[30:38]
                cbycord = line[38:46]
                cbzcord = line[46:54]
                cbxyz_line_edit = "%s%4d %s" %( line[:22], newnum_cb, line[27:] )
                pdb_cb_file.append( cbxyz_line_edit )
                newnum_cb = newnum_cb + 1
            if ((line[17:20] == "GLY") and (line[12:16] == " CA ")):

                cbxcord = line[30:38]
                cbycord = line[38:46]
                cbzcord = line[46:54]
                cbxyz_line_edit = "%s%4d %s" %( line[:22], newnum_cb, line[27:] )
                pdb_cb_file.append( cbxyz_line_edit )
                newnum_cb = newnum_cb + 1
            if line.startswith("ATOM") and (line[12:16].strip() not in ["C","N","O","CA","CB"] and line[77] in ['C','O','S','N']):
                thismass = masslst[line[77:78]] 
                cen_mass += thismass
                cen_xyz[0] += float(line[30:38]) * thismass
                cen_xyz[1] += float(line[38:46]) * thismass
                cen_xyz[2] += float(line[46:54]) * thismass
            
            if line.startswith("REMARK") and "torsion" not in line:
                if str(newnum_torsion) in missing_density_rsn_list:
                    newnum_torsion = newnum_torsion + 1
                    continue

                idealized_rsd = line.split()[4]
                if idealized_rsd != "X":
                    secstr        = line.split()[5]
                    phi           = line.split()[6]
                    psi           = line.split()[7]
                    omega         = line.split()[8]

                    torsion_line_edit = "%s %s %s %s %s %s\n" %( newnum_torsion, idealized_rsd, secstr, phi, psi, omega )
                    pdb_torsion_file.append( torsion_line_edit )
                    newnum_torsion = newnum_torsion + 1

            line = file.readline()
        file.close()

        if len(pdb_xyz_file) >= 1:
            if pdb_xyz_file[-1][17:20] == "GLY":
                cen_line_edit = "%s CEN%s" % (pdb_xyz_file[-1][:12], pdb_xyz_file[-1][16:])
                pdb_cen_file.append( cen_line_edit )
            elif pdb_xyz_file[-1][17:20] == "ALA":
                cen_line_edit = "%s CEN%s" % (pdb_cb_file[-1][:12], pdb_cb_file[-1][16:])
                pdb_cen_file.append( cen_line_edit )
            elif cen_mass != 0.0:
                coords = "%8.3f%8.3f%8.3f" % (cen_xyz[0]/cen_mass, cen_xyz[1]/cen_mass, cen_xyz[2]/cen_mass)
                cen_line_edit = "%s CEN%s%s%s" % (pdb_xyz_file[-1][:12], pdb_xyz_file[-1][16:30], coords, pdb_xyz_file[-1][54:] )
                pdb_cen_file.append( cen_line_edit )
    else:
        stderr.write("ERROR: there's no idealized and relaxed %s_0001_0001.pdb here!\n"  % pdb_fn )
        return 0
        exit()


    if request == "xyz":
        return pdb_xyz_file
    if request == "cb":
        return pdb_cb_file
    if request == "centroid":
        return pdb_cen_file
    if request == "torsion":
        return pdb_torsion_file

if __name__ == '__main__':
    if len( argv ) < 3:
        print
        print "USAGE: %s fasta_fn pdb_fn xyz/cb/centroid/torsion" % argv[0]
        print
        print "-"*75
        exit()

    print jump_over_missing_density( argv[1], argv[2], argv[3] )
    #print "".join(jump_over_missing_density( argv[1], argv[2], argv[3] ))


