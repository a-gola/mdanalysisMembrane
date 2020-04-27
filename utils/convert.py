#!/usr/bin/env python3

import os
import os.path

from MDAnalysis import Universe, Writer

from tqdm import tqdm

import numpy as np

import parmed as pmd

import subprocess


def xtc2gro(instr, intrj, timestep):
    """ Convert trajectory using MD Analysis IO ability

    Arguments:
        instr {str} -- input structure file name, e.g. .gro file
        intrj {str} -- input trajectroy file name, e.g. .xtc file
        timestep {int} -- closest time to look for in the trj

    """

    # root, oldext = os.path.basename(intrj).rsplit('.', 1)

    # outtrj = "{0}_every-{1}-frame{2}".format(root, n_every, ext)

    # if os.path.isfile(outtrj):
    #     outtrj = "{0}_every-{1}-frame_2{2}".format(root, n_every, ext)

    u = Universe(instr, intrj)

    # create a writer instance for the output trajectory
    # n_atoms = len(u.atoms)
    # w = Writer(outtrj, n_atoms)

    # loop through the trajectory and write a frame for every step
    for ts in tqdm(u.trajectory):
        if ts.time == int(timestep):
            u.atoms.write('extract_{}.gro'.format(timestep), reindex=False)
            quit()
        # w.write(ts)
    # w.close()
    # print("Converted {0!r} --> {1!r}".format(intrj, outtrj))


def trjconv(instr, intrj, ext=".ncdf", n_every=1):
    """ Convert trajectory using MD Analysis IO ability

    Arguments:
        instr {str} -- input structure file name, e.g. .gro file
        intrj {str} -- input trajectroy file name, e.g. .xtc file

    Keyword Arguments:
        ext {str} -- extension for the output trajectory (default: {".ncdf"})
    """

    root, oldext = os.path.basename(intrj).rsplit('.', 1)

    outtrj = "{0}_every-{1}-frame{2}".format(root, n_every, ext)

    if os.path.isfile(outtrj):
        outtrj = "{0}_every-{1}-frame_2{2}".format(root, n_every, ext)

    u = Universe(instr, intrj)

    # create a writer instance for the output trajectory
    n_atoms = len(u.atoms)
    w = Writer(outtrj, n_atoms)

    # loop through the trajectory and write a frame for every step
    for ts in tqdm(u.trajectory[::n_every]):
        w.write(ts)
    w.close()
    print("Converted {0!r} --> {1!r}".format(intrj, outtrj))


def gro2data(intop, instr):
    """ Convert Gromacs structure and topology file 'type' *.gro and *.top
        to lammps data file that can be open with OVITO.
        -1st stage is the conversion of the .top file to .psf file using Parmed
        -2nd stage is to open *.gro in VMD and add the topology information
             in VMD using the *.psf file
        -3rd stage is to write a lammps datafile (full type) using the
             topotools plugin in VMD

    Arguments:
        intop {str} -- input topology filename
        instr {str} -- input strucure filename
    """

    root, oldext = os.path.basename(instr).rsplit('.', 1)

    print("converting *.top to *.psf with ParmEd", flush=True)
    if not os.path.isfile(root+"_tmp.psf"):
        gmx_top = pmd.load_file(intop, xyz=instr)
        gmx_top.save(root+"_tmp.psf", vmd=True)
    else:
        print("*.psf found, continuing...")

    outdata = root + ".data"

    print("converting *.gro to *.data with VMD", flush=True)
    with open("vmd_tcl.tmp", 'w') as tmp_tcl:
        tmp_tcl.write('if {[catch {package require topotools 1.1} ver]} {\n\
                      vmdcon -error "$ver. This script requires at least \
                      TopoTools v1.1. Exiting..."\n   quit\n}\n')
        tmp_tcl.write("mol new {} autobonds no\n".format(instr))
        tmp_tcl.write("mol addfile {} \n".format(root+"_tmp.psf"))
        tmp_tcl.write("topo writelammpsdata {} full\n".format(outdata))
        tmp_tcl.write("quit")

    #os.system("vmd -dispdev text -e vmd_tcl.tmp")
    subprocess.call(["/bin/sh", "-i", "-c", "vmd -dispdev text -e vmd_tcl.tmp"]) # Note: MacOS users might have to replace /bin/sh with /bin/zsh

    # os.remove(root+"_tmp.psf")  # psf conversion is a time consuming step
                                  # might be useful to keep this file
    os.remove("vmd_tcl.tmp")

    print("Converted {0!r} --> {1!r}".format(instr, outdata), flush=True)


def changeMolIds(data):
    """ Change molecule IDs in LAMMPS data file to bypass Gromacs limitation
        of 100,000 molecule IDs per data file
    Arguments:
        data {str} -- input LAMMPS data filename
    """

    root, oldext = os.path.basename(data).rsplit('.', 1)

    print("Renumberring molecules in {}".format(data), flush=True)

    outdata = root + "_renum.data"

    infile = open(data, 'r').readlines()
    with open(outdata, 'w') as outfile:
        mol_id = 1
        atom_flag = False
        first_atom = True  # Flag to detect first atom being read
        for l in infile:
            if l.strip().startswith("Bonds"):
                # Exiting in Atoms section
                atom_flag = False
            elif l.strip().startswith("Atoms"):
                # Entering in Atoms section containing Mol Ids
                atom_flag = True
                outfile.write(l)
                continue

            if atom_flag:
                l_data = l.strip().split()
                if len(l_data) == 0:  # empty line
                    outfile.write(l)
                else:
                    if first_atom:
                        first_atom = False
                        former_id = int(l_data[1])  # keeping trace of former ID
                        l_data[1] = str(mol_id)
                    else:
                        if int(l_data[1]) != former_id: # Detecting mol id jump
                            mol_id += 1
                            former_id = int(l_data[1])
                        l_data[1] = str(mol_id)
                    new_line = " ".join(l_data)
                    outfile.write(new_line+'\n')
            else:
                outfile.write(l)

    print("Converted {0!r} --> {1!r}".format(data, outdata), flush=True)
