#!/bin/sh

list="cullpdblist"

# Enter the number of jobslots on each machine. 0 means as many as possible.
# Default is 100% which will run one job per CPU on each machine.
NP=0

# Enter your path
export PDB_DIR="/raid1/share/database/PDB_uncompressed/divided_dot_pdb/pdb"
export INET_HOST="your_fileserver_name_of_PDB_DIR"
export BLASTDB="/raid1/share/database/blast/db"
export BLASTMAT="/raid1/share/database/blast/matrices"

# Usage: python2.7 pdb2vall.py -p 7odcA
# For parallel generation, GNU parallel must be installed first. 'yum -y install parallel (CentOS 7)' or 'brew install parallel (macOS)'
# The argument '-p' should be like '7odcA'. Input '7ODCA' or '7odca' will cause an error.

cat ${list} | tr '[a-z]' '[A-Z]' | sed -e "s/\(.\{4\}\)/\L\1/g" | \
    parallel --noswap -j ${NP} -a - "(mkdir -p {} && cd {} && ../../pdb2vall/pdb2vall.py -d -p {} )" 
