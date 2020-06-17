
import numpy            as np
import pandas           as pd
import glob, linecache, os, sys
from Bio                import PDB
from Bio.PDB.PDBParser  import PDBParser
from Bio.PDB.DSSP       import dssp_dict_from_pdb_file
from Bio.PDB.PDBList    import PDBList
from io                 import StringIO, BytesIO
from sys                import argv
from urllib             import parse, request

def ask_ok(prompt, retries=4, complaint= 'Yes or no, please!'):
    while True:
        ok = input(prompt)
        if ok in ('y', 'ye', 'yes', 'Y'):
            return True
        if ok in ('n', 'no', 'nop', 'nope', 'N'):
            return False
        retries = retries - 1
        if retries < 0:
            raise IOError('refusenik user')
        print(complaint)

def checkgen(f):
    p = open(f, 'r')    # re-opening file
    data = linecache.getline(f, 19).split()[0]
    return float(data)

def checkpdb(f):
    p = open(f, 'r')    # re-opening file
    name = ''
    name = linecache.getline(f, 16).split()[1]
    return name

def genparse(f):

    p = open(f, 'r')        # re-opening file
    for i in range(18):     # read and ignore header lines
        header = p.readline()
    # Loop over lines and extract variables of interest
    final = []
    for line in p:
        line = line.strip()
        columns = line.split()
        if float(columns[0]) <= 260 and float(columns[0]) >= 190:
            if float(columns[0]).is_integer() == True:
                final.append(float(columns[1]))
    p.close()
    finale = np.transpose(final)
    if np.size(finale) == 71:
        return finale

def getpdbs(names):
    cwd = os.getcwd()
    pdbl = PDBList()
    pdbl.download_pdb_files(names, obsolete = False, file_format = "pdb", pdir=cwd)
    pdbl.download_pdb_files(names, obsolete = True, file_format = "pdb", pdir=cwd)
    # In module PDBList, I muted lines 300-301
    # and replaced line 293 with the line below:
    # final = {'pdb': '%s.pdb', 'mmCif': '%s.cif', 'xml': '%s.xml',

def nrmsdnew(spec, out, pdblist):

    m,n = spec.shape    # ".shape" returns a tuple of array dimensions;
                        # "spec.shape" gets the current shape of the passed array "spec".

    np.seterr(divide='ignore', invalid='ignore')
    #NUMERATOR WORK
    RMSD = []
    for i in range(1,n):
        if np.amax(spec[:,i]) == 0:
            print("error in array \"spec\": ", pdblist[i-1], "\t Index: ", pdblist.index(pdblist[i-1]))

        elif np.amax(out[:,i]) == 0:
            print("error in array \"out\": ", pdblist[i-1], "\t Index: ", pdblist.index(pdblist[i-1]))
        delta = spec[:,i] - out[:,i]
        nrm_max_spec = spec[:,i]/np.amax(np.abs(spec[:,i]))
        nrm_max_out = out[:,i]/np.amax(np.abs(out[:,i]))
        deltamax = nrm_max_spec - nrm_max_out
        RMSD.append(((np.sum(np.square(delta))) / (np.sum(np.square(deltamax))))**0.5)
    return dict(zip(pdblist, list(RMSD)))

def pdb2cd(name):
    f = name + ".pdb"
    dssp_tuple = dssp_dict_from_pdb_file(f)
    dssp_dict = dssp_tuple[0]
    p = PDBParser(QUIET = True).get_structure("file", f)



    # Initiates and fills array ("cc") with chains.
    cc = [chain.get_id()    for model in p    for chain in model]

    # Determines length of sequence, initiates an array ("ss") of same length.
    howLong = ss_out = 0
    for c in cc:
        howLong += len([_ for _ in p[0][c].get_residues() if PDB.is_aa(_)])
    if not howLong == len(dssp_tuple[1]): howLong = len(dssp_tuple[1])
    ss = np.arange(1, howLong+1)

    # Fills the array ("ss") with secondary structures.
    for i in ss:
        ss_lib = dssp_dict[dssp_tuple[1][i-3]]         # ss_lib = dssp_dict[(dssp_tuple[1][0][0], (' ', i-1, ' '))]
        dict_ss = ss_lib[1]
        if dict_ss == 'H':
            ss_out = 0
        if dict_ss == 'E':
            ss_out = 1
        if dict_ss == '-': # else:# dict_ss == '-':
            ss_out = 2
        ss[i-1] = ss_out
    # Returns the fractional composition of alpha helix, beta sheet or random coil.
    alpha = (ss == 0).sum() / ss.__len__()
    beta = (ss == 1).sum() / ss.__len__()
    coil = (ss == 2).sum() / ss.__len__()
    abc = [alpha, beta, coil]
    return abc

def queryuni(uni):
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from': 'ACC+ID',
    'to': 'PDB_ID',
    'format': 'tab',
    'query': uni
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
       c = pd.read_csv(
       io.StringIO(response.decode("utf-8")),
       sep='\s+',
       engine='python',
       usecols=[0,1]
       )

       return c
