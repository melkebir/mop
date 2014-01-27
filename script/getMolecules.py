#!/usr/bin/python
import urllib2
import sys

def getTopology(mol_id):
  l = urllib2.urlopen("http://compbio.chemistry.uq.edu.au/atb/molecule.py?molid={0}&outputType=top&atbVersion=v2Top&ffVersion=Gromos#files".format(mol_id)).read()
  qm0 = l.find("qm0") != -1
  # fetch mtb file
  f = open("{0}.mtb".format(mol_id), "w")
  if qm0:
    l = urllib2.urlopen("http://compbio.chemistry.uq.edu.au/atb/download.py?molid={0}&outputType=v2Top&file=itp_allatom_qm0".format(mol_id)).read()
  else:
    l = urllib2.urlopen("http://compbio.chemistry.uq.edu.au/atb/download.py?molid={0}&outputType=v2Top&file=itp_allatom".format(mol_id)).read()
  f.write(l)
  f.close()
  f = open("{0}.pdb".format(mol_id), "w")
  if qm0:
    l = urllib2.urlopen("http://compbio.chemistry.uq.edu.au/atb/download.py?molid={0}&outputType=v2Top&file=pdb_allatom_qm0_optimised".format(mol_id)).read()
  else:
    l = urllib2.urlopen("http://compbio.chemistry.uq.edu.au/atb/download.py?molid={0}&outputType=v2Top&file=pdb_allatom_optimised".format(mol_id)).read()
  f.write(l)
  f.close()
  return qm0

if __name__ == "__main__":
  for line in sys.stdin:
    m = line.rstrip("\n")
    print "Fetching {0}...".format(m),
    sys.stdout.flush()
    qm0 = getTopology(m)
    print "Done! QM0:", qm0
