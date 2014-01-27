#!/usr/bin/python
import sys
import os
import subprocess

def init_rep(rep_dir):
  rep = {}
  for filename in os.listdir(rep_dir):
    if filename.endswith(".lgf"):
      mol_id = filename.rstrip(".lgf")
      rep[mol_id] = rep_dir + "/" + filename 
  return rep

def run(input_file, rep_mol_id, rep_mol_filename, fragments_executable, shell_size):
  r = subprocess.check_output([" ".join([fragments_executable, " -s", shell_size, "-atb_id", rep_mol_id, input_file, rep_mol_filename])], shell=True)
  sys.stdout.write(r)

rep_dir = sys.argv[1]
input_file = sys.argv[2]
fragments_bin = sys.argv[3]
shell_size = sys.argv[4]

rep = init_rep(rep_dir)

first = True
sys.stdout.write("{\n")
sys.stdout.write("  molecules: [\n")
for mol_id in rep:
  if first:
    first = False
  else:
    sys.stdout.write(",\n")
  run(input_file, mol_id, rep[mol_id], fragments_bin, shell_size)
sys.stdout.write("\n")
sys.stdout.write("  ]\n")
sys.stdout.write("}\n")
