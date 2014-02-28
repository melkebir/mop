Extracting molecules
--------------------

    curl "http://compbio.chemistry.uq.edu.au/atb/index.py?molsPerPage=all&pageSelected=1&tab=existing_tab&display_moltype=lipid&display_curation=0" | grep molid | grep ">[0-9]\+<" | sed -e "s/.*>\([0-9]\+\).*/\1/g" | sort -n > lipids.txt
    ./getMolecules.py < lipids.txt

Generating repository
---------------------

    for f in *.pdb;do echo "Processing $f";../../../build/atb2lgf -mtb `basename $f .pdb`.mtb -pdb $f > `basename $f .pdb`.lgf;done
    grep '<a href="./molecule.py?molid=.*">Link</a>' atb.html | sed -e 's/.*molid=\(.*\)".*/\1/g' > ATB_ids.txt
