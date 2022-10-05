../../nma/nmanu.pl ../data/structure.ca.pdb hessian.dat 1 0 40
../../nma/diaghess/diaghess
../../nma/mc-eigen.pl eigenvec.dat > file.proj
../../nma/pca_anim_mc.pl -pdb ../data/structure.ca.pdb -evec eigenvec.dat -i file.proj -n 50 -pout traj.crd

