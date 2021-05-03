import numpy as np
import ase.io as aio
from pyscf import gto, dft
import argparse

parser = argparse.ArgumentParser(description='Process the input.')
parser.add_argument('--mol', type=str, nargs='+', dest='filename',
                    help='Molecular structure in xyz format')
parser.add_argument('--basis', type=str, nargs='+', dest='basis',
                    help='Basis set for DFT computation')
parser.add_argument('--auxbasis', type=str, nargs='+', dest='auxbasis',
                    help='Basis set for decomposition')
parser.add_argument('--func', type=str, nargs='+', dest='xc_func',
                    help='XC functional for computation')


args = parser.parse_args()


def read_xyz(filename):
    return aio.read(filename, format='xyz')


def convert_ASE_PySCF(molecules, basis, auxbasis):
    # Take an ASE atom and return  PySCF mol and auxmol objects.

    sym = molecules.get_chemical_symbols()
    pos = molecules.positions
    tmp_mol = []
    for s, (x, y, z) in zip(sym, pos):
        tmp_mol.append([s, (x, y, z)])

    pyscf_mol = gto.M(atom=tmp_mol, basis=basis)
    pyscf_auxmol = gto.M(atom=tmp_mol, basis=auxbasis)

    return pyscf_mol, pyscf_auxmol


def do_RKS(mol, xc):
    # Run an DFT computation for each molecule and return dm.
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.verbose = 1
    mf.run()
    # print("Convergence: ",mf.converged)
    # print("Energy: ",mf.e_tot)
    dm = mf.make_rdm1()

    return dm


def get_integrals(mol, auxmol):

    S = auxmol.intor('int1e_ovlp_sph')
    pmol = mol + auxmol

    eri2c = auxmol.intor('int2c2e_sph')
    eri3c = pmol.intor('int3c2e_sph', shls_slice=(
        0, mol.nbas, 0, mol.nbas, mol.nbas, mol.nbas+auxmol.nbas))
    eri3c = eri3c.reshape(mol.nao_nr(), mol.nao_nr(), -1)

    return S, eri2c, eri3c


def get_coeff(dm, eri2c, eri3c):
    rho = np.einsum('ijp,ij->p', eri3c, dm)
    c = np.linalg.solve(eri2c, rho)
    return c


def main():
    print('Reading molecule')
    compounds = aio.read(args.filename[0], format='xyz')
    mol, auxmol = convert_ASE_PySCF(compounds, args.basis[0], args.auxbasis[0])
    print('Performing DFT computation')
    dm = do_RKS(mol, args.xc_func[0])
    print('Computing integrals')
    S, eri2c, eri3c = get_integrals(mol, auxmol)
    print('Computing coefficients')
    c = get_coeff(dm, eri2c, eri3c)

    print('Saving coefficients')
    np.save('dm', dm)
    np.save('S', S)
    np.save('J', eri2c)
    np.save('coeff', c)


if __name__ == "__main__":
    main()
