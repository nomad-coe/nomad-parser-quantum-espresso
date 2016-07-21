import setup_paths
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import os
import sys
import json
import re
import logging
import nomadcore.unit_conversion.unit_conversion as unit_conversion
import math
import numpy as np
import QuantumEspressoCommon as QeC
from nomadcore.parser_backend import valueForStrValue
from QuantumEspressoCommon import RE_f, RE_i


LOGGER = logging.getLogger(__name__)


class QuantumEspressoParserPWSCF(QeC.ParserQuantumEspresso):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self):
        QeC.ParserQuantumEspresso.__init__(self)

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse
        different files"""
        self.secMethodIndex = None
        self.secSystemDescriptionIndex = None
        self.tmp = {}

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # reset values if same superContext is used to parse different files
        self.initialize_values()

    def onClose_section_basis_set_cell_dependent(
            self, backend, gIndex, section):
        pwc_ry = unit_conversion.convert_unit(
            section['basis_set_planewave_cutoff'][-1],
            'J', 'rydberg')
        backend.addValue('basis_set_cell_dependent_name', 'PW_%.1f' % (pwc_ry))
        backend.addValue('basis_set_cell_dependent_kind', 'plane_waves')

    def onOpen_section_method(
            self, backend, gIndex, section):
        self.secMethodIndex = gIndex
        self.cache_t_pseudopotential = {}

    def onOpen_section_system(
            self, backend, gIndex, section):
        self.secSystemIndex = gIndex

    def onClose_x_qe_t_section_pseudopotential(
            self, backend, gIndex, section):
        self.cache_t_pseudopotential[section['x_qe_t_pp_label'][0]] = section

    def onClose_section_method_atom_kind(
            self, backend, gIndex, section):
        pp_label = section['x_qe_pp_label'][0]
        pp = self.cache_t_pseudopotential[pp_label]
        for key, value in pp.simpleValues.items():
            if key == 'x_qe_t_pp_label' or key == 'x_qe_pp_valence':
                continue
            target = re.sub(r'^x_qe_t_',r'x_qe_',key)
            if target == key:
                raise RuntimeError('found non-temporary key in pseudopotential cache: "%s"' % (key))
            backend.addValue(target, value[-1])

    def onClose_section_method(
            self, backend, gIndex, section):
        # set flag if we deal with user-enforced XC functional
        if section['x_qe_t_xc_functional_shortname_enforced'] is not None:
            backend.addValue('x_qe_xc_functional_user_enforced', True)
        # translate XC functional to section_xc_functionals
        xc_functional_num = section['x_qe_t_xc_functional_num'][-1]

    def onClose_section_single_configuration_calculation(
            self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""
        backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodIndex)
        backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemIndex)

    def onOpen_section_eigenvalues(self, backend, gIndex, section):
        self.tmp['k_point'] = []
        self.tmp['k_energies'] = []

    def onClose_section_eigenvalues(self, backend, gIndex, section):
        if len(self.tmp['k_point']) > 0:
            backend.addArrayValues('eigenvalues_values', np.array([
                self.tmp['k_energies']
            ], dtype=np.float64))
            # k-points are in cartesian, but metaInfo specifies crystal
            k_point_cartesian = np.array(self.tmp['k_point'], dtype=np.float64)
            k_point_crystal = self.bmat_inv.dot(k_point_cartesian.T).T
            backend.addArrayValues('eigenvalues_kpoints', k_point_crystal)

    def onClose_section_system(self, backend, gIndex, section):
        # store direct lattice matrix for transformation crystal -> cartesian
        self.amat = np.array([
            section['x_qe_t_vec_a_x'], section['x_qe_t_vec_a_y'], section['x_qe_t_vec_a_z'],
        ], dtype=np.float64).T
        # store inverse for transformation cartesian -> crystal
        try:
            self.amat_inv = np.linalg.inv(self.amat)
        except np.linalg.linalg.LinAlgError:
            raise RuntimeError("error inverting bravais matrix " + str(self.amat))
        # store reciprocal lattice matrix for transformation crystal -> cartesian
        self.bmat = np.array([
            section['x_qe_t_vec_b_x'], section['x_qe_t_vec_b_y'], section['x_qe_t_vec_b_z'],
        ], dtype=np.float64).T
        # store inverse for transformation cartesian -> crystal
        try:
            self.bmat_inv = np.linalg.inv(self.bmat)
        except np.linalg.linalg.LinAlgError:
            raise RuntimeError("error inverting reciprocal cell matrix")
        backend.addArrayValues('simulation_cell', self.amat)
        backend.addArrayValues('x_qe_reciprocal_cell', self.bmat)
        # atom positions
        atpos_cart = np.array([
            section['x_qe_t_atpos_x'], section['x_qe_t_atpos_y'], section['x_qe_t_atpos_z']
        ], dtype=np.float64).T
        backend.addArrayValues('atom_positions',atpos_cart)
        backend.addArrayValues('atom_labels',np.asarray(section['x_qe_t_atom_labels']))
        backend.addArrayValues('x_qe_atom_idx',np.array(section['x_qe_t_atom_idx']))

    def onOpen_x_qe_t_section_kbands(self, backend, gIndex, section):
        self.tmp['this_k_energies'] = ''

    def onClose_x_qe_t_section_kbands(self, backend, gIndex, section):
        ek_split = []
        for energy in re.split(r'\s+', self.tmp['this_k_energies'].strip()):
            ek_split += [unit_conversion.convert_unit(valueForStrValue(energy, 'f'),'eV')]
        self.tmp['k_energies'].append(ek_split)
        self.tmp['k_point'].append([section['x_qe_t_k_x'][0], section['x_qe_t_k_y'][0], section['x_qe_t_k_z'][0]])

    def appendToTmp(self, tmpname, value):
        self.tmp[tmpname] += value

    def adHoc_alat(self, parser):
        alat = parser.lastMatch['x_qe_alat']
        unit_conversion.register_userdefined_quantity('usrAlat', 'm', alat)
        unit_conversion.register_userdefined_quantity('usrTpiba', '1/m', 2*math.pi/alat)

    def run_submatchers(self):
        """submatchers of section_run"""
        return [
            SM(name='header',
               startReStr=r".*(?:\s|\()lmax(?:x\))?\s*=",
               sections = ['section_basis_set_cell_dependent', 'section_method', 'section_system'],
               subMatchers=[
                   SM(name='enforced_XC',
                      startReStr=r"\s*IMPORTANT: XC functional enforced from input",
                      subMatchers=[
                          SM(name='xc_functional_enforced', required=True,
                             startReStr=r"\s*Exchange-correlation\s*=\s*(?P<x_qe_t_xc_functional_shortname_enforced>\S+)"
                          ),
                      ],
                   ),
                   SM(name='alat', required=True,
                      startReStr=r"\s*lattice parameter \((?:a_0|alat)\)\s*=\s*(?P<x_qe_alat__bohr>" + RE_f + r")\s*a\.u\.",
                      adHoc=self.adHoc_alat,
                   ),
                   SM(name='nat', required=True,
                      startReStr=r"\s*number of atoms/cell\s*=\s*(?P<number_of_atoms>\d+)",
                   ),
                   SM(name='nsp', required=True,
                      startReStr=r"\s*number of atomic types\s*=\s*(?P<x_qe_number_of_species>" + RE_i + r")",
                   ),
                   SM(name='nbnd', required=True,
                       startReStr=r"\s*number of Kohn-Sham states\s*=\s*(?P<x_qe_number_of_states>" + RE_i + r")"
                   ),
                   SM(name='ecutwfc', required=True,
                      startReStr=r"\s*kinetic-energy cutoff\s*=\s*(?P<basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry"
                   ),
                   SM(name='ecut_density', required=True,
                      startReStr=r"\s*charge density cutoff\s*=\s*(?P<x_qe_density_basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry"
                   ),
                   SM(name='xc_functional', required=True,
                      startReStr=r"\s*Exchange-correlation\s*=\s*(?P<x_qe_xc_functional_shortname>\S+)\s*\((?P<x_qe_t_xc_functional_num>[^\)]*)"
                   ),
                   SM(name='simulation_cell',
                      startReStr=r"\s*crystal axes: \(cart. coord.",
                      subMatchers=[
                          SM(name='cell_vec_a', repeats=True,
                             startReStr=r"\s*a\(\d\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_t_vec_a', 'usrAlat'),
                          ),
                      ],
                   ),
                   SM(name='reciprocal_cell',
                      startReStr=r"\s*reciprocal axes: \(cart. coord. in units 2 pi/(?:alat|a_0)\)",
                      subMatchers=[
                          SM(name='cell_vec_b', repeats=True,
                             startReStr=r"\s*b\(\d\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_t_vec_b', 'usrTpiba'),
                          ),
                      ],
                   ),
                   SM(name='pseudopotential', repeats=True,
                      startReStr=(r"\s*PseudoPot\.\s*#\s*(?P<x_qe_t_pp_idx>" + RE_i + r") for (?P<x_qe_t_pp_label>\S+) read from file" +
                                  r"(?::|\s*(?P<x_qe_t_pp_filename>\S+))\s*"),
                      sections=['x_qe_t_section_pseudopotential'],
                      subMatchers=[
                          SM(name='new_pp_filename',
                             startReStr=r"^\s*(?P<x_qe_t_pp_filename>\S+)\s*$",
                          ),
                          SM(name='pp_md5',
                             startReStr=r"\s*MD5 check sum:\s*(?P<x_qe_t_pp_md5sum>\S+)",
                          ),
                          SM(name='pp_type_val',
                             startReStr=r"\s*Pseudo is\s*(?P<x_qe_t_pp_type>.*?),\s*Zval\s*=\s*(?P<x_qe_t_pp_valence>" + RE_f + r")",
                          ),
                      ],
                   ),
                   SM(name='pp_atom_kind_map',
                      startReStr=r"\s*atomic species\s+valence\s+mass\s+pseudopotential",
                      subMatchers=[
                          SM(name='atom_kind', repeats=True,
                             startReStr=r"\s*(?P<method_atom_kind_label>\S+)\s+(?P<x_qe_pp_valence>" + RE_f + r")" +
                                        r"\s+(?P<x_qe_kind_mass>" + RE_f + r")\s+(?P<x_qe_pp_label>[^\s\(]+)" +
                                        r"\(\s*(?P<x_qe_pp_weight>" + RE_f + r")\s*\)",
                             sections=['section_method_atom_kind'],
                          ),
                      ],
                   ),
                   SM(name='nsymm',
                      startReStr=r"\s*(?P<x_qe_nsymm>\d+)\s*Sym\.\s*Ops\.",
                   ),
                   SM(name='atom_pos_cart_list',
                      startReStr=r"\s*site n.     atom                  positions \((?:a_0|alat) units\)",
                      subMatchers=[
                          SM(name='atom_pos_cart', repeats=True,
                             startReStr=(
                                 r"\s*(?P<x_qe_t_atom_idx>" + RE_i + r")" +
                                 r"\s+(?P<x_qe_t_atom_labels>\S+)\s+tau\(\s*" + RE_i + "\)\s*"
                                 r"=\s*\(\s*" + QeC.re_vec('x_qe_t_atpos', 'usrAlat')),
                          ),
                      ],
                   ),
               ],
            ), # header
            SM(name='self_consistent_calculation',
               startReStr=r"\s*Self-consistent Calculation",
               sections = ['section_single_configuration_calculation'],
               subMatchers=[
                   SM(name='iteration', repeats=True,
                      startReStr=r'\s*iteration\s*#',
                      sections=['section_scf_iteration'],
                      subMatchers=[
                          SM(name='e_total',
                             startReStr=r"\s*!?\s*total\s+energy\s*=\s*(?P<energy_total_scf_iteration>" + RE_f + r")",
                          ),
                      ],
                   ),
                   SM(name='scf_result', repeats=False,
                       startReStr=r'\s*End of self-consistent calculation',
                       sections=['section_eigenvalues'],
                       subMatchers=[
                          SM(name='bands', repeats=True,
                              sections=['x_qe_t_section_kbands'],
                              startReStr=r'\s*k\s*=\s*' + QeC.re_vec('x_qe_t_k', 'usrTpiba') + r'\s*\(\s*(?P<x_qe_t_k_pw>' + RE_i + ')\s*PWs',
                              subMatchers=[
                                  SM(name='kbnd', repeats=True,
                                      startReStr=r'\s*(?P<x_qe_t_k_point_energies>(?:\s*' + RE_f + ')+\s*$)',
                                      adHoc=lambda p: self.appendToTmp('this_k_energies', " " + p.lastMatch['x_qe_t_k_point_energies']),
                                  ),
                              ],
                          ),
                          SM(name='e_total',
                             startReStr=r'\s*!?\s*total\s+energy\s*=\s*(?P<energy_total>' + RE_f + ')',
                          ),
                       ],
                   ),
               ],
            ),
        ]

if __name__ == "__main__":
    parser = QuantumEspressoParserPWSCF()
    parser.parse()
