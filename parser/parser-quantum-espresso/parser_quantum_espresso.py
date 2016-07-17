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


LOGGER = logging.getLogger(__name__)


class QuantumEspressoParserPWSCF(QeC.ParserQuantumEspresso):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self):
        QeC.ParserQuantumEspresso.__init__(self)
        self.scfIterNr = 0

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse
        different files"""
        self.secMethodIndex = None
        self.secSystemDescriptionIndex = None

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # allows to reset values if the same superContext is used to parse
        # different files
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

    def onOpen_section_system(
            self, backend, gIndex, section):
        self.secSystemIndex = gIndex

    # just examples, you probably want to remove the following two triggers

    def onClose_section_single_configuration_calculation(
            self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""
        # backend.addValue("", self.scfIterNr)
        LOGGER.info(
            "closing section_single_configuration_calculation gIndex %d %s",
            gIndex, section.simpleValues)
        self.scfIterNr = 0
        backend.addValue('single_configuration_to_calculation_method_ref',
                         self.secMethodIndex)
        backend.addValue('single_configuration_calculation_to_system_ref',
                         self.secSystemIndex)

    def onClose_section_scf_iteration(self, backend, gIndex, section):
        """trigger called when section_scf_iteration is closed"""
        LOGGER.info(
            "closing section_scf_iteration bla gIndex %d %s",
            gIndex, section.simpleValues)
        self.scfIterNr += 1

    def onClose_section_system(self, backend, gIndex, section):
        backend.addArrayValues('simulation_cell', np.array([
            [section['x_qe_a1_x'], section['x_qe_a1_y'], section['x_qe_a1_z']],
            [section['x_qe_a2_x'], section['x_qe_a2_y'], section['x_qe_a2_z']],
            [section['x_qe_a3_x'], section['x_qe_a3_y'], section['x_qe_a3_z']],
        ]))
        backend.addArrayValues('x_qe_reciprocal_cell', np.array([
            [section['x_qe_b1_x'], section['x_qe_b1_y'], section['x_qe_b1_z']],
            [section['x_qe_b2_x'], section['x_qe_b2_y'], section['x_qe_b2_z']],
            [section['x_qe_b3_x'], section['x_qe_b3_y'], section['x_qe_b3_z']],
        ]))

    def adHoc_alat(self, parser):
        line = parser.fIn.readline()
        match = re.match(r"[^=]+=\s*(\S+)\s*a\.u\.", line)
        if match:
            alat_au = float(match.group(1))
        else:
            raise RuntimeError("This should never happen: %s", line)
        unit_conversion.register_userdefined_quantity('usrAlat', 'bohr', alat_au)
        unit_conversion.register_userdefined_quantity(
            'usrTpiba', '1/bohr', 2*math.pi/alat_au)

    def run_submatchers(self):
        return [
            SM(name='header',
               startReStr=r".*(?:\s|\()lmax(?:x\))?\s*=",
               sections = ['section_basis_set_cell_dependent', 'section_method', 'section_system'],
               subMatchers=[
                   SM(name='alat', required=True, forwardMatch=True,
                      startReStr=r"\s*lattice parameter \((?:a_0|alat)\)\s*=\s*(?P<x_qe_alat__bohr>\S+)\s*a\.u\.",
                      adHoc=self.adHoc_alat,
                   ),
                   SM(name='nat', required=True,
                      startReStr=r"\s*number of atoms/cell\s*=\s*(?P<x_qe_nat>\S+)",
                   ),
                   SM(name='nsp', required=True,
                      startReStr=r"\s*number of atomic types\s*=\s*(?P<x_qe_nsp>\S+)",
                   ),
                   SM(name='nbnd', required=True,
                       startReStr=r"\s*number of Kohn-Sham states\s*=\s*(?P<x_qe_nbnd>\S+)"
                   ),
                   SM(name='ecutwfc', required=True,
                      startReStr=r"\s*kinetic-energy cutoff\s*=\s*(?P<basis_set_planewave_cutoff__rydberg>\S+)\s*Ry"
                   ),
                   SM(name='ecut_density', required=True,
                      startReStr=r"\s*charge density cutoff\s*=\s*(?P<x_qe_density_basis_set_planewave_cutoff__rydberg>\S+)\s*Ry"
                   ),
                   SM(name='simulation_cell',
                      startReStr=r"\s*crystal axes: \(cart. coord.",
                      subMatchers=[
                          SM(name='cell_a1',
                             startReStr=r"\s*a\(1\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_a1', 'usrAlat'),
                          ),
                          SM(name='cell_a2',
                             startReStr=r"\s*a\(2\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_a2', 'usrAlat'),
                          ),
                          SM(name='cell_a3',
                             startReStr=r"\s*a\(3\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_a3', 'usrAlat'),
                          ),
                      ],
                   ),
                   SM(name='reciprocal_cell',
                      startReStr=r"\s*reciprocal axes: \(cart. coord. in units 2 pi/(?:alat|a_0)\)",
                      subMatchers=[
                          SM(name='cell_b1',
                             startReStr=r"\s*b\(1\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_b1', 'usrTpiba'),
                          ),
                          SM(name='cell_b2',
                             startReStr=r"\s*b\(2\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_b2', 'usrTpiba'),
                          ),
                          SM(name='cell_b3',
                             startReStr=r"\s*b\(3\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_b3', 'usrTpiba'),
                          ),
                      ],
                   ),
                   SM(name='pseudopotentials', forwardMatch=True,
                      startReStr=r"\s*PseudoPot\.\s*#\s*1",
                      # subMatchers=[
                      #     SM(
                      #         name='pseudopotential',
                      #         startReStr=r"\s*PseudoPot\.\s*#\s*\d+",
                      #         adHoc=adHoc
                   ),
                   SM(name='nsymm',
                      startReStr=r"\s*(?P<x_qe_nsymm>\d+)\s*Sym\.\s*Ops\.",
                   ),
               ],
            ), # header
            SM(name='self_consistent_calculation',
               startReStr=r"\s*Self-consistent Calculation",
               sections = ['section_single_configuration_calculation'],
               subMatchers=[
                   SM(name='iteration',
                      startReStr=r'\s*iteration\s*#',
                      sections=['section_scf_iteration'],
                      repeats=True,
                      subMatchers=[
                          SM(name='bands', repeats=True,
                              startReStr=r'\s*k\s*=(?:\s+(?P<x_qe_k>[^ \(]+))+',
                          ),
                          SM(name='e_total',
                             startReStr=r'\s*!?\s*total\s+energy\s*=\s*(?P<energy_total_scf_iteration>\S+)',
                          ),
                      ],
                   ),
               ],
            ),
        ]

if __name__ == "__main__":
    parser = QuantumEspressoParserPWSCF()
    parser.main()
