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


parserInfo = {
  "name": "parser_quantum_espresso",
  "version": "0.0.1"
}


def adHoc_alat(parser):
    line = parser.fIn.readline()
    match = re.match(r"[^=]+=\s*(" + QeC.RE_FORTRAN_FLOAT + r")\s*a\.u\.",
                     line)
    if match:
        alat_au = float(match.group(1))
    else:
        raise RuntimeError("This should never happen: %s", line)
    unit_conversion.register_userdefined_quantity('usrAlat', 'bohr', alat_au)
    unit_conversion.register_userdefined_quantity(
        'usrTpiba', '1/bohr', 2*math.pi/alat_au)


# description of the input
mainFileDescription = SM(
    name='root',
    weak=True,
    startReStr="",
    subMatchers=[
        SM(
            name='newRun',
            startReStr=QeC.RE_RUN_START,
            repeats=True,
            required=True,
            forwardMatch=False,
            fixedStartValues={'program_name': 'Quantum Espresso',
                              'program_basis_set_type': 'plane waves'},
            sections=['section_run', 'section_basis_set_cell_dependent', 'section_system'],
            subMatchers=[
                SM(
                    name='run_date',
                    startReStr=QeC.RE_RUN_DATE,
                    repeats=False,
                    required=False,
                    forwardMatch=False,
                ),
                SM(
                    name='alat',
                    startReStr=(
                        r"\s*lattice parameter \((?:a_0|alat)\)\s*=\s*" +
                        r"(?P<x_qe_alat__bohr>" +
                        QeC.RE_FORTRAN_FLOAT + r")\s*a\.u\."
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=True,
                    adHoc=adHoc_alat,
                ),
                SM(
                    name='nat',
                    startReStr=(
                        r"\s*number of atoms/cell\s*=\s*" +
                        r"(?P<x_qe_nat>" +
                        QeC.RE_FORTRAN_INT + r")\s*"
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=False,
                ),
                SM(
                    name='nsp',
                    startReStr=(
                        r"\s*number of atomic types\s*=\s*" +
                        r"(?P<x_qe_nsp>" +
                        QeC.RE_FORTRAN_INT + r")\s*"
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=False,
                ),
                SM(
                    name='nbnd',
                    startReStr=(
                        r"\s*number of Kohn-Sham states\s*=\s*" +
                        r"(?P<x_qe_nbnd>" +
                        QeC.RE_FORTRAN_INT + r")\s*"
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=False,
                ),
                SM(
                    name='ecutwfc',
                    startReStr=(
                        r"\s*kinetic-energy cutoff\s*=\s*" +
                        r"(?P<basis_set_planewave_cutoff__rydberg>" +
                        QeC.RE_FORTRAN_FLOAT + r")\s*Ry\s*"
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=False,
                ),
                SM(
                    name='ecut_density',
                    startReStr=(
                        r"\s*charge density cutoff\s*=\s*" +
                        r"(?P<x_qe_density_basis_set_planewave_cutoff__rydberg>" +
                        QeC.RE_FORTRAN_FLOAT + r")\s*Ry\s*"
                    ),
                    repeats=False,
                    required=True,
                    forwardMatch=False,
                ),
                SM(
                    name='simulation_cell',
                    startReStr=r"\s*crystal axes: \(cart. coord.",
                    endReStr=r"^\s*$", # empty line ends bravais matrix
                    subMatchers=[
                        SM(
                            name='cell_a1',
                            startReStr=r"\s*a\(1\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_a1', 'usrAlat'),
                        ),
                        SM(
                            name='cell_a2',
                            startReStr=r"\s*a\(2\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_a2', 'usrAlat'),
                        ),
                        SM(
                            name='cell_a3',
                            startReStr=r"\s*a\(3\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_a3', 'usrAlat'),
                        ),
                    ],
                ),
                SM(
                    name='reciprocal_cell',
                    startReStr=r"\s*reciprocal axes: \(cart. coord. in units 2 pi/(?:alat|a_0)\)",
                    endReStr=r"^\s*$", # empty line ends reciprocal matrix
                    subMatchers=[
                        SM(
                            name='cell_b1',
                            startReStr=r"\s*b\(1\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_b1', 'usrTpiba'),
                        ),
                        SM(
                            name='cell_b2',
                            startReStr=r"\s*b\(2\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_b2', 'usrTpiba'),
                        ),
                        SM(
                            name='cell_b3',
                            startReStr=r"\s*b\(3\)\s*=\s*\(\s*" +
                                QeC.re_vec('x_qe_b3', 'usrTpiba'),
                        ),
                    ],
                ),
                SM(
                    name='pseudopotentials',
                    startReStr=r"\s*PseudoPot\.\s*#\s*1",
                    endReStr=r"\s*\d+\s*Sym\.Ops\.",
                    # subMatchers=[
                    #     SM(
                    #         name='pseudopotential',
                    #         startReStr=r"\s*PseudoPot\.\s*#\s*\d+",
                    #         adHoc=adHoc
                    forwardMatch=True,
                ),
                SM(
                    name='nsymm',
                    startReStr=r"\s*(?P<x_qe_nsymm>\d+)\s*Sym\.\s*Ops\.",
                    endReStr=r"\s*Cartesian Axes",
                ),
            ],
        )
    ],
)


class QuantumEspressoParserContext(object):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self):
        self.scfIterNr = 0
        self.coverageIgnore = QeC.RE_COVERAGE_IGNORE

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse
        different files"""
        pass

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # allows to reset values if the same superContext is used to parse
        # different files
        self.initialize_values()

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

    def onClose_section_run(self, backend, gIndex, section):
        LOGGER.info("closing section run")
        string_run_start = section['x_qe_time_run_date_start'][-1]
        epoch = QeC.espresso_date_to_epoch(string_run_start)
        backend.addValue('time_run_date_start', epoch)

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



# which values to cache or forward (mapping meta name -> CachingLevel)
cachingLevelForMetaName = {}

if __name__ == "__main__":
    mainFunction(mainFileDescription, QeC.META_INFO_ENV, parserInfo,
                 cachingLevelForMetaName=cachingLevelForMetaName,
                 superContext=QuantumEspressoParserContext())
