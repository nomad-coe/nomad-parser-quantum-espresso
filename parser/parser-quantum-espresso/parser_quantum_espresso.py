import setup_paths
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import os
import sys
import json
import re
import logging
import calendar
import nomadcore.unit_conversion.unit_conversion as unit_conversion
import math
import numpy as np


LOGGER = logging.getLogger(__name__)

RE_FORTRAN_FLOAT = r"[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?"
RE_FORTRAN_INT = r"[+-]?\d+"


def re_vec(name, units=''):
    """generator for 3-component vector (space separated) regex"""
    if units:
        units = '__' + units
    res = (r'(?P<' + name + r'_x' + units + '>' + RE_FORTRAN_FLOAT + r')\s*' +
           r'(?P<' + name + r'_y' + units + '>' + RE_FORTRAN_FLOAT + r')\s*' +
           r'(?P<' + name + r'_z' + units + '>' + RE_FORTRAN_FLOAT + r')')
    return res


def adHoc_alat(parser):
    line = parser.fIn.readline()
    match = re.match(r"[^=]+=\s*(" + RE_FORTRAN_FLOAT + r")\s*a\.u\.", line)
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
            startReStr=(
                r"\s*Program\s+(?P<x_qe_program_name>\S+)\s+v\." +
                r"(?P<program_version>\S+(?:\s+\(svn\s+rev\.\s+\d+\s*\))?)" +
                r"\s+starts" +
                r"(?:(?:\s+on\s+(?P<x_qe_time_run_date_start>.+?)?)\s*$|"
                r"(?:\s*\.\.\.)\s*$)"),
            repeats=True,
            required=True,
            forwardMatch=False,
            fixedStartValues={'program_name': 'Quantum Espresso',
                              'program_basis_set_type': 'plane waves'},
            sections=['section_run', 'section_basis_set_cell_dependent', 'section_system'],
            subMatchers=[
                SM(
                    name='run_date',
                    startReStr=(
                        r"\s*Today is\s*(?P<x_qe_time_run_date_start>.+?)\s*$"
                    ),
                    repeats=False,
                    required=False,
                    forwardMatch=False,
                ),
                SM(
                    name='alat',
                    startReStr=(
                        r"\s*lattice parameter \((?:a_0|alat)\)\s*=\s*" +
                        r"(?P<x_qe_alat__bohr>" +
                        RE_FORTRAN_FLOAT + r")\s*a\.u\."
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
                        RE_FORTRAN_INT + r")\s*"
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
                        RE_FORTRAN_INT + r")\s*"
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
                        RE_FORTRAN_INT + r")\s*"
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
                        RE_FORTRAN_FLOAT + r")\s*Ry\s*"
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
                        RE_FORTRAN_FLOAT + r")\s*Ry\s*"
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
                                re_vec('x_qe_a1', 'usrAlat'),
                        ),
                        SM(
                            name='cell_a2',
                            startReStr=r"\s*a\(2\)\s*=\s*\(\s*" +
                                re_vec('x_qe_a2', 'usrAlat'),
                        ),
                        SM(
                            name='cell_a3',
                            startReStr=r"\s*a\(3\)\s*=\s*\(\s*" +
                                re_vec('x_qe_a3', 'usrAlat'),
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
                                re_vec('x_qe_b1', 'usrTpiba'),
                        ),
                        SM(
                            name='cell_b2',
                            startReStr=r"\s*b\(2\)\s*=\s*\(\s*" +
                                re_vec('x_qe_b2', 'usrTpiba'),
                        ),
                        SM(
                            name='cell_b3',
                            startReStr=r"\s*b\(3\)\s*=\s*\(\s*" +
                                re_vec('x_qe_b3', 'usrTpiba'),
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

# loading metadata from
# nomad-meta-info/meta_info/nomad_meta_info/quantum_espresso.nomadmetainfo.json
metaInfoPath = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../../../nomad-meta-info/meta_info/nomad_meta_info/" +
        "quantum_espresso.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(
    filePath=metaInfoPath, dependencyLoader=None,
    extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS, uri=None)

parserInfo = {
  "name": "parser_quantum_espresso",
  "version": "1.0"
}


class QuantumEspressoParserContext(object):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self):
        self.scfIterNr = 0
        coverageIgnoreList = [
           r"\s*",
           r"\s*Ultrasoft \(Vanderbilt\) Pseudopotentials\s*",
           r"\s*This program is part of the open-source Quantum ESPRESSO suite",
           r"\s*for quantum simulation of materials; please cite",
           r"\s*\"P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 \(2009\);\s*",
           r"\s*URL http://www.quantum-espresso.org\",\s*",
           r"\s*in publications or presentations arising from this work. More details at",
           r"\s*http://www.quantum-espresso.org/quote",
        ]
        self.coverageIgnore = re.compile(r"^(?:" +
            r"|".join(coverageIgnoreList) + r")$")

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
        epoch = espresso_date_to_epoch(string_run_start)
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

MONTHS = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
           'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
MONTH_NUMBER = { MONTHS[num]: num for num in range(0,12) }

def espresso_date_to_epoch(espresso_date):
    match = re.match(
        r"(\d+)\s*([A-Za-z]+)\s*(\d+)\s+at\s+(\d+):\s*(\d+):\s*(\d+)",
        espresso_date)
    if match:
        month = MONTH_NUMBER[match.group(2)]
        epoch = calendar.timegm(
            (int(match.group(3)), int(month), int(match.group(1)),
             int(match.group(4)), int(match.group(5)), int(match.group(6))))
        return(epoch)
    else:
        raise RuntimeError("unparsable date: %s", espresso_date)


# which values to cache or forward (mapping meta name -> CachingLevel)
cachingLevelForMetaName = {}

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo,
                 cachingLevelForMetaName=cachingLevelForMetaName,
                 superContext=QuantumEspressoParserContext())
