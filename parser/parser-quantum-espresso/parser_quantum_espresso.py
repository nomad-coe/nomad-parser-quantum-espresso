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


LOGGER = logging.getLogger(__name__)

RE_FORTRAN_FLOAT = r"[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?"
RE_FORTRAN_INT = r"[+-]?\d+"


def adHoc_alat(parser):
    line = parser.fIn.readline()
    match = re.match(r"[^=]+=\s*(" + RE_FORTRAN_FLOAT + r")\s*a\.u\.", line)
    if match:
        alat_au = float(match.group(1))
    else:
        raise RuntimeError("should not happen: %s", line)
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
                r"\s*Program\s+(?P<program_name>\S+)\s+v\." +
                r"(?P<program_version>\S+(?:\s+\(svn\s+rev\.\s+\d+\s*\))?)" +
                r"\s+starts" +
                r"(?:(?:\s+on\s+(?P<x_qe_time_run_date_start>.+?)?)\s*$|"
                r"(?:\s*\.\.\.)\s*$)"),
            repeats=True,
            required=True,
            forwardMatch=False,
            fixedStartValues={'program_name': 'Quantum Espresso',
                              'program_basis_set_type': 'plane waves'},
            sections=['section_run'],
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
                        r"\s*lattice parameter \(a_0\)\s*=\s*" +
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
