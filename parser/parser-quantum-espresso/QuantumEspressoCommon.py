import setup_paths
import calendar
import json
import os
import re
import numpy as np
import logging
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM, CachingLevel

############################################################
# This file contains functions and constants that are needed
# by more than one parser.
############################################################


LOGGER = logging.getLogger(__name__)


def re_vec(name, units='', split="\s*"):
    """generator for 3-component vector regex"""
    if units:
        units = '__' + units
    res = (
        r'(?P<' + name + r'_x' + units + '>\S+)' + split +
        r'(?P<' + name + r'_y' + units + '>\S+)' + split +
        r'(?P<' + name + r'_z' + units + '>\S+)'
        )
    return res


# loading metadata from
# nomad-meta-info/meta_info/nomad_meta_info/quantum_espresso.nomadmetainfo.json
META_INFO = loadJsonFile(
    filePath=os.path.normpath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../../../nomad-meta-info/meta_info/nomad_meta_info/quantum_espresso.nomadmetainfo.json")),
    dependencyLoader=None,
    extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
    uri=None)[0]

PARSER_INFO_DEFAULT = {
  "name": "parser_quantum_espresso",
  "version": "0.0.1"
}

# constants for date conversion
MONTHS = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
           'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
MONTH_NUMBER = { MONTHS[num]: num for num in range(0,12) }

class ParserQuantumEspresso(object):
    """Base class for all Quantum Espresso parsers"""
    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None):
        self.qe_program_name = None
        self.parserInfo = PARSER_INFO_DEFAULT.copy()
        self.cachingLevelForMetaName = {}
        for name in META_INFO.infoKinds:
            # set all temporaries to caching-only
            if name.startswith('x_qe_t_'):
                self.cachingLevelForMetaName[name] = CachingLevel.Cache
        # common prosa in espresso output
        self.coverageIgnoreList = [
            r"\s*",
            r"\s*Ultrasoft \(Vanderbilt\) Pseudopotentials\s*",
            r"\s*This program is part of the open-source Quantum ESPRESSO suite",
            r"\s*for quantum simulation of materials; please cite",
            r"\s*\"P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 \(2009\);\s*",
            r"\s*URL http://www.quantum-espresso.org\",\s*",
            r"\s*in publications or presentations arising from this work. More details at",
            r"\s*http://www.quantum-espresso.org/quote",
            r"\s*Current dimensions of program \S+ are:",
        ]
        self.coverageIgnore = None

    def parse(self):
        self.coverageIgnore = re.compile(r"^(?:" + r"|".join(self.coverageIgnoreList) + r")$")
        mainFunction(self.mainFileDescription(), META_INFO, self.parserInfo,
                    cachingLevelForMetaName=self.cachingLevelForMetaName,
                    superContext=self)

    def mainFileDescription(self):
        # assemble matchers and submatchers
        result = SM(
            name='root',
            weak=True,
            startReStr="",
            subMatchers=[
                SM(name='newRun', repeats=True, required=True,
                   startReStr=(
                       # program name, e.g. PWSCF, DOS
                       r"\s*Program\s+(?P<x_qe_program_name>\S+)\s+v\." +
                       # version
                       r"(?P<program_version>\S+(?:\s+\(svn\s+rev\.\s+\d+\s*\))?)" +
                       r"\s+starts" +
                       # newer espresso: "on $date"
                       r"(?:(?:\s+on\s+(?P<x_qe_time_run_date_start>.+?)?)\s*$|" +
                       # older espresso has just "..." and date on new line
                       r"(?:\s*\.\.\.)\s*$)"
                   ),
                   fixedStartValues={'program_name': 'Quantum Espresso',
                                     'program_basis_set_type': 'plane waves'},
                   sections=['section_run'],
                   subMatchers=([
                       # older espresso versions have start date on separate line
                       SM(name='run_date',
                          startReStr=r"\s*Today is\s*(?P<x_qe_time_run_date_start>.+?)\s*$"
                       ),
                   ] + self.run_submatchers()),
                )
            ]
        )
        return result

    def run_submatchers(self):
        return []

    def strValueTransform_strQeDate(self, espresso_date):
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
    strValueTransform_strQeDate.units = 's'

    def onClose_section_run(self, backend, gIndex, section):
        LOGGER.info("closing section run")
        string_run_start = section['x_qe_time_run_date_start'][-1]
        epoch = self.strValueTransform_strQeDate(string_run_start)
        backend.addValue('time_run_date_start', epoch)
