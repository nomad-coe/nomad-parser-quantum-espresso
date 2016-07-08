import setup_paths
from nomadcore.simple_parser import mainFunction, SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import os
import sys
import json

# description of the input
mainFileDescription = SM(
    name='root',
    weak=True,
    startReStr="",
    subMatchers=[
        SM(
            name='newRun',
            startReStr=r"\s*# SampleParser #\s*",
            repeats=True,
            required=True,
            forwardMatch=True,
            fixedStartValues={'program_name': 'Quantum Espresso',
                              'program_basis_set_type': 'plane waves'},
            sections=['section_run'],
            subMatchers=[
                SM(
                    name='header',
                    startReStr=r"\s*# SampleParser #\s*"
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
        logging.getLogger("nomadcore.parsing").info(
            "closing section_single_configuration_calculation gIndex %d %s",
            gIndex, section.simpleValues)
        self.scfIterNr = 0

    def onClose_section_scf_iteration(self, backend, gIndex, section):
        """trigger called when section_scf_iteration is closed"""
        logging.getLogger("nomadcore.parsing").info(
            "closing section_scf_iteration bla gIndex %d %s",
            gIndex, section.simpleValues)
        self.scfIterNr += 1

# which values to cache or forward (mapping meta name -> CachingLevel)
cachingLevelForMetaName = {}

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo,
                 cachingLevelForMetaName=cachingLevelForMetaName,
                 superContext=QuantumEspressoParserContext())
