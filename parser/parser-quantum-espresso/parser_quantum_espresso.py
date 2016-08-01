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
from QuantumEspressoXC import translate_qe_xc_num
from nomadcore.parser_backend import valueForStrValue


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
        self.alat = None

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser
        # reset values if same superContext is used to parse different files
        self.initialize_values()

    def onClose_section_basis_set_cell_dependent(
            self, backend, gIndex, section):
        if section['basis_set_planewave_cutoff'] is None:
            LOGGER.error("basis_set_planewave_cutoff is not set")
            return
        pwc_ry = unit_conversion.convert_unit(
            section['basis_set_planewave_cutoff'][-1],
            'J', 'rydberg')
        backend.addValue('basis_set_cell_dependent_name', 'PW_%.1f' % (pwc_ry))
        backend.addValue('basis_set_cell_dependent_kind', 'plane_waves')

    def onOpen_section_method(
            self, backend, gIndex, section):
        self.secMethodIndex = gIndex
        self.cache_t_pseudopotential = {}
        self.cache_t_pp_report = {}
        self.cache_t_pp_renorm_wfc = {}

    def onOpen_section_system(
            self, backend, gIndex, section):
        self.secSystemIndex = gIndex

    def onClose_x_qe_t_section_pseudopotential(
            self, backend, gIndex, section):
        if section['x_qe_t_pp_label'] is None:
            LOGGER.error("x_qe_t_pp_label is not set")
            return
        self.cache_t_pseudopotential[section['x_qe_t_pp_label'][0]] = section

    def onClose_section_method_atom_kind(
            self, backend, gIndex, section):
        pp_label = section['x_qe_pp_label'][0]
        pp = self.cache_t_pseudopotential[pp_label]
        for key, value in sorted(pp.simpleValues.items()):
            if key == 'x_qe_t_pp_label' or key == 'x_qe_pp_valence' or value is None:
                continue
            target = re.sub(r'^x_qe_t_',r'x_qe_',key)
            if target == key:
                raise Exception('found non-temporary key in pseudopotential cache: "%s"' % (key))
            backend.addValue(target, value[-1])
        if pp['x_qe_t_pp_idx'] is not None:
            pp_num = pp['x_qe_t_pp_idx'][-1]
            pp_report = self.cache_t_pp_report.get(pp_num, None)
            if pp_report is not None:
                backend.addValue('x_qe_pp_report_version', pp_report['x_qe_t_pp_report_version'][-1])
                backend.addValue('x_qe_pp_report_contents', "\n".join(pp_report['x_qe_t_pp_report_line']))
        else:
            LOGGER.error("x_qe_t_pp_idx is not set")
        if pp['x_qe_t_pp_filename'] is not None:
            fpath = pp['x_qe_t_pp_filename'][-1]
            basename = re.sub(r".*\/", r"", fpath)
            renorm_info = self.cache_t_pp_renorm_wfc.get(basename, None)
            if renorm_info is not None:
                backend.addValue('x_qe_pp_renormalized_wfc', renorm_info)

    def onClose_x_qe_t_section_pp_report(
            self, backend, gIndex, section):
        self.cache_t_pp_report[section['x_qe_t_pp_report_species'][0]] = section

    def onClose_x_qe_t_section_input_occupations(
            self, backend, gIndex, section):
        if section['x_qe_t_input_occupations'] is not None:
            iocc_joined = ' '.join(section['x_qe_t_input_occupations'])
            iocc_split = re.findall(RE_f, iocc_joined)
            LOGGER.error("implement proper output of x_qe_t_section_input_occupations")
        else:
            LOGGER.error("closing x_qe_t_section_input_occupations, but no input occupations parsed")

    def onClose_section_method(
            self, backend, gIndex, section):
        # set flag if we deal with user-enforced XC functional
        if section['x_qe_t_xc_functional_shortname_enforced'] is not None:
            backend.addValue('x_qe_xc_functional_user_enforced', True)
        # translate XC functional to section_xc_functionals
        xc_functionals = None
        if section['x_qe_xc_functional_num'] is not None:
            xc_functional_num = section['x_qe_xc_functional_num'][-1]
            xc_functionals = translate_qe_xc_num(xc_functional_num, section['x_qe_t_exact_exchange_fraction'])
        else:
            LOGGER.error("x_qe_xc_functional_num is not set")
        if xc_functionals is not None:
            for xc_functional in xc_functionals:
                self.addSectionDict(backend, 'section_XC_functionals', xc_functional)
        else:
            LOGGER.error("error getting xc_functionals")
        if section['x_qe_t_allocated_array_name'] is not None:
            backend.addArrayValues('x_qe_allocated_array_name', np.asarray(section['x_qe_t_allocated_array_name']))
        if section['x_qe_t_allocated_array_size'] is not None:
            backend.addArrayValues('x_qe_allocated_array_size', np.asarray(section['x_qe_t_allocated_array_size']))
        if section['x_qe_t_allocated_array_dimensions'] is not None:
            backend.addArrayValues('x_qe_allocated_array_dimensions', np.asarray(section['x_qe_t_allocated_array_dimensions']))
        if section['x_qe_t_temporary_array_name'] is not None:
            backend.addArrayValues('x_qe_temporary_array_name', np.asarray(section['x_qe_t_temporary_array_name']))
        if section['x_qe_t_temporary_array_size'] is not None:
            backend.addArrayValues('x_qe_temporary_array_size', np.asarray(section['x_qe_t_temporary_array_size']))
        if section['x_qe_t_temporary_array_dimensions'] is not None:
            backend.addArrayValues('x_qe_temporary_array_dimensions', np.asarray(section['x_qe_t_temporary_array_dimensions']))

    def onClose_section_single_configuration_calculation(
            self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""
        backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodIndex)
        backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemIndex)
        if section['x_qe_t_energy_decomposition_name'] is not None:
            backend.addArrayValues('x_qe_energy_decomposition_name', np.asarray(
                section['x_qe_t_energy_decomposition_name']))
        if section['x_qe_t_energy_decomposition_value'] is not None:
            backend.addArrayValues('x_qe_energy_decomposition_value', np.asarray(
                section['x_qe_t_energy_decomposition_value']))
        if section['x_qe_t_force_x'] is not None:
            backend.addArrayValues('atom_forces', np.array([
                section['x_qe_t_force_x'], section['x_qe_t_force_y'], section['x_qe_t_force_z']
                ]).T)
        if section['x_qe_t_stress_x'] is not None:
            backend.addArrayValues('stress_tensor', np.array([
                section['x_qe_t_stress_x'], section['x_qe_t_stress_y'], section['x_qe_t_stress_z']
                ]).T)

    def onClose_section_scf_iteration(
            self, backend, gIndex, section):
        """trigger called when section_scf_iteration is closed"""
        if section['x_qe_t_david_with_overlap'] is not None:
            backend.addValue('x_qe_diagonalization_scheme', 'davidson')
        elif section['x_qe_t_cg_diag']:
            backend.addValue('x_qe_diagonalization_scheme', 'conjugate_gradient')
        if section['x_qe_t_iteration_ethr'] is not None:
            backend.addValue('x_qe_iteration_ethr', section['x_qe_t_iteration_ethr'][-1])
        self.tmp['last_iteration'] = section['x_qe_iteration_number'][-1]

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
        if section['x_qe_t_vec_a_x'] is not None:
            self.amat = np.array([
                section['x_qe_t_vec_a_x'], section['x_qe_t_vec_a_y'], section['x_qe_t_vec_a_z'],
            ], dtype=np.float64).T
            # store inverse for transformation cartesian -> crystal
            try:
                self.amat_inv = np.linalg.inv(self.amat)
            except np.linalg.linalg.LinAlgError:
                raise Exception("error inverting bravais matrix " + str(self.amat))
            backend.addArrayValues('simulation_cell', self.amat)
        else:
            LOGGER.warning("No bravais matrix found in output")
        if section['x_qe_t_vec_b_x'] is not None:
            # store reciprocal lattice matrix for transformation crystal -> cartesian
            self.bmat = np.array([
                section['x_qe_t_vec_b_x'], section['x_qe_t_vec_b_y'], section['x_qe_t_vec_b_z'],
            ], dtype=np.float64).T
            # store inverse for transformation cartesian -> crystal
            try:
                self.bmat_inv = np.linalg.inv(self.bmat)
            except np.linalg.linalg.LinAlgError:
                raise Exception("error inverting reciprocal cell matrix")
            backend.addArrayValues('x_qe_reciprocal_cell', self.bmat)
        else:
            LOGGER.warning("No reciprocal cell matrix found in output")
        # atom positions
        if section['x_qe_t_atpos_x'] is not None:
            atpos_cart = np.array([
                section['x_qe_t_atpos_x'], section['x_qe_t_atpos_y'], section['x_qe_t_atpos_z']
            ], dtype=np.float64).T
            backend.addArrayValues('atom_positions',atpos_cart)
            backend.addArrayValues('atom_labels',np.asarray(section['x_qe_t_atom_labels']))
            backend.addArrayValues('x_qe_atom_idx',np.array(section['x_qe_t_atom_idx']))
        else:
            LOGGER.warning("No atom positions found in output")
        if section['x_qe_t_celldm'] is not None:
            celldm_joint = " ".join(section['x_qe_t_celldm'])
            celldm = [None, None, None, None, None, None]
            for match in re.findall(r"celldm\(\s*(\d+)\s*\)\s*=\s*(" + RE_f + r")", celldm_joint):
                celldm[int(match[0])-1] = valueForStrValue(match[1], 'f')
            celldm[0] = self.alat
            backend.addArrayValues('x_qe_celldm', np.array(celldm))
        else:
            LOGGER.warning("No QE cell dimensions found in output")
        if section['x_qe_t_k_info_vec_x'] is not None:
            backend.addArrayValues('x_qe_k_info_ik', np.array(section['x_qe_t_k_info_ik']))
            backend.addArrayValues('x_qe_k_info_wk', np.array(section['x_qe_t_k_info_wk']))
            backend.addArrayValues('x_qe_k_info_vec', np.array([
                section['x_qe_t_k_info_vec_x'], section['x_qe_t_k_info_vec_y'], section['x_qe_t_k_info_vec_z']
            ]).T)
        else:
            LOGGER.warning("No K-point info found in output")
        if section['x_qe_t_dense_FFT_grid_x'] is not None:
            backend.addArrayValues('x_qe_dense_FFT_grid', np.array([
                section['x_qe_t_dense_FFT_grid_x'], section['x_qe_t_dense_FFT_grid_y'], section['x_qe_t_dense_FFT_grid_z']
            ]).T)
        else:
            LOGGER.warning("No FFT grid info found in output")
        if section['x_qe_t_smooth_FFT_grid_x'] is not None:
            backend.addArrayValues('x_qe_smooth_FFT_grid', np.array([
                section['x_qe_t_smooth_FFT_grid_x'], section['x_qe_t_smooth_FFT_grid_y'], section['x_qe_t_smooth_FFT_grid_z']
            ]).T)
        if section['x_qe_t_vec_supercell_x'] is not None:
            backend.addArrayValues('x_qe_vec_supercell', np.array([
                section['x_qe_t_vec_supercell_x'], section['x_qe_t_vec_supercell_y'], section['x_qe_t_vec_supercell_z']
            ]).T)

    def onOpen_x_qe_t_section_kbands(self, backend, gIndex, section):
        self.tmp['this_k_energies'] = ''

    def onClose_x_qe_t_section_kbands(self, backend, gIndex, section):
        ek_split = []
        for energy in re.split(r'\s+', self.tmp['this_k_energies'].strip()):
            ek_split += [unit_conversion.convert_unit(valueForStrValue(energy, 'f'),'eV')]
        self.tmp['k_energies'].append(ek_split)
        self.tmp['k_point'].append([section['x_qe_t_k_x'][0], section['x_qe_t_k_y'][0], section['x_qe_t_k_z'][0]])

    def onClose_section_run(
            self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""
        if section['x_qe_t_profile_function'] is not None:
            backend.addArrayValues('x_qe_profile_function', np.asarray(
                section['x_qe_t_profile_function']))
            backend.addArrayValues('x_qe_profile_cputime', np.asarray(
                section['x_qe_t_profile_cputime']))
            backend.addArrayValues('x_qe_profile_walltime', np.asarray(
                section['x_qe_t_profile_walltime']))
            backend.addArrayValues('x_qe_profile_ncalls', np.asarray(
                section['x_qe_t_profile_ncalls']))

    def appendToTmp(self, tmpname, value):
        self.tmp[tmpname] += value

    def adHoc_alat(self, parser):
        alat = parser.lastMatch['x_qe_alat']
        self.alat = parser.lastMatch['x_qe_alat']
        unit_conversion.register_userdefined_quantity('usrAlat', 'm', alat)
        unit_conversion.register_userdefined_quantity('usrTpiba', '1/m', 2*math.pi/alat)

    def adHoc_pp_renorm(self, parser):
        self.cache_t_pp_renorm_wfc[parser.lastMatch[
            'x_qe_t_pp_renormalized_filename']] = parser.lastMatch['x_qe_t_pp_renormalized_wfc']

    def adHoc_profiling4(self, parser):
        parser.backend.addValue('x_qe_t_profile_ncalls', 1)
        parser.backend.addValue('x_qe_t_profile_walltime', parser.lastMatch['x_qe_t_profile_cputime'])

    def adHoc_profiling5(self, parser):
        if parser.lastMatch.get('x_qe_t_profile_ncalls', None) is None:
            parser.backend.addValue('x_qe_t_profile_ncalls', 1)
        parser.backend.addValue('x_qe_t_profile_walltime', parser.lastMatch['x_qe_t_profile_cputime'])

    def adHoc_profiling6(self, parser):
        if parser.lastMatch.get('x_qe_t_profile_ncalls', None) is None:
            parser.backend.addValue('x_qe_t_profile_ncalls', 1)
        parser.backend.addValue('x_qe_t_profile_walltime', parser.lastMatch['x_qe_t_profile_cputime'])

    def run_submatchers(self):
        """submatchers of section_run"""
        return [
            SM(name='header',
               startReStr=r"^\s*$",
               sections = ['section_basis_set_cell_dependent', 'section_method', 'section_system',
                           'x_qe_section_parallel', 'x_qe_section_compile_options'],
               subMatchers=[
                   SM(name='serial',
                      startReStr=r"\s*(?P<x_qe_compile_parallel_version>Serial) version\s*$",
                   ),
                   SM(name='serial_multithread',
                      startReStr=r"\s*(?P<x_qe_compile_parallel_version>Serial multi-threaded) version, running on\s*(?P<x_qe_nthreads>\d+)\s*processor cores\s*$",
                   ),
                   SM(name='parallel_mpi',
                      startReStr=r"\s*(?P<x_qe_compile_parallel_version>Parallel version \(MPI\)), running on\s*(?P<x_qe_nproc>\d+)\s*processors\s*$",
                      subMatchers=[
                          SM(name='npool',
                              startReStr=r"\s*K-points division:\s*npool\s*=\s*(?P<x_qe_npool>\d+)\s*$",
                          ),
                      ],
                   ),
                   SM(name='parallel_mpi_old',
                      startReStr=r"\s*(?P<x_qe_compile_parallel_version>Parallel version \(MPI\))\s*$",
                      subMatchers=[
                          SM(name='nproc',
                              startReStr=r"\s*Number of processors in use:\s*(?P<x_qe_nproc>\d+)\s*$",
                          ),
                      ],
                   ),
                   SM(name='qe_input_filename',
                      startReStr=r"\s*Reading input from\s*(?P<x_qe_input_filename>.*?)\s*$",
                   ),
                   SM(name='qe_dimensions',
                      startReStr=r"\s*Current dimensions of program\s*\S+\s*are:\s*$",
                      sections=['x_qe_section_compile_options'],
                      subMatchers=[
                          SM(name="qe_dim_ancient1",
                             startReStr=r"\s*ntypx\s*=\s*(?P<x_qe_ntypx>\d+)\s*npk\s*=\s*(?P<x_qe_npk>\d+)\s*lmax\s*=\s*(?P<x_qe_lmaxx>\d+)\s*$"
                          ),
                          SM(name="qe_dim_ancient2",
                             startReStr=r"\s*nchix\s*=\s*(?P<x_qe_nchix>\d+)\s*ndmx\s*=\s*(?P<x_qe_ndmx>\d+)\s*nbrx\s*=\s*(?P<x_qe_nbrx>\d+)\s*nqfx\s*=\s*(?P<x_qe_nqfx>\d+)\s*$",
                          ),
                          SM(name="qe_dim_old",
                             startReStr=r"\s*ntypx\s*=\s*(?P<x_qe_ntypx>\d+)\s*npk\s*=\s*(?P<x_qe_npk>\d+)\s*lmax\s*=\s*(?P<x_qe_lmaxx>\d+)\s*\s*ndmx\s*=\s*(?P<x_qe_ndmx>\d+)\s*$",
                          ),
                          SM(name="qe_dim_ntypx",
                             startReStr=r"\s*Max number of different atomic species \(ntypx\)\s*=\s*(?P<x_qe_ntypx>\d+)\s*$",
                          ),
                          SM(name="qe_dim_npk",
                             startReStr=r"\s*Max number of k-points \(npk\)\s*=\s*(?P<x_qe_npk>\d+)\s*$",
                          ),
                          SM(name="qe_dim_lmaxx",
                             startReStr=r"\s*Max angular momentum in pseudopotentials \(lmaxx?\)\s*=\s*(?P<x_qe_lmaxx>\d+)\s*$",
                          ),
                      ],
                   ),
                   SM(name='qe_input_filename2',
                      # second possible location for input filename
                      startReStr=r"\s*Reading input from\s*(?P<x_qe_input_filename>.*?)\s*$",
                   ),
                   SM(name='supercell1',
                      startReStr=r"\s*Found symmetry operation:\s*I\s*\+\s*\(\s*" + QeC.re_vec('x_qe_t_vec_supercell') + r"\s*\)\s*$",
                   ),
                   SM(name='supercell2',
                      startReStr=r"\s*This is a supercell, fractional translations? are disabled\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_supercell', True)
                   ),
                   SM(name='pseudopotential_report', repeats=True,
                      startReStr=r"\s*\|\s*pseudopotential report for atomic species\s*:\s*(?P<x_qe_t_pp_report_species>\d+)\s*\|\s*$",
                      endReStr=r"\s*={4,}\s*$",
                      sections=['x_qe_t_section_pp_report'],
                      subMatchers=[
                          SM(name='ppr_version',
                             startReStr=r"\s*\|\s*pseudo potential version(?P<x_qe_t_pp_report_version>(?:\s+\d+)+)\s*\|\s*$",
                          ),
                          SM(name='ppr_separator',
                             startReStr=r"\s*-{4,}\s*$",
                          ),
                          SM(name='ppr_line', repeats=True,
                             startReStr=r"\s*\|  (?P<x_qe_t_pp_report_line>.*?)\s*\|\s*$",
                          ),
                      ],
                   ),
                   SM(name='renormalized_pseudo_wavefunction', repeats=True,
                      startReStr=r"\s*file\s*(?P<x_qe_t_pp_renormalized_filename>.*?)\s*:\s*wavefunction\(s\)\s*(?P<x_qe_t_pp_renormalized_wfc>.*?)\s*renormalized",
                      adHoc=self.adHoc_pp_renorm,
                   ),
                   SM(name='enforced_XC',
                      startReStr=r"\s*IMPORTANT: XC functional enforced from input\s*:\s*",
                      subMatchers=[
                          SM(name='xc_functional_enforced', required=True, # details are parsed in xc_functional
                             startReStr=r"\s*Exchange-correlation\s*=\s*(?P<x_qe_t_xc_functional_shortname_enforced>\S+)\s*\([^\(]+\)\s*$"
                          ),
                          SM(name='exx_fraction_enforced',
                             startReStr=r"\s*EXX-fraction\s*=\s*" + RE_f + r"\s*$",
                          ),
                      ],
                   ),
                   SM(name='gamma_algorithms',
                      startReStr=r"\s*gamma-point specific algorithms are used\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_gamma_algorithms', True)
                   ),
                   SM(name='subspace_diagonalization',
                      startReStr=r"\s*Subspace diagonalization in iterative solution of the eigenvalue problem:\s*$",
                      subMatchers=[
                          SM(name='too_few_procs',
                             startReStr=r"\s*Too few procs for parallel algorithm\s*$",
                             subMatchers=[
                                 SM(name='min_procs',
                                    startReStr=r"\s*we need at least \d+ procs per pool\s*$",
                                 )
                             ],
                          ),
                          SM(name='serial_algorithm',
                             startReStr=r"\s*a serial algorithm will be used\s*$",
                             adHoc=lambda p: p.backend.addValue('x_qe_diagonalization_algorithm', 'serial')
                          ),
                      ],
                   ),
                   SM(name='mesh_sticks',
                      startReStr=r"\s*G-vector sticks info\s*$",
                      subMatchers=[
                          SM(name="sticks_separator",
                             startReStr=r"\s*-+\s*$",
                          ),
                          SM(name="sticks_heading", required=True,
                             startReStr=r"\s*sticks:\s*dense\s+smooth\s+PW\s+G-vecs:\s+dense\s+smooth\s+PW\s*$"
                          ),
                          SM(name='sticks_sum', required=True,
                             startReStr=(
                                 r"\s*Sum\s+(?P<x_qe_sticks_sum_dense>\d+)\s+(?P<x_qe_sticks_sum_smooth>\d+)\s+(?P<x_qe_sticks_sum_PW>\d+)" +
                                 r"\s+(?P<x_qe_sticks_sum_G_dense>\d+)\s+(?P<x_qe_sticks_sum_G_smooth>\d+)\s+(?P<x_qe_sticks_sum_G_PW>\d+)\s*$"
                             ),
                          ),
                          SM(name='sticks_tot',
                             startReStr=r"\s*Tot\s+(?P<x_qe_sticks_tot_dense>\d+)\s+(?P<x_qe_sticks_tot_smooth>\d+)\s+(?P<x_qe_sticks_tot_PW>\d+)\s*$",
                          ),
                      ],
                   ),
                   SM(name='mesh_sticks43',
                      startReStr=r"\s*Stick Mesh\s*$",
                      adHoc=lambda p: LOGGER.error("parse 4.3 sticks mesh properly"),
                      subMatchers=[
                          SM(name="sticks_separator43",
                             startReStr=r"\s*-+\s*$",
                          ),
                          SM(name="sticks_summary43", required=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>nst =.*)",
                          ),
                          SM(name="sticks_header43", required=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>n\.st\s+n\.stw\s+n\.sts\s+n\.g\s+n\.gw\s+n\.gs)",
                          ),
                          SM(name='sticks_line43', repeats=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>(?:min|max|)(?:\s+\d+){6})\s*$",
                          ),
                      ],
                   ),
                   SM(name='mesh_sticks40',
                      startReStr=r"\s*(?P<x_qe_sticks_old>Planes per process \(thick\)\s*:.*?)\s*$",
                      adHoc=lambda p: LOGGER.error("parse <= 4.0 sticks mesh properly"),
                      subMatchers=[
                          SM(name="sticks_smooth40",
                             startReStr=r"\s*(?P<x_qe_sticks_old>Planes per process \(smooth\)\s*:.*?)\s*$",
                          ),
                          SM(name="sticks_header40", required=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>Proc/\D+)\s*$",
                          ),
                          SM(name="sticks_header40_2", required=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>Pool\D+)\s*$",
                          ),
                          SM(name='sticks_line40', repeats=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>(?:\s+\d+){9})\s*$",
                          ),
                      ],
                   ),
                   SM(name='ibrav', required=True,
                      startReStr=r"\s*bravais-lattice index\s*=\s*(?P<x_qe_ibrav>" + RE_i + r")\s*$",
                   ),
                   SM(name='alat', required=True,
                      startReStr=r"\s*lattice parameter \((?:a_0|alat)\)\s*=\s*(?P<x_qe_alat__bohr>" + RE_f + r")\s*a\.u\.\s*$",
                      adHoc=self.adHoc_alat,
                   ),
                   SM(name='cell_volume', required=True,
                      startReStr=r"\s*unit-cell volume\s*=\s*(?P<x_qe_cell_volume__bohr3>" + RE_f + r")\s*\(a\.u\.\)\^3\s*$",
                   ),
                   SM(name='nat', required=True,
                      startReStr=r"\s*number of atoms/cell\s*=\s*(?P<number_of_atoms>\d+)\s*$",
                   ),
                   SM(name='nsp', required=True,
                      startReStr=r"\s*number of atomic types\s*=\s*(?P<x_qe_number_of_species>" + RE_i + r")\s*$",
                   ),
                   SM(name='nelec', required=True,
                      startReStr=r"\s*number of electrons\s*=\s*(?P<x_qe_number_of_electrons>" + RE_f + r")\s*$",
                   ),
                   SM(name='nbnd', required=True,
                       startReStr=r"\s*number of Kohn-Sham states\s*=\s*(?P<x_qe_number_of_states>" + RE_i + r")\s*$"
                   ),
                   SM(name='ecutwfc', required=True,
                      startReStr=r"\s*kinetic-energy cutoff\s*=\s*(?P<basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry\s*$"
                   ),
                   SM(name='ecut_density', required=True,
                      startReStr=r"\s*charge density cutoff\s*=\s*(?P<x_qe_density_basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry\s*$"
                   ),
                   SM(name='convergence_threshold',
                      startReStr=r"\s*convergence threshold\s*=\s*(?P<x_qe_potential_convergence_threshold>" + RE_f + r")\s*$",
                   ),
                   SM(name='mixing_beta',
                      startReStr=r"\s*mixing beta\s*=\s*(?P<x_qe_potential_mixing_beta>" + RE_f + r")\s*$",
                   ),
                   SM(name='mixing_scheme',
                      startReStr=r"\s*number of iterations used\s*=\s*(?P<x_qe_potential_mixing_iterations>\d+)\s*(?P<x_qe_potential_mixing_scheme>.*?)\s+mixing\s*$",
                   ),
                   SM(name='xc_functional', required=True,
                      startReStr=r"\s*Exchange-correlation\s*=\s*(?P<x_qe_xc_functional_shortname>.*?)\s*\((?P<x_qe_xc_functional_num>[^\)]*)\s*\)\s*$"
                   ),
                   SM(name='exx_fraction', required=True,
                      startReStr=r"\s*EXX-fraction\s*=\s*(?P<x_qe_t_exact_exchange_fraction>" + RE_f + r")\s*$",
                   ),
                   SM(name='celldm', repeats = True,
                      startReStr=r"(?P<x_qe_t_celldm>(?:\s*celldm\(\d+\)\s*=\s*" + RE_f + r")+)\s*$",
                   ),
                   SM(name='simulation_cell',
                      startReStr=r"\s*crystal axes: \(cart. coord.\s*in units of (?:a_0|alat)\s*\)\s*$",
                      subMatchers=[
                          SM(name='cell_vec_a', repeats=True,
                             startReStr=r"\s*a\(\d\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_t_vec_a', 'usrAlat') + r"\s*\)\s*$",
                          ),
                      ],
                   ),
                   SM(name='reciprocal_cell',
                      startReStr=r"\s*reciprocal axes: \(cart. coord. in units 2 pi/(?:alat|a_0)\)\s*$",
                      subMatchers=[
                          SM(name='cell_vec_b', repeats=True,
                             startReStr=r"\s*b\(\d\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_t_vec_b', 'usrTpiba') + r"\s*\)\s*$",
                          ),
                      ],
                   ),
                   SM(name='pseudopotential', repeats=True,
                      startReStr=(r"\s*PseudoPot\.\s*#\s*(?P<x_qe_t_pp_idx>" + RE_i + r") for\s+(?P<x_qe_t_pp_label>\S+)\s+read from file" +
                                  r"(?::|\s*(?P<x_qe_t_pp_filename>\S+))\s*"),
                      sections=['x_qe_t_section_pseudopotential'],
                      subMatchers=[
                          SM(name='new_pp_filename',
                             startReStr=r"^\s*(?P<x_qe_t_pp_filename>\S+)\s*$",
                          ),
                          SM(name='pp_md5',
                             startReStr=r"\s*MD5 check sum:\s*(?P<x_qe_t_pp_md5sum>\S+)\s*$",
                          ),
                          SM(name='pp_type_val',
                             startReStr=r"\s*Pseudo is\s*(?P<x_qe_t_pp_type>.*?),\s*Zval\s*=\s*(?P<x_qe_t_pp_valence>" + RE_f + r")\s*$",
                          ),
                          SM(name='pp_comment',
                             startReStr=r"\s*(?P<x_qe_t_pp_comment>.*?)\s*$",
                          ),
                          SM(name='pp_augmentation_shape',
                             startReStr=r"\s*Shape of augmentation charge:\s*(?P<x_qe_t_pp_augmentation_shape>.*?)\s*$",
                          ),
                          SM(name='pp_dimensions',
                             startReStr=r"\s*Using radial grid of\s*(?P<x_qe_t_pp_ndmx>\d+) points,\s*(?P<x_qe_t_pp_nbeta>\d+) beta functions with\s*:\s*$",
                          ),
                          SM(name='pp_beta', repeats=True,
                             startReStr=r"\s*l\(\s*(?P<x_qe_t_pp_l_idx>\d+)\s*\)\s*=\s*(?P<x_qe_t_pp_l>\d+)\s*$",
                          ),
                          SM(name='pp_coefficients',
                             startReStr=r"\s*Q\(r\) pseudized with\s*(?P<x_qe_t_pp_ncoefficients>\d+)\s*coefficients,\s*rinner\s*=\s*(?P<x_qe_t_rinner>(?:\s*" + RE_f + r")+)\s*$",
                             subMatchers=[
                                 SM(name='pp_coefficients2',
                                    startReStr=r"\s*(?P<x_qe_t_rinner>(?:\s*" + RE_f + r")+)\s*$",
                                 ),
                             ],
                          ),
                          SM(name='pp_coefficients3',
                             startReStr=r"\s*Q\(r\) pseudized with\s*(?P<x_qe_t_pp_ncoefficients>\d+)\s*coefficients\s*$",
                          ),
                       ],
                   ),
                   SM(name='pp_atom_kind_map',
                      startReStr=r"\s*atomic species\s+valence\s+mass\s+pseudopotential\s*$",
                      subMatchers=[
                          SM(name='atom_kind', repeats=True,
                             startReStr=r"\s*(?P<method_atom_kind_label>\S+)\s+(?P<x_qe_pp_valence>" + RE_f + r")" +
                                        r"\s+(?P<x_qe_kind_mass>" + RE_f + r")\s+(?P<x_qe_pp_label>[^\s\(]+)\s*" +
                                        r"\(\s*(?P<x_qe_pp_weight>" + RE_f + r")\s*\)\s*$",
                             sections=['section_method_atom_kind'],
                          ),
                      ],
                   ),
                   SM(name='nsymm',
                      startReStr=(r"\s*(?P<x_qe_nsymm>\d+)\s*Sym\.\s*Ops\.\s*[\(,]\s*(?P<x_qe_t_symm_inversion>\S+) inversion\s*[\),]\s*(?:found)?\s*"
                                  r"(?:\(\s*(?P<x_qe_nsymm_with_fractional_translation>\d+)\s*have fractional translation\s*\))?\s*$"),
                      adHoc=lambda p: p.backend.addValue('x_qe_symm_inversion', (p.lastMatch['x_qe_t_symm_inversion'] == 'with')),
                      subMatchers=[
                          SM(name='nsymm_ignored',
                             startReStr=r"\s*\(note:\s*(?P<x_qe_nsymm_ignored>\d+)\s*additional sym.ops. were found but ignored\s*$",
                          ),
                      ],
                   ),
                   SM(name='nosymm',
                      startReStr=r"\s*No symmetry\s*(?:found|!)\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_nsymm', 0)
                   ),
                   SM(name='atom_pos_cart_list',
                      startReStr=r"\s*Cartesian axes\s*$",
                      subMatchers=[
                          SM(name='cart_heading',
                             startReStr=r"\s*site n.     atom                  positions \((?:a_0|alat) units\)\s*$",
                          ),
                          SM(name='atom_pos_cart', repeats=True,
                             startReStr=(
                                 r"\s*(?P<x_qe_t_atom_idx>" + RE_i + r")" +
                                 r"\s+(?P<x_qe_t_atom_labels>\S+)\s+tau\(\s*" + RE_i + "\)\s*"
                                 r"=\s*\(\s*" + QeC.re_vec('x_qe_t_atpos', 'usrAlat') +
                                 r"\s*\)\s*$"),
                          ),
                          SM(name='kpoint_info',
                             startReStr=r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*$",
                             subMatchers=[
                                 SM(name="kpoint_heading",
                                    startReStr=r"\s*cart. coord. in units 2pi/(?:alat|a_0)\s*$",
                                 ),
                                 SM(name="kpoint_kpoints", repeats=True,
                                    startReStr=(r"\s*k\(\s*(?P<x_qe_t_k_info_ik>\d+)\s*\)\s*=\s*\(\s*" +
                                                QeC.re_vec('x_qe_t_k_info_vec', 'usrAlat') +
                                                r"\s*\),\s*wk\s*=\s*(?P<x_qe_t_k_info_wk>" + RE_f + r")\s*$"),

                                 ),
                             ],
                          ),
                          SM(name='kpoint_info_smearing_old',
                             startReStr=(r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*gaussian broad.\s*\(Ry\)\s*=" +
                                         r"\s*(?P<x_qe_smearing_width__rydberg>" + RE_f + r")\s*ngauss\s*=" +
                                         r"\s*(?P<x_qe_ngauss>" + RE_i + ")\s*$"),
                             adHoc=lambda p: LOGGER.error("translate ngauss to name of smearing"),
                             subMatchers=[
                                 SM(name="kpoint_heading",
                                    startReStr=r"\s*cart. coord. in units 2pi/(?:alat|a_0)\s*$",
                                 ),
                                 SM(name="kpoint_kpoints", repeats=True,
                                    startReStr=(r"\s*k\(\s*(?P<x_qe_t_k_info_ik>\d+)\s*\)\s*=\s*\(\s*" +
                                                QeC.re_vec('x_qe_t_k_info_vec', 'usrAlat') +
                                                r"\s*\),\s*wk\s*=\s*(?P<x_qe_t_k_info_wk>" + RE_f + r")\s*$"),

                                 ),
                             ],
                          ),
                      ],
                   ),
                   SM(name='dense_grid',
                      startReStr=(r"\s*Dense\s+grid:\s*(?P<x_qe_dense_g_vectors>\d+)\s*G-vectors\s*FFT\s+dimensions:\s*\(\s*" +
                                  QeC.re_vec("x_qe_t_dense_FFT_grid", split=r"\s*,\s*") + "\s*\)\s*$")
                   ),
                   SM(name='dense_grid_old',
                      startReStr=(r"\s*G\s+cutoff\s*=\s*(?P<x_qe_dense_g_cutoff>" + RE_f + r")\s*" +
                                  r"\(\s*(?P<x_qe_dense_g_vectors>\d+)\s*G-vectors\s*\)\s*FFT\s+grid:\s*\(\s*" +
                                  QeC.re_vec("x_qe_t_dense_FFT_grid", split=r"\s*,\s*") + "\s*\)\s*$"
                      ),
                   ),
                   SM(name='smooth_grid',
                      startReStr=(r"\s*Smooth\s+grid:\s*(?P<x_qe_smooth_g_vectors>\d+)\s*G-vectors\s*FFT\s+dimensions:\s*\(\s*" +
                                  QeC.re_vec("x_qe_t_smooth_FFT_grid", split=r"\s*,\s*") + "\s*\)\s*$")
                   ),
                   SM(name='smooth_grid_old',
                      startReStr=(r"\s*G\s+cutoff\s*=\s*(?P<x_qe_smooth_g_cutoff>" + RE_f + r")\s*" +
                                  r"\(\s*(?P<x_qe_smooth_g_vectors>\d+)\s*G-vectors\s*\)\s*smooth\s+grid:\s*\(\s*" +
                                  QeC.re_vec("x_qe_t_smooth_FFT_grid", split=r"\s*,\s*") + "\s*\)\s*$"
                      ),
                   ),
                   SM(name='core_charge_realspace',
                      startReStr=r"\s*Real space treatment of Q\(r\)\s*",
                      adHoc=lambda p: p.backend.addValue('x_qe_core_charge_realspace', True)
                   ),
                   SM(name='input_occupations',
                      startReStr=r"\s*Occupations\s*read\s*from\s*input\s*$",
                      sections=['x_qe_t_section_input_occupations'],
                      subMatchers=[
                          SM(name='input_occupations_occupations', repeats=True,
                              startReStr=r'\s*(?P<x_qe_t_input_occupations>(?:\s*' + RE_f + ')+\s*$)',
                          ),
                      ],
                   ),
                   SM(name='allocated_arrays',
                      startReStr=r"\s*Largest allocated arrays\s*est. size \(Mb\)\s*dimensions\s*$",
                      subMatchers=[
                          SM(name='allocated_array', repeats=True,
                             startReStr=(
                                 r"\s*(?P<x_qe_t_allocated_array_name>.*?)\s*(?P<x_qe_t_allocated_array_size__mebibyte>" +
                                 RE_f + r")\s*Mb\s*\(\s*(?P<x_qe_t_allocated_array_dimensions>(?:\s*\d+\s*,?)+)\s*\)\s*$"
                             ),
                          ),
                      ],
                   ),
                   SM(name='temporary_arrays',
                      startReStr=r"\s*Largest temporary arrays\s*est. size \(Mb\)\s*dimensions\s*$",
                      subMatchers=[
                          SM(name='temporary_array', repeats=True,
                             startReStr=(
                                 r"\s*(?P<x_qe_t_temporary_array_name>.*?)\s*(?P<x_qe_t_temporary_array_size__mebibyte>" +
                                 RE_f + r")\s*Mb\s*\(\s*(?P<x_qe_t_temporary_array_dimensions>(?:\s*\d+\s*,?)+)\s*\)\s*$"
                             ),
                          ),
                      ],
                   ),
                   SM(name='core_charge_check',
                      startReStr=(r"\s*Check: negative/imaginary core charge\s*=\s*(?P<x_qe_core_charge_negative>" +
                                  RE_f + r")\s*(?P<x_qe_core_charge_imaginary>" + RE_f + r")\s*$")
                   ),
                   SM(name='initial_potential',
                      startReStr=r"\s*Initial potential from\s*(?P<x_qe_starting_potential>.*?)\s*$",
                   ),
                   SM(name='starting_charge_negative',
                      startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_starting_charge_negative>" + RE_f +
                                  r")\s*$"),
                   ),
                   SM(name='initial_charge',
                      startReStr=(r"\s*starting charge\s*(?P<x_qe_starting_charge>" + RE_f +
                                  r")\s*,\s*renormalised to\s*(?P<x_qe_starting_charge_renormalized>" + RE_f +
                                  r")\s*$"),
                   ),
                   SM(name='starting_rho',
                      startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_starting_charge_negative_up>" + RE_f +
                                  r")\s*(?P<x_qe_starting_charge_negative_down>" + RE_f + r")\s*$"),
                   ),
                   SM(name='starting_wfc',
                      startReStr=r"\s*Starting wfc are\s*(?P<x_qe_starting_wfc>.*?)\s*$",
                   ),
                   SM(name='cputime_msg',
                      startReStr=(r"\s*total cpu time spent up to now is\s*(?P<x_qe_time_setup_cpu1_end>" + RE_f +
                                  r")\s*secs"),
                   ),
                   SM(name='per_process_mem',
                      startReStr=r"\s*per-process dynamical memory:\s*(?P<x_qe_per_process_mem__mebibyte>" + RE_f + ")\s*Mb",
                   ),
               ],
            ), # header
            SM(name='self_consistent_calculation',
               startReStr=r"\s*Self-consistent Calculation\s*$",
               sections = ['section_single_configuration_calculation'],
               subMatchers=[
                   SM(name='iteration', repeats=True,
                      startReStr=(r"\s*iteration\s*#\s*(?P<x_qe_iteration_number>\d+)\s*" +
                                  r"\s*ecut\s*=\s*(?P<x_qe_iteration_ecutwfc>" + RE_f +r")\s*Ry" +
                                  r"\s*beta\s*=\s*(?P<x_qe_iteration_beta>" + RE_f + r")\s*$"),
                      sections=['section_scf_iteration'],
                      subMatchers=[
                          SM(name='david_overlap',
                              startReStr=r"\s*(?P<x_qe_t_david_with_overlap>Davidson diagonalization with overlap.*)\s*$",
                          ),
                          SM(name='cg_style_diag',
                              startReStr=r"\s*(?P<x_qe_t_cg_diag>CG style diagonalization.*)\s*$",
                          ),
                          SM(name='ethr',
                             startReStr=(r"\s*ethr\s*=\s*(?P<x_qe_t_iteration_ethr>" + RE_f +
                                         r")\s*,\s*avg\s*#\s*of iterations\s*=\s*(?P<x_qe_t_iteration_avg>" + RE_f +
                                         r")\s*$")
                          ),
                          SM(name='david_overlap2',
                             startReStr=r"\s*(?P<x_qe_t_david_with_overlap>Davidson diagonalization with overlap.*)\s*$",
                          ),
                          SM(name='cg_style_diag2',
                              startReStr=r"\s*(?P<x_qe_t_cg_diag>CG style diagonalization.*)\s*$",
                          ),
                          SM(name='ethr2',
                             startReStr=(r"\s*ethr\s*=\s*(?P<x_qe_t_iteration_ethr>" + RE_f +
                                         r")\s*,\s*avg\s*#\s*of iterations\s*=\s*(?P<x_qe_t_iteration_avg>" + RE_f +
                                         r")\s*$")
                          ),
                          SM(name='iteration_rho',
                             startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_iteration_charge_negative_up>" + RE_f +
                                         r")\s*(?P<x_qe_iteration_charge_negative_down>" + RE_f + r")\s*$"),
                          ),
                          SM(name='cputime_iteration_msg',
                             startReStr=(r"\s*total cpu time spent up to now is\s*(?P<time_scf_iteration_cpu1_end>" + RE_f +
                                         r")\s*secs"),
                          ),
                          SM(name='e_total',
                             startReStr=(r"\s*!?\s*total\s+energy\s*=\s*(?P<energy_total_scf_iteration>" + RE_f + r")" +
                                         r"\s*Ry\s*$"),
                          ),
                          SM(name='harris',
                             startReStr=(r"\s*Harris-Foulkes estimate\s*=\s*(?P<x_qe_energy_total_harris_foulkes_estimate_iteration>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='estimate_accuracy',
                             startReStr=(r"\s*estimated scf accuracy\s*<\s*(?P<x_qe_energy_total_accuracy_estimate_iteration>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                      ],
                   ),
                   SM(name='scf_result', repeats=False,
                       startReStr=r'\s*End of self-consistent calculation',
                       sections=['section_eigenvalues'],
                       subMatchers=[
                          SM(name='bands', repeats=True,
                              sections=['x_qe_t_section_kbands'],
                              startReStr=r'\s*k\s*=\s*' + QeC.re_vec('x_qe_t_k', 'usrTpiba', '\s*') + r'\s*\(\s*(?P<x_qe_t_k_pw>' + RE_i + ')\s*PWs',
                              subMatchers=[
                                  SM(name='kbnd', repeats=True,
                                      startReStr=r'\s*(?P<x_qe_t_k_point_energies>(?:\s*' + RE_f + ')+\s*$)',
                                      adHoc=lambda p: self.appendToTmp('this_k_energies', " " + p.lastMatch['x_qe_t_k_point_energies']),
                                  ),
                              ],
                          ),
                          SM(name='highest_occupied',
                             startReStr=(r"\s*highest occupied level \(ev\):\s*(?P<x_qe_highest_occupied__eV>" + RE_f +
                                         r")\s*$"),
                             adHoc=lambda p: LOGGER.error("apply highest occupied"),
                          ),
                          SM(name='highest_occupied_lowest_unoccupied',
                             startReStr=(r"\s*highest occupied, lowest unoccupied level \(ev\):\s*" +
                                         r"(?P<x_qe_highest_occupied__eV>" + RE_f + r")\s*" +
                                         r"(?P<x_qe_lowest_unoccupied__eV>" + RE_f + r")\s*$"),
                             adHoc=lambda p: LOGGER.error("apply highest occupied/lowest unoccupied"),
                          ),
                          SM(name='e_fermi',
                             startReStr=(r"\s*the Fermi energy is\s*(?P<x_qe_energy_fermi>" + RE_f + ")\s*ev\s*$"),
                             adHoc=lambda p: LOGGER.error("apply fermi energy"),
                          ),
                          SM(name='e_total',
                             startReStr=r'\s*!?\s*total\s+energy\s*=\s*(?P<energy_total>' + RE_f + ')\s*Ry\s*$',
                          ),
                          SM(name='harris',
                             startReStr=(r"\s*Harris-Foulkes estimate\s*=\s*(?P<x_qe_energy_total_harris_foulkes_estimate>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='estimate_accuracy',
                             startReStr=(r"\s*estimated scf accuracy\s*<\s*(?P<x_qe_energy_total_accuracy_estimate>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='total_AE_energy',
                             startReStr=(r"\s*total all-electron energy\s*=\s*(?P<x_qe_energy_total_paw_all_electron>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='energy_decomposition',
                             startReStr=r"\s*The total energy is the sum of the following terms\s*:\s*$",
                             subMatchers=[
                                 SM(name='energy_decomposition_contribution', repeats=True,
                                    startReStr=(r"\s*(?P<x_qe_t_energy_decomposition_name>.*?)\s*=\s*" +
                                                r"(?P<x_qe_t_energy_decomposition_value__rydberg>" + RE_f + r")\s*Ry"),
                                 ),
                             ],
                          ),
                          SM(name="convergence_iterations",
                             startReStr=r"\s*convergence has been achieved in\s*(?P<x_qe_convergence_iterations>\d+)\s*iterations\s*",
                          ),
                          SM(name="convergence_achieved",
                             startReStr=r"\s*convergence has been achieved\s*$",
                             adHoc=lambda p: p.backend.addValue('x_qe_convergence_iterations', p.superContext.tmp['last_iteration'])
                          ),
                          SM(name="warning_save_mgga",
                             startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                          ),
                          SM(name="atom_forces",
                             startReStr=r"\s*Forces acting on atoms\s*\(Ry/au\):\s*$",
                             subMatchers=[
                                 SM(name="atom_force", repeats=True,
                                    startReStr=(r"\s*atom\s*(?P<x_qe_t_force_atom_idx>\d+)"+
                                                r"\s*type\s*\d+\s*force\s*=\s*"+
                                                QeC.re_vec('x_qe_t_force', 'rydberg_bohr_1') + "\s*$")
                                 ),
                                 SM(name="total_force",
                                    startReStr=(r"\s*Total\s+force\s*=\s*(?P<x_qe_force_total__rydberg_bohr_1>" + RE_f +
                                                r")\s*Total SCF correction =\s*" +
                                                r"(?P<x_qe_force_total_scf_correction__rydberg_bohr_1>" + RE_f +
                                                r")\s*$"),
                                 ),
                             ],
                          ),
                          SM(name="stress_tensor",
                             startReStr=r"\s*entering subroutine stress \.\.\.\s*$",
                             subMatchers=[
                                 SM(name="stress_header",
                                    startReStr=(r"\s*total\s*stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*" +
                                                r"(?P<x_qe_pressure__kilobar>" + RE_f + r")\s*$"),
                                 ),
                                 SM(name="stress_components", repeats=True,
                                    startReStr=(r"\s*" + QeC.re_vec('x_qe_t_stress', 'rydberg_bohr_3') + "\s*" +
                                                RE_f + r"\s*" + RE_f + r"\s*" + RE_f + r"\s*$")
                                 ),
                             ],
                          ),
                          SM(name="stress_message",
                             startReStr=r"\s*Message from routine stress\s*",
                             subMatchers=[
                                 SM(name="no_mgga",
                                    startReStr=r"\s*Meta-GGA and stress not implemented\s*$",
                                    adHoc=lambda p: p.backend.addValue("x_qe_stress_unimplemented", "Meta-GGA")
                                 ),
                             ],
                          ),
                          SM(name="write_datafile",
                             startReStr=r"\s*Writing output data file\s*(?P<x_qe_output_datafile>.*?)\s*$",
                          ),
                          SM(name="warning_save_mgga2",
                             startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                          ),
                          SM(name="profiling0", repeats=True,
                             startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" + RE_f +
                                         r")\s*s\s*CPU\s*(?P<x_qe_t_profile_walltime__s>" + RE_f + r")\s*s\s*WALL" +
                                         r"\s*\(\s*(?P<x_qe_t_profile_ncalls>\d+)\s*calls\s*\)\s*$"),
                          ),
                          SM(name="profiling_caller", repeats=True,
                             startReStr=r"\s*(?:Called by)?\s*(?P<x_qe_t_profile_caller>\S+?)(?:\s*routines)?:?\s*$",
                             subMatchers=[
                                 SM(name="profiling", repeats=True,
                                    startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" + RE_f +
                                                r")\s*s\s*CPU\s*(?P<x_qe_t_profile_walltime__s>" + RE_f + r")\s*s\s*WALL" +
                                                r"\s*\(\s*(?P<x_qe_t_profile_ncalls>\d+)\s*calls\s*\)\s*$"),
                                 ),
                             ],
                          ),
                          SM(name="profiling2", repeats=True,
                             startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*" +
                                         r"(?P<x_qe_t_profile_cputime__strQeTimespan>.*)CPU" +
                                         r"\s*(?P<x_qe_t_profile_walltime__strQeTimespan>.*)\s*WALL\s*$"),
                             adHoc=lambda p: p.backend.addValue('x_qe_t_profile_ncalls', 1),
                          ),
                          # profiling of 'old' espresso
                          SM(name="profiling3",
                             startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" + RE_f +
                                         r")\s*s\s*CPU\s*time,\s*(?P<x_qe_t_profile_walltime__s>" + RE_f + r")\s*s\s*wall\s*time\s*$"),
                             adHoc=lambda p: p.backend.addValue('x_qe_t_profile_ncalls', 1),
                          ),
                          SM(name="profiling4", repeats=True,
                             startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" + RE_f +
                                         r")\s*s\s*CPU\s*$"),
                             adHoc=self.adHoc_profiling4,
                          ),
                          SM(name="profiling5", repeats=True,
                             startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" +
                                         RE_f + r")\s*s\s*CPU" +
                                         r"\s*(?:\(\s*(?P<x_qe_t_profile_ncalls>\d+)\s*calls,\s*" + RE_f +
                                         r"\s*s\s*avg\s*\))?\s*$"),
                             adHoc=self.adHoc_profiling5,
                          ),
                          SM(name="profiling6_caller_general", repeats=True,
                             startReStr=r"\s*(?:Called by)?\s*(?P<x_qe_t_profile_caller>\S+?)(?:\s*routines)?:?\s*$",
                             subMatchers=[
                                 SM(name="profiling6", repeats=True,
                                    startReStr=(r"\s*(?P<x_qe_t_profile_function>\S+)\s*:\s*(?P<x_qe_t_profile_cputime__s>" +
                                                RE_f + r")\s*s\s*CPU" +
                                                r"(?:\s*\(\s*(?P<x_qe_t_profile_ncalls>\d+)\s*calls,\s*" + RE_f +
                                                r"\s*s\s*avg\s*\))?\s*$"),
                                    adHoc=self.adHoc_profiling6
                                 ),
                             ],
                          ),
                       ],
                   ),
               ],
            ),
        ]

if __name__ == "__main__":
    parser = QuantumEspressoParserPWSCF()
    parser.parse()
