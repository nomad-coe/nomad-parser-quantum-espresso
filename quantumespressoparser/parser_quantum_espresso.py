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
import quantumespressoparser.QuantumEspressoCommon as QeC
from nomadcore.parser_backend import valueForStrValue
from quantumespressoparser.QuantumEspressoCommon import RE_f, RE_i, cRE_f, cRE_i
from quantumespressoparser.QuantumEspressoXC import translate_qe_xc_num
from nomadcore.parser_backend import valueForStrValue


LOGGER = logging.getLogger(__name__)


# Lookup table mapping string to bool flag
QE_SPIN_NONCOLLINEAR = {
    'Noncollinear': True,
    'Non magnetic': False,
}

# Lookup table mapping diagonalization scheme match to string
QE_DIAGONALIZATION = {
    'Davidson diagonalization with overlap': 'davidson',
    'CG style diagonalization': 'conjugate_gradient',
}

# Lookup table for VCSMD notation of atom position units
QE_VCSMD_ATPOS_UNITS = {
    'cryst coord': 'crystal',
    'cart coord (alat unit)': 'alat',
}

QE_MD_RELAX_SAMPLING_METHOD = {
    'langevin_overdamped_dynamics': 'langevin_dynamics',
    'BFGS': 'geometry_optimization',
    'damped_dynamics': 'geometry_optimization',
    'molecular_dynamics': 'molecular_dynamics',
    'vcsmd': 'molecular_dynamics',
    'vcsmd_wentzcovitch_damped_minization': 'geometry_optimization',
}

class QuantumEspressoParserPWSCF(QeC.ParserQuantumEspresso):
    """main place to keep the parser status, open ancillary files,..."""
    def __init__(self, metainfo_to_keep=None, backend=None, default_units=None,
        metainfo_units=None, debug=True, log_level=logging.ERROR, store=True):
        QeC.ParserQuantumEspresso.__init__(
            self, re_program_name=re.compile(r"^PWSCF$"), metainfo_to_keep=metainfo_to_keep,
            backend=backend, default_units=default_units, metainfo_units=metainfo_units,
            debug=debug, log_level=log_level, store=store)

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse
        different files"""
        self.sectionIdx = {}
        self.openSectionIdx = {}
        self.tmp = {}
        self.alat = None
        self.section = {}

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
        method_basis_set_gIndex = backend.openSection('section_method_basis_set')
        backend.addValue('mapping_section_method_basis_set_cell_associated', gIndex)
        backend.addValue('method_basis_set_kind', self.tmp['method_basis_set_kind'])
        backend.closeSection('section_method_basis_set', method_basis_set_gIndex)

    def onOpen_section_method(
            self, backend, gIndex, section):
        self.sectionIdx['section_method'] = gIndex
        self.cache_t_pseudopotential = {}
        self.cache_t_pp_report = {}
        self.cache_t_pp_renorm_wfc = {}
        self.cache_t_method = section
        self.atom_kind_idx = -1

    def onClose_x_qe_t_section_pseudopotential(
            self, backend, gIndex, section):
        if section['x_qe_t_pp_label'] is None:
            LOGGER.error("x_qe_t_pp_label is not set")
            return
        self.cache_t_pseudopotential[section['x_qe_t_pp_label'][0]] = section

    def onClose_section_method_atom_kind(
            self, backend, gIndex, section):
        self.atom_kind_idx += 1
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
            if len(basename) > 0:
                backend.addValue('method_atom_kind_pseudopotential_name',
                                 basename)
            renorm_info = self.cache_t_pp_renorm_wfc.get(basename, None)
            if renorm_info is not None:
                backend.addValue('x_qe_pp_renormalized_wfc', renorm_info)
        # set charge/magnetization integration radius if applicable
        integration_radius=self.cache_t_method['x_qe_t_species_integration_radius']
        if integration_radius:
            backend.addValue('x_qe_species_integration_radius',
                             unit_conversion.convert_unit(
                                 integration_radius[self.atom_kind_idx],
                                 'usrAlat'))
        # DFT-D data
        cache_dft_d = self.tmp.get('dispersion_correction', {})
        dft_d = cache_dft_d.get(pp_label, None)
        if dft_d is not None:
            for k, v in dft_d.items():
                backend.addValue(k, v)

    def onClose_x_qe_t_section_pp_report(
            self, backend, gIndex, section):
        self.cache_t_pp_report[section['x_qe_t_pp_report_species'][0]] = section

    def onClose_x_qe_t_section_pp_warning(
            self, backend, gIndex, section):
        self.cache_t_pp_renorm_wfc[
            section['x_qe_t_pp_warning_filename'][-1]] = section['x_qe_t_pp_warning_wfclabel'][-1]

    def onOpen_x_qe_t_section_input_occupations(
            self, backend, gIndex, section):
        self.tmp['occ_spin'] = 'none'
        self.tmp['occ_vals_spin'] = {}
        self.tmp['occ_vals_spin']['none'] = []
        self.tmp['occ_vals'] = self.tmp['occ_vals_spin']['none']

    def adHoc_input_occupations_spin(self, parser):
        self.tmp['occ_spin'] = parser.lastMatch['x_qe_t_input_occupations_spin']
        self.tmp['occ_vals_spin'][self.tmp['occ_spin']] = []
        self.tmp['occ_vals'] = self.tmp['occ_vals_spin'][self.tmp['occ_spin']]

    def onClose_x_qe_t_section_input_occupations(
            self, backend, gIndex, section):
        if len(self.tmp['occ_vals_spin']['none']) > 0:
            # spin-unpolarized case
            input_occ = np.array([ self.tmp['occ_vals_spin']['none'] ])
        else:
            input_occ = np.array([
                self.tmp['occ_vals_spin']['up'],
                self.tmp['occ_vals_spin']['down']
            ])
        if input_occ.size < 1:
            LOGGER.error("closing x_qe_t_section_input_occupations, but no input occupations parsed")
        else:
            LOGGER.error("FIXME: implement proper output of x_qe_t_section_input_occupations (shape %s)", str(input_occ.shape))

    def onClose_section_method(
            self, backend, gIndex, section):
        # set flag if we deal with user-enforced XC functional
        if section['x_qe_t_xc_functional_shortname_enforced'] is not None:
            backend.addValue('x_qe_xc_functional_user_enforced', True)
        # translate XC functional to section_xc_functionals
        method_xc_functionals = None
        xc_functionals = None
        if section['x_qe_xc_functional_num'] is not None:
            xc_functional_num = section['x_qe_xc_functional_num'][-1]
            exx_fraction = None
            if section['x_qe_t_exact_exchange_fraction']:
                # first SCF in EXX is done without HF-X, cache for later
                self.tmp['exx_fraction'] = section['x_qe_t_exact_exchange_fraction'][-1]
                self.tmp['xc_functional_num'] = xc_functional_num
                exx_fraction = 0.0
            elif section['x_qe_exact_exchange_fraction']:
                exx_fraction = section['x_qe_exact_exchange_fraction'][-1]
            (method_xc_functionals, xc_functionals) = translate_qe_xc_num(xc_functional_num, exx_fraction)
        else:
            LOGGER.error("x_qe_xc_functional_num is not set")
        if method_xc_functionals is not None:
            # NOTE: value of XC_functional generated by translate_qe_xc_num
            #       does not fully respect the metaInfo definition
            #       when XC_functional_parameters are involved.
            # Therefore, remove it here
            method_xc_functionals.pop('XC_functional', None)
            self.addDict(backend, method_xc_functionals)
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
        if section['x_qe_t_spin_orbit_magn'] is not None:
            backend.addValue('x_qe_spin_orbit', (section['x_qe_t_spin_orbit_mode'][-1] == 'with'))
            noncollinear = QE_SPIN_NONCOLLINEAR.get(section['x_qe_t_spin_orbit_magn'][-1], None)
            if noncollinear is None:
                LOGGER.error("unimplemented value for 'x_qe_t_spin_orbit_magn': '%s'",
                             str(section['x_qe_t_spin_orbit_magn'][-1]))
            else:
                backend.addValue('x_qe_spin_noncollinear', noncollinear)
        # TODO check for LDA+U and switch to 'DFT+U' in that clase
        backend.addValue('electronic_structure_method', 'DFT')

    def onOpen_section_single_configuration_calculation(
            self, backend, gIndex, section):
        exx_refine = self.tmp.pop('exx_refine', None)
        md_relax = self.tmp.get('md_relax', None)
        have_new_system = self.tmp.pop('have_new_system', None)
        old_scc = self.section.get('single_configuration_calculation', None)
        if old_scc is not None:
            if md_relax is not None:
                LOGGER.info('new scc due to md_relax=="%s"', str(md_relax))
                if self.tmp.get('frames', None) is None:
                    # add previous scc
                    self.tmp['frames'] = [self.sectionIdx['single_configuration_calculation']]
            elif exx_refine is not None:
                LOGGER.info('new scc due to exx_refine=="%s"', str(exx_refine))
            else:
                del(self.section['single_configuration_calculation'])
                del(self.sectionIdx['single_configuration_calculation'])
                LOGGER.warn("encountered new section_single_configuration_calculation without knowing why")
        if md_relax:
            if have_new_system:
                self.tmp['frames'].append(gIndex)
            elif exx_refine:
                LOGGER.info("EXX refinement in MD/relax run")
            else:
                raise Exception('running MD/relax calculation, but no new cell or atom coordinates found')

        if exx_refine:
            backend.addValue('x_qe_exx_refine', True)
            exx_fraction = self.tmp.pop('exx_fraction', None)
            if exx_fraction:
                # we need to add a section_method including exx
                method_gIndex = backend.openSection('section_method')
                backend.addValue('x_qe_xc_functional_num', self.tmp.pop('xc_functional_num'))
                backend.addValue('x_qe_exact_exchange_fraction',  exx_fraction)
                self.addSectionDict(
                    backend, 'section_method_to_method_refs', {
                        'method_to_method_ref': (method_gIndex-1),
                        'method_to_method_kind': 'starting_point',
                    }
                )
                backend.closeSection('section_method', method_gIndex)
            self.addSectionDict(
                backend, 'section_calculation_to_calculation_refs', {
                    'calculation_to_calculation_ref': (gIndex-1),
                    'calculation_to_calculation_kind': 'starting_point',
                }
            )
        self.close_header_sections(backend)
        # reset temporary storage for band structures
        self.tmp['k_energies'] = []
        self.tmp['k_occupations'] = []
        self.tmp['kspin'] = {}
        self.section['single_configuration_calculation'] = section
        self.sectionIdx['single_configuration_calculation'] = gIndex

    def onClose_section_single_configuration_calculation(
            self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""

        if 'extra_SCF' in self.tmp:
            backend.addValue('x_qe_extra_SCF', True)
            del(self.tmp['extra_SCF'])

        backend.addValue('single_configuration_to_calculation_method_ref', self.sectionIdx['section_method'])
        backend.addValue('single_configuration_calculation_to_system_ref', self.sectionIdx['section_system'])
        # extract k band structure data if available
        self.create_section_eigenvalues(backend, section)
        if section['x_qe_t_energy_decomposition_name'] is not None:
            backend.addArrayValues('x_qe_energy_decomposition_name', np.asarray(
                section['x_qe_t_energy_decomposition_name']))
        if section['x_qe_t_energy_decomposition_value'] is not None:
            backend.addArrayValues('x_qe_energy_decomposition_value', np.asarray(
                section['x_qe_t_energy_decomposition_value']))
        if section['x_qe_t_force_x'] is not None:
            # constraints etc. not part of the reported forces, so correct metaInfo is 'atom_forces_raw'
            backend.addArrayValues('atom_forces_raw', np.array([
                section['x_qe_t_force_x'], section['x_qe_t_force_y'], section['x_qe_t_force_z']
                ]).T)
        if section['x_qe_t_dispersion_force_x'] is not None:
            backend.addArrayValues('x_qe_atom_dispersion_force', np.array([
                section['x_qe_t_dispersion_force_x'], section['x_qe_t_dispersion_force_y'], section['x_qe_t_dispersion_force_z']
                ]).T)
        if section['x_qe_t_stress_x'] is not None:
            backend.addArrayValues('stress_tensor', np.array([
                section['x_qe_t_stress_x'], section['x_qe_t_stress_y'], section['x_qe_t_stress_z']
                ]).T)
        had_energy_reference = (self.section['section_system']['number_of_electrons'] is not None)
        HOMO = section['x_qe_t_energy_reference_highest_occupied']
        if HOMO is not None:
            if len(HOMO)>1:
                LOGGER.error('more than one value for HOMO: %s', str(HOMO))
            had_energy_reference = True
            backend.addArrayValues('energy_reference_highest_occupied', np.asarray([HOMO[0]]))
        LUMO = section['x_qe_t_energy_reference_lowest_unoccupied']
        if LUMO is not None:
            if len(LUMO)>1:
                LOGGER.error('more than one value for LUMO: %s', str(LUMO))
            backend.addArrayValues('energy_reference_lowest_unoccupied', np.asarray([LUMO[0]]))
        E_Fermi_up = section['x_qe_t_energy_reference_fermi_up']
        E_Fermi = section['x_qe_t_energy_reference_fermi']
        if E_Fermi_up is not None:
            if len(E_Fermi_up)>1:
                LOGGER.error('more than one value for E_Fermi_up: %s', str(E_Fermi_up))
            backend.addArrayValues('energy_reference_fermi', np.asarray([
                E_Fermi_up[0], section['x_qe_t_energy_reference_fermi_down'][0]
            ]))
            had_energy_reference = True
        elif E_Fermi is not None:
            if len(E_Fermi)>1:
                LOGGER.error('more than one value for E_Fermi: %s', str(E_Fermi))
            backend.addArrayValues('energy_reference_fermi', np.asarray([E_Fermi[0]]))
            had_energy_reference = True
        if not (had_energy_reference):
            LOGGER.error("Neither HOMO, Fermi energy nor number of electrons are defined")
        # setup section_system for next scf
        if section['x_qe_t_md_atom_labels'] or section['x_qe_t_md_vec_a_x']:
            next_system_gIndex = backend.openSection('section_system')
            # we cannot simply do
            #   backend.addValue('x_qe_t_vec_a_x', section['x_qe_t_md_vec_a_x'])
            # as this adds an outer list with one element (WTF)
            new_system = {}
            if section['x_qe_t_md_vec_a_units']:
                # we got new cell vectors
                new_system['x_qe_t_vec_a_units'] = section['x_qe_t_md_vec_a_units']
                new_system['x_qe_t_vec_a_x'] = section['x_qe_t_md_vec_a_x']
                new_system['x_qe_t_vec_a_y'] = section['x_qe_t_md_vec_a_y']
                new_system['x_qe_t_vec_a_z'] = section['x_qe_t_md_vec_a_z']
            if section['x_qe_t_md_atom_labels']:
                # we got new atom positions and labels
                new_system['x_qe_t_atom_labels'] = section['x_qe_t_md_atom_labels']
                new_system['x_qe_t_atpos_units'] = section['x_qe_t_md_atom_positions_units']
                new_system['x_qe_t_atpos_x'] = section['x_qe_t_md_atom_positions_x']
                new_system['x_qe_t_atpos_y'] = section['x_qe_t_md_atom_positions_y']
                new_system['x_qe_t_atpos_z'] = section['x_qe_t_md_atom_positions_z']
            if section['x_qe_t_md_k_info_vec_x']:
                # we got new atom positions and labels
                new_system['x_qe_t_k_info_wk'] = section['x_qe_t_md_k_info_wk']
                new_system['x_qe_t_k_info_vec_x'] = section['x_qe_t_md_k_info_vec_x']
                new_system['x_qe_t_k_info_vec_y'] = section['x_qe_t_md_k_info_vec_y']
                new_system['x_qe_t_k_info_vec_z'] = section['x_qe_t_md_k_info_vec_z']
            if section['x_qe_t_md_k_info_ik']:
                new_system['x_qe_t_k_info_ik'] = section['x_qe_t_md_k_info_ik']
            for target, data in new_system.items():
                for val in data:
                    backend.addValue(target, val)
            backend.closeSection('section_system', next_system_gIndex)

    def onClose_section_scf_iteration(
            self, backend, gIndex, section):
        """trigger called when section_scf_iteration is closed"""
        if section["x_qe_t_iter_mpersite_idx"] is not None:
            backend.addArrayValues("x_qe_iter_mpersite_idx", np.asarray(section["x_qe_t_iter_mpersite_idx"]))
            backend.addArrayValues("x_qe_iter_mpersite_charge", np.asarray(section["x_qe_t_iter_mpersite_charge"]))
            backend.addArrayValues("x_qe_iter_mpersite_magn", np.asarray(section["x_qe_t_iter_mpersite_magn"]))
            backend.addArrayValues("x_qe_iter_mpersite_constr", np.asarray(section["x_qe_t_iter_mpersite_constr"]))
        self.tmp['last_iteration'] = section['x_qe_iteration_number'][-1]

    def create_section_eigenvalues(self, backend, src_sec):
        if src_sec['x_qe_t_k_x'] is None or len(src_sec['x_qe_t_k_x']) < 1:
            LOGGER.error("no k-points, not creating section_eigenvalues!")
            return
        sec_eigenvalues_gIndex = backend.openSection('section_eigenvalues')
        # prepare numpy arrays
        k_energies = np.array([self.tmp['k_energies']], dtype=np.float64)
        k_energies = unit_conversion.convert_unit(k_energies, 'eV')
        k_occupations = np.array([self.tmp['k_occupations']], dtype=np.float64)
        npw = None
        if src_sec['x_qe_t_k_pw'] is not None:
            npw = np.array(src_sec['x_qe_t_k_pw'])
        k_point_cartesian = np.array([
            src_sec['x_qe_t_k_x'], src_sec['x_qe_t_k_y'], src_sec['x_qe_t_k_z']
        ], dtype=np.float64).T
        # check if we are dealing with spin-polarized data
        #   QE represents this as 2*k-points, with repeating coordinates
        nk = len(k_point_cartesian)
        if (nk & 1):
            # odd number of k points cannot describe spin-polarized case
            nspin = 1
        else:
            kd = (k_point_cartesian[0:nk//2,:] - k_point_cartesian[nk//2:,:])
            kd_len = np.sqrt(np.einsum('ij,ij->i', kd, kd))
            # LOGGER.error("kdl: %s", str(kd_len))
            # difference in k coordinates
            # values in 1/m are in the order of 1e+10, so threshold is sufficient
            if np.all(kd_len<0.1):
                # first half of k coodinates same as second half -> spin
                # polarized case
                nspin = 2
            else:
                nspin = 1
        # sanity checks
        len_kspin = len(self.tmp['kspin'])
        if (nspin == 2) and (len_kspin!=2):
            raise Exception("nspin=2 and kspin!=2")
        elif (nspin == 1) and (len_kspin>1):
            raise Exception("nspin=1 and kspin>1")
        elif (len_kspin>2):
            raise Exception("total fuckup: kspin=%d" % (len_kspin))
        # transform spin-polarized case
        if nspin == 2:
            k_point_cartesian[0:nk//2,:]
            if npw is not None:
                npw = npw[0:nk//2]
            # put spin channel into first dimension
            k_energies = np.concatenate((
                k_energies[:,0:nk//2,:],
                k_energies[:,nk//2:,:]), axis=0)
            k_occupations = np.concatenate((
                k_occupations[:,0:nk//2,:],
                k_occupations[:,nk//2:,:]), axis=0)
        # k-points are in cartesian, but metaInfo specifies crystal
        k_point_crystal = self.bmat_inv.dot(k_point_cartesian.T).T
        # emit data
        if npw is not None:
            backend.addArrayValues('x_qe_eigenvalues_number_of_planewaves', npw)
        if k_occupations.shape[2] > 0:
            backend.addArrayValues('eigenvalues_occupation', k_occupations)
        backend.addArrayValues('eigenvalues_kpoints', k_point_crystal)
        backend.addArrayValues('eigenvalues_values', k_energies)
        backend.closeSection('section_eigenvalues', sec_eigenvalues_gIndex)

    def onClose_section_system(self, backend, gIndex, section):
        old_system = self.section.get('section_system', None)
        if old_system is not None and self.tmp.get('md_relax', None) is None:
            raise Exception('encountered new section_system without knowing why')
        if section['number_of_atoms']:
            self.number_of_atoms = section['number_of_atoms'][-1]

        # store direct lattice matrix and inverse for transformation crystal <-> cartesian
        if section['x_qe_t_vec_a_x'] is not None:
            # there may have been more than one instance of the bravais matrix
            # in different units (seen in relax/MD)
            self.amat = np.array([
                section['x_qe_t_vec_a_x'][-3:], section['x_qe_t_vec_a_y'][-3:], section['x_qe_t_vec_a_z'][-3:],
            ], dtype=np.float64).T
            # explicit unit conversion
            amat_units = section['x_qe_t_vec_a_units'][-1]
            if amat_units == 'a_0' or amat_units == 'alat':
                self.amat = unit_conversion.convert_unit(self.amat, 'usrAlat')
            elif amat_units == 'bohr':
                self.amat = unit_conversion.convert_unit(self.amat, 'bohr')
            elif amat_units == 'angstrom':
                self.amat = unit_conversion.convert_unit(self.amat, 'angstrom')
            else:
                raise RuntimeError("unknown amat_units: %s" % (amat_units))

            # store inverse for transformation cartesian -> crystal
            try:
                self.amat_inv = np.linalg.inv(self.amat)
            except np.linalg.linalg.LinAlgError:
                raise Exception("error inverting bravais matrix " + str(self.amat))
            LOGGER.info('NewCell')
        elif old_system is not None:
            # we did not get new cell vectors, recycle old ones
            LOGGER.info('OldCell')
        else:
            raise Exception("missing bravais vectors")
        backend.addArrayValues('simulation_cell', self.amat)
        self.tmp['have_new_system'] = True

        # store reciprocal lattice matrix and inverse for transformation crystal <-> cartesian
        if section['x_qe_t_vec_b_x'] is not None:
            self.bmat = np.array([
                section['x_qe_t_vec_b_x'], section['x_qe_t_vec_b_y'], section['x_qe_t_vec_b_z'],
            ], dtype=np.float64).T
            # store inverse for transformation cartesian -> crystal
            try:
                self.bmat_inv = np.linalg.inv(self.bmat)
            except np.linalg.linalg.LinAlgError:
                raise Exception("error inverting reciprocal cell matrix")
        elif section['x_qe_t_vec_a_x'] is not None:
            # we got new lattice vectors, but no reciprocal ones, calculate
            # on-the-fly
            LOGGER.info('calculating bmat on the fly from amat')
            abmat = np.zeros((3,3), dtype=np.float64)
            abmat[0] = np.cross(self.amat[1],self.amat[2])
            abmat[1] = np.cross(self.amat[2],self.amat[0])
            abmat[2] = np.cross(self.amat[0],self.amat[1])
            abmat *= 2*math.pi / np.dot(abmat[0],self.amat[0])
            self.bmat = abmat
            # store inverse for transformation cartesian -> crystal
            try:
                self.bmat_inv = np.linalg.inv(self.bmat)
            except np.linalg.linalg.LinAlgError:
                raise Exception("error inverting reciprocal cell matrix")
        elif old_system is not None:
            # keep what we had
            pass
        else:
            raise Exception("missing reciprocal cell vectors")
        backend.addArrayValues('x_qe_reciprocal_cell', self.bmat)
        # atom positions
        if section['x_qe_t_atpos_x'] is not None:
            # there may have been more than one instance of the atom positions
            # in different units (seen in relax/MD)
            atpos = np.array([
                section['x_qe_t_atpos_x'][-self.number_of_atoms:],
                section['x_qe_t_atpos_y'][-self.number_of_atoms:],
                section['x_qe_t_atpos_z'][-self.number_of_atoms:]
            ], dtype=np.float64).T
            atpos_units = section['x_qe_t_atpos_units'][-1]
            if atpos_units == 'a_0' or atpos_units == 'alat':
                atpos_cart = unit_conversion.convert_unit(atpos, 'usrAlat')
            elif atpos_units == 'bohr':
                atpos_cart = unit_conversion.convert_unit(atpos, 'bohr')
            elif atpos_units == 'angstrom':
                atpos_cart = unit_conversion.convert_unit(atpos, 'angstrom')
            elif atpos_units == 'crystal' or atpos_units == 'cryst. coord.':
                atpos_cart = self.amat.dot(atpos.T).T
            else:
                raise RuntimeError("unknown atpos_units: %s" % (atpos_units))
            LOGGER.info('NewAtpos')
        elif old_system is not None and old_system['atom_positions'] is not None:
            atpos_cart = old_system['atom_positions'][-1]
            LOGGER.info('OldAtpos')
        else:
            raise Exception("missing atom positions")
        backend.addArrayValues('atom_positions',atpos_cart)

        if section['x_qe_t_atom_labels'] is not None:
            backend.addArrayValues('atom_labels',np.asarray(section['x_qe_t_atom_labels'][-self.number_of_atoms:]))
        elif old_system is not None:
            backend.addArrayValues('atom_labels',old_system['atom_labels'][-1])
        else:
            raise Exception("missing atom labels")

        if section['x_qe_t_atom_idx'] is not None:
            backend.addArrayValues('x_qe_atom_idx', np.asarray(section['x_qe_t_atom_idx'][-self.number_of_atoms:]))
        elif old_system is not None:
            backend.addArrayValues('x_qe_atom_idx', old_system['x_qe_atom_idx'][-1])
        else:
            raise Exception("missing x_qe_atom_idx")

        if section['x_qe_t_starting_magnetization_species'] is not None:
            # build dict with per-species magnetization
            sp_magn = {}
            for (label, magn) in zip(
                    section['x_qe_t_starting_magnetization_species'],
                    section['x_qe_t_starting_magnetization_value']):
                sp_magn[label] = magn
            at_magn = []
            # transform to per-atom magnetization
            for label in section['atom_labels'][-1]:
                at_magn.append(sp_magn[label])
            backend.addArrayValues('x_qe_atom_starting_magnetization',np.array(at_magn))

        if section['x_qe_t_celldm'] is not None:
            celldm_joint = " ".join(section['x_qe_t_celldm'])
            celldm = [None, None, None, None, None, None]
            for match in re.findall(r"celldm\(\s*(\d+)\s*\)\s*=\s*(" + RE_f + r")", celldm_joint):
                celldm[int(match[0])-1] = valueForStrValue(match[1], 'f')
            celldm[0] = self.alat
            backend.addArrayValues('x_qe_celldm', np.array(celldm))

        if section['x_qe_t_k_info_vec_x'] is not None:
            backend.addArrayValues('x_qe_k_info_vec', np.array([
                section['x_qe_t_k_info_vec_x'], section['x_qe_t_k_info_vec_y'], section['x_qe_t_k_info_vec_z']
            ]).T)
        elif old_system is not None and old_system['x_qe_k_info_vec'] is not None:
            # unless espresso explicitly writes new k-points, sampling is kept fixed
            backend.addArrayValues('x_qe_k_info_vec', old_system['x_qe_k_info_vec'][-1])
        else:
            LOGGER.error("No K-point info found in output")

        if section['x_qe_t_k_info_ik'] is not None:
            backend.addArrayValues('x_qe_k_info_ik', np.array(section['x_qe_t_k_info_ik']))
        elif old_system is not None and old_system['x_qe_k_info_ik'] is not None:
            # unless espresso explicitly writes new k-points, sampling is kept fixed
            backend.addArrayValues('x_qe_k_info_ik', old_system['x_qe_k_info_ik'][-1])
        else:
            LOGGER.error("No K-point index info found in output")

        if section['x_qe_t_k_info_wk'] is not None:
            backend.addArrayValues('x_qe_k_info_wk', np.array(section['x_qe_t_k_info_wk']))
        elif old_system is not None and old_system['x_qe_k_info_wk'] is not None:
            # unless espresso explicitly writes new k-points, sampling is kept fixed
            backend.addArrayValues('x_qe_k_info_wk', old_system['x_qe_k_info_wk'][-1])
        else:
            LOGGER.error("No K-point weight info found in output")


        if section['x_qe_t_dense_FFT_grid_x'] is not None:
            backend.addArrayValues('x_qe_dense_FFT_grid', np.array([
                section['x_qe_t_dense_FFT_grid_x'], section['x_qe_t_dense_FFT_grid_y'], section['x_qe_t_dense_FFT_grid_z']
            ]).T)
        elif old_system is not None and old_system['x_qe_dense_FFT_grid'] is not None:
            # unless espresso explicitly writes new FFT grid info, sampling is kept fixed
            backend.addArrayValues('x_qe_dense_FFT_grid', old_system['x_qe_dense_FFT_grid'][-1])
        else:
            LOGGER.warning("No dense FFT grid info found in output")

        if section['x_qe_t_smooth_FFT_grid_x'] is not None:
            backend.addArrayValues('x_qe_smooth_FFT_grid', np.array([
                section['x_qe_t_smooth_FFT_grid_x'], section['x_qe_t_smooth_FFT_grid_y'], section['x_qe_t_smooth_FFT_grid_z']
            ]).T)
        elif old_system is not None and old_system['x_qe_smooth_FFT_grid'] is not None:
            backend.addArrayValues('x_qe_smooth_FFT_grid', old_system['x_qe_smooth_FFT_grid'][-1])

        if section['x_qe_t_vec_supercell_x'] is not None:
            backend.addArrayValues('x_qe_vec_supercell', np.array([
                section['x_qe_t_vec_supercell_x'], section['x_qe_t_vec_supercell_y'], section['x_qe_t_vec_supercell_z']
            ]).T)

        if section['x_qe_t_number_of_electrons_up'] is not None:
            # spin polarized case, with explicit up/down electrons
            if len(section['x_qe_t_number_of_electrons_up'])>1:
                LOGGER.error("got multiple nelec_up: %s", str(
                    section['x_qe_t_number_of_electrons_up']))
            backend.addArrayValues('number_of_electrons', np.array([
                section['x_qe_t_number_of_electrons_up'][0],
                section['x_qe_t_number_of_electrons_down'][0]
            ]))
        elif section['x_qe_t_number_of_electrons'] is not None:
            if len(section['x_qe_t_number_of_electrons'])!=1:
                LOGGER.error("got wrong nelec: %s", str(section['x_qe_t_number_of_electrons']))
            backend.addArrayValues('number_of_electrons', np.array([
                section['x_qe_t_number_of_electrons'][-1]
            ]))
        elif old_system is not None and old_system['number_of_electrons'] is not None:
            backend.addArrayValues('number_of_electrons', old_system['number_of_electrons'][-1])
        else:
            raise Exception("missing info about number of electrons in system")

        backend.addArrayValues('configuration_periodic_dimensions', np.asarray([True, True, True]))
        self.sectionIdx['section_system'] = gIndex
        self.section['section_system'] = section

    def onOpen_section_run(self, backend, gIndex, section):
        """trigger called when section_single_configuration_calculation
        is closed"""
        self.initialize_values()
        self.tmp.pop('x_qe_t_profile_caller', None)
        self.tmp.pop('x_qe_t_profile_category', None)
        # manually open header sections, closed at the beginning of scf
        for sec in self.header_sections():
            gIndex = backend.openSection(sec)
            self.openSectionIdx[sec] = gIndex

    def adHoc_final_scf_MD(self, parser):
        """final SCF calculation in VC-relax runs needs open header sections"""
        # manually open header sections, closed at the beginning of scf
        for sec in self.header_sections():
            gIndex = parser.backend.openSection(sec)
            self.openSectionIdx[sec] = gIndex

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
            backend.addArrayValues('x_qe_profile_category', np.asarray(
                section['x_qe_t_profile_category_list']))
            backend.addArrayValues('x_qe_profile_caller', np.asarray(
                section['x_qe_t_profile_caller_list']))

        frames = self.tmp.get('frames', None)
        if frames:
            sampling_method_gIndex = backend.openSection('section_sampling_method')
            backend.addValue('sampling_method', QE_MD_RELAX_SAMPLING_METHOD[self.tmp['md_relax']])
            backend.closeSection('section_sampling_method', sampling_method_gIndex)
            frame_sequence_gIndex = backend.openSection('section_frame_sequence')
            backend.addValue('frame_sequence_to_sampling_ref', sampling_method_gIndex)
            backend.addArrayValues('frame_sequence_local_frames_ref', np.array(self.tmp['frames']))
            backend.closeSection('section_frame_sequence', frame_sequence_gIndex)

    def appendToTmp(self, tmpname, value):
        self.tmp[tmpname] += value

    def setTmp(self, tmpname, value):
        self.tmp[tmpname] = value

    def setTmpUnlessExists(self, tmpname, value):
        if self.tmp.get(tmpname,None) is None:
            self.tmp[tmpname] = value

    def popTmp(self, tmpname, fallback=None):
        return self.tmp.pop(tmpname, fallback)

    def adHoc_dispersion_correction_values(self, parser):
        if 'dispersion_correction' not in self.tmp:
            self.tmp['dispersion_correction'] = {}
        self.tmp['dispersion_correction'][parser.lastMatch['x_qe_t_species_dispersion_correction_label']] = {
            'x_qe_dispersion_correction_vdw_radius': parser.lastMatch['x_qe_t_species_dispersion_correction_vdw_radius'],
            'x_qe_dispersion_correction_C6': parser.lastMatch['x_qe_t_species_dispersion_correction_C6'],
        }

    def adHoc_alat(self, parser):
        alat = parser.lastMatch['x_qe_alat']
        self.alat = parser.lastMatch['x_qe_alat']
        unit_conversion.register_userdefined_quantity('usrAlat', 'm', alat)
        unit_conversion.register_userdefined_quantity('usrTpiba', '1/m', 2*math.pi/alat)

    def adHoc_pp_renorm(self, parser):
        self.cache_t_pp_renorm_wfc[parser.lastMatch[
            'x_qe_t_pp_renormalized_filename']] = parser.lastMatch['x_qe_t_pp_renormalized_wfc']

    def adHoc_bands_spin(self, parser):
        k_x = self.section['single_configuration_calculation']['x_qe_t_k_x']
        if k_x is None:
            nk_current = 0
        else:
            nk_current = len(k_x)
        self.tmp['kspin'][parser.lastMatch['x_qe_t_spin_channel'].lower()] = nk_current

    def adHoc_profiling_complete(self, parser):
        # we have effectively 3 SMs in one, split between the 3 cases
        if parser.lastMatch.get('x_qe_t_profile_caller', None) is not None:
            # caller info
            self.setTmp('x_qe_t_profile_caller', parser.lastMatch['x_qe_t_profile_caller'])
        elif parser.lastMatch.get('x_qe_t_profile_category', None) is not None:
            # category info
            self.setTmp('x_qe_t_profile_category', parser.lastMatch['x_qe_t_profile_category'])
            self.tmp.pop('x_qe_t_profile_caller', None)
        else:
            # actual profiling data
            if parser.lastMatch.get('x_qe_t_profile_cputime', None) is None:
                parser.backend.addValue('x_qe_t_profile_cputime', QeC.NAN)
            if parser.lastMatch.get('x_qe_t_profile_walltime', None) is None:
                parser.backend.addValue('x_qe_t_profile_walltime', QeC.NAN)
            if parser.lastMatch.get('x_qe_t_profile_ncalls', None) is None:
                parser.backend.addValue('x_qe_t_profile_ncalls', QeC.NAN)
            parser.backend.addValue('x_qe_t_profile_caller_list', self.tmp.get('x_qe_t_profile_caller', ''))
            parser.backend.addValue('x_qe_t_profile_category_list', self.tmp.get('x_qe_t_profile_category', ''))

    def SMs_summaryf90(self):
        return [
            SM(name='ibrav', required=True,
               # first line printed by summary.f90
               startReStr=r"\s*bravais-lattice index\s*=\s*(?P<x_qe_ibrav>" + RE_i + r")\s*$",
               subMatchers=[
                   # other info coming from/triggered by summary.f90
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
                      startReStr=(r"\s*number of electrons\s*=\s*(?P<x_qe_t_number_of_electrons>" + RE_f +
                                  r")(?:\s*\(up:\s*(?P<x_qe_t_number_of_electrons_up>" + RE_f +
                                  r")\s*,\s*down:\s*(?P<x_qe_t_number_of_electrons_down>" + RE_f + ")\s*\))?\s*$"),
                   ),
                   SM(name='nbnd', required=True,
                       startReStr=r"\s*number of Kohn-Sham states\s*=\s*(?P<x_qe_number_of_states>" + RE_i + r")\s*$"
                   ),
                   SM(name='ecutwfc', required=True, sections=['section_basis_set_cell_dependent'],
                      startReStr=r"\s*kinetic-energy cutoff\s*=\s*(?P<basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry\s*$",
                      adHoc=lambda p: self.setTmp('method_basis_set_kind', 'wavefunction'),
                   ),
                   SM(name='ecut_density', required=True, sections=['section_basis_set_cell_dependent'],
                      startReStr=r"\s*charge density cutoff\s*=\s*(?P<basis_set_planewave_cutoff__rydberg>" + RE_f + r")\s*Ry\s*$",
                      adHoc=lambda p: self.setTmp('method_basis_set_kind', 'density'),
                   ),
                   SM(name='ecutfock',
                      startReStr=r"\s*cutoff for Fock operator\s*=\s*(?P<x_qe_fock_operator_cutoff__rydberg>" + RE_f + r")\s*Ry\s*$"
                   ),
                   SM(name='convergence_threshold',
                      startReStr=r"\s*convergence threshold\s*=\s*(?P<scf_threshold_energy_change__rydberg>" + RE_f + r")\s*$",
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
                   SM(name='exx_fraction',
                      startReStr=r"\s*EXX-fraction\s*=\s*(?P<x_qe_t_exact_exchange_fraction>" + RE_f + r")\s*$",
                   ),
                   SM(name='nstep',
                      startReStr=r"\s*nstep\s*=\s*(?P<x_qe_md_max_steps>" + RE_i + r")\s*$",
                   ),
                   SM(name='spin_orbit_mode',
                      startReStr=r"\s*(?P<x_qe_t_spin_orbit_magn>.*?)\s*calculation\s*(?P<x_qe_t_spin_orbit_mode>with(?:out)?)\s*spin-orbit\s*$",
                   ),
                   SM(name='berry_efield',
                      startReStr=r"\s*Using Berry phase electric field\s*$",
                      subMatchers=[
                          SM(name='berry_efield_direction',
                             startReStr=r"\s*Direction\s*:\s*(?P<x_qe_berry_efield_direction>\d+)\s*$",
                          ),
                          SM(name='berry_efield_intensity',
                             # Ry unit is not printed in 4.0
                             startReStr=(r"\s*Intensity \((?:Ry\s*)?a.u.\)\s*:\s*(?P<x_qe_berry_efield_intensity__rydberg>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name='berry_efield_strings',
                             startReStr=(r"\s*Strings composed by\s*:\s*(?P<x_qe_berry_efield_strings_nk>" + RE_i +
                                         r")\s*k-points\s*$"),
                          ),
                          SM(name='berry_efield_niter',
                             startReStr=(r"\s*Number of iterative cycles:\s*(?P<x_qe_berry_efield_niter>" + RE_i +
                                         r")\s*$"),
                          ),
                      ],
                      adHoc=lambda p: p.backend.addValue('x_qe_berry_efield', True),
                   ),
                   SM(name='assume_isolated',
                      startReStr=r"\s*Assuming isolated system,\s*(?P<x_qe_isolated_system_method>.*?)\s*method",
                   ),
                   SM(name='celldm', repeats = True,
                      startReStr=r"(?P<x_qe_t_celldm>(?:\s*celldm\(\d+\)\s*=\s*" + RE_f + r")+)\s*$",
                   ),
                   SM(name='simulation_cell',
                      startReStr=r"\s*crystal axes: \(cart. coord.\s*in units of (?P<x_qe_t_vec_a_units>a_0|alat)\s*\)\s*$",
                      subMatchers=[
                          SM(name='cell_vec_a', repeats=True,
                             startReStr=r"\s*a\(\d\)\s*=\s*\(\s*" + QeC.re_vec('x_qe_t_vec_a') + r"\s*\)\s*$",
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
                          SM(name='pp_integral_directions',
                             startReStr=(r"\s*Setup to integrate on\s*" +
                                         r"(?P<x_qe_t_pp_integral_ndirections>\d+)\s+directions:\s*" +
                                         r"integral exact up to l =\s*(?P<x_qe_t_pp_integral_lmax_exact>\d+)\s*$"),
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
                   SM(name='vdw_kernel_table_file',
                      startReStr=r"\s*vdW kernel table read from file .*?\s*$",
                      adHoc=lambda p: LOGGER.info("do sth with vdW kernel table file"),
                      subMatchers=[
                          SM(name='vdw_kernel_table_md5sum',
                             startReStr=r"\s*MD5 check sum:\s*\S+\s*$",
                          ),
                      ],
                   ),
                   SM(name='pp_atom_kind_map',
                      startReStr=r"\s*atomic species\s+valence\s+mass\s+pseudopotential\s*$",
                      subMatchers=[
                          SM(name='atom_kind', repeats=True,
                             startReStr=r"\s*(?P<method_atom_kind_label>\S+)\s+(?P<method_atom_kind_explicit_electrons>" + RE_f + r")" +
                                        r"\s+(?P<x_qe_kind_mass>" + RE_f + r")\s+(?P<x_qe_pp_label>[^\s\(]+)\s*" +
                                        r"\(\s*(?P<x_qe_pp_weight>" + RE_f + r")\s*\)\s*$",
                             sections=['section_method_atom_kind'],
                          ),
                      ],
                   ),
                   SM(name='starting_magnetization',
                      startReStr=r"\s*Starting magnetic structure\s*$",
                      subMatchers=[
                          SM(name='starting_magnetization_header',
                             startReStr=r"\s*atomic species\s+magnetization\s*$",
                             subMatchers=[
                                 SM(name='starting_magnetization_data', repeats=True,
                                    startReStr=(r"\s*(?P<x_qe_t_starting_magnetization_species>\S+)\s*" +
                                                r"(?P<x_qe_t_starting_magnetization_value>" + RE_f + r")\s*$"),
                                 ),
                             ],
                          ),
                      ],
                   ),
                   SM(name='cell_mass',
                      startReStr=r"\s*cell mass =\s*(?P<x_qe_md_cell_mass>" + RE_f + r")\s*AMU/\(a\.u\.\)\^2",
                   ),
                   SM(name='nsymm',
                      startReStr=(r"\s*(?P<x_qe_nsymm>\d+)\s*Sym\.\s*Ops\.\s*[\(,]\s*(?P<x_qe_t_symm_inversion>\S+) inversion\s*[\),]\s*(?:found)?\s*"
                                  r"(?:\(\s*(?P<x_qe_nsymm_with_fractional_translation>\d+)\s*have fractional translation\s*\))?\s*$"),
                      adHoc=lambda p: p.backend.addValue('x_qe_symm_inversion', (p.lastMatch['x_qe_t_symm_inversion'] == 'with')),
                      subMatchers=[
                          SM(name='nsymm_ignored',
                             startReStr=r"\s*\(note:\s*(?P<x_qe_nsymm_ignored>\d+)\s*additional sym.ops. were found but ignored\s*$",
                             subMatchers=[
                                 SM(name='nsymm_ignored_cIgnore', coverageIgnore=True,
                                    startReStr=r"\s*their fractional transl?ations are incommensurate with FFT grid\)\s*$",
                                 ),
                             ],
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
                             startReStr=r"\s*site n.     atom                  positions \((?P<x_qe_t_atpos_units>a_0|alat) units\)\s*$",
                             subMatchers=[
                                 SM(name='atom_pos_cart', repeats=True,
                                    startReStr=(
                                        r"\s*(?P<x_qe_t_atom_idx>" + RE_i + r")" +
                                        r"\s+(?P<x_qe_t_atom_labels>\S+)\s+tau\(\s*" + RE_i + "\)\s*"
                                        r"=\s*\(\s*" + QeC.re_vec('x_qe_t_atpos') +
                                        r"\s*\)\s*$"),
                                 ),
                             ],
                          ),
                      ],
                   ),
                   SM(name='atom_pos_cryst_list',
                      # happens in verbose mode
                      startReStr=r"\s*Crystallographic axes\s*$",
                      subMatchers=[
                          SM(name='cryst_heading',
                             startReStr=r"\s*site n.     atom                  positions \((?P<x_qe_t_atpos_units>cryst\. coord\.)\)\s*$",
                             subMatchers=[
                                 SM(name='atom_pos_cyst', repeats=True,
                                    startReStr=(
                                        r"\s*(?P<x_qe_t_atom_idx>" + RE_i + r")" +
                                        r"\s+(?P<x_qe_t_atom_labels>\S+)\s+tau\(\s*" + RE_i + "\)\s*"
                                        r"=\s*\(\s*" + QeC.re_vec('x_qe_t_atpos') +
                                        r"\s*\)\s*$"),
                                 ),
                             ],
                          ),
                      ],
                   ),
                   SM(name='kpoint_info',
                      startReStr=r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*$",
                      subMatchers=self.SMs_kpoints(),
                   ),
                   SM(name='kpoint_info_smearing_old',
                      startReStr=(r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*gaussian broad.\s*\(Ry\)\s*=" +
                                  r"\s*(?P<smearing_width__rydberg>" + RE_f + r")\s*ngauss\s*=" +
                                  r"\s*(?P<x_qe_smearing_ngauss>" + RE_i + ")\s*$"),
                      adHoc=lambda p: p.backend.addValue('smearing_kind', QeC.QE_SMEARING_KIND.get(str(p.lastMatch['x_qe_smearing_ngauss']))),
                      subMatchers=self.SMs_kpoints(),
                   ),
                   SM(name='kpoint_info_smearing_new',
                      startReStr=(r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*(?P<x_qe_smearing_kind>.+?)\s*,\s*width\s*\(Ry\)\s*=" +
                                  r"\s*(?P<smearing_width__rydberg>" + RE_f + r")\s*$"),
                      adHoc=lambda p: p.backend.addValue('smearing_kind', QeC.QE_SMEARING_KIND[p.lastMatch['x_qe_smearing_kind']]),
                      subMatchers=self.SMs_kpoints(),
                   ),
                   SM(name='kpoint_info_tetrahedra',
                      startReStr=r"\s*number of k points=\s*(?P<x_qe_nk>\d+)\s*\((?P<x_qe_smearing_kind>tetrahedron method)\)\s*$",
                      adHoc=lambda p: p.backend.addValue('smearing_kind', QeC.QE_SMEARING_KIND[p.lastMatch['x_qe_smearing_kind']]),
                      subMatchers=self.SMs_kpoints(),
                   ),
                   SM(name='kpoint_info_nokpoints',
                      startReStr=r"\s*(?P<x_qe_warning>Number of k-points >= 100: set verbosity='high' to print them\.)\s*$",
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
                          SM(name='input_occupations_spin', repeats=True,
                             startReStr=r"\s*Spin-(?P<x_qe_t_input_occupations_spin>up|down)\s*$",
                             adHoc=self.adHoc_input_occupations_spin,
                             subMatchers=[
                                 SM(name='input_occupations_occupations', repeats=True,
                                     startReStr=r'\s*(?P<x_qe_t_input_occupations>(?:\s*' + RE_f + ')+\s*$)',
                                     adHoc=lambda p: self.tmp['occ_vals'].extend(cRE_f.findall(p.lastMatch['x_qe_t_input_occupations'])),
                                 ),
                             ],
                          ),
                          SM(name='input_occupations_occupations', repeats=True,
                              startReStr=r'\s*(?P<x_qe_t_input_occupations>(?:\s*' + RE_f + ')+\s*$)',
                              adHoc=lambda p: self.tmp['occ_vals'].extend(cRE_f.findall(p.lastMatch['x_qe_t_input_occupations']))
                          ),
                      ],
                   ),
               ],
            ),
        ]

    def adHoc_add_kbnd_occ_storage(self, parser):
        self.tmp['k_energies'].append([])
        self.tmp['k_occupations'].append([])

    def SMs_bands(self, suffix=''):
        return [
            SM(name='bands' + suffix, repeats=True,
                startReStr=(r'\s*k\s*=\s*' + QeC.re_vec('x_qe_t_k', 'usrTpiba', '\s*') +
                            r'(?:\s*\(\s*(?P<x_qe_t_k_pw>' + RE_i + r')\s*PWs\s*\))?' +
                            r'\s*band(?:s|\s+energies)\s*\(\s*[eE][vV]\s*\)\s*:?\s*$'),
                # create new empty list for this k point's eigenvalues
                adHoc=self.adHoc_add_kbnd_occ_storage,
                subMatchers=[
                    SM(name='kbnd' + suffix, repeats=True,
                        startReStr=r'\s*(?P<x_qe_t_k_point_energies>(?:\s*' + RE_f + ')+\s*$)',
                        # extend list by eigenvalues in this line
                        adHoc=lambda p: self.tmp['k_energies'][-1].extend(cRE_f.findall(p.lastMatch['x_qe_t_k_point_energies'])),
                    ),
                    SM(name='occhead' + suffix,
                       # only printed in verbose mode
                       startReStr=r"\s*occupation\s*numbers\s*$",
                       subMatchers=[
                           SM(name='kbndocc' + suffix, repeats=True,
                               startReStr=r'\s*(?P<x_qe_t_k_point_energies>(?:\s*' + RE_f + ')+\s*$)',
                               # extend list by eigenvalues in this line
                               adHoc=lambda p: self.tmp['k_occupations'][-1].extend(cRE_f.findall(p.lastMatch['x_qe_t_k_point_energies'])),
                           ),
                       ],
                    ),
                ],
            ),
        ]

    def SMs_read_input_file(self, suffix=''):
        return [
            SM(name='read_input_file' + suffix,
               # first case: input is read from stdin...
               startReStr=r"\s*(?:Waiting for input\.\.\.|Reading input from\s*(?P<x_qe_input_filename>.*?))\s*$",
               subMatchers=[
                   # debug msg appearing in some old QE versions
                   SM(name="ignore_ierr_debug", repeats=True, coverageIgnore=True,
                      startReStr=r"\s*ierr\s*=\s*\d+",
                   ),
                   SM(name="warning_input_card_ignored" + suffix, repeats=True,
                       startReStr=r"\s*Warning:\s*(?P<x_qe_warning>card.*ignored)\s*$",
                   ),
                   SM(name="msg_read_cards" + suffix, repeats=True,
                      startReStr=r"\s*Message from routine read_cards\s*:\s*$",
                      subMatchers=[
                          SM(name="warning_read_cards" + suffix,
                             startReStr=r"\s*(?P<x_qe_warning>DEPRECATED: no units specified in CELL_PARAMETERS card)\s*$"
                          ),
                      ],
                   ),
               ],
            ),
        ]

    def SMs_kpoints(self, suffix='',
                    coordinates=r"cart. coord. in units 2pi/(?:alat|a_0)",
                    units='usrTpiba', target='x_qe_t_k_info'):
        return [
            SM(name="kpoint_heading" + suffix,
               startReStr=r"\s*" + coordinates + r"\s*$",
               # stop on the first empty line. keeps k-in-fractional-coords (only printed in high-verbosity mode) from interfering
               endReStr=r"\s*$",
               subMatchers=[
                   SM(name="kpoint_kpoints" + suffix, repeats=True,
                      startReStr=(
                          r"\s*k\(\s*(?P<" + target + "_ik>\d+)\s*\)\s*=\s*\(\s*" + QeC.re_vec(target + "_vec", units) +
                          r"\s*\),\s*wk\s*=\s*(?P<" + target + "_wk>" + RE_f + r")\s*$"),
                   ),
               ],
            ),
        ]

    def SMs_md_system_new(self, suffix='', eatEndFinalSubMatchers=[]):
        return [
            SM(name="eat_start_final" + suffix,
               startReStr=r"\s*Begin final coordinates\s*$"
            ),
            SM(name="new_cell_volume" + suffix,
               startReStr=(r"\s*new unit-cell volume =\s*" +
                           r"(?P<x_qe_t_md_new_volume__bohr3>" + RE_f + r") a\.u\.\^3" +
                           r"\s*\(\s*" + RE_f + "\s*Ang\^3\s*\)\s*$"),
            ),
            SM(name="cellparam" + suffix,
               startReStr=(r"CELL_PARAMETERS\s*\(\s*(?P<x_qe_t_md_vec_a_units>\S+?)\s*" +
                           r"(?:=\s*(?P<x_qe_t_md_vec_a_alat>" + RE_f + ")\s*)?\)$"),
               # QE use syntax of its _input_ files here
               subMatchers=[
                   SM(name='md_cell_vec_a' + suffix, repeats=True,
                      startReStr=r"\s*" + QeC.re_vec('x_qe_t_md_vec_a') + r"\s*$",
                   ),
               ],
            ),
            SM(name="atpos" + suffix,
               startReStr="ATOMIC_POSITIONS\s*\(\s*(?P<x_qe_t_md_atom_positions_units>\S+)\s*\)$",
               # QE use syntax of its _input_ files here
               subMatchers=[
                   SM(name='atpos_data' + suffix, repeats=True,
                      startReStr=(r"\s*(?P<x_qe_t_md_atom_labels>\S+)\s+" +
                                  QeC.re_vec('x_qe_t_md_atom_positions') +
                                  r"(?:\s+(?P<x_qe_t_md_atom_free_x>\d)\s+(?P<x_qe_t_md_atom_free_y>\d)\s+(?P<x_qe_t_md_atom_free_z>\d))?" +
                                  r"\s*$"),
                      fixedStartValues={
                          'x_qe_t_md_atom_free_x': True,
                          'x_qe_t_md_atom_free_y': True,
                          'x_qe_t_md_atom_free_z': True,
                      },
                   ),
               ],
            ),
            SM(name="eat_end_final" + suffix,
               startReStr=r"\s*End final coordinates\s*$",
               subMatchers=eatEndFinalSubMatchers,
            ),
        ]

    def SMs_vcsmd_system_new(self, suffix='', eatEndFinalSubMatchers=[]):
        return [
            # SM(name="eat_start_final" + suffix,
            #    startReStr=r"\s*Begin final coordinates\s*$"
            # ),
            SM(name="cellparam" + suffix,
               startReStr=r"\s*(?:Final estimate of|new) lattice vectors\s*\(\s*(?P<x_qe_t_md_vec_a_units>.+?)\s*units?\s*\)\s*:?\s*$",
               subMatchers=[
                   SM(name='md_cell_vec_a' + suffix, repeats=True,
                      startReStr=r"\s*" + QeC.re_vec('x_qe_t_md_vec_a') + r"\s*$",
                   ),
                   SM(name="new_cell_volume" + suffix,
                       startReStr=(r"\s*(?:final|new) unit-cell volume =\s*" +
                                  r"(?P<x_qe_t_md_new_volume__bohr3>" + RE_f + r") \(a\.u\.\)\^3"),
                   ),
                   SM(name="new_cell_alat" + suffix,
                      startReStr=r"\s*input alat =\s*" + RE_f + "\s*\(a\.u\.\)",
                   ),
               ],
            ),
            SM(name="atpos" + suffix, repeats=True,
               startReStr="\s*new positions in\s*(?P<x_qe_t_md_atom_positions_units_vcsmd>.*?)\s*$",
               adHoc=lambda p: p.backend.addValue('x_qe_t_md_atom_positions_units', QE_VCSMD_ATPOS_UNITS[p.lastMatch['x_qe_t_md_atom_positions_units_vcsmd']]),
               # QE use syntax of its _input_ files here
               subMatchers=[
                   SM(name='atpos_data' + suffix, repeats=True,
                      startReStr=(r"\s*(?P<x_qe_t_md_atom_labels>\S+)\s+" +
                                  QeC.re_vec('x_qe_t_md_atom_positions') +
                                  r"(?:\s+(?P<x_qe_t_md_atom_free_x>\d)\s+(?P<x_qe_t_md_atom_free_y>\d)\s+(?P<x_qe_t_md_atom_free_z>\d))?" +
                                  r"\s*$"),
                      fixedStartValues={
                          'x_qe_t_md_atom_free_x': True,
                          'x_qe_t_md_atom_free_y': True,
                          'x_qe_t_md_atom_free_z': True,
                      },
                   ),
               ],
            ),
            # SM(name="eat_end_final" + suffix,
            #    startReStr=r"\s*End final coordinates\s*$",
            #    subMatchers=eatEndFinalSubMatchers,
            # ),
        ]

    def SMs_relax_bfgs(self):
        return self.SMs_md_system_new(suffix='extra_SCF') + [
            SM(name="bfgs_info",
               startReStr=r"\s*BFGS Geometry Optimization\s*$",
               adHoc = lambda p: self.setTmp('md_relax', 'BFGS')
            ),
            SM(name="lsda_relax_nonmag",
               startReStr=r"\s*lsda relaxation :  a final configuration with zero\s*$",
               adHoc = lambda p: self.setTmp('extra_SCF', True)
            ),
            SM(name="bfgs_scf_cycles",
               startReStr=r"\s*number of scf cycles\s*=\s*(?P<x_qe_t_md_bfgs_scf_cycles>\d+)\s*$",
               subMatchers=[
                   SM(name="bfgs_steps",
                      startReStr=r"\s*number of bfgs steps\s*=\s*(?P<x_qe_t_md_bfgs_steps>\d+)\s*$",
                   ),
                   SM(name="bfgs_energy_old",
                      startReStr=r"\s*energy\s+old\s*=\s*(?P<x_qe_t_md_bfgs_energy_old__rydberg>" + RE_f + r")\s*Ry\s*$",
                   ),
                   SM(name="bfgs_energy_new",
                      startReStr=r"\s*energy\s+new\s*=\s*(?P<x_qe_t_md_bfgs_energy_new__rydberg>" + RE_f + r")\s*Ry\s*$",
                   ),
                   SM(name="bfgs_enthalpy_old",
                      startReStr=r"\s*enthalpy\s+old\s*=\s*(?P<x_qe_t_md_bfgs_enthalpy_old__rydberg>" + RE_f + r")\s*Ry\s*$",
                   ),
                   SM(name="bfgs_enthalpy_new",
                      startReStr=r"\s*enthalpy\s+new\s*=\s*(?P<x_qe_t_md_bfgs_enthalpy_new__rydberg>" + RE_f + r")\s*Ry\s*$",
                   ),
                   SM(name="bfgs_case",
                      startReStr=r"\s*CASE:\s*(?P<x_qe_t_md_bfgs_case>en.*_new\s*\S+\s*en.*_old)\s*$",
                   ),
                   SM(name="bfgs_warning",
                      startReStr=r"\s*WARNING:\s*(?P<x_qe_warning>bfgs curvature condition failed, Theta=\s*\S+)\s*$",
                   ),
                   SM(name="bfgs_reset",
                      startReStr=r"\s*uphill step: resetting bfgs history\s*$",
                      fixedStartValues={'x_qe_t_md_bfgs_reset': 'uphill'},
                   ),
                   SM(name="bfgs_trust_new",
                      startReStr=r"\s*new trust radius\s*=\s*(?P<x_qe_t_md_bfgs_trust_new__bohr>" + RE_f + r")\s*bohr\s*$",
                   ),
                   SM(name="bfgs_conv_thr_new",
                      startReStr=r"\s*new conv_thr\s*=\s*(?P<x_qe_t_md_bfgs_conv_thr_new__rydberg>" + RE_f + r")\s*Ry\s*$",
                   ),
               ] + self.SMs_md_system_new(suffix='BFGS1'),
            ),
            SM(name="bfgs_converged",
               # we know the actual numbers already
               startReStr=r"\s*bfgs converged in\s+\d+\s+scf cycles and\s+\d+\s+bfgs steps",
               fixedStartValues={ 'x_qe_t_md_bfgs_converged': True},
               subMatchers=[
                   SM(name="bfgs_converged_criteria",
                      startReStr=r"\s*\(criteria:\s*(?P<x_qe_t_md_bfgs_converged_criteria>.*?)\)\s*$",
                   ),
               ],
            ),
            SM(name="bfgs_end",
               startReStr=r"\s*End of BFGS Geometry Optimization",
               subMatchers=[
                   SM(name="bfgs_final_energy",
                      startReStr=(r"\s*Final energy\s*=\s*(?P<x_qe_t_md_bfgs_final_energy__rydberg>" +
                                  RE_f + r")\s*Ry\s*$"),
                   ),
                   SM(name="bfgs_final_enthalpy",
                      startReStr=(r"\s*Final enthalpy\s*=\s*(?P<x_qe_t_md_bfgs_final_enthalpy__rydberg>" +
                                  RE_f + r")\s*Ry\s*$"),
                   ),
               ] + self.SMs_md_system_new(suffix='BFGS2'),
            ),
        ]

    def SMs_molecular_dynamics_end(self, suffix=''):
        return [
            SM(name="md_end" + suffix,
               startReStr=r"\s*End of molecular dynamics calculation\s*$",
               fixedStartValues={ 'x_qe_t_md_end': True },
               subMatchers=[
                   SM(name="md_diffusion_coefficients_header" + suffix,
                      startReStr=r"\s*diffusion coefficients\s*:\s*$",
                      subMatchers=[
                          SM(name="md_diffusion_coeffients" + suffix, repeats=True,
                             startReStr=(r"\s*atom\s*(?P<x_qe_t_md_diffusion_atomidx>\d+)\s*D\s*=\s*" +
                                         r"(?P<x_qe_t_md_diffusion_coefficient__cm2_s>" + RE_f +
                                         r")\s*cm\^2/s\s*$"),
                          ),
                          SM(name="md_diffusion_mean" + suffix,
                             startReStr=(r"\s*<\s*D\s*>\s*=\s*(?P<x_qe_t_md_diffusion_coefficient_mean__cm2_s>" +
                                         RE_f + r")\s*cm\^2/s\s*$"),
                          ),
                      ],
                   ),
               ],
            ),
        ]

    def SMs_molecular_dynamics(self):
        return [
            SM(name="md_info",
               # This is only shown after first SCF
               startReStr=r"\s*Molecular Dynamics Calculation\s*$",
               subMatchers=[
                   SM(name="md_atom", repeats=True,
                      startReStr=(r"\s*mass\s*(?P<x_qe_t_md_atom_mass_label>\S+)\s*=\s*" +
                                  r"(?P<x_qe_t_md_atom_mass_value>" + RE_f + r")\s*$"),
                   ),
                   SM(name="md_timestep",
                      startReStr=(r"\s*Time step\s*=\s*" + RE_f + r"\s*a\.u\.,\s*" +
                                  r"(?P<x_qe_t_md_timestep_size__femtoseconds>" + RE_f +
                                  r")\s*femto-seconds\s*$"),
                   ),
               ],
               adHoc = lambda p: self.setTmp('md_relax', 'molecular_dynamics')
            ),
            SM(name="md_info_relax_damped",
               startReStr=r"\s*Damped Dynamics Calculation\s*$",
               adHoc = lambda p: self.setTmp('md_relax', 'damped_dynamics'),
            ),
            SM(name="md_info_relax_damped",
               startReStr=r"\s*Over-damped Langevin Dynamics Calculation\s*$",
               adHoc = lambda p: self.setTmp('md_relax', 'langevin_overdamped_dynamics'),
            ),
            # older espresso writes this _before_ 'entering dynamics'
            SM(name="md_maxSteps1",
               startReStr=r"\s*The maximum number of steps has been reached.\s*$",
               fixedStartValues={ 'x_qe_t_md_max_steps_reached': True },
               subMatchers=self.SMs_molecular_dynamics_end(suffix='1')
            ),
            SM(name="md_step",
               startReStr=(r"\s*Entering Dynamics:\s*iteration\s*=\s*(?P<x_qe_t_md_iteration>" + RE_i +
                           r")\s*$"),
               subMatchers=[
                   SM(name="md_time",
                      startReStr=(r"\s*time\s*=\s*(?P<x_qe_t_md_time__picoseconds>" + RE_f +
                                  r")\s*pico-seconds"),
                   ),
                   SM(name="md_nat2_distance",
                      startReStr=r"\s*DISTANCE\s*=\s*(?P<x_qe_t_new_nat2_distance__bohr>" + RE_f + r")\s*$"
                   ),
                   SM(name="md_projected_velocity",
                       startReStr=r"\s*<vel\(dt\)\|acc\(dt\)>\s*=\s*(?P<x_qe_t_projected_velocity>" + RE_f + r")\s*$",
                   ),
               ] + self.SMs_md_system_new(suffix='MD') + [
                   SM(name="md_ekin",
                      startReStr=(r"\s*kinetic energy\s*\(Ekin\)\s*=\s*(?P<x_qe_t_md_kinetic_energy__rydberg>" + RE_f +
                                  r")\s*Ry\s*$"),
                   ),
                   SM(name="md_temperature",
                      startReStr=(r"\s*temperature\s*=\s*(?P<x_qe_t_md_temperature__kelvin>" + RE_f +
                                  r")\s*K\s*$"),
                   ),
                   SM(name="md_ekin_etot",
                      startReStr=(r"\s*Ekin\s*\+\s*Etot\s*\(const\)\s*=\s*(?P<x_qe_t_md_ekin_etot__rydberg>" + RE_f +
                                  r")\s*Ry\s*$"),
                   ),
                   SM(name="md_linear_momentum",
                      startReStr=(r"\s*Linear momentum\s*:\s*" + QeC.re_vec('x_qe_t_md_linear_momentum') +
                                  r"\s*$"),
                   ),
                   # newer espresso writes this _after_ 'entering dynamics'
                   SM(name="md_maxSteps2",
                      startReStr=r"\s*The maximum number of steps has been reached.\s*$",
                      fixedStartValues={ 'x_qe_t_md_max_steps_reached': True },
                      subMatchers=self.SMs_molecular_dynamics_end(suffix='2'),
                   ),
               ],
            ),
            SM(name="damped_converged",
               startReStr=r"\s*Damped Dynamics: convergence achieved in\s*(?P<x_qe_t_relax_converged_steps>\d+)\s+steps\s*$",
               subMatchers=[
                   SM(name="damped_end",
                      startReStr=r"\s*End of damped dynamics calculation\s*$",
                      adHoc=lambda p: LOGGER.error("TODO: do sth with end-of-damped-dynamics data"),
                      subMatchers=[
                          SM(name="damped_final_energy",
                              startReStr=r"\s*Final energy =\s*(?P<x_qe_t_relax_final_energy__rydberg>" + RE_f + r")\s*Ry\s*$",
                          ),
                      ] + self.SMs_md_system_new(
                              suffix="DMDFINAL",
                              eatEndFinalSubMatchers=[
                                  SM(name="DMDFi1", coverageIgnore=True,
                                     startReStr=r"\s*Entering Dynamics:\s*iteration.*$",
                                  ),
                                  SM(name="DMDFi2", coverageIgnore=True,
                                     startReStr=r"\s*<vel\(dt\)\|acc\(dt\)>\s*=.*$",
                                  ),
                                  SM(name="DMDFi3", coverageIgnore=True,
                                     startReStr=r"\s*ATOMIC_POSITIONS \(.*\)\s*$",
                                     subMatchers=[
                                         SM(name="DMDFi4", repeats=True, coverageIgnore=True,
                                            startReStr=r"\S+(?:\s+" + RE_f + r"){3}(?:(?:\s+\d+){3})?\s*$",
                                         ),
                                     ],
                                  ),
                              ],
                          ),
                   ),
               ],
            ),
        ]

    def SMs_vcs_molecular_dynamics(self):
        return [
            SM(name='vcs_wentzcovitch_damped_minimization',
                startReStr=r"\s*Wentzcovitch Damped Cell(?:-| )Dynamics Minimization:?\s*$",
               adHoc=lambda p: self.setTmp('md_relax', 'vcsmd_wentzcovitch_damped_minization'),
               subMatchers=[
                   SM(name='vcs_went_dmin_thresholds',
                      startReStr=(
                          r"\s*convergence thresholds:? EPSE =\s*(?P<x_qe_t_relax_threshold_energy__rydberg>" + RE_f + r")" +
                          r"\s*EPSF =\s*(?P<x_qe_t_relax_threshold_force__rydberg_bohr_1>" + RE_f + r")" +
                          r"(?:\s*EPSP =\s*(?P<x_qe_t_relax_threshold_pressure__kilobar>" + RE_f + r"))?" +
                          r"\s*$"
                      ),
                   ),
                   SM(name='vcs_went_dmin_converged',
                      startReStr=(
                          r"\s*convergence achieved, Efinal=\s*(?P<x_qe_t_relax_final_energy__rydberg>" + RE_f + r")\s*$"
                      ),
                      subMatchers=[
                      ] + self.SMs_vcsmd_system_new(suffix='VCS') + [
                      ] + self.SMs_md_system_new(suffix='VCS') + [
                      ],
                   ),
               ],
            ),
            SM(name='vcs_enter',
               startReStr=(r"\s*Entering Dynamics;\s*it\s*=\s*(?P<x_qe_t_md_iteration>" + RE_i + r")\s*" +
                           r"time\s*=\s*(?P<x_qe_t_md_time__picoseconds>" + RE_f +
                           r")\s*pico-seconds\s*$"),
               # there are several algorithms in vcsmd, but apparently only
               # printed when used for minimization and at the end of the calculation ?!
               adHoc = lambda p: self.setTmpUnlessExists('md_relax', 'vcsmd'),
               subMatchers=[
               ] + self.SMs_vcsmd_system_new(suffix='VCSMD1') + [
                   SM(name='vcs_energies',
                      startReStr=(r"\s*Ekin\s*=\s*(?P<x_qe_t_md_kinetic_energy__rydberg>" + RE_f + r")\s*Ry" +
                                  r"\s*T\s*=\s*(?P<x_qe_t_md_temperature__kelvin>" + RE_f + r")\s*K" +
                                  r"\s*Etot\s*=\s*(?P<x_qe_t_md_total_energy>" + RE_f + r")\s*$"),
                   ),
               ] + self.SMs_md_system_new(suffix='VCS1') + [
                   SM(name='vcs_maxiter',
                      startReStr=r"\s*Maximum number of iterations reached, stopping\s*$",
                      fixedStartValues={ 'x_qe_t_md_max_steps_reached': True },
                   ),
               ],
            ),
            SM(name='vcs_went_converged',
               startReStr=(
                   r"\s*Wentzcovitch Damped Dynamics: convergence achieved, Efinal=\s*(?P<x_qe_t_relax_final_energy__rydberg>" + RE_f + r")\s*$"
               ),
               adHoc=lambda p: self.setTmp('md_relax', 'vcsmd_wentzcovitch_damped_minization'),
               subMatchers=[
               ] + self.SMs_vcsmd_system_new(suffix='VCS2') + [
               ] + self.SMs_md_system_new(suffix='VCS2') + [
               ],
            ),
        ]

    def SMs_diagonalization(self, namespace=None):
        return [
            SM(name='david_or_cg', repeats=True,
               startReStr=r"\s*(?P<x_qe_t_" + namespace + r"_algorithm>Davidson diagonalization with overlap|CG style diagonalization)\s*$",
               adHoc=lambda p: p.backend.addValue('x_qe_' + namespace + '_algorithm', QE_DIAGONALIZATION[p.lastMatch['x_qe_t_' + namespace + '_algorithm']]),
               sections=['x_qe_section_' + namespace],
               subMatchers=[
                   SM(name='warn_not_converged', repeats=True, floating=True,
                      startReStr=(r"\s*WARNING:\s*(?P<x_qe_" + namespace + r"_warn_n_unconverged_eigenvalues>" + RE_i +
                                  r")\s*eigenvalues not converged(?:\s+in\s+regterg)?\s*$"),
                   ),
                   SM(name='c_bands_not_converged', repeats=True,
                      startReStr=(r"\s*c_bands:\s*(?P<x_qe_" + namespace + r"_c_bands_n_unconverged_eigenvalues>" + RE_i +
                                  r")\s*eigenvalues not converged\s*$"),
                   ),
                   SM(name='ethr', repeats=True,
                      startReStr=(r"\s*ethr\s*=\s*(?P<x_qe_" + namespace + "_ethr>" + RE_f +
                                  r")\s*,\s*avg\s*#\s*of iterations\s*=\s*(?P<x_qe_" + namespace + "_iteration_avg>" + RE_f +
                                  r")\s*$"),
                      subMatchers=[
                          SM(name='redo_with_lower_ethr_cIgn1', coverageIgnore=True,
                             startReStr=r"\s*Threshold \(ethr\) on eigenvalues was too large:\s*$",
                             subMatchers=[
                                 SM(name='redo_with_lower_ethr_cIgn2', coverageIgnore=True,
                                    startReStr=r"\s*Diagonalizing with lowered threshold\s*$",
                                 ),
                             ],
                          ),
                          SM(name='warn_ef_above_band', repeats=True,
                             startReStr=r"\s*Warning: (?P<x_qe_t_warning>ef =\s*" + RE_f + "\s*is above the highest band at k-point\s*\d+)\s*$",
                             adHoc=lambda p: self.setTmp('x_qe_t_warning', p.lastMatch['x_qe_t_warning']),
                             subMatchers=[
                                 SM(name="warn_ef_above_band1",
                                    startReStr=r"\s*(?P<x_qe_t_warning>e\s*=\s*" + RE_f + r")\s*$",
                                    adHoc=lambda p: p.backend.addValue('x_qe_warning', (self.popTmp('x_qe_t_warning') + "\n" + p.lastMatch['x_qe_t_warning'])),
                                 )
                             ],
                          ),
                      ],
                   ),
               ],
            ),
        ]

    def header_sections(self):
        return ['section_method',
                'section_system', 'x_qe_section_parallel',
                'x_qe_section_compile_options']

    def close_header_sections(self, backend):
        # close header sections if they are open
        for sec in self.header_sections():
            sec_gIndex = self.openSectionIdx.pop(sec,None)
            if sec_gIndex is not None:
                backend.closeSection(sec, sec_gIndex)

    def run_submatchers(self):
        """submatchers of section_run"""
        return [
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
            # pure informational msg about how code was compiled
            SM(name='cIgn_pp_compiletype', coverageIgnore=True,
               startReStr=(
                   r"\s*(?:For Norm-Conserving or\s*)?Ultrasoft \(Vanderbilt\) Pseudopotentials\s*(?:(?:and|or) PAW)?\s*"
               ),
            ),
        ] + self.SMs_read_input_file() + [
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
                   SM(name="qe_warn_efield_symm",
                      startReStr=r"\s*(?P<x_qe_warning>Presently no symmetry can be used with electric field)\s*",
                   ),
               ],
            ),
        ] + self.SMs_read_input_file(suffix='2') + [
            SM(name='input_positions_cell_dir_header1',
               startReStr=r"\s*Atomic positions and unit cell read from directory:\s*$",
               subMatchers=[
                   SM(name='input_positions_cell_dir1',
                      startReStr=r"\s*(?P<x_qe_input_positions_cell_dirname>.+?)\s*$",
                   ),
               ],
            ),
            SM(name='supercell1',
               startReStr=r"\s*Found symmetry operation:\s*I\s*\+\s*\(\s*" + QeC.re_vec('x_qe_t_vec_supercell') + r"\s*\)\s*$",
            ),
            SM(name='supercell2',
               startReStr=r"\s*This is a supercell, fractional translations? are disabled\s*$",
               adHoc=lambda p: p.backend.addValue('x_qe_supercell', True)
            ),
            SM(name='supercell3', repeats=True,
               # observed in PWSCF 4.0
               # "     Found additional translation:   -0.5000   -0.5000    0.0000"
               startReStr=r"\s*Found additional translation:\s*" + QeC.re_vec('x_qe_t_vec_supercell') + r"\s*$",
               adHoc=lambda p: p.backend.addValue('x_qe_supercell', True),
            ),
            SM(name='supercell41', repeats=True,
               # observed in PWSCF 4.1
               # "     Fractionary translation:   -0.5000   -0.5000    0.0000is a symmetry operation:"
               startReStr=r"\s*Fractionary\s+translation:\s+" + QeC.re_vec('x_qe_t_vec_supercell') + r"\s*is a symmetry operation:\s*$",
               subMatchers=[
                   SM(name="supercell41_ft_disabled",
                      startReStr=r"\s*This is a supercell, fractionary translation are disabled:\s*$",
                   ),
               ],
               adHoc=lambda p: p.backend.addValue('x_qe_supercell', True),
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
            SM(name='enforced_XC',
               startReStr=r"\s*IMPORTANT: XC functional enforced from input\s*:\s*",
               subMatchers=[
                   SM(name='xc_functional_enforced', required=True, # details are parsed in xc_functional
                      startReStr=r"\s*Exchange-correlation\s*=\s*(?P<x_qe_t_xc_functional_shortname_enforced>\S+)\s*\([^\(]+\)\s*$"
                   ),
                   SM(name='exx_fraction_enforced',
                      startReStr=r"\s*EXX-fraction\s*=\s*" + RE_f + r"\s*$",
                   ),
                   SM(name='xc_functional_enforced_txt1', coverageIgnore=True,
                      startReStr=r"\s*Any further DFT definition will be discarded",
                   ),
                   SM(name='xc_functional_enforced_txt1', coverageIgnore=True,
                      startReStr=r"\s*Please, verify this is what you really want",
                   ),
               ],
            ),
            SM(name='renormalized_pseudo_wavefunction', repeats=True,
               startReStr=r"\s*file\s*(?P<x_qe_t_pp_renormalized_filename>.*?)\s*:\s*wavefunction\(s\)\s*(?P<x_qe_t_pp_renormalized_wfc>.*?)\s*renormalized\s*$",
               adHoc=self.adHoc_pp_renorm,
            ),
            SM(name='pseudopotential_warning', repeats=True,
               # same as 'renormalized_pseudo_wavefunction', but seen in QE 4.1
               startReStr=(r"\s*WARNING: Pseudopotential #\s*(?P<x_qe_t_pp_warning_idx>\d+)\s*"+
                           r"file\s*:\s*(?P<x_qe_t_pp_warning_filename>.+?)\s*$"),
               subMatchers=[
                   SM(name='pp_warning_not_normalized', repeats=True,
                      startReStr=(r"\s*WARNING: WFC #\s*(?P<x_qe_t_pp_warning_wfcidx>\d+)\s*" +
                                  r"\((?P<x_qe_t_pp_warning_wfclabel>[^\)]+)\) IS NOT CORRECTLY NORMALIZED:\s*" +
                                  r"norm=\s*(?P<x_qe_t_pp_warning_wfcnorm>" + RE_f + r")\s*$"),
                      subMatchers=[
                          SM(name='pp_warning_normalization_msg', repeats=True,
                             startReStr=r"\s*WARNING: WFC HAS BEEN NOW RENORMALIZED\s*!?\s*$",
                          ),
                      ],
                   ),
               ],
               sections=['x_qe_t_section_pp_warning'],
            ),
            SM(name='input_positions_cell_dir_header2',
               startReStr=r"\s*Atomic positions and unit cell read from directory:\s*$",
               subMatchers=[
                   SM(name='input_positions_cell_dir2',
                      startReStr=r"\s*(?P<x_qe_input_positions_cell_dirname>.+?)\s*$",
                   ),
               ],
            ),
            SM(name='dispersion_correction_obsolete_iosys',
               startReStr=r"\s*Message from routine iosys:\s*$",
               subMatchers=[
                   SM(name='dispersion_correction_obsolete',
                      startReStr=r"\s*(?P<x_qe_warning>london is obsolete, use \"vdw_corr='grimme-d2'\" instead)\s*$",
                   )
               ],
            ),
            SM(name='dispersion_correction',
               startReStr=r"\s*Parameters for Dispersion Correction:\s*$",
               subMatchers=[
                   SM(name='dispersion_correction_header',
                      startReStr=r"\s*atom\s*VdW radius\s*C_6\s*$",
                      subMatchers=[
                          SM(name='dispersion_correction_values', repeats=True,
                             startReStr=(r"\s*(?P<x_qe_t_species_dispersion_correction_label>.+?)\s+" +
                                         r"(?P<x_qe_t_species_dispersion_correction_vdw_radius>" + RE_f + r")" +
                                         r"\s*(?P<x_qe_t_species_dispersion_correction_C6>" + RE_f + r")\s*$"),
                             adHoc=self.adHoc_dispersion_correction_values,
                          ),
                      ],
                   ),
               ],
               adHoc=lambda p: p.backend.addValue('x_qe_dispersion_correction', True)
            ),
            SM(name='gamma_algorithms',
               startReStr=r"\s*gamma-point specific algorithms are used\s*$",
               adHoc=lambda p: p.backend.addValue('x_qe_gamma_algorithms', True)
            ),
            SM(name='warning_setup',
               startReStr=r"\s*Message from routine setup:\s*$",
               subMatchers=[
                  SM(name='warning_metallic',
                     startReStr=r"\s*(?P<x_qe_warning>the system is metallic, specify occupations)\s*$",
                  ),
                  SM(name='warning_dynamics_symm',
                     startReStr=r"\s*(?P<x_qe_warning>Dynamics, you should have no symmetries)\s*$",
                  ),
               ],
            ),
            SM(name='symm_incompatible_FFT', repeats=True,
               startReStr=r"\s*(?P<x_qe_t_warning>warning: symmetry operation #\s*\d+\s*not compatible with FFT grid\.)\s*$",
               adHoc=lambda p: self.setTmp('x_qe_t_warning', p.lastMatch['x_qe_t_warning']),
               subMatchers=[
                   SM(name='symm_incompatible_FFT_row1',
                      startReStr=r"\s*(?P<x_qe_t_warning>" + RE_f + r"\s*" + RE_f + r"\s*" + RE_f + r")\s*$",
                      adHoc=lambda p: self.appendToTmp('x_qe_t_warning', "\n" + p.lastMatch['x_qe_t_warning']),
                   ),
                   SM(name='symm_incompatible_FFT_row2',
                      startReStr=r"\s*(?P<x_qe_t_warning>" + RE_f + r"\s*" + RE_f + r"\s*" + RE_f + r")\s*$",
                      adHoc=lambda p: self.appendToTmp('x_qe_t_warning', "\n" + p.lastMatch['x_qe_t_warning']),
                   ),
                   SM(name='symm_incompatible_FFT_row3',
                      startReStr=r"\s*(?P<x_qe_t_warning>" + RE_f + r"\s*" + RE_f + r"\s*" + RE_f + r")\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_warning', (self.popTmp('x_qe_t_warning') + "\n" + p.lastMatch['x_qe_t_warning'])),
                   ),
               ],
            ),
            SM(name='exx_grid_same',
               startReStr=r"\s*EXX: grid of k\+q points same as grid of k-points",
               fixedStartValues={ "x_qe_exx_grid_same_as_k_grid": True },
            ),
            SM(name='exx_k_plus_q_grid',
               startReStr=r"\s*EXX: setup a grid of \d+ q-points centered on each k-point\s*$",
               adHoc=lambda p: LOGGER.info("do something with EXX k+q grid"),
               subMatchers=[
                   SM(name='exx_k_plus_q_grid_header', coverageIgnore=True,
                      startReStr=r"\s*\(k\+q\)-points:\s*$",
                   ),
                   SM(name='exx_k_plus_q_grid_line', repeats=True,
                      startReStr=r"\s*" + RE_f + r"\s+" + RE_f + r"\s+" + RE_f + r"\s+" + RE_i + r"\s+" + RE_i + "\s*$",
                   ),
               ],
            ),
            SM(name='forbidden_symm1', repeats=True,
               startReStr=r"\s*warning: (?P<x_qe_t_warning>symmetry operation #\s*\d+\s*not allowed.\s*fractional translation:)\s*$",
               adHoc=lambda p: self.setTmp('x_qe_t_warning', p.lastMatch['x_qe_t_warning']),
               subMatchers=[
                   SM(name="forbidden_symm_ft1",
                      startReStr=r"\s*(?P<x_qe_t_warning>" + RE_f + r"\s*" + RE_f +  r"\s*" + RE_f + r"\s*in crystal coordinates)\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_warning', (self.popTmp('x_qe_t_warning') + "\n" + p.lastMatch['x_qe_t_warning'])),
                   )
               ],
            ),
            SM(name='subspace_diagonalization',
               startReStr=r"\s*(?:Subspace diagonalization in iterative solution of the eigenvalue problem:|Iterative solution of the eigenvalue problem)\s*$",
               subMatchers=[
                   SM(name='too_few_procs',
                      startReStr=r"\s*Too few procs for parallel algorithm\s*$",
                      subMatchers=[
                          SM(name='min_procs',
                             startReStr=r"\s*we need at least \d+ procs per pool\s*$",
                          )
                      ],
                   ),
                   SM(name='too_few_proc_42',
                      # fallback for 4.2, msg on one line
                      startReStr=r"\s*Too few procs for parallel algorithm: we need at least \d* procs per pool\s*$",
                   ),
                   SM(name='serial_algorithm',
                      startReStr=r"\s*a serial algorithm will be used\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_diagonalization_algorithm', 'serial')
                   ),
               ],
            ),
            SM(name='final_scf_MD', floating=True,
               # entry point for final scf calculation after VC-relax
               # used to reset position (a.k.a. as a goto) in our SM hierarchy
               startReStr=r"\s*A final scf calculation at the relaxed structure\.\s*$",
               adHoc=self.adHoc_final_scf_MD,
               subMatchers=[
                   SM(name='final_scf_MD2',
                      startReStr=r"\s*The G-vectors are recalculated(?:\.|\s*for the final unit cell)\s*$",
                   ),
                   SM(name='final_scf_MD3',
                      startReStr=r"\s*Results may differ from those at the preceding step\.\s*$",
                   ),
               ],
            ),
            SM(name='forbidden_symm2', repeats=True,
               startReStr=r"\s*warning: (?P<x_qe_t_warning>symmetry operation #\s*\d+\s*not allowed.\s*fractional translation:)\s*$",
               adHoc=lambda p: self.setTmp('x_qe_t_warning', p.lastMatch['x_qe_t_warning']),
               subMatchers=[
                   SM(name="forbidden_symm_ft2",
                      startReStr=r"\s*(?P<x_qe_t_warning>" + RE_f + r"\s*" + RE_f +  r"\s*" + RE_f + r"\s*in crystal coordinates)\s*$",
                      adHoc=lambda p: p.backend.addValue('x_qe_warning', (self.popTmp('x_qe_t_warning') + "\n" + p.lastMatch['x_qe_t_warning'])),
                   )
               ],
            ),
            SM(name='mesh_sticks',
               startReStr=r"\s*G-vector sticks info\s*$",
               subMatchers=[
                   SM(name="sticks_heading", required=True,
                      startReStr=r"\s*sticks:\s*dense\s+smooth\s+PW\s+G-vecs:\s+dense\s+smooth\s+PW\s*$",
                      subMatchers=[
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
               ],
            ),
            SM(name='mesh_sticks43',
               startReStr=r"\s*Stick Mesh\s*$",
               adHoc=lambda p: LOGGER.info("parse 4.3 sticks mesh properly"),
               subMatchers=[
                   SM(name="sticks_summary43",
                      startReStr=r"\s*(?P<x_qe_sticks_old>nst =.*)\s*$",
                   ),
                   SM(name="sticks_header43",
                      startReStr=r"\s*(?P<x_qe_sticks_old>n\.st\s+n\.stw\s+n\.sts\s+n\.g\s+n\.gw\s+n\.gs)\s*$",
                      subMatchers=[
                          SM(name='sticks_line43', repeats=True,
                             startReStr=r"\s*(?P<x_qe_sticks_old>(?:min|max|)(?:\s+\d+){6})\s*$",
                          ),
                      ],
                   ),
               ],
            ),
            SM(name='mesh_sticks40',
               startReStr=r"\s*(?P<x_qe_sticks_old>Planes per process \(thick\)\s*:.*?)\s*$",
               adHoc=lambda p: LOGGER.info("parse <= 4.0 sticks mesh properly"),
               subMatchers=[
                   SM(name="sticks_smooth40",
                      startReStr=r"\s*(?P<x_qe_sticks_old>Planes per process \(smooth\)\s*:.*?)\s*$",
                   ),
                   SM(name="sticks_header40", required=True,
                      startReStr=r"\s*Proc/\s+planes\s+cols\s+G\s+planes\s+cols\s+G\s+columns\s+G\s*$",
                      subMatchers=[
                          SM(name="sticks_header40_2", required=True,
                             startReStr=r"\s*Pool\s+\(dense grid\)\s+\(smooth grid\)\s+\(wavefct grid\)\s*$",
                             subMatchers=[
                                 SM(name='sticks_line40', repeats=True,
                                    startReStr=r"\s*(?P<x_qe_sticks_old>(?:\s+\d+){9})\s*$",
                                 ),
                             ],
                          ),
                      ],
                   ),
               ],
            ),
            SM(name='atom_radii',
               startReStr=r"\s*Generating pointlists\s*\.\.\.\s*$",
               subMatchers=[
                   SM(name='new_r_m', repeats=True,
                      # radius is in alat units, but they are not yet defined.
                      # convert manually in onClose hook...
                      startReStr=(r"\s*new\s+r_m\s*:\s*(?P<x_qe_t_species_integration_radius>" + RE_f + r")\s*" +
                                  r"(?:\((?:alat|a_0)\s*units\)\s*" + RE_f + r"\s*\(a\.u\.\)\s*for type\s*" +
                                  r"(?P<x_qe_t_species_integration_radius_idx>" + RE_i + r"))?\s*$"
                      ),
                   ),
               ],
            ),
        ] + self.SMs_summaryf90() + [
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
            SM(name='martyna_tuckerman_parameters',
               startReStr=(r"\s*alpha, beta MT =\s*(?P<x_qe_isolated_system_method_martyna_tuckerman_alpha>" + RE_f + r")\s*" +
                           r"(?P<x_qe_isolated_system_method_martyna_tuckerman_beta>" + RE_f + r")\s*$"),
               fixedStartValues={ 'x_qe_isolated_system_method': 'Martyna-Tuckerman' },
            ),
            SM(name='core_charge_check',
               startReStr=(r"\s*Check: negative/imaginary core charge\s*=\s*(?P<x_qe_core_charge_negative>" +
                           RE_f + r")\s*(?P<x_qe_core_charge_imaginary>" + RE_f + r")\s*$")
            ),
            SM(name='init_aug_dense_cIgn', coverageIgnore=True,
               startReStr=r"\s*Initializing real-space augmentation for DENSE grid\s*$",
               subMatchers=[
                   SM(name='init_aug_smooth_dense_cIgn', coverageIgnore=True,
                      startReStr=r"\s*SMOOTH grid -> DENSE grid\s*$",
                   ),
               ],
            ),
            SM(name='input_potential_recalculated_file_header',
               startReStr=r"\s*The potential is recalculated from file :\s*$",
               subMatchers=[
                   SM(name='input_potential_recalculated_file',
                      startReStr=r"\s*(?P<x_qe_input_potential_recalculated_file>.+?)\s*$",
                   ),
               ],
            ),
            SM(name='initial_density_from_file',
               startReStr=r"\s*The initial density is read from file\s*:\s*$",
               subMatchers=[
                   SM(name='initial_density_file',
                      startReStr=r"\s*(?P<x_qe_starting_density_file>.+\.(?:dat|xml))\s*$",
                   ),
               ],
            ),
            SM(name='initial_potential',
               startReStr=r"\s*Initial potential from\s*(?P<x_qe_starting_potential>.*?)\s*$",
            ),
            SM(name='starting_charge_negative',
               startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_starting_charge_negative>" + RE_f +
                           r")\s*$"),
            ),
            SM(name='starting_charge_negative_spin_ignore', repeats=True,
               startReStr=(r"\s*Check: negative starting charge=\(component\d\):?\s*" + RE_f + r"\s*$"),
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
            SM(name='vdw_df-citation', coverageIgnore=True,
               # xc_functionals does reflect this, so simply ignore
               startReStr=r"\s*% You are using vdW-DF, which was implemented by the Thonhauser group. %\s*",
               subMatchers=[
                   SM(name='vdw_df-citation1', coverageIgnore=True,
                      startReStr=r"\s*% Please cite the following two papers that made this development      %\s*",
                   ),
                   SM(name='vdw_df-citation2', coverageIgnore=True,
                      startReStr=r"\s*% possible and the two reviews that describe the various versions:     %\s*",
                   ),
                   SM(name='vdw_df-citation3', coverageIgnore=True,
                      startReStr=r"\s*%   T. Thonhauser et al., PRL 115, 136402 \(2015\).                      %\s*",
                   ),
                   SM(name='vdw_df-citation4', coverageIgnore=True,
                      startReStr=r"\s*%   T. Thonhauser et al., PRB 76, 125112 \(2007\).                       %\s*",
                   ),
                   SM(name='vdw_df-citation5', coverageIgnore=True,
                      startReStr=r"\s*%   K. Berland et al., Rep. Prog. Phys. 78, 066501 \(2015\).             %\s*",
                   ),
                   SM(name='vdw_df-citation6', coverageIgnore=True,
                      startReStr=r"\s*%   D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 \(2009\). %\s*",
                   ),
                   SM(name='vdw_df-citation7', coverageIgnore=True,
                      startReStr=r"\s*% If you are calculating the stress with vdW-DF, please also cite:     %\s*",
                   ),
                   SM(name='vdw_df-citation8', coverageIgnore=True,
                      startReStr=r"\s*%   R. Sabatini et al., J. Phys.: Condens. Matter 24, 424209 \(2012\).   %\s*",
                   ),
               ],
            ),
            SM(name='starting_wfc',
               startReStr=r"\s*Starting wfc\s*(?P<x_qe_starting_wfc>.*?)\s*$",
            ),
            SM(name="write_datafileINITcIgn", coverageIgnore=True,
               startReStr=r"\s*Writing output data file\s*.*?\s*$", # (?P<x_qe_t_output_datafile>.*?)\s*$",
               subMatchers=[
                  SM(name="warning_save_mgga2",
                     startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                  ),
                  SM(name='starting_charge_negative2',
                     startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_starting_charge_negative>" + RE_f +
                                 r")\s*$"),
                  ),
                  SM(name='starting_charge_negative_spin_ignore2', repeats=True,
                     startReStr=(r"\s*Check: negative starting charge=\(component\d\):?\s*" + RE_f + r"\s*$"),
                  ),
               ],
            ),
            SM(name="md_new_old_atomic_charge_approximation",
               startReStr=r"\s*NEW-OLD atomic charge density approx\. for the potential\s*$",
               fixedStartValues={ 'x_qe_extrapolation_charge': 'atomic' },
               subMatchers=[
                   SM(name='core_charge_check',
                      startReStr=(r"\s*Check: negative/imaginary core charge\s*=\s*(?P<x_qe_core_charge_negative>" +
                                  RE_f + r")\s*(?P<x_qe_core_charge_imaginary>" + RE_f + r")\s*$")
                   ),
                   SM(name='starting_charge_negative2',
                      startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_starting_charge_negative>" + RE_f +
                                  r")\s*$"),
                   ),
                   SM(name='starting_charge_negative_spin_ignore2', repeats=True,
                      startReStr=(r"\s*Check: negative starting charge=\(component\d\):?\s*" + RE_f + r"\s*$"),
                   ),
                   SM(name='starting_rho2',
                      startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_starting_charge_negative_up>" + RE_f +
                                  r")\s*(?P<x_qe_starting_charge_negative_down>" + RE_f + r")\s*$"),
                   ),
               ],
            ),
            SM(name='paw_dealloc_cIgn', coverageIgnore=True,
               startReStr=r"\s*Checking if some PAW data can be deallocated\.\.\.\s*$",
            ),
            SM(name='cputime_msg',
               startReStr=(r"\s*total cpu time spent up to now is\s*(?P<x_qe_time_setup_cpu1_end>" + RE_f +
                           r")\s*secs\s*$"),
            ),
            SM(name='per_process_mem',
               startReStr=r"\s*per-process dynamical memory:\s*(?P<x_qe_per_process_mem__mebibyte>" + RE_f + ")\s*Mb\s*$",
            ),
            # ------------------------------------------
            # end of header
            # ------------------------------------------
            SM(name='self_consistent_calculation', repeats=True,
               startReStr=r"\s*Self-consistent Calculation\s*$",
               sections = ['section_single_configuration_calculation'],
               subMatchers=[
                   SM(name='iteration', repeats=True,
                      startReStr=(r"\s*iteration\s*#\s*(?P<x_qe_iteration_number>\d+)\s*" +
                                  r"\s*ecut\s*=\s*(?P<x_qe_iteration_ecutwfc__rydberg>" + RE_f +r")\s*Ry" +
                                  r"\s*beta\s*=\s*(?P<x_qe_iteration_beta>" + RE_f + r")\s*$"),
                      sections=['section_scf_iteration'],
                      subMatchers=[
                      ] + self.SMs_diagonalization(namespace='scf_diagonalization') + [
                          SM(name="iter_warning_save_mggaI",
                             startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                          ),
                          SM(name='iteration_rho',
                             startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_iteration_charge_negative_up>" + RE_f +
                                         r")\s*(?P<x_qe_iteration_charge_negative_down>" + RE_f + r")\s*$"),
                          ),
                          SM(name='iteration_per_site_magnetization_header',
                             startReStr=(r"\s*Magnetic\s*moment\s*per\s*site:?\s*$"),
                             subMatchers=[
                                 SM(name='iteration_per_site_magnetization', repeats=True,
                                    startReStr=(r"\s*atom:\s*(?P<x_qe_t_iter_mpersite_idx>" + RE_i + r")\s*" +
                                                r"charge:\s*(?P<x_qe_t_iter_mpersite_charge>" + RE_f + r")\s*" +
                                                r"magn:\s*(?P<x_qe_t_iter_mpersite_magn>" + RE_f + r")\s*" +
                                                r"constr:\s*(?P<x_qe_t_iter_mpersite_constr>" + RE_f + r")\s*$"),
                                 ),
                             ],
                          ),
                          SM(name='iteration_Efield',
                             startReStr=(r"\s*Expectation value of exp\(iGx\):\s*\(\s*" +
                                         r"(?P<x_qe_iteration_efield_eeigx_re>" + RE_f + r")\s*,\s*" +
                                         r"(?P<x_qe_iteration_efield_eeigx_im>" + RE_f + r")\s*\)(?:\s*" + RE_f + r")?\s*$"),
                             subMatchers=[
                                 SM(name='iteration_Efield_dkfact',
                                    startReStr=(r"\s*" + RE_f + r"\s*$"),
                                 ),
                                 SM(name='iteration_Efield_dipole_electric',
                                    startReStr=(r"\s*Electronic Dipole per cell \((?:Ry\s*)?a.u.\)\s*" +
                                                r"(?P<x_qe_iteration_efield_dipole_electronic__rydberg>" + RE_f + r")\s*$"),
                                 ),
                                 SM(name='iteration_Efield_dipole_ionic',
                                    startReStr=(r"\s*Ionic Dipole per cell \((?:Ry\s*)?a.u.\)\s*" +
                                                r"(?P<x_qe_iteration_efield_dipole_ionic__rydberg>" + RE_f + r")\s*$"),
                                 ),
                             ],
                          ),
                          SM(name='cputime_iteration_msg',
                             startReStr=(r"\s*total cpu time spent up to now is\s*(?P<time_scf_iteration_cpu1_end>" + RE_f +
                                         r")\s*secs\s*$"),
                          ),
                          SM(name='e_total',
                             startReStr=(r"\s*!?\s*total\s+energy\s*=\s*(?P<energy_total_scf_iteration__rydberg>" + RE_f + r")" +
                                         r"\s*Ry\s*$"),
                          ),
                          SM(name='harris',
                             startReStr=(r"\s*Harris-Foulkes estimate\s*=\s*(?P<x_qe_energy_total_harris_foulkes_estimate_iteration__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='estimate_accuracy',
                             startReStr=(r"\s*estimated scf accuracy\s*<\s*(?P<x_qe_energy_total_accuracy_estimate_iteration__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='total_magnetization',
                             startReStr=(r"\s*total magnetization\s*=\s*(?P<x_qe_magnetization_total_iteration__mubohr>" + RE_f +
                                         r")\s*Bohr mag/cell\s*$"),
                          ),
                          SM(name='absolute_magnetization',
                             startReStr=(r"\s*absolute magnetization\s*=\s*(?P<x_qe_magnetization_absolute_iteration__mubohr>" + RE_f +
                                         r")\s*Bohr mag/cell\s*$"),
                          ),
                      ],
                   ),
                   SM(name='scf_result', repeats=False,
                       startReStr=r'\s*End of self-consistent calculation\s*$',
                       subMatchers=self.SMs_bands() + [
                          SM(name='bands_spin', repeats=True,
                              startReStr=r"\s*-+\s*SPIN\s+(?P<x_qe_t_spin_channel>UP|DOWN)\s*-+\s*$",
                              adHoc=self.adHoc_bands_spin,
                              subMatchers=self.SMs_bands(suffix='_spin'),
                          ),
                          SM(name='highest_occupied',
                             startReStr=(r"\s*highest occupied level \(ev\):\s*(?P<x_qe_t_energy_reference_highest_occupied__eV>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name='highest_occupied_lowest_unoccupied',
                             startReStr=(r"\s*highest occupied, lowest unoccupied level \(ev\):\s*" +
                                         r"(?P<x_qe_t_energy_reference_highest_occupied__eV>" + RE_f + r")\s*" +
                                         r"(?P<x_qe_t_energy_reference_lowest_unoccupied__eV>" + RE_f + r")\s*$"),
                          ),
                          SM(name='e_fermi',
                             startReStr=(r"\s*the Fermi energy is\s*(?P<x_qe_t_energy_reference_fermi__eV>" + RE_f + ")\s*ev\s*$"),
                          ),
                          SM(name='e_fermi_spin',
                             startReStr=(r"\s*the spin up/dw Fermi energies are\s*" +
                                         r"(?P<x_qe_t_energy_reference_fermi_up__eV>" + RE_f + ")\s*" +
                                         r"(?P<x_qe_t_energy_reference_fermi_down__eV>" + RE_f + ")\s*ev\s*$"),
                          ),
                          SM(name='e_total',
                             startReStr=r'\s*!?\s*total\s+energy\s*=\s*(?P<energy_total__rydberg>' + RE_f + ')\s*Ry\s*$',
                          ),
                          SM(name='harris',
                             startReStr=(r"\s*Harris-Foulkes estimate\s*=\s*(?P<x_qe_energy_total_harris_foulkes_estimate__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='estimate_accuracy',
                             startReStr=(r"\s*estimated scf accuracy\s*<\s*(?P<x_qe_energy_total_accuracy_estimate__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='total_AE_energy',
                             startReStr=(r"\s*total all-electron energy\s*=\s*(?P<x_qe_energy_total_paw_all_electron__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='energy_decomposition',
                             startReStr=r"\s*The total energy is the sum of the following terms\s*:\s*$",
                             subMatchers=[
                                 SM(name='energy_decomposition_contribution', repeats=True,
                                    startReStr=(r"\s*(?P<x_qe_t_energy_decomposition_name>.*?)\s*=\s*" +
                                                r"(?P<x_qe_t_energy_decomposition_value__rydberg>" + RE_f + r")\s*Ry\s*$"),
                                 ),
                             ],
                          ),
                          SM(name='total_magnetization',
                             startReStr=(r"\s*total magnetization\s*=\s*(?P<x_qe_magnetization_total__mubohr>" + RE_f +
                                         r")\s*Bohr mag/cell\s*$"),
                          ),
                          SM(name='absolute_magnetization',
                             startReStr=(r"\s*absolute magnetization\s*=\s*(?P<x_qe_magnetization_absolute__mubohr>" + RE_f +
                                         r")\s*Bohr mag/cell\s*$"),
                          ),
                          SM(name="convergence_iterations",
                             startReStr=r"\s*convergence has been achieved in\s*(?P<x_qe_convergence_iterations>\d+)\s*iterations\s*",
                          ),
                          SM(name="convergence_achieved",
                             startReStr=r"\s*convergence has been achieved\s*$",
                             adHoc=lambda p: p.backend.addValue('x_qe_convergence_iterations', p.superContext.tmp['last_iteration'])
                          ),
                          SM(name="warning_save_mgga1",
                             startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                          ),
                          # EXX writes them here
                          SM(name='e_totalEXX',
                             startReStr=r'\s*!?\s*total\s+energy\s*=\s*(?P<energy_total__rydberg>' + RE_f + ')\s*Ry\s*$',
                          ),
                          SM(name='harrisEXX',
                             startReStr=(r"\s*Harris-Foulkes estimate\s*=\s*(?P<x_qe_energy_total_harris_foulkes_estimate__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='exchange_errorEXX',
                             startReStr=(r"\s*est. exchange err \(dexx\)\s*=\s*(?P<x_qe_energy_exchange_error_estimate__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='fock_averageEXX',
                             startReStr=(r"\s*- averaged Fock potential\s*=\s*(?P<x_qe_energy_exchange_average_fock_potential__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name='fock_energyEXX',
                             startReStr=(r"\s*\+ Fock energy\s*=\s*(?P<x_qe_energy_fock__rydberg>" +
                                         RE_f + r")\s*Ry\s*$"),
                          ),
                          SM(name="exx_self_consistency",
                             startReStr=r"\s*EXX self-consistency reached\s*$",
                             fixedStartValues={ "x_qe_exx_self_consistency": True },
                          ),
                          SM(name='md_starting_rho_new2',
                             startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_t_md_starting_charge_negative_new_up>" + RE_f +
                                         r")\s*(?P<x_qe_t_md_starting_charge_negative_new_down>" + RE_f + r")\s*$"),
                          ),
                          SM(name="atom_forces",
                             startReStr=r"\s*Forces acting on atoms\s*\(Ry/au\):\s*$",
                             subMatchers=[
                                 SM(name="negative_rho_atforces",
                                    # we had this already earlier
                                    startReStr=r"\s*negative rho \(up, down\):\s*" + RE_f + r"\s*" + RE_f + r"\s*$",
                                 ),
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
                          SM(name="atoms_dispersion_forces",
                             startReStr=(r"\s*Total Dispersion Force =\s*(?P<x_qe_dispersion_force_total__rydberg_bohr_1>" +
                                         RE_f + r")\s*"),
                             subMatchers=[
                                 SM(name="atom_dispersion_force_header",
                                    startReStr=r"\s*Dispersion forces acting on atoms:\s*$",
                                    subMatchers=[
                                        SM(name="atom_dispersion_force", repeats=True,
                                           startReStr=(r"\s*atom\s*(?P<x_qe_t_dispersion_force_atom_idx>\d+)"+
                                                       r"\s*type\s*\d+\s*force\s*=\s*"+
                                                       QeC.re_vec('x_qe_t_dispersion_force', 'rydberg_bohr_1') + "\s*$")
                                        ),
                                    ],
                                 ),
                             ],
                          ),
                          SM(name="force_ethr_warning",
                             startReStr=r"\s*(?P<x_qe_warning>SCF correction compared to forces is(?:\s+too)?\s+large[:,] reduce conv_thr(?:\s+to get better values)?)\s*$",
                          ),
                          SM(name="stress_tensor",
                             startReStr=r"\s*entering subroutine stress \.\.\.\s*$",
                             subMatchers=[
                                 SM(name="negative_rho_stress",
                                    # we had this already earlier
                                    startReStr=r"\s*negative rho \(up, down\):\s*" + RE_f + r"\s*" + RE_f + r"\s*$",
                                 ),
                                 SM(name="stress_header",
                                    startReStr=(r"\s*total\s*stress\s*\(Ry/bohr\*\*3\)\s*\(kbar\)\s*P=\s*" +
                                                r"(?P<x_qe_pressure__kilobar>" + RE_f + r")\s*$"),
                                    subMatchers=[
                                        SM(name="stress_components", repeats=True,
                                           startReStr=(r"\s*" + QeC.re_vec('x_qe_t_stress', 'rydberg_bohr_3') + "\s*" +
                                                       RE_f + r"\s*" + RE_f + r"\s*" + RE_f + r"\s*$")
                                        ),
                                    ],
                                 ),
                             ],
                          ),
                          SM(name="stress_message",
                             startReStr=r"\s*Message from routine stress\s*:?\s*$",
                             subMatchers=[
                                 SM(name="no_mgga",
                                    startReStr=r"\s*Meta-GGA and stress not implemented\s*$",
                                    adHoc=lambda p: p.backend.addValue("x_qe_stress_unimplemented", "Meta-GGA")
                                 ),
                             ],
                          ),
                      ] + self.SMs_relax_bfgs() + [
                      ] + self.SMs_molecular_dynamics() + [
                      ] + self.SMs_vcs_molecular_dynamics() + [
                          SM(name="write_datafile",
                             startReStr=r"\s*Writing output data file\s*(?P<x_qe_output_datafile>.*?)\s*$",
                             subMatchers=[
                                SM(name="warning_save_mgga2",
                                   startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                                ),
                             ],
                          ),
                          # The following belongs actually to the 'next' scf calculation!
                          SM(name='md_starting_charge_negative_old',
                             startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_t_md_starting_charge_negative_old>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name="md_wave_function_extrapolation",
                             startReStr=r"\s*(?P<x_qe_t_md_extrapolation_wfc>\S+ order)\s*wave-functions extrapolation\s*$",
                             subMatchers=[
                                 SM(name='warning_header',
                                    startReStr=r"\s*Message from extrapolate_wfcs:?\s*$",
                                    subMatchers=[
                                        SM(name='warning_msg',
                                           startReStr=r"\s*(?P<x_qe_warning>\S+.*?)\s*$",
                                        ),
                                    ],
                                 ),
                             ],
                          ),
                          SM(name='md_starting_charge_negative_spin_ignore', repeats=True,
                             startReStr=(r"\s*Check: negative starting charge=\(component\d\):?\s*" + RE_f + r"\s*$"),
                          ),
                          SM(name="md_new_old_atomic_charge_approximation",
                             startReStr=r"\s*NEW-OLD atomic charge density approx\. for the potential\s*$",
                             fixedStartValues={ 'x_qe_t_md_extrapolation_charge': 'atomic' }
                          ),
                          SM(name="md_first_order_charge_approximation",
                             startReStr=r"\s*first order charge density extrapolation\s*$",
                             fixedStartValues={ 'x_qe_t_md_extrapolation_charge': 'first-order' }
                          ),
                          SM(name="md_second_order_charge_approximation",
                             startReStr=r"\s*second order charge density extrapolation\s*$",
                             fixedStartValues={ 'x_qe_t_md_extrapolation_charge': 'second-order' }
                          ),
                          SM(name="md_new_kpoints_info_old",
                             # seen in old espresso
                             startReStr=r"\s*NEW K-POINTS\s*$",
                             subMatchers=[
                                 SM(name="md_new_k_kpoints", repeats=True,
                                    startReStr=(
                                        r"\s*" + QeC.re_vec("x_qe_t_md_k_info_vec", 'usrTpiba') +
                                        r"\s+(?P<x_qe_t_md_k_info_wk>" + RE_f + r")\s*$"),
                                 ),
                             ],
                          ),
                          SM(name="md_new_kpoints_info",
                             startReStr=r"\s*NEW k-points:\s*$",
                             subMatchers=[
                                 SM(name="md_new_kpoints_info_kpoints", repeats=True,
                                    startReStr=(
                                        r"\s*k\(\s*(?P<x_qe_t_md_k_info_ik>\d+)\s*\)\s*=\s*\(\s*" +
                                        QeC.re_vec("x_qe_t_md_k_info_vec", "usrTpiba") +
                                        r"\s*\),\s*wk\s*=\s*(?P<x_qe_t_md_k_info_wk>" + RE_f + r")\s*$"),
                                 ),
                             ],
                          ),
                          SM(name='md_core_charge_check',
                             startReStr=(r"\s*Check: negative/imaginary core charge\s*=\s*(?P<x_qe_t_md_core_charge_negative>" +
                                         RE_f + r")\s*(?P<x_qe_t_md_core_charge_imaginary>" + RE_f + r")\s*$")
                          ),
                          SM(name='md_starting_charge_negative_new',
                             startReStr=(r"\s*Check: negative starting charge=\s*(?P<x_qe_t_md_starting_charge_negative_new>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name='md_starting_charge_negative_spin_ignore', repeats=True,
                             startReStr=(r"\s*Check: negative starting charge=\(component\d\):?\s*" + RE_f + r"\s*$"),
                          ),
                          SM(name='md_starting_rho_new',
                             startReStr=(r"\s*negative rho \(up, down\):\s*(?P<x_qe_t_md_starting_charge_negative_new_up>" + RE_f +
                                         r")\s*(?P<x_qe_t_md_starting_charge_negative_new_down>" + RE_f + r")\s*$"),
                          ),
                          SM(name='extrapolated_charge',
                             startReStr=(r"\s*extrapolated charge\s*(?P<x_qe_t_md_starting_charge>" + RE_f +
                                         r")\s*,\s*renormalised to\s*(?P<x_qe_t_md_starting_charge_renormalized>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name="md_wave_function_extrapolation2",
                             startReStr=r"\s*(?P<x_qe_t_md_extrapolation_wfc>\S+ order)\s*wave-functions extrapolation\s*$",
                             subMatchers=[
                                 SM(name='warning_header',
                                    startReStr=r"\s*Message from extrapolate_wfcs:?\s*$",
                                    subMatchers=[
                                        SM(name='warning_msg',
                                           startReStr=r"\s*(?P<x_qe_warning>\S+.*?)\s*$",
                                        ),
                                    ],
                                 ),
                             ],
                          ),
                          # this does not belong to next scf!
                          SM(name='exx_refine',
                             startReStr=r"\s*EXX: now go back to refine exchange calculation\s*$",
                             # this applies to the _next_ single_config_calc, need to write it there
                             adHoc=lambda p: self.setTmp('exx_refine', True)
                          ),
                          SM(name="md_write_datafile_cputime",
                             startReStr=(r"\s*total cpu time spent up to now is\s*(?P<x_qe_t_md_write_datafile_cputime>" + RE_f +
                                         r")\s*secs\s*$"),
                          ),
                          SM(name="md_write_datafile_mem_dynamical",
                             startReStr=(
                                 r"\s*per-process dynamical memory:\s*(?P<x_qe_t_md_write_datafile_mem_dynamical__megabyte>" +
                                 RE_f + r")\s*Mb\s*$"),
                          ),
                          SM(name='md_martyna_tuckerman_parameters',
                             startReStr=(r"\s*alpha, beta MT =\s*(?P<x_qe_t_md_isolated_system_method_martyna_tuckerman_alpha>" + RE_f + r")\s*" +
                                         r"(?P<x_qe_t_md_isolated_system_method_martyna_tuckerman_beta>" + RE_f + r")\s*$"),
                          ),
                       ],
                   ),
               ],
            ),
            SM(name="band_structure_calculation",
               startReStr=r"\s*Band Structure Calculation\s*$",
               sections = ['section_single_configuration_calculation'],
               subMatchers = [
               ] + self.SMs_diagonalization(namespace='bands_diagonalization') + [
                   SM(name='cputime_nscf_msg',
                      startReStr=(r"\s*total cpu time spent up to now is\s*(?P<time_single_configuration_calculation_cpu1_end>" + RE_f +
                                  r")\s*secs\s*$"),
                   ),
                   SM(name='end_band_structure_calculation',
                      startReStr=r"\s*End of band structure calculation\s*$",
                      subMatchers=self.SMs_bands() + [
                          SM(name='bands_spinBS', repeats=True,
                              startReStr=r"\s*-+\s*SPIN\s+(?P<x_qe_t_spin_channel>UP|DOWN)\s*-+\s*$",
                              adHoc=self.adHoc_bands_spin,
                              subMatchers=self.SMs_bands(suffix='_spinBS'),
                          ),
                          SM(name='end_bands_nobands',
                             startReStr=r"\s*(?P<x_qe_warning>Number of k-points >= \d+: set verbosity='high' to print the bands\.)",
                          ),
                          SM(name='highest_occupiedBS',
                             startReStr=(r"\s*highest occupied level \(ev\):\s*(?P<x_qe_t_energy_reference_highest_occupied__eV>" + RE_f +
                                         r")\s*$"),
                          ),
                          SM(name='highest_occupied_lowest_unoccupiedBS',
                             startReStr=(r"\s*highest occupied, lowest unoccupied level \(ev\):\s*" +
                                         r"(?P<x_qe_t_energy_reference_highest_occupied__eV>" + RE_f + r")\s*" +
                                         r"(?P<x_qe_t_energy_reference_lowest_unoccupied__eV>" + RE_f + r")\s*$"),
                          ),
                          SM(name='e_fermiBS',
                             startReStr=(r"\s*the Fermi energy is\s*(?P<x_qe_t_energy_reference_fermi__eV>" + RE_f + ")\s*ev\s*$"),
                          ),
                          SM(name='e_fermi_spinBS',
                             startReStr=(r"\s*the spin up/dw Fermi energies are\s*" +
                                         r"(?P<x_qe_t_energy_reference_fermi_up__eV>" + RE_f + ")\s*" +
                                         r"(?P<x_qe_t_energy_reference_fermi_down__eV>" + RE_f + ")\s*ev\s*$"),
                          ),
                          SM(name="write_datafileBS",
                             startReStr=r"\s*Writing output data file\s*(?P<x_qe_output_datafile>.*?)\s*$",
                             subMatchers=[
                                SM(name="warning_save_mgga2",
                                   startReStr=r"\s*Warning:\s*(?P<x_qe_warning>cannot save meta-gga kinetic terms: not implemented\.)\s*$",
                                ),
                             ],
                          ),
                      ],
                   ),
               ]
            ),
            SM(name="profiling_protector",
               startReStr=r"\s*(?:PWSCF|init_run)\s*:\s*\d+.*$",
               forwardMatch=True,
               subMatchers=[
                   SM(name="profiling", repeats=True,
                      # ugly: 3 SMs in one...
                      startReStr=(r"\s*(?:" +
                                  r"Called by\s*(?P<x_qe_t_profile_caller>\S+?):?" + r'|' +
                                  r"(?P<x_qe_t_profile_category>.*?)\s*routines:?" + r'|' +
                                  r"(?P<x_qe_t_profile_function>\S+)\s*:\s*" +
                                  r"(?:(?P<x_qe_t_profile_cputime__strQeTimespan>.*)\s*(?:CPU\s*time\s*,|CPU)\s*)?" +
                                  r"(?:(?P<x_qe_t_profile_walltime__strQeTimespan>.*)\s*[wW][aA][lL][lL](?:\s*[tT][iI][mM][eE])?\s*)?"
                                  r"(?:\(\s*(?P<x_qe_t_profile_ncalls>\d+)\s*calls\s*(?:,\s*\S+\s*s\s*avg\s*)?\)\s*)?" +
                                  r")\s*$"),
                      adHoc=self.adHoc_profiling_complete,
                   ),
               ],
            ),
        ]

if __name__ == "__main__":
    parser = QuantumEspressoParserPWSCF()
    parser.parse()
