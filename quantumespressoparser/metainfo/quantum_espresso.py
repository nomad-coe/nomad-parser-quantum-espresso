#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public

m_package = Package(
    name='quantum_espresso_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='quantum_espresso.nomadmetainfo.json'))


class x_qe_section_parallel(MSection):
    '''
    section for run-time parallization options of Quantum Espresso
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_section_parallel'))

    x_qe_nthreads = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of OpenMP threads
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nthreads'))

    x_qe_nproc = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of MPI ranks
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nproc'))

    x_qe_npool = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of K-Point pools
        ''',
        a_legacy=LegacyDefinition(name='x_qe_npool'))


class x_qe_section_compile_options(MSection):
    '''
    section for compile-time options of Quantum Espresso
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_section_compile_options'))

    x_qe_ntypx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of different atom species
        ''',
        a_legacy=LegacyDefinition(name='x_qe_ntypx'))

    x_qe_ndmx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum dimension of radial grid (Pseudopotential)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_ndmx'))

    x_qe_npk = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of k-points
        ''',
        a_legacy=LegacyDefinition(name='x_qe_npk'))

    x_qe_lmaxx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum non local angular momentum (Pseudopotential)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_lmaxx'))

    x_qe_nbrx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of beta functions (Pseudopotential)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nbrx'))

    x_qe_nqfx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of coefficients in Q smoothing (Pseudopotential)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nqfx'))

    x_qe_nchix = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of atomic wavefunctions per Pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nchix'))

    x_qe_compile_parallel_version = Quantity(
        type=str,
        shape=[],
        description='''
        Parallelization compile-time options
        ''',
        a_legacy=LegacyDefinition(name='x_qe_compile_parallel_version'))


class x_qe_t_section_pp_report(MSection):
    '''
    section to collect 'pseudopotential report' information in new QE, used only for
    'old', non-UPF pseudopotentials
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_t_section_pp_report'))

    x_qe_t_pp_report_species = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: PP report: species number
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_report_species'))

    x_qe_t_pp_report_version = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: PP report: pp version
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_report_version'))

    x_qe_t_pp_report_line = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: PP report: parsed line
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_report_line'))


class x_qe_t_section_pp_warning(MSection):
    '''
    section to collect 'pseudopotential warning' information in old QE
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_t_section_pp_warning'))

    x_qe_t_pp_warning_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: pp index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_warning_idx'))

    x_qe_t_pp_warning_filename = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: filename
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_warning_filename'))

    x_qe_t_pp_warning_wfcidx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: pseudo-wavefunction index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_warning_wfcidx'))

    x_qe_t_pp_warning_wfclabel = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: pseudo-wavefunction label
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_warning_wfclabel'))

    x_qe_t_pp_warning_wfcnorm = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: pseudo-wavefunction original norm
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_warning_wfcnorm'))


class x_qe_t_section_pseudopotential(MSection):
    '''
    pseudo-section for collecting pseudopotential data (atomic number lookup requires
    table printed later in output)
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_t_section_pseudopotential'))

    x_qe_t_pp_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Index of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_idx'))

    x_qe_t_pp_ndmx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Radial grid of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_ndmx'))

    x_qe_t_pp_label = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Label of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_label'))

    x_qe_t_pp_filename = Quantity(
        type=str,
        shape=[],
        description='''
        Filename of pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_filename'))

    x_qe_t_pp_type = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Type of pseudopotential, e.g. 'Norm-conserving' or 'Ultrasoft'
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_type'))

    x_qe_t_pp_md5sum = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: MD5 checksum of pseudopotential file
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_md5sum'))

    x_qe_t_pp_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: comment about pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_comment'))

    x_qe_t_pp_integral_ndirections = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: number of integration directions (PAW)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_integral_ndirections'))

    x_qe_t_pp_integral_lmax_exact = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: maximum l for which integration is exact (PAW)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_integral_lmax_exact'))

    x_qe_t_pp_augmentation_shape = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: shape of augmentation charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_augmentation_shape'))

    x_qe_t_pp_valence = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Number of Valence electrons in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_valence'))

    x_qe_t_pp_nbeta = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Number of beta functions in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_nbeta'))

    x_qe_t_pp_l_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: beta function l index in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_l_idx'))

    x_qe_t_pp_l = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: beta function l in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_l'))

    x_qe_t_pp_ncoefficients = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Number of coefficients in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_ncoefficients'))

    x_qe_t_rinner = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Inner Radii of pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_rinner'))


class x_qe_section_scf_diagonalization(MSection):
    '''
    section for diagonalization info in QE scf iterations
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_section_scf_diagonalization'))

    x_qe_t_scf_diagonalization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Diagonalization algorithm
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_scf_diagonalization_algorithm'))

    x_qe_scf_diagonalization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Diagonalization algorithm
        ''',
        a_legacy=LegacyDefinition(name='x_qe_scf_diagonalization_algorithm'))

    x_qe_scf_diagonalization_warn_n_unconverged_eigenvalues = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of uncoverged eigenvalues (Warning)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_scf_diagonalization_warn_n_unconverged_eigenvalues'))

    x_qe_scf_diagonalization_c_bands_n_unconverged_eigenvalues = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of uncoverged eigenvalues (Warning from function c_bands)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_scf_diagonalization_c_bands_n_unconverged_eigenvalues'))

    x_qe_scf_diagonalization_ethr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Convergence Threshold in scf diagonalization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_scf_diagonalization_ethr'))

    x_qe_scf_diagonalization_iteration_avg = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Average of iterations in scf diagonalization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_scf_diagonalization_iteration_avg'))


class x_qe_section_bands_diagonalization(MSection):
    '''
    section for diagonalization info in QE band structure calculation
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_section_bands_diagonalization'))

    x_qe_t_bands_diagonalization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Diagonalization algorithm
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_bands_diagonalization_algorithm'))

    x_qe_bands_diagonalization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Diagonalization algorithm
        ''',
        a_legacy=LegacyDefinition(name='x_qe_bands_diagonalization_algorithm'))

    x_qe_bands_diagonalization_warn_n_unconverged_eigenvalues = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of uncoverged eigenvalues (Warning)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_bands_diagonalization_warn_n_unconverged_eigenvalues'))

    x_qe_bands_diagonalization_c_bands_n_unconverged_eigenvalues = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of uncoverged eigenvalues (Warning from function c_bands)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_bands_diagonalization_c_bands_n_unconverged_eigenvalues'))

    x_qe_bands_diagonalization_ethr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Convergence Threshold in bands diagonalization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_bands_diagonalization_ethr'))

    x_qe_bands_diagonalization_iteration_avg = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Average of iterations in bands diagonalization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_bands_diagonalization_iteration_avg'))


class x_qe_t_section_input_occupations(MSection):
    '''
    Temporary: Section for User-specified band occupations
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_qe_t_section_input_occupations'))

    x_qe_t_input_occupations_spin = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: User-specified band occupations, spin channel
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_input_occupations_spin'))

    x_qe_t_input_occupations = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: User-specified band occupations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_input_occupations'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_qe_extra_SCF = Quantity(
        type=bool,
        shape=[],
        description='''
        Extra SCF without electronic history at the end of relaxation. Triggered in
        magnetic simulations when relax converges to non-magnetic solution
        ''',
        a_legacy=LegacyDefinition(name='x_qe_extra_SCF'))

    x_qe_t_spin_channel = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for spin channel
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_spin_channel'))

    x_qe_t_k_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_x'))

    x_qe_t_k_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_y'))

    x_qe_t_k_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_z'))

    x_qe_t_k_pw = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: number of plane waves for this k-point
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_pw'))

    x_qe_t_k_point_energies = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: k-point band energies
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_point_energies'))

    x_qe_energy_total_harris_foulkes_estimate = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Harris-Foulkes estimate of total energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_total_harris_foulkes_estimate'))

    x_qe_energy_total_accuracy_estimate = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Accuracy estimate of total energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_total_accuracy_estimate'))

    x_qe_energy_exchange_error_estimate = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Estimated error on exchange
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_exchange_error_estimate'))

    x_qe_energy_exchange_average_fock_potential = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Averaged Fock potential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_exchange_average_fock_potential'))

    x_qe_energy_fock = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Fock energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_fock'))

    x_qe_energy_total_paw_all_electron = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        All-electron total energy from PAW
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_total_paw_all_electron'))

    x_qe_t_energy_reference_highest_occupied = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Energy of highest occupied state
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_reference_highest_occupied'))

    x_qe_t_energy_reference_lowest_unoccupied = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Energy of lowest unoccupied state
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_reference_lowest_unoccupied'))

    x_qe_t_energy_reference_fermi = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Fermi Energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_reference_fermi'))

    x_qe_t_energy_reference_fermi_up = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Fermi Energy (spin up)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_reference_fermi_up'))

    x_qe_t_energy_reference_fermi_down = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Fermi Energy (spin down)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_reference_fermi_down'))

    x_qe_t_energy_decomposition_name = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Total energy decomposition: contribution name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_decomposition_name'))

    x_qe_t_energy_decomposition_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Total energy decomposition: contribution value
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_energy_decomposition_value'))

    x_qe_energy_decomposition_name = Quantity(
        type=str,
        shape=['x_qe_number_of_energy_components'],
        description='''
        Total energy decomposition: contribution name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_decomposition_name'))

    x_qe_energy_decomposition_value = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_number_of_energy_components'],
        description='''
        Total energy decomposition: contribution value
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_decomposition_value'))

    x_qe_magnetization_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Total per-cell magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_magnetization_total'))

    x_qe_magnetization_absolute = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Absolute per-cell magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_magnetization_absolute'))

    x_qe_convergence_iterations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of iterations after which self-consistency has been achieved
        ''',
        a_legacy=LegacyDefinition(name='x_qe_convergence_iterations'))

    x_qe_exx_refine = Quantity(
        type=bool,
        shape=[],
        description='''
        Flag: Exact-exchange refinement is active
        ''',
        a_legacy=LegacyDefinition(name='x_qe_exx_refine'))

    x_qe_exx_self_consistency = Quantity(
        type=bool,
        shape=[],
        description='''
        Exact-exchange has been reached (flag)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_exx_self_consistency'))

    x_qe_output_datafile = Quantity(
        type=str,
        shape=[],
        description='''
        Output datafile
        ''',
        a_legacy=LegacyDefinition(name='x_qe_output_datafile'))

    x_qe_t_force_atom_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_force_atom_idx'))

    x_qe_t_force_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_force_x'))

    x_qe_t_force_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_force_y'))

    x_qe_t_force_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_force_z'))

    x_qe_t_dispersion_force_atom_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dispersion_force_atom_idx'))

    x_qe_atom_dispersion_force = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_atom_dispersion_force'))

    x_qe_t_dispersion_force_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dispersion_force_x'))

    x_qe_t_dispersion_force_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dispersion_force_y'))

    x_qe_t_dispersion_force_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dispersion_force_z'))

    x_qe_dispersion_force_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dispersion_force_total'))

    x_qe_force_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_force_total'))

    x_qe_force_total_scf_correction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_force_total_scf_correction'))

    x_qe_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pressure'))

    x_qe_t_stress_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_stress_x'))

    x_qe_t_stress_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_stress_y'))

    x_qe_t_stress_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_stress_z'))

    x_qe_stress_unimplemented = Quantity(
        type=str,
        shape=[],
        description='''
        Reason why stress tensor is not implemented
        ''',
        a_legacy=LegacyDefinition(name='x_qe_stress_unimplemented'))

    x_qe_t_md_iteration = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: MD step: iteration number
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_iteration'))

    x_qe_t_projected_velocity = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: MD step: projected velocity
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_projected_velocity'))

    x_qe_t_md_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        MD step: time
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_time'))

    x_qe_t_md_vec_a_units = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for new direct lattice vectors (vc-relax), units
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_vec_a_units'))

    x_qe_t_md_vec_a_alat = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new direct lattice vectors (vc-relax), lattice parameter a
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_vec_a_alat'))

    x_qe_t_md_vec_a_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new direct lattice vectors (vc-relax), x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_vec_a_x'))

    x_qe_t_md_vec_a_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new direct lattice vectors (vc-relax), y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_vec_a_y'))

    x_qe_t_md_vec_a_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new direct lattice vectors (vc-relax), z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_vec_a_z'))

    x_qe_t_md_atom_positions_units = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax), units
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_positions_units'))

    x_qe_t_md_atom_positions_units_vcsmd = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax via VCSMD), units
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_positions_units_vcsmd'))

    x_qe_t_md_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax), atom labels
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_labels'))

    x_qe_t_md_atom_positions_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax), x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_positions_x'))

    x_qe_t_md_atom_positions_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax), y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_positions_y'))

    x_qe_t_md_atom_positions_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new atom positions (MD, (vc-)relax), z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_positions_z'))

    x_qe_t_md_atom_free_x = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for new atom fixed flag (MD, (vc-)relax), x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_free_x'))

    x_qe_t_md_atom_free_y = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for new atom fixed flag (MD, (vc-)relax), y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_free_y'))

    x_qe_t_md_atom_free_z = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for new atom fixed flag (MD, (vc-)relax), z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_free_z'))

    x_qe_t_new_nat2_distance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new 2-atom distance (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_new_nat2_distance'))

    x_qe_t_md_atom_mass_label = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for MD setup, atom mass, labels
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_mass_label'))

    x_qe_t_md_atom_mass_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD setup, atom mass, value
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_atom_mass_value'))

    x_qe_t_md_timestep_size = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD setup, timestep size
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_timestep_size'))

    x_qe_t_md_kinetic_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD, kinetic energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_kinetic_energy'))

    x_qe_t_md_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD, temperature
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_temperature'))

    x_qe_t_md_total_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD, total energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_total_energy'))

    x_qe_t_md_ekin_etot = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for MD, sum of energies
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_ekin_etot'))

    x_qe_t_md_linear_momentum_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for linear momentum (MD, (vc-)relax), x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_linear_momentum_x'))

    x_qe_t_md_linear_momentum_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for linear momentum (MD, (vc-)relax), y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_linear_momentum_y'))

    x_qe_t_md_linear_momentum_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for linear momentum (MD, (vc-)relax), z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_linear_momentum_z'))

    x_qe_t_md_write_datafile_cputime = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for cpu time after write-datafile (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_write_datafile_cputime'))

    x_qe_t_md_write_datafile_mem_dynamical = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for dynamical memory after write-datafile (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_write_datafile_mem_dynamical'))

    x_qe_t_md_extrapolation_charge = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for charge extrapolation scheme (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_extrapolation_charge'))

    x_qe_t_md_extrapolation_wfc = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for wave function extrapolation scheme (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_extrapolation_wfc'))

    x_qe_t_md_starting_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for extrapolated starting charge (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge'))

    x_qe_t_md_starting_charge_renormalized = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for extrapolated starting charge, renormalized (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge_renormalized'))

    x_qe_t_md_max_steps_reached = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for max_steps-reached flag (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_max_steps_reached'))

    x_qe_t_md_end = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for end-of-md flag (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_end'))

    x_qe_t_md_diffusion_atomidx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary storage for diffusion coeffients (MD), atom index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_diffusion_atomidx'))

    x_qe_t_md_diffusion_coefficient = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for diffusion coeffients (MD), atom coeffient
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_diffusion_coefficient'))

    x_qe_t_md_diffusion_coefficient_mean = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for diffusion coeffients (MD), mean coeffient
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_diffusion_coefficient_mean'))

    x_qe_t_md_bfgs_scf_cycles = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary storage for number of scf cycles (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_scf_cycles'))

    x_qe_t_md_bfgs_steps = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary storage for number of steps (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_steps'))

    x_qe_t_md_bfgs_energy_old = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for 'old' energy (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_energy_old'))

    x_qe_t_md_bfgs_energy_new = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for 'new' energy (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_energy_new'))

    x_qe_t_md_bfgs_enthalpy_old = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for 'old' enthalpy (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_enthalpy_old'))

    x_qe_t_md_bfgs_enthalpy_new = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for 'new' enthalpy (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_enthalpy_new'))

    x_qe_t_md_bfgs_case = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for BFGS case, energy comparison (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_case'))

    x_qe_t_md_bfgs_reset = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for BFGS history reset reason (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_reset'))

    x_qe_t_md_bfgs_trust_new = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new trust radius (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_trust_new'))

    x_qe_t_md_bfgs_conv_thr_new = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new electronic convergence threshold (relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_conv_thr_new'))

    x_qe_t_md_starting_charge_negative_old = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for old negative starting charge (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge_negative_old'))

    x_qe_t_md_starting_charge_negative_new = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new negative starting charge (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge_negative_new'))

    x_qe_t_md_starting_charge_negative_new_up = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new negative starting charge (MD, (vc-)relax), spin up
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge_negative_new_up'))

    x_qe_t_md_starting_charge_negative_new_down = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new negative starting charge (MD, (vc-)relax), spin down
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_starting_charge_negative_new_down'))

    x_qe_t_md_bfgs_converged = Quantity(
        type=bool,
        shape=[],
        description='''
        Temporary storage for 'converged' flag ((vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_converged'))

    x_qe_t_md_bfgs_converged_criteria = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for converged criteria ((vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_converged_criteria'))

    x_qe_t_md_bfgs_final_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for final energy ((vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_final_energy'))

    x_qe_t_md_bfgs_final_enthalpy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for final enthalpy ((vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_bfgs_final_enthalpy'))

    x_qe_t_md_new_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for new cell volume ((vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_new_volume'))

    x_qe_t_md_isolated_system_method_martyna_tuckerman_alpha = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD: Isolated system with Martyna-Tuckerman method, parameter alpha
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_isolated_system_method_martyna_tuckerman_alpha'))

    x_qe_t_md_isolated_system_method_martyna_tuckerman_beta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD: Isolated system with Martyna-Tuckerman method, parameter beta
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_isolated_system_method_martyna_tuckerman_beta'))

    x_qe_t_md_core_charge_negative = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD: QE check: negative core charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_core_charge_negative'))

    x_qe_t_md_core_charge_imaginary = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD: QE check: imaginary core charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_core_charge_imaginary'))

    x_qe_t_relax_converged_steps = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary Relax: number of steps after which structure relaxation converged
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_relax_converged_steps'))

    x_qe_t_relax_final_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary Relax: final energy in relaxation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_relax_final_energy'))

    x_qe_t_relax_threshold_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary Relax: convergence threshold on energy in relaxation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_relax_threshold_energy'))

    x_qe_t_relax_threshold_force = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary Relax: convergence threshold on force components in relaxation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_relax_threshold_force'))

    x_qe_t_relax_threshold_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary Relax: convergence threshold on pressure in relaxation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_relax_threshold_pressure'))

    x_qe_t_md_k_info_ik = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary MD storage for k-point info, k-index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_k_info_ik'))

    x_qe_t_md_k_info_vec_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD storage for k-point info, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_k_info_vec_x'))

    x_qe_t_md_k_info_vec_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD storage for k-point info, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_k_info_vec_y'))

    x_qe_t_md_k_info_vec_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD storage for k-point info, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_k_info_vec_z'))

    x_qe_t_md_k_info_wk = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary MD storage for k-point info, weight
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_md_k_info_wk'))

    x_qe_section_bands_diagonalization = SubSection(
        sub_section=SectionProxy('x_qe_section_bands_diagonalization'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_section_bands_diagonalization'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_qe_program_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of program from Quantum Espresso suite
        ''',
        a_legacy=LegacyDefinition(name='x_qe_program_name'))

    x_qe_input_filename = Quantity(
        type=str,
        shape=[],
        description='''
        Filename input was read from
        ''',
        a_legacy=LegacyDefinition(name='x_qe_input_filename'))

    x_qe_t_warning = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Warnings from Quantum Espresso
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_warning'))

    x_qe_warning = Quantity(
        type=str,
        shape=[],
        description='''
        Warnings from Quantum Espresso
        ''',
        a_legacy=LegacyDefinition(name='x_qe_warning'))

    x_qe_profile_caller = Quantity(
        type=str,
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: caller name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_caller'))

    x_qe_profile_category = Quantity(
        type=str,
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: category
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_category'))

    x_qe_profile_function = Quantity(
        type=str,
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: function name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_function'))

    x_qe_profile_cputime = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: cputime spent in function
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_cputime'))

    x_qe_profile_walltime = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: wallclock time spent in function
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_walltime'))

    x_qe_profile_ncalls = Quantity(
        type=np.dtype(np.int32),
        shape=['x_qe_number_of_profiling_entries'],
        description='''
        QE profiling: how often was function called
        ''',
        a_legacy=LegacyDefinition(name='x_qe_profile_ncalls'))

    x_qe_t_profile_function = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: QE profiling: function name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_function'))

    x_qe_t_profile_cputime = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: QE profiling: cputime spent in function
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_cputime'))

    x_qe_t_profile_walltime = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: QE profiling: wallclock time spent in function
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_walltime'))

    x_qe_t_profile_ncalls = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: QE profiling: how often was function called
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_ncalls'))

    x_qe_t_profile_caller = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: QE profiling: who was the caller
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_caller'))

    x_qe_t_profile_caller_list = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: QE profiling: who was the caller (list for each function)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_caller_list'))

    x_qe_t_profile_category = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: QE profiling: category
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_category'))

    x_qe_t_profile_category_list = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: QE profiling: category (list for each function)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_profile_category_list'))

    x_qe_input_positions_cell_dirname = Quantity(
        type=str,
        shape=[],
        description='''
        Directory where initial atom_positions and simulation_cell were read from
        ''',
        a_legacy=LegacyDefinition(name='x_qe_input_positions_cell_dirname'))

    x_qe_input_potential_recalculated_file = Quantity(
        type=str,
        shape=[],
        description='''
        File that was used to recalculate initial potential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_input_potential_recalculated_file'))

    x_qe_section_parallel = SubSection(
        sub_section=SectionProxy('x_qe_section_parallel'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_section_parallel'))

    x_qe_section_compile_options = SubSection(
        sub_section=SectionProxy('x_qe_section_compile_options'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_section_compile_options'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_qe_t_species_dispersion_correction_label = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: DFT-D species label
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_species_dispersion_correction_label'))

    x_qe_t_species_dispersion_correction_vdw_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: DFT-D species vdW radius
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_species_dispersion_correction_vdw_radius'))

    x_qe_t_species_dispersion_correction_C6 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: DFT-D species C6 coefficient
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_species_dispersion_correction_C6'))

    x_qe_dispersion_correction = Quantity(
        type=bool,
        shape=[],
        description='''
        Calculation includes semi-empirical DFT-D dispersion correction
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dispersion_correction'))

    x_qe_gamma_algorithms = Quantity(
        type=bool,
        shape=[],
        description='''
        Usage of gamma-only optimized algorithms
        ''',
        a_legacy=LegacyDefinition(name='x_qe_gamma_algorithms'))

    x_qe_exx_grid_same_as_k_grid = Quantity(
        type=bool,
        shape=[],
        description='''
        Exact-exchange k+q grid is the same as k grid (flag)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_exx_grid_same_as_k_grid'))

    x_qe_diagonalization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Algorithm used in subspace diagonalization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_diagonalization_algorithm'))

    x_qe_sticks_sum_dense = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_dense'))

    x_qe_sticks_sum_smooth = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_smooth'))

    x_qe_sticks_sum_PW = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_PW'))

    x_qe_sticks_sum_G_dense = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_G_dense'))

    x_qe_sticks_sum_G_smooth = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_G_smooth'))

    x_qe_sticks_sum_G_PW = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_sum_G_PW'))

    x_qe_sticks_tot_dense = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_tot_dense'))

    x_qe_sticks_tot_smooth = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_tot_smooth'))

    x_qe_sticks_tot_PW = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_tot_PW'))

    x_qe_sticks_old = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_sticks_old'))

    x_qe_t_species_integration_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: radius used to integrate charge/magnetization over (per species)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_species_integration_radius'))

    x_qe_t_species_integration_radius_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: radius used to integrate charge/magnetization over (per species),
        species index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_species_integration_radius_idx'))

    x_qe_fock_operator_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Cutoff for defining the direct-space grid used to compute Fock exchange in EXX
        ''',
        a_legacy=LegacyDefinition(name='x_qe_fock_operator_cutoff'))

    x_qe_t_xc_functional_shortname_enforced = Quantity(
        type=str,
        shape=[],
        description='''
        Short name of User-enforced XC functional; overrides the setting implied by the
        pseudopotentials
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_xc_functional_shortname_enforced'))

    x_qe_potential_convergence_threshold = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Convergence threshold for potentials
        ''',
        a_legacy=LegacyDefinition(name='x_qe_potential_convergence_threshold'))

    x_qe_potential_mixing_beta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Mixing scheme: parameter beta
        ''',
        a_legacy=LegacyDefinition(name='x_qe_potential_mixing_beta'))

    x_qe_potential_mixing_iterations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Mixing scheme: number of previous iterations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_potential_mixing_iterations'))

    x_qe_potential_mixing_scheme = Quantity(
        type=str,
        shape=[],
        description='''
        Mixing scheme: type of mixing
        ''',
        a_legacy=LegacyDefinition(name='x_qe_potential_mixing_scheme'))

    x_qe_xc_functional_user_enforced = Quantity(
        type=bool,
        shape=[],
        description='''
        True if the user enforced setting the XC functional; overrides the setting implied
        by the pseudopotentials
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_functional_user_enforced'))

    x_qe_xc_functional_shortname = Quantity(
        type=str,
        shape=[],
        description='''
        Short name of XC functional used in calculation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_functional_shortname'))

    x_qe_xc_functional_num = Quantity(
        type=str,
        shape=[],
        description='''
        QE Index number representation of XC functional
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_functional_num'))

    x_qe_xc_iexch_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (density exchange component) in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_iexch_name'))

    x_qe_xc_icorr_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (density correlation component) in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_icorr_name'))

    x_qe_xc_igcx_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (gradient exchange component) in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcx_name'))

    x_qe_xc_igcc_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (gradient correlation component) in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcc_name'))

    x_qe_xc_imeta_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (meta-gga component) in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_imeta_name'))

    x_qe_xc_inlc_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional (Van-der-Waals non-local component) in Quantum Espresso
        context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_inlc_name'))

    x_qe_xc_iexch_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (density exchange component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_iexch_comment'))

    x_qe_xc_icorr_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (density correlation component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_icorr_comment'))

    x_qe_xc_igcx_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (gradient exchange component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcx_comment'))

    x_qe_xc_igcc_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (gradient correlation component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcc_comment'))

    x_qe_xc_imeta_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (meta-gga component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_imeta_comment'))

    x_qe_xc_inlc_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about XC functional (Van-der-Waals non-local component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_inlc_comment'))

    x_qe_xc_iexch = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (density exchange
        component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_iexch'))

    x_qe_xc_icorr = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (density
        correlation component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_icorr'))

    x_qe_xc_igcx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (gradient exchange
        component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcx'))

    x_qe_xc_igcc = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (gradient
        correlation component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_igcc'))

    x_qe_xc_imeta = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (meta-gga
        component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_imeta'))

    x_qe_xc_inlc = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Quantum Espresso internal code-specific index of XC functional (Van-der-Waals non-
        local component)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_inlc'))

    x_qe_t_exact_exchange_fraction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: store fraction of exact-exchange before defining section_xc_functionals
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_exact_exchange_fraction'))

    x_qe_exact_exchange_fraction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Fraction of exact-exchange in EXX-refinement
        ''',
        a_legacy=LegacyDefinition(name='x_qe_exact_exchange_fraction'))

    x_qe_md_max_steps = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum number of ionic+electronic steps in dynamics (MD/relax) calculation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_md_max_steps'))

    x_qe_berry_efield_direction = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Finite E-field: direction
        ''',
        a_legacy=LegacyDefinition(name='x_qe_berry_efield_direction'))

    x_qe_berry_efield_intensity = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Berry phase with E-field: intensity
        ''',
        a_legacy=LegacyDefinition(name='x_qe_berry_efield_intensity'))

    x_qe_berry_efield_strings_nk = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Berry phase with E-field: number of k-points in string
        ''',
        a_legacy=LegacyDefinition(name='x_qe_berry_efield_strings_nk'))

    x_qe_berry_efield_niter = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Berry phase with E-field: number of iterative cycles
        ''',
        a_legacy=LegacyDefinition(name='x_qe_berry_efield_niter'))

    x_qe_berry_efield = Quantity(
        type=bool,
        shape=[],
        description='''
        Berry phase with E-field: flag if berry-efield-calc was done
        ''',
        a_legacy=LegacyDefinition(name='x_qe_berry_efield'))

    x_qe_t_spin_orbit_magn = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: spin-orbit msg: magnetic mode (non-)collinear / (non-)magnetic
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_spin_orbit_magn'))

    x_qe_t_spin_orbit_mode = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: spin-orbit msg: with/without spin-orbit
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_spin_orbit_mode'))

    x_qe_spin_orbit = Quantity(
        type=bool,
        shape=[],
        description='''
        Spin-orbit coupling flag: with/without spin-orbit
        ''',
        a_legacy=LegacyDefinition(name='x_qe_spin_orbit'))

    x_qe_spin_noncollinear = Quantity(
        type=bool,
        shape=[],
        description='''
        Noncollinear spin mode
        ''',
        a_legacy=LegacyDefinition(name='x_qe_spin_noncollinear'))

    x_qe_t_pp_renormalized_filename = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential: filename
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_renormalized_filename'))

    x_qe_t_pp_renormalized_wfc = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_pp_renormalized_wfc'))

    x_qe_t_allocated_array_name = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: allocated arrays, name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_allocated_array_name'))

    x_qe_t_allocated_array_size = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: allocated arrays, size
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_allocated_array_size'))

    x_qe_t_allocated_array_dimensions = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: allocated arrays, dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_allocated_array_dimensions'))

    x_qe_allocated_array_name = Quantity(
        type=str,
        shape=['x_qe_allocated_arrays'],
        description='''
        Allocated arrays, name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_allocated_array_name'))

    x_qe_allocated_array_size = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_allocated_arrays'],
        description='''
        Allocated arrays, size
        ''',
        a_legacy=LegacyDefinition(name='x_qe_allocated_array_size'))

    x_qe_allocated_array_dimensions = Quantity(
        type=str,
        shape=['x_qe_allocated_arrays'],
        description='''
        Allocated arrays, dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_allocated_array_dimensions'))

    x_qe_t_temporary_array_name = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: temporary arrays, name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_temporary_array_name'))

    x_qe_t_temporary_array_size = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: temporary arrays, size
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_temporary_array_size'))

    x_qe_t_temporary_array_dimensions = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: temporary arrays, dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_temporary_array_dimensions'))

    x_qe_temporary_array_name = Quantity(
        type=str,
        shape=['x_qe_temporary_arrays'],
        description='''
        Temporary arrays, name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_temporary_array_name'))

    x_qe_temporary_array_size = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_temporary_arrays'],
        description='''
        Temporary arrays, size
        ''',
        a_legacy=LegacyDefinition(name='x_qe_temporary_array_size'))

    x_qe_temporary_array_dimensions = Quantity(
        type=str,
        shape=['x_qe_temporary_arrays'],
        description='''
        Temporary arrays, dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_temporary_array_dimensions'))

    x_qe_core_charge_negative = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        QE check: negative core charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_core_charge_negative'))

    x_qe_core_charge_imaginary = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        QE check: imaginary core charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_core_charge_imaginary'))

    x_qe_core_charge_realspace = Quantity(
        type=bool,
        shape=[],
        description='''
        QE flag: core charge treated in real space
        ''',
        a_legacy=LegacyDefinition(name='x_qe_core_charge_realspace'))

    x_qe_starting_density_file = Quantity(
        type=str,
        shape=[],
        description='''
        Starting density from file
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_density_file'))

    x_qe_starting_potential = Quantity(
        type=str,
        shape=[],
        description='''
        Starting potential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_potential'))

    x_qe_starting_charge_negative = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Starting charge (warning about negative starting charge)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_charge_negative'))

    x_qe_starting_charge_negative_up = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Starting charge up (warning about negative starting charge)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_charge_negative_up'))

    x_qe_starting_charge_negative_down = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Starting charge down (warning about negative starting charge)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_charge_negative_down'))

    x_qe_starting_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Starting charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_charge'))

    x_qe_starting_charge_renormalized = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Starting charge, renormalized
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_charge_renormalized'))

    x_qe_starting_wfc = Quantity(
        type=str,
        shape=[],
        description='''
        Starting Wave functions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_starting_wfc'))

    x_qe_time_setup_cpu1_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CPU time, setup up until first iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_time_setup_cpu1_end'))

    x_qe_per_process_mem = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Per-process dynamical memory
        ''',
        a_legacy=LegacyDefinition(name='x_qe_per_process_mem'))

    x_qe_isolated_system_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used if system is assumed to be isolated
        ''',
        a_legacy=LegacyDefinition(name='x_qe_isolated_system_method'))

    x_qe_isolated_system_method_martyna_tuckerman_alpha = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Isolated system with Martyna-Tuckerman method, parameter alpha
        ''',
        a_legacy=LegacyDefinition(name='x_qe_isolated_system_method_martyna_tuckerman_alpha'))

    x_qe_isolated_system_method_martyna_tuckerman_beta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Isolated system with Martyna-Tuckerman method, parameter beta
        ''',
        a_legacy=LegacyDefinition(name='x_qe_isolated_system_method_martyna_tuckerman_beta'))

    x_qe_input_occupations = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_k_points', 'number_of_eigen_values'],
        description='''
        User-specified band occupations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_input_occupations'))

    x_qe_extrapolation_charge = Quantity(
        type=str,
        shape=[],
        description='''
        Charge extrapolation scheme (MD, (vc-)relax)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_extrapolation_charge'))

    x_qe_t_section_pp_report = SubSection(
        sub_section=SectionProxy('x_qe_t_section_pp_report'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_t_section_pp_report'))

    x_qe_t_section_pp_warning = SubSection(
        sub_section=SectionProxy('x_qe_t_section_pp_warning'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_t_section_pp_warning'))

    x_qe_t_section_pseudopotential = SubSection(
        sub_section=SectionProxy('x_qe_t_section_pseudopotential'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_t_section_pseudopotential'))

    x_qe_t_section_input_occupations = SubSection(
        sub_section=SectionProxy('x_qe_t_section_input_occupations'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_t_section_input_occupations'))


class section_method_atom_kind(public.section_method_atom_kind):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method_atom_kind'))

    x_qe_dispersion_correction_vdw_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        DFT-D species vdW radius
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dispersion_correction_vdw_radius'))

    x_qe_dispersion_correction_C6 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        DFT-D species C6 coefficient
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dispersion_correction_C6'))

    x_qe_species_integration_radius = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Radius used to integrate charge/magnetization over (per species)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_species_integration_radius'))

    x_qe_pp_renormalized_wfc = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: renormalized WFCs in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_renormalized_wfc'))

    x_qe_pp_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Index of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_idx'))

    x_qe_pp_label = Quantity(
        type=str,
        shape=[],
        description='''
        Label of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_label'))

    x_qe_pp_filename = Quantity(
        type=str,
        shape=[],
        description='''
        Filename of pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_filename'))

    x_qe_pp_type = Quantity(
        type=str,
        shape=[],
        description='''
        Type of pseudopotential, e.g. 'Norm-conserving' or 'Ultrasoft'
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_type'))

    x_qe_pp_md5sum = Quantity(
        type=str,
        shape=[],
        description='''
        MD5 checksum of pseudopotential file
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_md5sum'))

    x_qe_pp_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Comment about pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_comment'))

    x_qe_pp_integral_ndirections = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of integration directions (PAW)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_integral_ndirections'))

    x_qe_pp_integral_lmax_exact = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum l for which integration is exact (PAW)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_integral_lmax_exact'))

    x_qe_pp_augmentation_shape = Quantity(
        type=str,
        shape=[],
        description='''
        Shape of augmentation charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_augmentation_shape'))

    x_qe_pp_report_version = Quantity(
        type=str,
        shape=[],
        description='''
        Pseudopotential report: version of PP
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_report_version'))

    x_qe_pp_report_contents = Quantity(
        type=str,
        shape=[],
        description='''
        Pseudopotential report: contents of PP report
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_report_contents'))

    x_qe_pp_valence = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of Valence electrons in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_valence'))

    x_qe_pp_weight = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_weight'))

    x_qe_pp_ncoefficients = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of coefficients in pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_ncoefficients'))

    x_qe_rinner = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Inner Radii of pseudopotential
        ''',
        a_legacy=LegacyDefinition(name='x_qe_rinner'))

    x_qe_kind_mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Atomic mass of species
        ''',
        a_legacy=LegacyDefinition(name='x_qe_kind_mass'))

    x_qe_pp_ndmx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Radial grid of Pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_ndmx'))

    x_qe_pp_nbeta = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of beta functions in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_nbeta'))

    x_qe_pp_l_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Beta function l index in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_l_idx'))

    x_qe_pp_l = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Beta function l in pseudopotential on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_pp_l'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_qe_ibrav = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Bravais lattice index, constant during a run
        ''',
        a_legacy=LegacyDefinition(name='x_qe_ibrav'))

    x_qe_alat = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Lattice Parameter 'a', constant during a run and used as unit in other quantities
        ''',
        a_legacy=LegacyDefinition(name='x_qe_alat'))

    x_qe_cell_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Volume of unit cell
        ''',
        a_legacy=LegacyDefinition(name='x_qe_cell_volume'))

    x_qe_number_of_species = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of Atom species, a.k.a. unique Atom labels; a label may include symmetry-
        breaking suffices, e.g. 'Fe1' and 'Fe2', as some quantities can only prescribed
        per species and not per site
        ''',
        a_legacy=LegacyDefinition(name='x_qe_number_of_species'))

    x_qe_t_number_of_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Number of electrons in system
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_qe_t_number_of_electrons'))

    x_qe_t_number_of_electrons_up = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Number of electrons in system (spin up)
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_qe_t_number_of_electrons_up'))

    x_qe_t_number_of_electrons_down = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Number of electrons in system (spin down)
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_qe_t_number_of_electrons_down'))

    x_qe_number_of_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of Kohn-Sham states/bands
        ''',
        a_legacy=LegacyDefinition(name='x_qe_number_of_states'))

    x_qe_md_cell_mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Mass of cell in MD/relax calculation
        ''',
        a_legacy=LegacyDefinition(name='x_qe_md_cell_mass'))

    x_qe_t_celldm = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for QE cell dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_celldm'))

    x_qe_celldm = Quantity(
        type=np.dtype(np.float64),
        shape=[6],
        description='''
        QE cell dimensions
        ''',
        a_legacy=LegacyDefinition(name='x_qe_celldm'))

    x_qe_t_vec_supercell_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for supercell translation vector in fractional coordinates,
        x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_supercell_x'))

    x_qe_t_vec_supercell_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for supercell translation vector in fractional coordinates,
        y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_supercell_y'))

    x_qe_t_vec_supercell_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for supercell translation vector in fractional coordinates,
        z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_supercell_z'))

    x_qe_vec_supercell = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_number_of_supercell_translations', 3],
        description='''
        Supercell translation vector(s) in fractional coordinates
        ''',
        a_legacy=LegacyDefinition(name='x_qe_vec_supercell'))

    x_qe_supercell = Quantity(
        type=bool,
        shape=[],
        description='''
        Supercell flag
        ''',
        a_legacy=LegacyDefinition(name='x_qe_supercell'))

    x_qe_t_vec_a_units = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storage for direct lattice vectors, units
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_a_units'))

    x_qe_t_vec_a_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for direct lattice vectors, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_a_x'))

    x_qe_t_vec_a_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for direct lattice vectors, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_a_y'))

    x_qe_t_vec_a_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for direct lattice vectors, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_a_z'))

    x_qe_t_vec_b_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for reciprocal lattice vectors, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_b_x'))

    x_qe_t_vec_b_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for reciprocal lattice vectors, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_b_y'))

    x_qe_t_vec_b_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for reciprocal lattice vectors, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_vec_b_z'))

    x_qe_reciprocal_cell = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='1 / meter',
        description='''
        Reciprocal Lattice vectors (in Cartesian coordinates). The first index runs over
        the $x,y,z$ Cartesian coordinates, and the second index runs over the 3 lattice
        vectors.
        ''',
        categories=[public.configuration_core],
        a_legacy=LegacyDefinition(name='x_qe_reciprocal_cell'))

    x_qe_t_starting_magnetization_species = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Starting magnetic configuration: Species name
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_starting_magnetization_species'))

    x_qe_t_starting_magnetization_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: Starting magnetic configuration: Species magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_starting_magnetization_value'))

    x_qe_atom_starting_magnetization = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Starting magnetic configuration: Atom magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_atom_starting_magnetization'))

    x_qe_nsymm = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of detected symmetry operations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nsymm'))

    x_qe_nsymm_with_fractional_translation = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of detected symmetry operations including fractional translations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nsymm_with_fractional_translation'))

    x_qe_nsymm_ignored = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of ignored symmetry operations, due to uncommensurable fractional
        translations
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nsymm_ignored'))

    x_qe_t_symm_inversion = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Inversion symmetry
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_symm_inversion'))

    x_qe_symm_inversion = Quantity(
        type=bool,
        shape=[],
        description='''
        Inversion symmetry
        ''',
        a_legacy=LegacyDefinition(name='x_qe_symm_inversion'))

    x_qe_atom_idx = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms'],
        description='''
        Index of atom on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_atom_idx'))

    x_qe_t_atpos_units = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Units for atom position
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atpos_units'))

    x_qe_t_atom_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Index of atom on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atom_idx'))

    x_qe_t_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary: Label of atom on Espresso side
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atom_labels'))

    x_qe_t_atpos_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for atom position, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atpos_x'))

    x_qe_t_atpos_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for atom position, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atpos_y'))

    x_qe_t_atpos_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for atom position, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_atpos_z'))

    x_qe_nk = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        K-point info, number of k-points
        ''',
        a_legacy=LegacyDefinition(name='x_qe_nk'))

    x_qe_smearing_ngauss = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        K-point info, QE number represenation of smearing technique
        ''',
        a_legacy=LegacyDefinition(name='x_qe_smearing_ngauss'))

    x_qe_smearing_kind = Quantity(
        type=str,
        shape=[],
        description='''
        K-point info, QE string represenation of smearing technique
        ''',
        a_legacy=LegacyDefinition(name='x_qe_smearing_kind'))

    x_qe_t_k_info_ik = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary storage for k-point info, k-index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_info_ik'))

    x_qe_t_k_info_vec_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point info, x-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_info_vec_x'))

    x_qe_t_k_info_vec_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point info, y-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_info_vec_y'))

    x_qe_t_k_info_vec_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point info, z-component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_info_vec_z'))

    x_qe_t_k_info_wk = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storage for k-point info, weight
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_k_info_wk'))

    x_qe_k_info_ik = Quantity(
        type=np.dtype(np.int32),
        shape=['x_qe_nk'],
        description='''
        K-point info, k-index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_k_info_ik'))

    x_qe_k_info_vec = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_nk', 3],
        description='''
        K-point info, cartesian coordinate
        ''',
        a_legacy=LegacyDefinition(name='x_qe_k_info_vec'))

    x_qe_k_info_wk = Quantity(
        type=np.dtype(np.float64),
        shape=['x_qe_nk'],
        description='''
        K-point info, weight
        ''',
        a_legacy=LegacyDefinition(name='x_qe_k_info_wk'))

    x_qe_dense_g_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Dense-grid info, G cutoff
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dense_g_cutoff'))

    x_qe_dense_g_vectors = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Dense-grid info, number of G vectors
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dense_g_vectors'))

    x_qe_t_dense_FFT_grid_x = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Dense-grid info, FFT grid x
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dense_FFT_grid_x'))

    x_qe_t_dense_FFT_grid_y = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Dense-grid info, FFT grid y
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dense_FFT_grid_y'))

    x_qe_t_dense_FFT_grid_z = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Dense-grid info, FFT grid z
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_dense_FFT_grid_z'))

    x_qe_dense_FFT_grid = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        Dense-grid info, FFT grid
        ''',
        a_legacy=LegacyDefinition(name='x_qe_dense_FFT_grid'))

    x_qe_smooth_g_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smooth-grid info, G cutoff
        ''',
        a_legacy=LegacyDefinition(name='x_qe_smooth_g_cutoff'))

    x_qe_smooth_g_vectors = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Smooth-grid info, number of G vectors
        ''',
        a_legacy=LegacyDefinition(name='x_qe_smooth_g_vectors'))

    x_qe_t_smooth_FFT_grid_x = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Smooth-grid info, FFT grid x
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_smooth_FFT_grid_x'))

    x_qe_t_smooth_FFT_grid_y = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Smooth-grid info, FFT grid y
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_smooth_FFT_grid_y'))

    x_qe_t_smooth_FFT_grid_z = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: Smooth-grid info, FFT grid z
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_smooth_FFT_grid_z'))

    x_qe_smooth_FFT_grid = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        Smooth-grid info, FFT grid
        ''',
        a_legacy=LegacyDefinition(name='x_qe_smooth_FFT_grid'))


class section_XC_functionals(public.section_XC_functionals):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_XC_functionals'))

    x_qe_xc_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of XC functional component in Quantum Espresso context
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_name'))

    x_qe_xc_comment = Quantity(
        type=str,
        shape=[],
        description='''
        Quantum Espresso comment about meaning of XC functional component
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_comment'))

    x_qe_xc_index_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of Index within Quantum Espresso where XC functional component was set from
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_index_name'))

    x_qe_xc_index = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Index value within Quantum Espresso where XC functional component was set from
        ''',
        a_legacy=LegacyDefinition(name='x_qe_xc_index'))


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    x_qe_t_iter_mpersite_idx = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Temporary: iteration per-site magnetization data, atom index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_iter_mpersite_idx'))

    x_qe_t_iter_mpersite_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: iteration per-site magnetization data, atom charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_iter_mpersite_charge'))

    x_qe_t_iter_mpersite_magn = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: iteration per-site magnetization data, atom magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_iter_mpersite_magn'))

    x_qe_t_iter_mpersite_constr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary: iteration per-site magnetization data, constraints
        ''',
        a_legacy=LegacyDefinition(name='x_qe_t_iter_mpersite_constr'))

    x_qe_iter_mpersite_idx = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms'],
        description='''
        iteration per-site magnetization data, atom index
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iter_mpersite_idx'))

    x_qe_iter_mpersite_charge = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        iteration per-site magnetization data, atom charge
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iter_mpersite_charge'))

    x_qe_iter_mpersite_magn = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        iteration per-site magnetization data, atom magnetization
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iter_mpersite_magn'))

    x_qe_iter_mpersite_constr = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        iteration per-site magnetization data, constraints
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iter_mpersite_constr'))

    x_qe_iteration_efield_eeigx_re = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        E-field: expectation value of exp(iGx), real part, in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_efield_eeigx_re'))

    x_qe_iteration_efield_eeigx_im = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        E-field: expectation value of exp(iGx), imaginary part, in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_efield_eeigx_im'))

    x_qe_iteration_efield_dipole_electronic = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        E-field: Electronic dipole, in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_efield_dipole_electronic'))

    x_qe_iteration_efield_dipole_ionic = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        E-field: Ionic dipole, in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_efield_dipole_ionic'))

    x_qe_iteration_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Iteration number
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_number'))

    x_qe_iteration_ecutwfc = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        PW cutoff used during iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_ecutwfc'))

    x_qe_iteration_beta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Mixing parameter Beta during iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_beta'))

    x_qe_iteration_charge_negative_up = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Charge in iteration (up)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_charge_negative_up'))

    x_qe_iteration_charge_negative_down = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Charge in iteration (down)
        ''',
        a_legacy=LegacyDefinition(name='x_qe_iteration_charge_negative_down'))

    x_qe_energy_total_harris_foulkes_estimate_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Harris-Foulkes estimate of total energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_total_harris_foulkes_estimate_iteration'))

    x_qe_energy_total_accuracy_estimate_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Accuracy estimate of total energy
        ''',
        a_legacy=LegacyDefinition(name='x_qe_energy_total_accuracy_estimate_iteration'))

    x_qe_magnetization_total_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Total per-cell magnetization in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_magnetization_total_iteration'))

    x_qe_magnetization_absolute_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Absolute per-cell magnetization in iteration
        ''',
        a_legacy=LegacyDefinition(name='x_qe_magnetization_absolute_iteration'))

    x_qe_section_scf_diagonalization = SubSection(
        sub_section=SectionProxy('x_qe_section_scf_diagonalization'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_qe_section_scf_diagonalization'))


class section_eigenvalues(public.section_eigenvalues):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_eigenvalues'))

    x_qe_eigenvalues_number_of_planewaves = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_eigenvalues_kpoints'],
        description='''
        Number of plane waves for each k-point
        ''',
        a_legacy=LegacyDefinition(name='x_qe_eigenvalues_number_of_planewaves'))


m_package.__init_metainfo__()
