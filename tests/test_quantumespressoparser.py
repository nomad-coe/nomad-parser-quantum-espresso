#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from quantumespressoparser import QuantumEspressoParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return QuantumEspressoParser()


def test_scf(parser):
    archive = EntryArchive()
    parser.parse('tests/data/HO_scf/benchmark2.out', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == '5.2.1 (svn rev. 11920)'
    assert sec_run.x_qe_input_filename == 'uspp1.in'
    assert sec_run.time_run_date_start.magnitude == 1451140876.0
    assert sec_run.x_qe_section_compile_options[0].x_qe_lmaxx == 3
    assert sec_run.x_qe_section_parallel[0].x_qe_nproc == 4
    assert sec_run.section_basis_set_cell_dependent[0].basis_set_cell_dependent_name == 'PW_25.0'
    assert sec_run.section_basis_set_cell_dependent[1].basis_set_planewave_cutoff.magnitude == approx(2.17987236e-16,)
    assert sec_run.section_sampling_method[0].sampling_method == 'geometry_optimization'
    assert 'rdiaghg' in sec_run.x_qe_profile_function
    assert sec_run.time_run_date_end.magnitude == 1451140881.0
    assert sec_run.run_clean_end

    sec_method = sec_run.section_method[0]
    assert sec_method.x_qe_sticks_sum_G_smooth == 135043
    assert 'NL pseudopotentials' in sec_method.x_qe_allocated_array_name
    assert sec_method.x_qe_allocated_array_size[2] == 33554432.
    assert sec_method.x_qe_temporary_array_dimensions[3] == '262144,    8'
    assert sec_method.x_qe_per_process_mem == approx(2.84373811e+08)
    assert sec_method.x_qe_potential_mixing_scheme == 'plain'
    assert sec_method.x_qe_starting_charge == 7.99998
    assert len(sec_method.section_XC_functionals) == 2
    assert sec_method.x_qe_xc_igcc_name == 'pbc'
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_X_PBE'
    sec_atoms = sec_method.section_method_atom_kind
    assert len(sec_atoms) == 2
    assert sec_atoms[1].method_atom_kind_label == 'H'
    assert sec_atoms[0].x_qe_pp_md5sum == '7e325307d184e51bd80757047dcf04f9'
    assert sec_atoms[1].x_qe_pp_ncoefficients == 8
    assert sec_atoms[0].x_qe_kind_mass == 16.0

    sec_system = sec_run.section_system[0]
    assert sec_system.atom_labels == ['O', 'H', 'H']
    assert sec_system.atom_positions[2][0].magnitude == approx(5.12015994e-10)
    assert False not in sec_system.configuration_periodic_dimensions
    assert sec_system.x_qe_reciprocal_cell[2][2].magnitude == approx(5.93674971e+09)
    assert len(sec_system.x_qe_k_info_vec) == 1
    assert sec_system.x_qe_cell_volume == approx(1.18547769e-27)
    assert sec_system.x_qe_nsymm == 4
    assert sec_system.x_qe_dense_FFT_grid[1] == 64
    assert sec_system.number_of_electrons[0] == 8

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == approx(-7.49748432e-17)
    assert 'ewald contribution' in sec_scc.x_qe_energy_decomposition_name
    assert sec_scc.x_qe_energy_decomposition_value[1] == approx(7.42289975e-17)
    assert sec_scc.atom_forces_raw[1][1].magnitude == approx(-3.57815176e-10)
    assert sec_scc.stress_tensor[2][2].magnitude == approx(-1.68e+08)
    assert np.shape(sec_scc.section_eigenvalues[0].eigenvalues_kpoints) == (1, 3)
    assert np.shape(sec_scc.section_eigenvalues[0].eigenvalues_values) == (1, 1, 4)
    assert sec_scc.section_eigenvalues[0].eigenvalues_values[0][0][2].magnitude == approx(-1.42427094e-18)
    assert sec_scc.energy_reference_highest_occupied.magnitude == approx(-1.15444837e-18)
    assert sec_scc.x_qe_output_datafile == 'pwscf.save'
    sec_scfs = sec_scc.section_scf_iteration
    assert len(sec_scfs) == 8
    assert sec_scfs[4].energy_total_scf_iteration.magnitude == approx(-7.49748038e-17)
    assert sec_scfs[1].x_qe_energy_total_accuracy_estimate_iteration == approx(1.15623477e-18)
    assert sec_scfs[6].x_qe_iteration_ecutwfc == approx(5.4496809027589626e-17)
    assert sec_scfs[0].time_scf_iteration_cpu1_end.magnitude == 1.2
    assert sec_scfs[3].x_qe_iteration_charge_negative_up == 0.06614


def test_multirun(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Mn_multirun/9064627814752884918776106151027.log', archive, None)

    sec_runs = archive.section_run
    assert len(sec_runs) == 3
    assert sec_runs[0].section_method[0].smearing_width == approx(2.3978595972139434e-20)
    assert len(sec_runs[1].section_single_configuration_calculation[0].section_scf_iteration) == 111
    assert sec_runs[2].section_single_configuration_calculation[0].section_scf_iteration[45].x_qe_iter_mpersite_magn[6] == -0.3325
    assert sec_runs[0].section_system[0].x_qe_atom_starting_magnetization[1] == 0.133
    assert np.shape(sec_runs[0].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_values) == (2, 20, 100)
    assert np.shape(sec_runs[1].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_values) == (2, 20, 100)
    assert np.shape(sec_runs[2].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_values) == (2, 20, 100)
    assert len(sec_runs[0].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_kpoints) == 20
    assert len(sec_runs[1].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_kpoints) == 20
    assert len(sec_runs[2].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_kpoints) == 20
    assert sec_runs[0].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_kpoints[10][1] == approx(-0.1667096)
    assert sec_runs[1].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_values[0][3][-5].magnitude == approx(1.42385437e-19)
    assert sec_runs[2].section_single_configuration_calculation[0].section_eigenvalues[0].eigenvalues_values.magnitude[1][-10][35] == approx(-7.25180392e-18)


def test_md(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si_md/out.out', archive, None)

    assert archive.section_run[0].section_sampling_method[0].sampling_method == 'molecular_dynamics'
    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 50
    assert archive.section_run[0].section_system[6].atom_positions[1][2].magnitude == approx(6.66987013e-11)
    assert sec_sccs[-3].atom_forces_raw[1][1].magnitude == approx(9.55685747e-10)
    assert len(sec_sccs[22].section_scf_iteration) == 3


def test_dos(parser):
    archive = EntryArchive()
    parser.parse('tests/data/W_dos/w.dos.out', archive, None)

    sec_dos = archive.section_run[0].section_single_configuration_calculation[0].section_dos[0]
    assert np.shape(sec_dos.dos_values) == (1, 1801)
    assert len(sec_dos.dos_energies) == 1801
    assert sec_dos.dos_energies[269].magnitude == approx(1.23207383e-18)
    assert sec_dos.dos_values[0][150] == approx(2.8991809650870246e+17)
    assert sec_dos.dos_integrated_values[0][1316] == 8.582
