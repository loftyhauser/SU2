#!/usr/bin/env python

## \file parallel_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    test_list = []

    ##########################
    ### Compressible Euler ###
    ##########################

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 20
    channel.test_vals = [-2.647975, 2.818090, 0.022280, 0.004644]
    channel.su2_exec  = "parallel_computation.py -f"
    channel.timeout   = 1600
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.014140, -3.537888, 0.333403, 0.021227]
    naca0012.su2_exec  = "parallel_computation.py -f"
    naca0012.timeout   = 1600
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.944740, 4.782451, -0.208522, 0.036742]
    wedge.su2_exec  = "parallel_computation.py -f"
    wedge.timeout   = 1600
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-9.277150, -8.694005, 0.281703, 0.011821]
    oneram6.su2_exec  = "parallel_computation.py -f"
    oneram6.timeout   = 3200
    oneram6.tol       = 0.00001
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-7.379831, -1.886302, 0.300000, 0.019471]
    fixedCL_naca0012.su2_exec  = "parallel_computation.py -f"
    fixedCL_naca0012.timeout   = 1600
    fixedCL_naca0012.tol       = 0.00001
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.217981, 4.256386, 0.009084, 0.016823]
    polar_naca0012.su2_exec  = "compute_polar.py -i 11"
    polar_naca0012.timeout   = 1600
    polar_naca0012.tol       = 0.00001
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.540009, 6.916653, 0.000000, 1.868976] #last 4 columns
    bluntbody.su2_exec  = "parallel_computation.py -f"
    bluntbody.timeout   = 1600
    bluntbody.tol       = 0.00001
    test_list.append(bluntbody)

    # Equivalent area NACA64-206
    ea_naca64206           = TestCase('ea_naca64206')
    ea_naca64206.cfg_dir   = "optimization_euler/equivalentarea_naca64206"
    ea_naca64206.cfg_file  = "NACA64206.cfg"
    ea_naca64206.test_iter = 10
    ea_naca64206.test_vals = [-1.076215, -0.391987, -0.000701, 67775.0]
    ea_naca64206.su2_exec  = "mpirun -n 2 SU2_CFD"
    ea_naca64206.timeout   = 1600
    ea_naca64206.tol       = 0.00001
    test_list.append(ea_naca64206)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-9.336395, -3.849426, 0.001112, 0.036276, 2.361500, -2.325300, -2.279700, -2.279700]
    flatplate.su2_exec  = "mpirun -n 2 SU2_CFD"
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Custom objective function
    flatplate_udobj           = TestCase('flatplate_udobj')
    flatplate_udobj.cfg_dir   = "user_defined_functions"
    flatplate_udobj.cfg_file  = "lam_flatplate.cfg"
    flatplate_udobj.test_iter = 20
    flatplate_udobj.test_vals = [-6.653802, -1.18143, -0.794887, 0.000611, -3.6850e-04, 7.3568e-04, -1.1042e-03, 5.9669e+02, 2.9980e+02, 2.9689e+02, 1.7147]
    flatplate_udobj.su2_exec  = "mpirun -n 2 SU2_CFD"
    flatplate_udobj.timeout   = 1600
    flatplate_udobj.tol       = 0.00001
    test_list.append(flatplate_udobj)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.754517, -1.286785, -0.213640, 0.706519]
    cylinder.su2_exec  = "parallel_computation.py -f"
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.858484, -1.396528, -1.854558, 110.033249]
    cylinder_lowmach.su2_exec  = "parallel_computation.py -f"
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050847, 0.648238, 0.000200, 13.639839]
    poiseuille.su2_exec  = "parallel_computation.py -f"
    poiseuille.timeout   = 1600
    poiseuille.tol       = 0.001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.492859, -7.672756, -0.000000, 2.085796]
    poiseuille_profile.su2_exec  = "parallel_computation.py -f"
    poiseuille_profile.timeout   = 1600
    poiseuille_profile.tol       = 0.00001
    test_list.append(poiseuille_profile)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.004689, -5.265793, 0.809463, 0.062016]
    rae2822_sa.su2_exec  = "parallel_computation.py -f"
    rae2822_sa.timeout   = 1600
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510641, 4.870022, 0.813722, 0.062439]
    rae2822_sst.su2_exec  = "parallel_computation.py -f"
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.435890, 4.870022, 0.813722, 0.062439]
    rae2822_sst_sust.su2_exec  = "parallel_computation.py -f"
    rae2822_sst_sust.timeout   = 1600
    rae2822_sst_sust.tol       = 0.00001
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.147548, -6.729213, -0.176227, 0.057731]
    turb_flatplate.su2_exec  = "parallel_computation.py -f"
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.388839, -6.689413, 0.230321, 0.157640] #last 4 columns
    turb_oneram6.su2_exec  = "parallel_computation.py -f"
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # ONERA M6 Wing - Newton-Krylov
    turb_oneram6_nk           = TestCase('turb_oneram6_nk')
    turb_oneram6_nk.cfg_dir   = "rans/oneram6"
    turb_oneram6_nk.cfg_file  = "turb_ONERAM6_nk.cfg"
    turb_oneram6_nk.test_iter = 20
    turb_oneram6_nk.test_vals = [-4.892257, -4.514011, -11.432312, 0.221025, 0.045570, 2, -0.899459, 3.1384e+01]
    turb_oneram6_nk.su2_exec  = "mpirun -n 2 SU2_CFD"
    turb_oneram6_nk.timeout   = 600
    turb_oneram6_nk.tol       = 0.0001
    test_list.append(turb_oneram6_nk)

    # NACA0012 (SA, FUN3D finest grid results: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-8.621456, -10.378269, 1.064502, 0.019710, 20.000000, -1.811700, 20.000000, -5.171326]
    turb_naca0012_sa.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-11.456387, -12.800055, -5.865784, 1.049989, 0.019163, -1.838252]
    turb_naca0012_sst.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D finest grid results: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-11.370785, -12.641676, -5.748419, 1.005233, 0.019017, -2.057144]
    turb_naca0012_sst_sust.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_sust.timeout   = 3200
    turb_naca0012_sst_sust.tol       = 0.00001
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [ -5.216685, -9.562448, -1.565778, 1.022393, 0.040542, -3.729648]
    turb_naca0012_sst_fixedvalues.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    turb_naca0012_sst_fixedvalues.tol       = 0.00001
    test_list.append(turb_naca0012_sst_fixedvalues)

    # NACA0012 (SST, explicit Euler for flow and turbulence equations)
    turb_naca0012_sst_expliciteuler           = TestCase('turb_naca0012_sst_expliciteuler')
    turb_naca0012_sst_expliciteuler.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_expliciteuler.cfg_file  = "turb_NACA0012_sst_expliciteuler.cfg"
    turb_naca0012_sst_expliciteuler.test_iter = 10
    turb_naca0012_sst_expliciteuler.test_vals = [-3.532228, -3.157766, 3.364025, 1.124824, 0.501717, -float("inf")]
    turb_naca0012_sst_expliciteuler.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_expliciteuler.timeout   = 3200
    turb_naca0012_sst_expliciteuler.tol       = 0.00001
    test_list.append(turb_naca0012_sst_expliciteuler)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389575, -8.409529, 0.000048, 0.056329] #last 4 columns
    propeller.su2_exec  = "parallel_computation.py -f"
    propeller.timeout   = 3200
    propeller.tol       = 0.00001
    test_list.append(propeller)

    #######################################
    ### Axisymmetric Compressible RANS  ###
    #######################################
    
    # Axisymmetric air nozzle (transonic)
    axi_rans_air_nozzle           = TestCase('axi_rans_air_nozzle')
    axi_rans_air_nozzle.cfg_dir   = "axisymmetric_rans/air_nozzle"
    axi_rans_air_nozzle.cfg_file  = "air_nozzle.cfg"
    axi_rans_air_nozzle.test_iter = 10
    axi_rans_air_nozzle.test_vals = [-12.096377, -6.636625, -8.786639, -2.399099]
    axi_rans_air_nozzle.su2_exec  = "mpirun -n 2 SU2_CFD"
    axi_rans_air_nozzle.timeout   = 1600
    axi_rans_air_nozzle.tol       = 0.0001
    test_list.append(axi_rans_air_nozzle)

    #################################
    ## Compressible RANS Restart  ###
    #################################

    # NACA0012 SST Multigrid restart
    turb_naca0012_sst_restart_mg           = TestCase('turb_naca0012_sst_restart_mg')
    turb_naca0012_sst_restart_mg.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_restart_mg.cfg_file  = "turb_NACA0012_sst_multigrid_restart.cfg"
    turb_naca0012_sst_restart_mg.test_iter = 20
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-7.619889, -7.729499, -1.981039, -0.000016, 0.079062]
    turb_naca0012_sst_restart_mg.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)

    ###########################
    ### Turbulence modeling ###
    ###########################

    # SA Baseline (Identical to RANS SA RAE2822)
    turbmod_sa_bsl_rae2822           = TestCase('turbmod_sa_bsl_rae2822')
    turbmod_sa_bsl_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_bsl_rae2822.cfg_file  = "turb_SA_BSL_RAE2822.cfg"
    turbmod_sa_bsl_rae2822.test_iter = 20
    turbmod_sa_bsl_rae2822.test_vals = [-2.004689, 0.742306, 0.497308, -5.265793, 0.809463, 0.062016]
    turbmod_sa_bsl_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_bsl_rae2822.timeout   = 1600
    turbmod_sa_bsl_rae2822.new_output = True
    turbmod_sa_bsl_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_bsl_rae2822)

    # SA Negative
    turbmod_sa_neg_rae2822           = TestCase('turbmod_sa_neg_rae2822')
    turbmod_sa_neg_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_neg_rae2822.cfg_file  = "turb_SA_NEG_RAE2822.cfg"
    turbmod_sa_neg_rae2822.test_iter = 10
    turbmod_sa_neg_rae2822.test_vals = [-1.094542, 3.161741, 2.333560, 2.864805, 1.494301, 0.530135]
    turbmod_sa_neg_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_neg_rae2822.timeout   = 1600
    turbmod_sa_neg_rae2822.new_output = True
    turbmod_sa_neg_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_neg_rae2822)

    # SA Compressibility Correction
    turbmod_sa_comp_rae2822           = TestCase('turbmod_sa_comp_rae2822')
    turbmod_sa_comp_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_rae2822.cfg_file  = "turb_SA_COMP_RAE2822.cfg"
    turbmod_sa_comp_rae2822.test_iter = 20
    turbmod_sa_comp_rae2822.test_vals = [-2.004687, 0.742304, 0.497309, -5.266081, 0.809467, 0.062029]
    turbmod_sa_comp_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_comp_rae2822.timeout   = 1600
    turbmod_sa_comp_rae2822.new_output = True
    turbmod_sa_comp_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_comp_rae2822)

    # SA Edwards
    turbmod_sa_edw_rae2822           = TestCase('turbmod_sa_edw_rae2822')
    turbmod_sa_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_edw_rae2822.cfg_file  = "turb_SA_EDW_RAE2822.cfg"
    turbmod_sa_edw_rae2822.test_iter = 20
    turbmod_sa_edw_rae2822.test_vals = [-2.004687, 0.742306, 0.497310, -5.290769, 0.809485, 0.062036]
    turbmod_sa_edw_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_edw_rae2822.timeout   = 1600
    turbmod_sa_edw_rae2822.new_output = True
    turbmod_sa_edw_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_edw_rae2822)

    # SA Compressibility and Edwards
    turbmod_sa_comp_edw_rae2822           = TestCase('turbmod_sa_comp_edw_rae2822')
    turbmod_sa_comp_edw_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_comp_edw_rae2822.cfg_file  = "turb_SA_COMP_EDW_RAE2822.cfg"
    turbmod_sa_comp_edw_rae2822.test_iter = 20
    turbmod_sa_comp_edw_rae2822.test_vals = [-2.004685, 0.742307, 0.497311, -5.290750, 0.809487, 0.062045]
    turbmod_sa_comp_edw_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_comp_edw_rae2822.timeout   = 1600
    turbmod_sa_comp_edw_rae2822.new_output = True
    turbmod_sa_comp_edw_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_comp_edw_rae2822)

    # SA QCR
    turbmod_sa_qcr_rae2822           = TestCase('turbmod_sa_qcr_rae2822')
    turbmod_sa_qcr_rae2822.cfg_dir   = "turbulence_models/sa/rae2822"
    turbmod_sa_qcr_rae2822.cfg_file  = "turb_SA_QCR_RAE2822.cfg"
    turbmod_sa_qcr_rae2822.test_iter = 20
    turbmod_sa_qcr_rae2822.test_vals = [-2.004793, 0.742353, 0.497315, -5.265974, 0.807841, 0.062027]
    turbmod_sa_qcr_rae2822.su2_exec  = "mpirun -n 2 SU2_CFD"
    turbmod_sa_qcr_rae2822.timeout   = 1600
    turbmod_sa_qcr_rae2822.new_output = True
    turbmod_sa_qcr_rae2822.tol       = 0.00001
    test_list.append(turbmod_sa_qcr_rae2822)

    ############################
    ###      Transition      ###
    ############################

    # Schubauer-Klebanoff Natural Transition Case
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.test_vals    = [-7.994740, -13.240225, 0.000046, 0.007987]
    schubauer_klebanoff_transition.su2_exec     = "parallel_computation.py -f"
    schubauer_klebanoff_transition.timeout      = 1600
    schubauer_klebanoff_transition.tol          = 0.00001
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.300815, -14.587362, 0.300920, 0.019552]
    contadj_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.00001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.130993, -12.702085, 0.685900, 0.007594]
    contadj_oneram6.su2_exec  = "parallel_computation.py -f"
    contadj_oneram6.timeout   = 1600
    contadj_oneram6.tol       = 0.00001
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872691, -2.755572, 853010.000000, 0.000000] #last 4 columns
    contadj_wedge.su2_exec  = "parallel_computation.py -f"
    contadj_wedge.timeout   = 1600
    contadj_wedge.tol       = 0.00001
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixed_CL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixed_CL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixed_CL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixed_CL_naca0012.test_iter = 100
    contadj_fixed_CL_naca0012.test_vals = [0.275856, -5.200511, 0.342710, 0.000105]
    contadj_fixed_CL_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_fixed_CL_naca0012.timeout   = 1600
    contadj_fixed_CL_naca0012.tol       = 0.00001
    test_list.append(contadj_fixed_CL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [-3.651435, -9.113086, 2.056700, -0.000000]
    contadj_ns_cylinder.su2_exec  = "parallel_computation.py -f"
    contadj_ns_cylinder.timeout   = 1600
    contadj_ns_cylinder.tol       = 0.00001
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210] #last 4 columns
    contadj_ns_naca0012_sub.su2_exec  = "parallel_computation.py -f"
    contadj_ns_naca0012_sub.timeout   = 1600
    contadj_ns_naca0012_sub.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_sub)

    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [ -1.039664, -6.575019, 1.772300, 0.012495] #last 4 columns
    contadj_ns_naca0012_trans.su2_exec  = "parallel_computation.py -f"
    contadj_ns_naca0012_trans.timeout   = 1600
    contadj_ns_naca0012_trans.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_trans)

    #######################################################
    ### Cont. adj. compressible RANS (frozen viscosity) ###
    #######################################################

    # Adjoint turbulent NACA0012
    contadj_rans_naca0012           = TestCase('contadj_rans_naca0012')
    contadj_rans_naca0012.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
    contadj_rans_naca0012.test_iter = 20
    contadj_rans_naca0012.test_vals = [-0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)

    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 18
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012_bin.su2_exec  = "parallel_computation.py -f"
    contadj_rans_naca0012_bin.timeout   = 1600
    contadj_rans_naca0012_bin.tol       = 0.00001
    test_list.append(contadj_rans_naca0012_bin)

    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.372407, -10.874841, -0.212470, 0.005448]
    contadj_rans_rae2822.su2_exec  = "parallel_computation.py -f"
    contadj_rans_rae2822.timeout   = 1600
    contadj_rans_rae2822.tol       = 0.00001
    test_list.append(contadj_rans_rae2822)

    #############################
    ### Compressibele RANS UQ ###
    #############################

    # NACA0012 1c
    turb_naca0012_1c           = TestCase('turb_naca0012_1c')
    turb_naca0012_1c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_1c.cfg_file  = "turb_NACA0012_uq_1c.cfg"
    turb_naca0012_1c.test_iter = 10
    turb_naca0012_1c.test_vals = [-4.976520, 1.141381, 0.460006, -0.078852]
    turb_naca0012_1c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_1c.timeout   = 1600
    turb_naca0012_1c.tol       = 0.00001
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.483311, 0.968826, 0.303278, -0.113945]
    turb_naca0012_2c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_2c.timeout   = 1600
    turb_naca0012_2c.tol       = 0.00001
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.584310, 0.931348, 0.279056, -0.113209]
    turb_naca0012_3c.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_3c.timeout   = 1600
    turb_naca0012_3c.tol       = 0.00001
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.128788, 1.077306, 0.584228, -0.048152]
    turb_naca0012_p1c1.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_p1c1.timeout   = 1600
    turb_naca0012_p1c1.tol       = 0.00001
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.554534, 0.943774, 0.400400, -0.095514]
    turb_naca0012_p1c2.su2_exec  = "parallel_computation.py -f"
    turb_naca0012_p1c2.timeout   = 1600
    turb_naca0012_p1c2.tol       = 0.00001
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.162664, 0.066378, 1.399789, 2.220404]
    square_cylinder.su2_exec  = "parallel_computation.py -f"
    square_cylinder.timeout   = 1600
    square_cylinder.tol       = 0.00001
    square_cylinder.unsteady  = True
    test_list.append(square_cylinder)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714758, -5.883004, -0.215005, 0.023783] #last 4 columns
    ddes_flatplate.su2_exec  = "parallel_computation.py -f"
    ddes_flatplate.timeout   = 1600
    ddes_flatplate.tol       = 0.00001
    ddes_flatplate.unsteady  = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.004011, -6.876230, 1.487888, 0.421869]
    unst_inc_turb_naca0015_sa.su2_exec  = "parallel_computation.py -f"
    unst_inc_turb_naca0015_sa.timeout   = 1600
    unst_inc_turb_naca0015_sa.tol       = 0.00001
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # Flat plate
    flatplate_unsteady           = TestCase('flatplate_unsteady')
    flatplate_unsteady.cfg_dir   = "navierstokes/flatplate"
    flatplate_unsteady.cfg_file  = "lam_flatplate_unst.cfg"
    flatplate_unsteady.test_iter = 3
    flatplate_unsteady.test_vals = [7.9509e-06, -8.868859, -8.231652, -6.283262, -5.466675, -3.391163, 0.002078, -0.343642]
    flatplate_unsteady.su2_exec  = "mpirun -n 2 SU2_CFD"
    flatplate_unsteady.timeout   = 1600
    flatplate_unsteady.tol       = 0.00001
    flatplate_unsteady.unsteady  = True
    test_list.append(flatplate_unsteady)

    ######################################
    ### NICFD                          ###
    ######################################

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 100
    edge_VW.test_vals = [-5.048044, 1.115667, -0.000009, 0.000000]
    edge_VW.su2_exec  = "parallel_computation.py -f"
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 100
    edge_PPR.test_vals = [-5.400790, 0.739723, -0.000035, 0.000000]
    edge_PPR.su2_exec  = "parallel_computation.py -f"
    edge_PPR.timeout   = 1600
    edge_PPR.tol       = 0.00001
    test_list.append(edge_PPR)

    ##########################
    ###   Python wrapper   ###
    ##########################

    # NACA0012
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 100
    pywrapper_naca0012.test_vals = [-6.747210, -6.149915, 0.333445, 0.021241] #last 4 columns
    pywrapper_naca0012.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    test_list.append(pywrapper_naca0012)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-11.456387, -12.800055, -5.865784, 1.049989, 0.019163, -1.838252]
    pywrapper_turb_naca0012_sst.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    test_list.append(pywrapper_turb_naca0012_sst)

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-1.162664, 0.066378, 1.399789, 2.220404] #last 4 columns
    pywrapper_square_cylinder.su2_exec  = "mpirun -np 2 SU2_CFD.py --parallel -f"
    pywrapper_square_cylinder.timeout   = 1600
    pywrapper_square_cylinder.tol       = 0.00001
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.851428, 2.192348, 0.000000, 0.000000] #last 4 columns
    mms_fvm_ns.su2_exec  = "mpirun -n 2 SU2_CFD"
    mms_fvm_ns.timeout   = 1600
    mms_fvm_ns.tol       = 0.0001
    test_list.append(mms_fvm_ns)

    # DG, compressible, euler
    ringleb_dg_euler           = TestCase('ringleb_dg_euler')
    ringleb_dg_euler.cfg_dir   = "mms/dg_ringleb"
    ringleb_dg_euler.cfg_file  = "ringleb_dg.cfg"
    ringleb_dg_euler.test_iter = 100
    ringleb_dg_euler.test_vals = [-5.136652, -4.724941, 0.000000, 0.000000] #last 4 columns
    ringleb_dg_euler.su2_exec  = "SU2_CFD"
    ringleb_dg_euler.timeout   = 1600
    ringleb_dg_euler.tol       = 0.0001
    test_list.append(ringleb_dg_euler)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    pass_list = [ test.run_test() for test in test_list ]

    ######################################
    ### RUN SU2_SOL TESTS              ###
    ######################################

    # parallel STL output using
    stl_writer_test                = TestCase('stl_writer_test')
    stl_writer_test.cfg_dir        = "rans/oneram6"
    stl_writer_test.cfg_file       = "turb_ONERAM6.cfg"
    stl_writer_test.test_iter      = 1
    stl_writer_test.su2_exec       = "mpirun -n 2 SU2_SOL"
    stl_writer_test.timeout        = 1600
    stl_writer_test.reference_file = "surface_flow.stl.ref"
    stl_writer_test.test_file      = "surface_flow.stl"
    pass_list.append(stl_writer_test.run_filediff())
    test_list.append(stl_writer_test)

    pass_list.append(sphere_ffd_def_bspline.run_def())
    test_list.append(sphere_ffd_def_bspline)

    # Tests summary
    print('==================================================================')
    print('Summary of the parallel tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print('  passed - %s'%test.tag)
        else:
            print('* FAILED - %s'%test.tag)

    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()
