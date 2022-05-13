#!/usr/bin/env python

## \file serial_regression.py
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

from __future__ import print_function, division, absolute_import
import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values.
       This will be used to do checks when code is pushed to github
       to make sure nothing is broken. '''

    test_list = []

    #########################
    ## Compressible Euler ###
    #########################

    # Dry run test Euler
    channel_d           = TestCase('dry run Euler')
    channel_d.cfg_dir   = "euler/channel"
    channel_d.cfg_file  = "inv_channel_RK.cfg"
    channel_d.su2_exec  = "SU2_CFD -d"
    channel_d.timeout   = 1600
    test_list.append(channel_d)

    # Channel
    channel           = TestCase('channel')
    channel.cfg_dir   = "euler/channel"
    channel.cfg_file  = "inv_channel_RK.cfg"
    channel.test_iter = 10
    channel.test_vals = [-2.475872, 3.046370, -0.203974, 0.036018]
    channel.su2_exec  = "SU2_CFD"
    channel.timeout   = 1600
    channel.new_output = True
    channel.tol       = 0.00001
    test_list.append(channel)

    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.023999, -3.515034, 0.339426, 0.022217] #last 4 columns
    naca0012.su2_exec  = "SU2_CFD"
    naca0012.timeout   = 1600
    naca0012.new_output= True
    naca0012.tol       = 0.00001
    test_list.append(naca0012)

    # Supersonic wedge
    wedge           = TestCase('wedge')
    wedge.cfg_dir   = "euler/wedge"
    wedge.cfg_file  = "inv_wedge_HLLC.cfg"
    wedge.test_iter = 20
    wedge.test_vals = [-0.942862, 4.784581, -0.208106, 0.036665] #last 4 columns
    wedge.su2_exec  = "SU2_CFD"
    wedge.timeout   = 1600
    wedge.new_output= True
    wedge.tol       = 0.00001
    test_list.append(wedge)

    # ONERA M6 Wing
    oneram6           = TestCase('oneram6')
    oneram6.cfg_dir   = "euler/oneram6"
    oneram6.cfg_file  = "inv_ONERAM6.cfg"
    oneram6.test_iter = 10
    oneram6.test_vals = [-9.279396, -8.697739, 0.281703, 0.011821]
    oneram6.su2_exec  = "SU2_CFD"
    oneram6.timeout   = 9600
    oneram6.new_output = True
    oneram6.tol       = 0.00001
    test_list.append(oneram6)

    # Fixed CL NACA0012
    fixedCL_naca0012           = TestCase('fixedcl_naca0012')
    fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    fixedCL_naca0012.cfg_file  = "inv_NACA0012.cfg"
    fixedCL_naca0012.test_iter = 10
    fixedCL_naca0012.test_vals = [-7.382410, -1.879887, 0.300000, 0.019471]
    fixedCL_naca0012.su2_exec  = "SU2_CFD"
    fixedCL_naca0012.new_output = True
    fixedCL_naca0012.timeout   = 1600
    fixedCL_naca0012.tol       = 0.00001
    test_list.append(fixedCL_naca0012)

    # Polar sweep of the inviscid NACA0012
    polar_naca0012           = TestCase('polar_naca0012')
    polar_naca0012.cfg_dir   = "polar/naca0012"
    polar_naca0012.cfg_file  = "inv_NACA0012.cfg"
    polar_naca0012.polar     = True
    polar_naca0012.new_output= True
    polar_naca0012.test_iter = 10
    polar_naca0012.test_vals = [-1.243326, 4.224483, 0.016432, 0.016145]
    polar_naca0012.su2_exec  = "compute_polar.py -n 1 -i 11"
    polar_naca0012.timeout   = 1600
    polar_naca0012.tol       = 0.00001
    test_list.append(polar_naca0012)

    # HYPERSONIC FLOW PAST BLUNT BODY
    bluntbody           = TestCase('bluntbody')
    bluntbody.cfg_dir   = "euler/bluntbody"
    bluntbody.cfg_file  = "blunt.cfg"
    bluntbody.new_output = True
    bluntbody.test_iter = 20
    bluntbody.test_vals = [0.540009, 6.916653, -0.000000, 1.868975] #last 4 columns
    bluntbody.su2_exec  = "SU2_CFD"
    bluntbody.timeout   = 1600
    bluntbody.tol       = 0.00001
    test_list.append(bluntbody)

    ##########################
    ###  Compressible N-S  ###
    ##########################

    # Dry run test NS
    flatplate_d           = TestCase('dry run NS')
    flatplate_d.cfg_dir   = "navierstokes/flatplate"
    flatplate_d.cfg_file  = "lam_flatplate.cfg"
    flatplate_d.su2_exec  = "SU2_CFD -d"
    flatplate_d.timeout   = 1600
    test_list.append(flatplate_d)

    # Laminar flat plate
    flatplate           = TestCase('flatplate')
    flatplate.cfg_dir   = "navierstokes/flatplate"
    flatplate.cfg_file  = "lam_flatplate.cfg"
    flatplate.test_iter = 100
    flatplate.test_vals = [-9.856258, -4.371983, 0.001112, 0.036277, 2.361500, -2.325300, -2.279500, -2.279500]
    flatplate.su2_exec  = "SU2_CFD"
    flatplate.new_output = True
    flatplate.timeout   = 1600
    flatplate.tol       = 0.00001
    test_list.append(flatplate)

    # Laminar cylinder (steady)
    cylinder           = TestCase('cylinder')
    cylinder.cfg_dir   = "navierstokes/cylinder"
    cylinder.cfg_file  = "lam_cylinder.cfg"
    cylinder.test_iter = 25
    cylinder.test_vals = [-6.765430, -1.297426, 0.019508, 0.310015] #last 4 columns
    cylinder.su2_exec  = "SU2_CFD"
    cylinder.new_output = True
    cylinder.timeout   = 1600
    cylinder.tol       = 0.00001
    test_list.append(cylinder)

    # Laminar cylinder (low Mach correction)
    cylinder_lowmach           = TestCase('cylinder_lowmach')
    cylinder_lowmach.cfg_dir   = "navierstokes/cylinder"
    cylinder_lowmach.cfg_file  = "cylinder_lowmach.cfg"
    cylinder_lowmach.test_iter = 25
    cylinder_lowmach.test_vals = [-6.850123, -1.388088, -0.056090, 108.140176]
    cylinder_lowmach.su2_exec  = "SU2_CFD"
    cylinder_lowmach.new_output = True
    cylinder_lowmach.timeout   = 1600
    cylinder_lowmach.tol       = 0.00001
    test_list.append(cylinder_lowmach)

    # 2D Poiseuille flow (body force driven with periodic inlet / outlet)
    poiseuille           = TestCase('poiseuille')
    poiseuille.cfg_dir   = "navierstokes/poiseuille"
    poiseuille.cfg_file  = "lam_poiseuille.cfg"
    poiseuille.test_iter = 10
    poiseuille.test_vals = [-5.050732, 0.648355, 0.012273, 13.643219] #last 4 columns
    poiseuille.su2_exec  = "SU2_CFD"
    poiseuille.new_output = True
    poiseuille.timeout   = 1600
    poiseuille.tol       = 0.00001
    test_list.append(poiseuille)

    # 2D Poiseuille flow (inlet profile file)
    poiseuille_profile           = TestCase('poiseuille_profile')
    poiseuille_profile.cfg_dir   = "navierstokes/poiseuille"
    poiseuille_profile.cfg_file  = "profile_poiseuille.cfg"
    poiseuille_profile.test_iter = 10
    poiseuille_profile.test_vals = [-12.494720, -7.711373, -0.000000, 2.085796] #last 4 columns
    poiseuille_profile.su2_exec  = "SU2_CFD"
    poiseuille_profile.new_output = True
    poiseuille_profile.timeout   = 1600
    poiseuille_profile.tol       = 0.00001
    test_list.append(poiseuille_profile)

    ##########################
    ### Compressible RANS  ###
    ##########################

    # Dry run RANS
    rae2822_sa_d           = TestCase('dry run RANS')
    rae2822_sa_d .cfg_dir   = "rans/rae2822"
    rae2822_sa_d .cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa_d .su2_exec  = "SU2_CFD -d"
    rae2822_sa_d .timeout   = 1600
    test_list.append(rae2822_sa_d)

    # RAE2822 SA
    rae2822_sa           = TestCase('rae2822_sa')
    rae2822_sa.cfg_dir   = "rans/rae2822"
    rae2822_sa.cfg_file  = "turb_SA_RAE2822.cfg"
    rae2822_sa.test_iter = 20
    rae2822_sa.test_vals = [-2.020123, -5.269330, 0.807147, 0.060499]
    rae2822_sa.su2_exec  = "SU2_CFD"
    rae2822_sa.timeout   = 1600
    rae2822_sa.new_output = True
    rae2822_sa.tol       = 0.00001
    test_list.append(rae2822_sa)

    # RAE2822 SST
    rae2822_sst           = TestCase('rae2822_sst')
    rae2822_sst.cfg_dir   = "rans/rae2822"
    rae2822_sst.cfg_file  = "turb_SST_RAE2822.cfg"
    rae2822_sst.test_iter = 20
    rae2822_sst.test_vals = [-0.510639, 4.872266, 0.812659, 0.061095]
    rae2822_sst.su2_exec  = "SU2_CFD"
    rae2822_sst.new_output = True
    rae2822_sst.timeout   = 1600
    rae2822_sst.tol       = 0.00001
    test_list.append(rae2822_sst)

    # RAE2822 SST_SUST
    rae2822_sst_sust           = TestCase('rae2822_sst_sust')
    rae2822_sst_sust.cfg_dir   = "rans/rae2822"
    rae2822_sst_sust.cfg_file  = "turb_SST_SUST_RAE2822.cfg"
    rae2822_sst_sust.test_iter = 20
    rae2822_sst_sust.test_vals = [-2.431741, 4.872266, 0.812658, 0.061095]
    rae2822_sst_sust.su2_exec  = "SU2_CFD"
    rae2822_sst_sust.timeout   = 1600
    rae2822_sst_sust.tol       = 0.00001
    test_list.append(rae2822_sst_sust)

    # Flat plate
    turb_flatplate           = TestCase('turb_flatplate')
    turb_flatplate.cfg_dir   = "rans/flatplate"
    turb_flatplate.cfg_file  = "turb_SA_flatplate.cfg"
    turb_flatplate.test_iter = 20
    turb_flatplate.test_vals = [-4.157169, -6.737133, -0.176253, 0.057446] #last 4 columns
    turb_flatplate.su2_exec  = "SU2_CFD"
    turb_flatplate.new_output  = True
    turb_flatplate.timeout   = 1600
    turb_flatplate.tol       = 0.00001
    test_list.append(turb_flatplate)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SST
    turb_wallfunction_flatplate_sst           = TestCase('turb_sst_wallfunction_flatplate')
    turb_wallfunction_flatplate_sst.cfg_dir   = "wallfunctions/flatplate/compressible_SST"
    turb_wallfunction_flatplate_sst.cfg_file  = "turb_SST_flatplate.cfg"
    turb_wallfunction_flatplate_sst.test_iter = 10
    turb_wallfunction_flatplate_sst.test_vals = [-4.229955, -1.930560, -1.998477, 1.250383, -1.635663, 1.462396, 10.000000, -2.151959, 0.072873, 0.002514] #last 10 columns
    turb_wallfunction_flatplate_sst.su2_exec  = "SU2_CFD"
    turb_wallfunction_flatplate_sst.timeout   = 1600
    turb_wallfunction_flatplate_sst.tol       = 0.00001
    test_list.append(turb_wallfunction_flatplate_sst)

    # FLAT PLATE, WALL FUNCTIONS, COMPRESSIBLE SA
    turb_wallfunction_flatplate_sa           = TestCase('turb_sa_wallfunction_flatplate')
    turb_wallfunction_flatplate_sa.cfg_dir   = "wallfunctions/flatplate/compressible_SA"
    turb_wallfunction_flatplate_sa.cfg_file  = "turb_SA_flatplate.cfg"
    turb_wallfunction_flatplate_sa.test_iter = 10
    turb_wallfunction_flatplate_sa.test_vals = [-4.436048, -2.044706, -2.114644, 0.979771, -5.393729, 10, -1.589465, 0.069744, 0.002686] #last 9 columns
    turb_wallfunction_flatplate_sa.su2_exec  = "SU2_CFD"
    turb_wallfunction_flatplate_sa.timeout   = 1600
    turb_wallfunction_flatplate_sa.tol       = 0.00001
    test_list.append(turb_wallfunction_flatplate_sa)

    # ONERA M6 Wing
    turb_oneram6           = TestCase('turb_oneram6')
    turb_oneram6.cfg_dir   = "rans/oneram6"
    turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
    turb_oneram6.test_iter = 10
    turb_oneram6.test_vals = [-2.388841, -6.689414, 0.230321, 0.157640]#last 4 columns
    turb_oneram6.su2_exec  = "SU2_CFD"
    turb_oneram6.new_output = True
    turb_oneram6.timeout   = 3200
    turb_oneram6.tol       = 0.00001
    test_list.append(turb_oneram6)

    # NACA0012 (SA, FUN3D results for finest grid: CL=1.0983, CD=0.01242)
    turb_naca0012_sa           = TestCase('turb_naca0012_sa')
    turb_naca0012_sa.cfg_dir   = "rans/naca0012"
    turb_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    turb_naca0012_sa.test_iter = 10
    turb_naca0012_sa.test_vals = [-8.629583, -10.377793, 1.064488, 0.019711, 20.000000, -2.173971, 20.000000, -5.213344]
    turb_naca0012_sa.su2_exec  = "SU2_CFD"
    turb_naca0012_sa.new_output = True
    turb_naca0012_sa.timeout   = 3200
    turb_naca0012_sa.tol       = 0.00001
    test_list.append(turb_naca0012_sa)

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst           = TestCase('turb_naca0012_sst')
    turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    turb_naca0012_sst.test_iter = 10
    turb_naca0012_sst.test_vals = [-11.451010, -12.798258, -5.863895, 1.049989, 0.019163, -1.925018]
    turb_naca0012_sst.su2_exec  = "SU2_CFD"
    turb_naca0012_sst.new_output  = True
    turb_naca0012_sst.timeout   = 3200
    turb_naca0012_sst.tol       = 0.00001
    test_list.append(turb_naca0012_sst)

    # NACA0012 (SST_SUST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    turb_naca0012_sst_sust           = TestCase('turb_naca0012_sst_sust')
    turb_naca0012_sst_sust.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_sust.cfg_file  = "turb_NACA0012_sst_sust.cfg"
    turb_naca0012_sst_sust.test_iter = 10
    turb_naca0012_sst_sust.test_vals = [-11.367386, -12.640857, -5.747260, 1.005233, 0.019017, -1.985871]
    turb_naca0012_sst_sust.su2_exec  = "SU2_CFD"
    turb_naca0012_sst_sust.timeout   = 3200
    turb_naca0012_sst_sust.tol       = 0.00001
    test_list.append(turb_naca0012_sst_sust)

    # NACA0012 (SST, fixed values for turbulence quantities)
    turb_naca0012_sst_fixedvalues           = TestCase('turb_naca0012_sst_fixedvalues')
    turb_naca0012_sst_fixedvalues.cfg_dir   = "rans/naca0012"
    turb_naca0012_sst_fixedvalues.cfg_file  = "turb_NACA0012_sst_fixedvalues.cfg"
    turb_naca0012_sst_fixedvalues.test_iter = 10
    turb_naca0012_sst_fixedvalues.test_vals = [-5.206744, -9.562435, -1.566603, 1.022029, 0.040549, -3.483576]
    turb_naca0012_sst_fixedvalues.su2_exec  = "SU2_CFD"
    turb_naca0012_sst_fixedvalues.timeout   = 3200
    turb_naca0012_sst_fixedvalues.tol       = 0.00001
    test_list.append(turb_naca0012_sst_fixedvalues)

    # PROPELLER
    propeller           = TestCase('propeller')
    propeller.cfg_dir   = "rans/propeller"
    propeller.cfg_file  = "propeller.cfg"
    propeller.test_iter = 10
    propeller.test_vals = [-3.389575, -8.409529, 0.000048, 0.056329] #last 4 columns
    propeller.su2_exec  = "SU2_CFD"
    propeller.new_output = True
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
    axi_rans_air_nozzle.test_vals = [ -12.092891, -6.630495, -8.784840, -2.399099]
    axi_rans_air_nozzle.su2_exec  = "SU2_CFD"
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
    turb_naca0012_sst_restart_mg.test_iter = 50
    turb_naca0012_sst_restart_mg.ntest_vals = 5
    turb_naca0012_sst_restart_mg.test_vals = [-7.653235, -7.729550, -1.981855, -0.000015, 0.079061]
    turb_naca0012_sst_restart_mg.su2_exec  = "SU2_CFD"
    turb_naca0012_sst_restart_mg.new_output  = True
    turb_naca0012_sst_restart_mg.timeout   = 3200
    turb_naca0012_sst_restart_mg.tol       = 0.000001
    test_list.append(turb_naca0012_sst_restart_mg)

    #########################
    ###    Transition     ###
    #########################

    # Schubauer-Klebanoff Natural Transition
    schubauer_klebanoff_transition              = TestCase('Schubauer_Klebanoff')
    schubauer_klebanoff_transition.cfg_dir      = "transition/Schubauer_Klebanoff"
    schubauer_klebanoff_transition.cfg_file     = "transitional_BC_model_ConfigFile.cfg"
    schubauer_klebanoff_transition.test_iter    = 10
    schubauer_klebanoff_transition.new_output   = True
    schubauer_klebanoff_transition.test_vals    = [-8.029786, -13.240213, 0.000053, 0.007986] #last 4 columns
    schubauer_klebanoff_transition.su2_exec     = "SU2_CFD"
    schubauer_klebanoff_transition.timeout      = 1600
    schubauer_klebanoff_transition.tol          = 0.00001
    test_list.append(schubauer_klebanoff_transition)

    #####################################
    ### Cont. adj. compressible Euler ###
    #####################################

    # Dry run Cont. Adj. Euler
    contadj_naca0012_d           = TestCase('dry run Cont. Adj. Euler')
    contadj_naca0012_d.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012_d.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012_d.su2_exec  = "SU2_CFD -d"
    contadj_naca0012_d.timeout   = 1600
    test_list.append(contadj_naca0012_d)

    # Inviscid NACA0012
    contadj_naca0012           = TestCase('contadj_naca0012')
    contadj_naca0012.cfg_dir   = "cont_adj_euler/naca0012"
    contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
    contadj_naca0012.test_iter = 5
    contadj_naca0012.test_vals = [-9.289565, -14.563861, 0.300920, 0.019552] #last 4 columns
    contadj_naca0012.su2_exec  = "SU2_CFD"
    contadj_naca0012.new_output = True
    contadj_naca0012.timeout   = 1600
    contadj_naca0012.tol       = 0.001
    test_list.append(contadj_naca0012)

    # Inviscid ONERA M6
    contadj_oneram6           = TestCase('contadj_oneram6')
    contadj_oneram6.cfg_dir   = "cont_adj_euler/oneram6"
    contadj_oneram6.cfg_file  = "inv_ONERAM6.cfg"
    contadj_oneram6.test_iter = 10
    contadj_oneram6.test_vals = [-12.133160, -12.706697, 0.685900, 0.007594] #last 4 columns
    contadj_oneram6.su2_exec  = "SU2_CFD"
    contadj_oneram6.new_output = True
    contadj_oneram6.timeout   = 1600
    contadj_oneram6.tol       = 0.00001
    test_list.append(contadj_oneram6)

    # Inviscid WEDGE: tests averaged outflow total pressure adjoint
    contadj_wedge             = TestCase('contadj_wedge')
    contadj_wedge.cfg_dir   = "cont_adj_euler/wedge"
    contadj_wedge.cfg_file  = "inv_wedge_ROE.cfg"
    contadj_wedge.test_iter = 10
    contadj_wedge.test_vals = [2.872691, -2.755572, 853000.000000, 0.000000] #last 4 columns
    contadj_wedge.su2_exec  = "SU2_CFD"
    contadj_wedge.new_output = True
    contadj_wedge.timeout   = 1600
    contadj_wedge.tol       = 0.00001
    test_list.append(contadj_wedge)

    # Inviscid fixed CL NACA0012
    contadj_fixedCL_naca0012           = TestCase('contadj_fixedcl_naca0012')
    contadj_fixedCL_naca0012.cfg_dir   = "fixed_cl/naca0012"
    contadj_fixedCL_naca0012.cfg_file  = "inv_NACA0012_ContAdj.cfg"
    contadj_fixedCL_naca0012.test_iter = 100
    contadj_fixedCL_naca0012.test_vals = [0.293213, -5.201710, 0.360590, -0.000022]
    contadj_fixedCL_naca0012.su2_exec  = "SU2_CFD"
    contadj_fixedCL_naca0012.new_output= True
    contadj_fixedCL_naca0012.timeout   = 1600
    contadj_fixedCL_naca0012.tol       = 0.00001
    test_list.append(contadj_fixedCL_naca0012)

    ###################################
    ### Cont. adj. compressible N-S ###
    ###################################

    # Dry run Cont. Adj. NS
    contadj_ns_cylinder_d           = TestCase('dry run Cont. Adj. NS')
    contadj_ns_cylinder_d.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder_d.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder_d.su2_exec  = "SU2_CFD -d"
    contadj_ns_cylinder_d.timeout   = 1600
    test_list.append(contadj_ns_cylinder_d)

    # Adjoint laminar cylinder
    contadj_ns_cylinder           = TestCase('contadj_ns_cylinder')
    contadj_ns_cylinder.cfg_dir   = "cont_adj_navierstokes/cylinder"
    contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
    contadj_ns_cylinder.test_iter = 20
    contadj_ns_cylinder.test_vals = [ -3.665848, -9.132055, 2.056700, -0.000000] #last 4 columns
    contadj_ns_cylinder.su2_exec  = "SU2_CFD"
    contadj_ns_cylinder.new_output = True
    contadj_ns_cylinder.timeout   = 1600
    contadj_ns_cylinder.tol       = 0.00001
    test_list.append(contadj_ns_cylinder)

    # Adjoint laminar naca0012 subsonic
    contadj_ns_naca0012_sub           = TestCase('contadj_ns_naca0012_sub')
    contadj_ns_naca0012_sub.cfg_dir   = "cont_adj_navierstokes/naca0012_sub"
    contadj_ns_naca0012_sub.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_sub.test_iter = 20
    contadj_ns_naca0012_sub.test_vals = [-2.743268, -8.215193, 0.518810, 0.001210] #last 4 columns
    contadj_ns_naca0012_sub.su2_exec  = "SU2_CFD"
    contadj_ns_naca0012_sub.new_output =True
    contadj_ns_naca0012_sub.timeout   = 1600
    contadj_ns_naca0012_sub.tol       = 0.00001
    test_list.append(contadj_ns_naca0012_sub)

    # Adjoint laminar naca0012 transonic
    contadj_ns_naca0012_trans           = TestCase('contadj_ns_naca0012_trans')
    contadj_ns_naca0012_trans.cfg_dir   = "cont_adj_navierstokes/naca0012_trans"
    contadj_ns_naca0012_trans.cfg_file  = "lam_NACA0012.cfg"
    contadj_ns_naca0012_trans.test_iter = 20
    contadj_ns_naca0012_trans.test_vals = [-1.039664, -6.575019, 1.772300, 0.012495] #last 4 columns
    contadj_ns_naca0012_trans.su2_exec  = "SU2_CFD"
    contadj_ns_naca0012_trans.new_output = True
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
    contadj_rans_naca0012.test_vals = [ -0.794162, -5.761722, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012.su2_exec  = "SU2_CFD"
    contadj_rans_naca0012.new_output = True
    contadj_rans_naca0012.timeout   = 1600
    contadj_rans_naca0012.tol       = 0.00001
    test_list.append(contadj_rans_naca0012)

    # Adjoint turbulent NACA0012 with binary restarts
    contadj_rans_naca0012_bin           = TestCase('contadj_rans_naca0012_bin')
    contadj_rans_naca0012_bin.cfg_dir   = "cont_adj_rans/naca0012"
    contadj_rans_naca0012_bin.cfg_file  = "turb_nasa_binary.cfg"
    contadj_rans_naca0012_bin.test_iter = 18
    contadj_rans_naca0012_bin.test_vals = [-0.794169, -5.761671, 19.214000, -0.000000] #last 4 columns
    contadj_rans_naca0012_bin.su2_exec  = "SU2_CFD"
    contadj_rans_naca0012_bin.new_output = True
    contadj_rans_naca0012_bin.timeout   = 1600
    contadj_rans_naca0012_bin.tol       = 0.00001
    test_list.append(contadj_rans_naca0012_bin)

    # Adjoint turbulent RAE2822
    contadj_rans_rae2822           = TestCase('contadj_rans_rae2822')
    contadj_rans_rae2822.cfg_dir   = "cont_adj_rans/rae2822"
    contadj_rans_rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
    contadj_rans_rae2822.test_iter = 20
    contadj_rans_rae2822.test_vals = [-5.369688, -10.872211, -0.212470, 0.005448]
    contadj_rans_rae2822.su2_exec  = "SU2_CFD"
    contadj_rans_rae2822.new_output = True
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
    turb_naca0012_1c.test_vals = [-4.985899, 1.137425, 0.375998, -0.078790]
    turb_naca0012_1c.su2_exec  = "SU2_CFD"
    turb_naca0012_1c.new_output = True
    turb_naca0012_1c.timeout   = 1600
    turb_naca0012_1c.tol       = 0.00001
    test_list.append(turb_naca0012_1c)

    # NACA0012 2c
    turb_naca0012_2c           = TestCase('turb_naca0012_2c')
    turb_naca0012_2c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_2c.cfg_file  = "turb_NACA0012_uq_2c.cfg"
    turb_naca0012_2c.test_iter = 10
    turb_naca0012_2c.test_vals = [-5.483374, 0.968842, 0.255857, -0.109404] #last 4 columns
    turb_naca0012_2c.su2_exec  = "SU2_CFD"
    turb_naca0012_2c.new_output = True
    turb_naca0012_2c.timeout   = 1600
    turb_naca0012_2c.tol       = 0.00001
    test_list.append(turb_naca0012_2c)

    # NACA0012 3c
    turb_naca0012_3c           = TestCase('turb_naca0012_3c')
    turb_naca0012_3c.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_3c.cfg_file  = "turb_NACA0012_uq_3c.cfg"
    turb_naca0012_3c.test_iter = 10
    turb_naca0012_3c.test_vals = [-5.584300, 0.931352, 0.248990, -0.109405] #last 4 columns
    turb_naca0012_3c.su2_exec  = "SU2_CFD"
    turb_naca0012_3c.new_output = True
    turb_naca0012_3c.timeout   = 1600
    turb_naca0012_3c.tol       = 0.00001
    test_list.append(turb_naca0012_3c)

    # NACA0012 p1c1
    turb_naca0012_p1c1           = TestCase('turb_naca0012_p1c1')
    turb_naca0012_p1c1.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c1.cfg_file  = "turb_NACA0012_uq_p1c1.cfg"
    turb_naca0012_p1c1.test_iter = 10
    turb_naca0012_p1c1.test_vals = [-5.127531, 1.077321, 0.543822, -0.018663] #last 4 columns
    turb_naca0012_p1c1.su2_exec  = "SU2_CFD"
    turb_naca0012_p1c1.new_output = True
    turb_naca0012_p1c1.timeout   = 1600
    turb_naca0012_p1c1.tol       = 0.00001
    test_list.append(turb_naca0012_p1c1)

    # NACA0012 p1c2
    turb_naca0012_p1c2           = TestCase('turb_naca0012_p1c2')
    turb_naca0012_p1c2.cfg_dir   = "rans_uq/naca0012"
    turb_naca0012_p1c2.cfg_file  = "turb_NACA0012_uq_p1c2.cfg"
    turb_naca0012_p1c2.test_iter = 10
    turb_naca0012_p1c2.test_vals = [-5.554491, 0.943648, 0.320568, -0.091089] #last 4 columns
    turb_naca0012_p1c2.su2_exec  = "SU2_CFD"
    turb_naca0012_p1c2.new_output = True
    turb_naca0012_p1c2.timeout   = 1600
    turb_naca0012_p1c2.tol       = 0.00001
    test_list.append(turb_naca0012_p1c2)

    ######################################
    ### Moving Wall                    ###
    ######################################

    # Rotating NACA 0012
    rot_naca0012           = TestCase('rot_naca0012')
    rot_naca0012.cfg_dir   = "rotating/naca0012"
    rot_naca0012.cfg_file  = "rot_NACA0012.cfg"
    rot_naca0012.test_iter = 25
    rot_naca0012.test_vals = [-2.688979, 2.857521, -0.079219, 0.002135]
    rot_naca0012.su2_exec  = "SU2_CFD"
    rot_naca0012.timeout   = 1600
    rot_naca0012.tol       = 0.00001
    test_list.append(rot_naca0012)

    # Lid-driven cavity
    cavity           = TestCase('cavity')
    cavity.cfg_dir   = "moving_wall/cavity"
    cavity.cfg_file  = "lam_cavity.cfg"
    cavity.test_iter = 25
    cavity.test_vals = [-5.627934, -0.164470, 0.051972, 2.547039]
    cavity.su2_exec  = "SU2_CFD"
    cavity.new_output = True
    cavity.timeout   = 1600
    cavity.tol       = 0.00001
    test_list.append(cavity)

    # Spinning cylinder
    spinning_cylinder           = TestCase('spinning_cylinder')
    spinning_cylinder.cfg_dir   = "moving_wall/spinning_cylinder"
    spinning_cylinder.cfg_file  = "spinning_cylinder.cfg"
    spinning_cylinder.test_iter = 25
    spinning_cylinder.test_vals = [-7.889994, -2.469385, 1.708162, 1.670039] #last 4 columns
    spinning_cylinder.su2_exec  = "SU2_CFD"
    spinning_cylinder.new_output = True
    spinning_cylinder.timeout   = 1600
    spinning_cylinder.tol       = 0.00001
    test_list.append(spinning_cylinder)

    ######################################
    ### Unsteady                       ###
    ######################################

    # Square cylinder
    square_cylinder           = TestCase('square_cylinder')
    square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    square_cylinder.cfg_file  = "turb_square.cfg"
    square_cylinder.test_iter = 3
    square_cylinder.test_vals = [-1.162561, 0.066414, 1.399788, 2.220411] #last 4 columns
    square_cylinder.su2_exec  = "SU2_CFD"
    square_cylinder.timeout   = 1600
    square_cylinder.tol       = 0.00001
    square_cylinder.unsteady  = True
    square_cylinder.new_output = True
    test_list.append(square_cylinder)

    # Delayed Detached Eddy Simulation
    ddes_flatplate        = TestCase('ddes_flatplate')
    ddes_flatplate.cfg_dir   = "ddes/flatplate"
    ddes_flatplate.cfg_file  = "ddes_flatplate.cfg"
    ddes_flatplate.test_iter = 10
    ddes_flatplate.test_vals = [-2.714758, -5.883004, -0.215005, 0.023783] #last 4 columns
    ddes_flatplate.su2_exec  = "SU2_CFD"
    ddes_flatplate.timeout   = 1600
    ddes_flatplate.tol       = 0.00001
    ddes_flatplate.unsteady  = True
    ddes_flatplate.new_output = True
    test_list.append(ddes_flatplate)

    # unsteady pitching NACA0015, SA
    unst_inc_turb_naca0015_sa           = TestCase('unst_inc_turb_naca0015_sa')
    unst_inc_turb_naca0015_sa.cfg_dir   = "unsteady/pitching_naca0015_rans_inc"
    unst_inc_turb_naca0015_sa.cfg_file  = "config_incomp_turb_sa.cfg"
    unst_inc_turb_naca0015_sa.test_iter = 1
    unst_inc_turb_naca0015_sa.test_vals = [-3.007635, -6.879789, 1.445300, 0.419281] #last 4 columns
    unst_inc_turb_naca0015_sa.su2_exec  = "SU2_CFD"
    unst_inc_turb_naca0015_sa.timeout   = 1600
    unst_inc_turb_naca0015_sa.tol       = 0.00001
    unst_inc_turb_naca0015_sa.unsteady  = True
    test_list.append(unst_inc_turb_naca0015_sa)

    # unsteady pitching NACA0012, Euler, Deforming
    unst_deforming_naca0012           = TestCase('unst_deforming_naca0012')
    unst_deforming_naca0012.cfg_dir   = "disc_adj_euler/naca0012_pitching_def"
    unst_deforming_naca0012.cfg_file  = "inv_NACA0012_pitching_deform.cfg"
    unst_deforming_naca0012.test_iter = 5
    unst_deforming_naca0012.test_vals = [-3.665128, -3.793593, -3.716506, -3.148308] #last 4 columns
    unst_deforming_naca0012.su2_exec  = "SU2_CFD"
    unst_deforming_naca0012.timeout   = 1600
    unst_deforming_naca0012.tol       = 0.00001
    unst_deforming_naca0012.unsteady  = True
    test_list.append(unst_deforming_naca0012)

    ######################################
    ### NICFD                          ###
    ######################################

    # ls89_sa
    ls89_sa           = TestCase('ls89_sa')
    ls89_sa.cfg_dir   = "nicf/LS89"
    ls89_sa.cfg_file  = "turb_SA_PR.cfg"
    ls89_sa.test_iter = 20
    ls89_sa.test_vals = [-5.050483, -13.389547, 0.174939, 0.430757] #last 4 columns
    ls89_sa.su2_exec  = "SU2_CFD"
    ls89_sa.new_output= True
    ls89_sa.timeout   = 1600
    ls89_sa.tol       = 0.00001
    test_list.append(ls89_sa)

    # Rarefaction shock wave edge_VW
    edge_VW           = TestCase('edge_VW')
    edge_VW.cfg_dir   = "nicf/edge"
    edge_VW.cfg_file  = "edge_VW.cfg"
    edge_VW.test_iter = 20
    edge_VW.test_vals = [-0.711552, 5.490479, -0.000975, 0.000000] #last 4 columns
    edge_VW.su2_exec  = "SU2_CFD"
    edge_VW.new_output = True
    edge_VW.timeout   = 1600
    edge_VW.tol       = 0.00001
    test_list.append(edge_VW)

    # Rarefaction shock wave edge_PPR
    edge_PPR           = TestCase('edge_PPR')
    edge_PPR.cfg_dir   = "nicf/edge"
    edge_PPR.cfg_file  = "edge_PPR.cfg"
    edge_PPR.test_iter = 20
    edge_PPR.test_vals = [-1.670439, 4.522842, 0.001027, 0.000000] #last 4 columns
    edge_PPR.su2_exec  = "SU2_CFD"
    edge_PPR.new_output = True
    edge_PPR.timeout   = 1600
    edge_PPR.tol       = 0.00001
    test_list.append(edge_PPR)

    ######################################
    ### Sliding Mesh                   ###
    ######################################

    # Dry run Multizone
    uniform_flow_d         = TestCase('dry run Multizone')
    uniform_flow_d.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow_d.cfg_file  = "uniform_NN.cfg"
    uniform_flow_d.su2_exec  = "SU2_CFD -d"
    uniform_flow_d.timeout   = 1600
    test_list.append(uniform_flow_d)

    # Uniform flow
    uniform_flow         = TestCase('uniform_flow')
    uniform_flow.cfg_dir   = "sliding_interface/uniform_flow"
    uniform_flow.cfg_file  = "uniform_NN.cfg"
    uniform_flow.test_iter = 2
    uniform_flow.test_vals = [2.000000, 0.000000, -0.205134, -13.250256] #last 4 columns
    uniform_flow.su2_exec  = "SU2_CFD"
    uniform_flow.timeout   = 1600
    uniform_flow.tol       = 0.000001
    uniform_flow.unsteady  = True
    uniform_flow.multizone = True
    test_list.append(uniform_flow)

   # Channel_2D
    channel_2D           = TestCase('channel_2D')
    channel_2D.cfg_dir   = "sliding_interface/channel_2D"
    channel_2D.cfg_file  = "channel_2D_WA.cfg"
    channel_2D.test_iter = 2
    channel_2D.test_vals = [2.000000, 0.000000, 0.398005, 0.352783, 0.405475] #last 5 columns
    channel_2D.su2_exec  = "SU2_CFD"
    channel_2D.timeout   = 100
    channel_2D.tol       = 0.00001
    channel_2D.unsteady  = True
    channel_2D.multizone = True
    test_list.append(channel_2D)

    # Channel_3D
    channel_3D           = TestCase('channel_3D')
    channel_3D.cfg_dir   = "sliding_interface/channel_3D"
    channel_3D.cfg_file  = "channel_3D_WA.cfg"
    channel_3D.test_iter = 1
    channel_3D.test_vals = [1.000000, 0.000000, 0.661408, 0.769988, 0.696033] #last 5 columns
    channel_3D.su2_exec  = "SU2_CFD"
    channel_3D.timeout   = 1600
    channel_3D.tol       = 0.00001
    channel_3D.unsteady  = True
    channel_3D.multizone = True
    test_list.append(channel_3D)

    # Pipe
    pipe           = TestCase('pipe')
    pipe.cfg_dir   = "sliding_interface/pipe"
    pipe.cfg_file  = "pipe_NN.cfg"
    pipe.test_iter = 2
    pipe.test_vals = [0.491954, 0.677756, 0.963981, 1.006936] #last 4 columns
    pipe.su2_exec  = "SU2_CFD"
    pipe.timeout   = 1600
    pipe.tol       = 0.00001
    pipe.unsteady  = True
    pipe.multizone = True
    test_list.append(pipe)

    # Rotating cylinders
    rotating_cylinders           = TestCase('rotating_cylinders')
    rotating_cylinders.cfg_dir   = "sliding_interface/rotating_cylinders"
    rotating_cylinders.cfg_file  = "rot_cylinders_WA.cfg"
    rotating_cylinders.test_iter = 3
    rotating_cylinders.test_vals = [3.000000, 0.000000, 0.777574, 1.134794, 1.224127] #last 4 columns
    rotating_cylinders.su2_exec  = "SU2_CFD"
    rotating_cylinders.timeout   = 1600
    rotating_cylinders.tol       = 0.00001
    rotating_cylinders.unsteady  = True
    rotating_cylinders.multizone = True
    test_list.append(rotating_cylinders)

    # Supersonic vortex shedding
    supersonic_vortex_shedding           = TestCase('supersonic_vortex_shedding')
    supersonic_vortex_shedding.cfg_dir   = "sliding_interface/supersonic_vortex_shedding"
    supersonic_vortex_shedding.cfg_file  = "sup_vor_shed_WA.cfg"
    supersonic_vortex_shedding.test_iter = 5
    supersonic_vortex_shedding.test_vals = [5.000000, 0.000000, 1.227921, 1.638901] #last 4 columns
    supersonic_vortex_shedding.su2_exec  = "SU2_CFD"
    supersonic_vortex_shedding.timeout   = 1600
    supersonic_vortex_shedding.tol       = 0.00001
    supersonic_vortex_shedding.unsteady  = True
    supersonic_vortex_shedding.multizone = True
    test_list.append(supersonic_vortex_shedding)

    # Bars_SST_2D
    bars_SST_2D           = TestCase('bars_SST_2D')
    bars_SST_2D.cfg_dir   = "sliding_interface/bars_SST_2D"
    bars_SST_2D.cfg_file  = "bars.cfg"
    bars_SST_2D.test_iter = 13
    bars_SST_2D.test_vals = [13.000000, -0.619686, -1.564594]
    bars_SST_2D.su2_exec  = "SU2_CFD"
    bars_SST_2D.timeout   = 1600
    bars_SST_2D.tol       = 0.00001
    bars_SST_2D.multizone = True
    test_list.append(bars_SST_2D)

    # Sliding mesh with incompressible flows (steady)
    slinc_steady           = TestCase('slinc_steady')
    slinc_steady.cfg_dir   = "sliding_interface/incompressible_steady"
    slinc_steady.cfg_file  = "config.cfg"
    slinc_steady.test_iter = 19
    slinc_steady.test_vals = [19.000000, -1.803326, -2.097400] #last 3 columns
    slinc_steady.su2_exec  = "SU2_CFD"
    slinc_steady.timeout   = 100
    slinc_steady.tol       = 0.00001
    slinc_steady.multizone = True
    test_list.append(slinc_steady)

    # Sliding mesh with incompressible flows (unsteady)
    # slinc_unsteady           = TestCase('slinc_unsteady')
    # slinc_unsteady.cfg_dir   = "sliding_interface/incompressible_unsteady"
    # slinc_unsteady.cfg_file  = "config.cfg"
    # slinc_unsteady.test_iter = 19
    # slinc_unsteady.test_vals = [-3.515218,1.930028,0.000000,0.000000] #last 4 columns
    # slinc_unsteady.su2_exec  = "SU2_CFD"
    # slinc_unsteady.timeout   = 100
    # slinc_unsteady.tol       = 0.00001
    # slinc_unsteady.unsteady  = True
    # test_list.append(slinc_unsteady)

    ##############################################
    ### Method of Manufactured Solutions (MMS) ###
    ##############################################

    # FVM, compressible, laminar N-S
    mms_fvm_ns           = TestCase('mms_fvm_ns')
    mms_fvm_ns.cfg_dir   = "mms/fvm_navierstokes"
    mms_fvm_ns.cfg_file  = "lam_mms_roe.cfg"
    mms_fvm_ns.test_iter = 20
    mms_fvm_ns.test_vals = [-2.851428, 2.192348, 0.000000, 0.000000] #last 4 columns
    mms_fvm_ns.su2_exec  = "SU2_CFD"
    mms_fvm_ns.timeout   = 1600
    mms_fvm_ns.tol       = 0.0001
    test_list.append(mms_fvm_ns)

    ######################################
    ### RUN TESTS                      ###
    ######################################

    pass_list = [ test.run_test() for test in test_list ]


    ######################################
    ### RUN SU2_GEO TESTS              ###
    ######################################

    # NACA0012
    naca0012_geo           = TestCase('naca0012_geo')
    naca0012_geo.cfg_dir   = "optimization_euler/steady_naca0012"
    naca0012_geo.cfg_file  = "inv_NACA0012_adv.cfg"
    naca0012_geo.test_vals = [1.0000, 62.0455, 0.120011, 0.0000] #chord, LE radius, ToC, Alpha
    naca0012_geo.su2_exec  = "SU2_GEO"
    naca0012_geo.timeout   = 1600
    naca0012_geo.tol       = 0.00001
    pass_list.append(naca0012_geo.run_geo())
    test_list.append(naca0012_geo)

    ######################################
    ### RUN SU2_DEF TESTS              ###
    ######################################

    # intersection prevention
    intersect_def            = TestCase('intersectionprevention')
    intersect_def.cfg_dir   = "deformation/intersection_prevention"
    intersect_def.cfg_file  = "def_intersect.cfg"
    intersect_def.test_iter = 10
    intersect_def.test_vals = [0.000112] #residual
    intersect_def.su2_exec  = "SU2_DEF"
    intersect_def.timeout   = 1600
    intersect_def.tol       = 1e-04

    pass_list.append(intersect_def.run_def())
    test_list.append(intersect_def)

    # Inviscid NACA0012 (triangles)
    naca0012_def            = TestCase('naca0012_def')
    naca0012_def.cfg_dir   = "deformation/naca0012"
    naca0012_def.cfg_file  = "def_NACA0012.cfg"
    naca0012_def.test_iter = 10
    naca0012_def.test_vals = [0.00344658] #residual
    naca0012_def.su2_exec  = "SU2_DEF"
    naca0012_def.timeout   = 1600
    naca0012_def.tol       = 1e-08

    pass_list.append(naca0012_def.run_def())
    test_list.append(naca0012_def)

    # Inviscid NACA0012 based on SURFACE_FILE input (surface_bump.dat)
    naca0012_def_file            = TestCase('naca0012_def_file')
    naca0012_def_file.cfg_dir   = "deformation/naca0012"
    naca0012_def_file.cfg_file  = "surface_file_NACA0012.cfg"
    naca0012_def_file.test_iter = 10
    naca0012_def_file.test_vals = [0.00344658] #residual
    naca0012_def_file.su2_exec  = "SU2_DEF"
    naca0012_def_file.timeout   = 1600
    naca0012_def_file.tol       = 1e-8

    pass_list.append(naca0012_def_file.run_def())
    test_list.append(naca0012_def_file)

    # RAE2822 (mixed tris + quads)
    rae2822_def            = TestCase('rae2822_def')
    rae2822_def.cfg_dir   = "deformation/rae2822"
    rae2822_def.cfg_file  = "def_RAE2822.cfg"
    rae2822_def.test_iter = 10
    rae2822_def.test_vals = [7.94218e-09] #residual
    rae2822_def.su2_exec  = "SU2_DEF"
    rae2822_def.timeout   = 1600
    rae2822_def.tol       = 1e-13

    pass_list.append(rae2822_def.run_def())
    test_list.append(rae2822_def)

    # Turb NACA4412 (quads, wall distance)
    naca4412_def            = TestCase('naca4412_def')
    naca4412_def.cfg_dir   = "deformation/naca4412"
    naca4412_def.cfg_file  = "def_NACA4412.cfg"
    naca4412_def.test_iter = 10
    naca4412_def.test_vals = [8.855370e-13] #residual
    naca4412_def.su2_exec  = "SU2_DEF"
    naca4412_def.timeout   = 1600
    naca4412_def.tol       = 1e-12

    pass_list.append(naca4412_def.run_def())
    test_list.append(naca4412_def)

    # Brick of tets (inverse volume)
    brick_tets_def            = TestCase('brick_tets_def')
    brick_tets_def.cfg_dir   = "deformation/brick_tets"
    brick_tets_def.cfg_file  = "def_brick_tets.cfg"
    brick_tets_def.test_iter = 10
    brick_tets_def.test_vals = [8.973010e-04] #residual
    brick_tets_def.su2_exec  = "SU2_DEF"
    brick_tets_def.timeout   = 1600
    brick_tets_def.tol       = 1e-09

    pass_list.append(brick_tets_def.run_def())
    test_list.append(brick_tets_def)

    # Brick of isotropic hexas (inverse volume)
    brick_hex_def           = TestCase('brick_hex_def')
    brick_hex_def.cfg_dir   = "deformation/brick_hex"
    brick_hex_def.cfg_file  = "def_brick_hex.cfg"
    brick_hex_def.test_iter = 10
    brick_hex_def.test_vals = [2.082100e-04] #residual
    brick_hex_def.su2_exec  = "SU2_DEF"
    brick_hex_def.timeout   = 1600
    brick_hex_def.tol       = 1e-09

    pass_list.append(brick_hex_def.run_def())
    test_list.append(brick_hex_def)

    # Brick with a pyramid layer (inverse volume)
    brick_pyra_def           = TestCase('brick_pyra_def')
    brick_pyra_def.cfg_dir   = "deformation/brick_pyra"
    brick_pyra_def.cfg_file  = "def_brick_pyra.cfg"
    brick_pyra_def.test_iter = 10
    brick_pyra_def.test_vals = [0.00150063] #residual
    brick_pyra_def.su2_exec  = "SU2_DEF"
    brick_pyra_def.timeout   = 1600
    brick_pyra_def.tol       = 1e-08

    pass_list.append(brick_pyra_def.run_def())
    test_list.append(brick_pyra_def)

    # Brick of isotropic prisms (inverse volume)
    brick_prism_def           = TestCase('brick_prism_def')
    brick_prism_def.cfg_dir   = "deformation/brick_prism"
    brick_prism_def.cfg_file  = "def_brick_prism.cfg"
    brick_prism_def.test_iter = 10
    brick_prism_def.test_vals = [0.00212069] #residual
    brick_prism_def.su2_exec  = "SU2_DEF"
    brick_prism_def.timeout   = 1600
    brick_prism_def.tol       = 1e-08

    pass_list.append(brick_prism_def.run_def())
    test_list.append(brick_prism_def)

    # Brick of prisms with high aspect ratio cells near the wall (wall distance)
    brick_prism_rans_def           = TestCase('brick_prism_rans_def')
    brick_prism_rans_def.cfg_dir   = "deformation/brick_prism_rans"
    brick_prism_rans_def.cfg_file  = "def_brick_prism_rans.cfg"
    brick_prism_rans_def.test_iter = 10
    brick_prism_rans_def.test_vals = [4.8066e-08] #residual
    brick_prism_rans_def.su2_exec  = "SU2_DEF"
    brick_prism_rans_def.timeout   = 1600
    brick_prism_rans_def.tol       = 1e-12

    pass_list.append(brick_prism_rans_def.run_def())
    test_list.append(brick_prism_rans_def)

    # Brick of hexas with high aspect ratio cells near the wall (inverse volume)
    brick_hex_rans_def           = TestCase('brick_hex_rans_def')
    brick_hex_rans_def.cfg_dir   = "deformation/brick_hex_rans"
    brick_hex_rans_def.cfg_file  = "def_brick_hex_rans.cfg"
    brick_hex_rans_def.test_iter = 10
    brick_hex_rans_def.test_vals = [2.260750e-07] #residual
    brick_hex_rans_def.su2_exec  = "SU2_DEF"
    brick_hex_rans_def.timeout   = 1600
    brick_hex_rans_def.tol       = 1e-12

    pass_list.append(brick_hex_rans_def.run_def())
    test_list.append(brick_hex_rans_def)

    # Cylindrical FFD test
    cylinder_ffd_def           = TestCase('cylinder_ffd_def')
    cylinder_ffd_def.cfg_dir   = "deformation/cylindrical_ffd"
    cylinder_ffd_def.cfg_file  = "def_cylindrical.cfg"
    cylinder_ffd_def.test_iter = 10
    cylinder_ffd_def.test_vals = [0.000470133] #residual
    cylinder_ffd_def.su2_exec  = "SU2_DEF"
    cylinder_ffd_def.timeout   = 1600
    cylinder_ffd_def.tol       = 1e-09

    pass_list.append(cylinder_ffd_def.run_def())
    test_list.append(cylinder_ffd_def)

    # Spherical FFD test
    sphere_ffd_def           = TestCase('sphere_ffd_def')
    sphere_ffd_def.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def.cfg_file  = "def_spherical.cfg"
    sphere_ffd_def.test_iter = 10
    sphere_ffd_def.test_vals = [0.00356699] #residual
    sphere_ffd_def.su2_exec  = "SU2_DEF"
    sphere_ffd_def.timeout   = 1600
    sphere_ffd_def.tol       = 1e-08

    pass_list.append(sphere_ffd_def.run_def())
    test_list.append(sphere_ffd_def)

    # Spherical FFD test using BSplines
    sphere_ffd_def_bspline           = TestCase('sphere_ffd_def_bspline')
    sphere_ffd_def_bspline.cfg_dir   = "deformation/spherical_ffd"
    sphere_ffd_def_bspline.cfg_file  = "def_spherical_bspline.cfg"
    sphere_ffd_def_bspline.test_iter = 10
    sphere_ffd_def_bspline.test_vals = [0.00206808] #residual
    sphere_ffd_def_bspline.su2_exec  = "SU2_DEF"
    sphere_ffd_def_bspline.timeout   = 1600
    sphere_ffd_def_bspline.tol       = 1e-08

    pass_list.append(sphere_ffd_def_bspline.run_def())
    test_list.append(sphere_ffd_def_bspline)

    ######################################
    ### RUN PYTHON TESTS               ###
    ######################################

    # test continuous_adjoint.py
    contadj_euler_py = TestCase('contadj_euler_py')
    contadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    contadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    contadj_euler_py.test_iter = 10
    contadj_euler_py.su2_exec  = "continuous_adjoint.py -f"
    contadj_euler_py.timeout   = 1600
    contadj_euler_py.reference_file = "of_grad_cd.dat.ref"
    contadj_euler_py.test_file = "of_grad_cd.dat"
    contadj_euler_py.new_output = True
    pass_list.append(contadj_euler_py.run_filediff())
    test_list.append(contadj_euler_py)

    # test shape_optimization.py
    shape_opt_euler_py           = TestCase('shape_opt_euler_py')
    shape_opt_euler_py.cfg_dir   = "optimization_euler/steady_naca0012"
    shape_opt_euler_py.cfg_file  = "inv_NACA0012_adv.cfg"
    shape_opt_euler_py.test_iter = 1
    shape_opt_euler_py.test_vals = [1, 1, 2.134974E-05, 0.003847] #last 4 columns
    shape_opt_euler_py.su2_exec  = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
    shape_opt_euler_py.timeout   = 1600
    shape_opt_euler_py.new_output = True
    shape_opt_euler_py.tol       = 0.00001
    pass_list.append(shape_opt_euler_py.run_opt())
    test_list.append(shape_opt_euler_py)

    # Multiple functionals with the continuous adjoint
    contadj_multi_py            = TestCase('contadj_multi_py')
    contadj_multi_py.cfg_dir    = "cont_adj_euler/wedge"
    contadj_multi_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
    contadj_multi_py.test_iter  = 10
    contadj_multi_py.su2_exec   = "continuous_adjoint.py -f"
    contadj_multi_py.timeout    = 1600
    contadj_multi_py.new_output = True
    contadj_multi_py.reference_file = "of_grad_combo.dat.ref"
    contadj_multi_py.test_file  = "of_grad_combo.dat"
    pass_list.append(contadj_multi_py.run_filediff())
    test_list.append(contadj_multi_py)

    # Optimization with multiple objectives, with gradients evaluated individually
    # the difference in gradient value relative to combined case
    # is due to lack of solution file for the adjoint and small number of iterations
#    opt_multiobj_py            = TestCase('opt_multiobj_py')
#    opt_multiobj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
#    opt_multiobj_py.cfg_file   = "inv_wedge_ROE_multiobj.cfg"
#    opt_multiobj_py.test_iter  = 1
#    opt_multiobj_py.test_vals = [1, 1, 1.084701E+02, 3.799222E+00] #last 4 columns
#    opt_multiobj_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
#    opt_multiobj_py.timeout    = 1600
#    opt_multiobj_py.tol       = 0.00001
#    pass_list.append(opt_multiobj_py.run_opt())
#    test_list.append(opt_multiobj_py)
#
#    # test optimization, with multiple objectives and gradient evaluated as 'combo'
#    opt_multiobjcombo_py            = TestCase('opt_multiobjcombo_py')
#    opt_multiobjcombo_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
#    opt_multiobjcombo_py.cfg_file   = "inv_wedge_ROE_multiobj_combo.cfg"
#    opt_multiobjcombo_py.test_iter  = 1
#    opt_multiobjcombo_py.test_vals = [1, 1, 1.084701E+02, 3.789322E+00] #last 4 columns
#    opt_multiobjcombo_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT -f"
#    opt_multiobjcombo_py.timeout    = 1600
#    opt_multiobjcombo_py.tol       = 0.00001
#    pass_list.append(opt_multiobjcombo_py.run_opt())
#    test_list.append(opt_multiobjcombo_py)

    # test optimization, with multiple objectives evaluated on a single surface
    opt_multiobj1surf_py            = TestCase('opt_multiobj1surf_py')
    opt_multiobj1surf_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_multiobj1surf_py.cfg_file   = "inv_wedge_ROE_multiobj_1surf.cfg"
    opt_multiobj1surf_py.test_iter  = 1
    opt_multiobj1surf_py.test_vals = [1.000000, 1.000000, 30.428280, 2.039416] #last 4 columns
    opt_multiobj1surf_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT -f "
    opt_multiobj1surf_py.timeout    = 1600
    opt_multiobj1surf_py.tol       = 0.00001
    pass_list.append(opt_multiobj1surf_py.run_opt())
    test_list.append(opt_multiobj1surf_py)

    # test optimization, with a single objective evaluated on multiple surfaces
    opt_2surf1obj_py            = TestCase('opt_2surf1obj_py')
    opt_2surf1obj_py.cfg_dir    = "optimization_euler/multiobjective_wedge"
    opt_2surf1obj_py.cfg_file   = "inv_wedge_ROE_2surf_1obj.cfg"
    opt_2surf1obj_py.test_iter  = 1
    opt_2surf1obj_py.test_vals = [1.000000, 1.000000, 2.005694, 0.000185] #last 4 columns
    opt_2surf1obj_py.su2_exec   = "shape_optimization.py -g CONTINUOUS_ADJOINT  -f"
    opt_2surf1obj_py.timeout    = 1600
    opt_2surf1obj_py.tol       = 0.00001
    pass_list.append(opt_2surf1obj_py.run_opt())
    test_list.append(opt_2surf1obj_py)

    ##########################
    ###   Python wrapper   ###
    ##########################

    # NACA0012
    pywrapper_naca0012           = TestCase('pywrapper_naca0012')
    pywrapper_naca0012.cfg_dir   = "euler/naca0012"
    pywrapper_naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    pywrapper_naca0012.test_iter = 20
    pywrapper_naca0012.test_vals = [-4.023999, -3.515034, 0.339426, 0.022217] #last 4 columns
    pywrapper_naca0012.su2_exec  = "SU2_CFD.py -f"
    pywrapper_naca0012.new_output  = True
    pywrapper_naca0012.timeout   = 1600
    pywrapper_naca0012.tol       = 0.00001
    test_list.append(pywrapper_naca0012)
    pass_list.append(pywrapper_naca0012.run_test())

    # NACA0012 (SST, FUN3D results for finest grid: CL=1.0840, CD=0.01253)
    pywrapper_turb_naca0012_sst           = TestCase('pywrapper_turb_naca0012_sst')
    pywrapper_turb_naca0012_sst.cfg_dir   = "rans/naca0012"
    pywrapper_turb_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    pywrapper_turb_naca0012_sst.test_iter = 10
    pywrapper_turb_naca0012_sst.test_vals = [-11.451010, -12.798258, -5.863895, 1.049989, 0.019163, -1.925018]
    pywrapper_turb_naca0012_sst.su2_exec  = "SU2_CFD.py -f"
    pywrapper_turb_naca0012_sst.new_output = True
    pywrapper_turb_naca0012_sst.timeout   = 3200
    pywrapper_turb_naca0012_sst.tol       = 0.00001
    test_list.append(pywrapper_turb_naca0012_sst)
    pass_list.append(pywrapper_turb_naca0012_sst.run_test())

    # Square cylinder
    pywrapper_square_cylinder           = TestCase('pywrapper_square_cylinder')
    pywrapper_square_cylinder.cfg_dir   = "unsteady/square_cylinder"
    pywrapper_square_cylinder.cfg_file  = "turb_square.cfg"
    pywrapper_square_cylinder.test_iter = 3
    pywrapper_square_cylinder.test_vals = [-1.162560, 0.066414, 1.399788, 2.220411]
    pywrapper_square_cylinder.su2_exec  = "SU2_CFD.py -f"
    pywrapper_square_cylinder.timeout   = 1600
    pywrapper_square_cylinder.tol       = 0.00001
    pywrapper_square_cylinder.unsteady  = True
    test_list.append(pywrapper_square_cylinder)
    pass_list.append(pywrapper_square_cylinder.run_test())

    # Rigid motion
    pywrapper_rigidMotion               = TestCase('pywrapper_rigidMotion')
    pywrapper_rigidMotion.cfg_dir       = "py_wrapper/flatPlate_rigidMotion"
    pywrapper_rigidMotion.cfg_file      = "flatPlate_rigidMotion_Conf.cfg"
    pywrapper_rigidMotion.test_iter     = 5
    pywrapper_rigidMotion.test_vals     = [-1.614170, 2.242953, 0.350050, 0.093137]
    pywrapper_rigidMotion.su2_exec      = "python launch_flatPlate_rigidMotion.py -f"
    pywrapper_rigidMotion.new_output      = True
    pywrapper_rigidMotion.timeout       = 1600
    pywrapper_rigidMotion.tol           = 0.00001
    pywrapper_rigidMotion.unsteady      = True
    test_list.append(pywrapper_rigidMotion)
    pass_list.append(pywrapper_rigidMotion.run_test())

    # Tests summary
    print('==================================================================')
    print('Summary of the serial tests')
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
