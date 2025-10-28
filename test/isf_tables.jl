#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Version 0.1, 2025.
#
# If you found this piece of software useful, please cite the paper:
# J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain,
# "Primal and Dual Approaches for the Chamber Enumeration
# of real hyperplane arrangements", 2025.
#
# Authors:
# - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
# - Jean Charles Gilbert (INRIA, France),
# - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
#
# Copyright 2025, INRIA (France) and Univ. of Sherbrooke (Canada).
#
# ISF is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program. If not, see
# <https://doc.qt.io/archives/3.3/license.html>.
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

### File used to run some benchmark and comparisons - does not contain any function.

# using DelimitedFiles
# using Printf
# using HiGHS

#include("../src/isf_main.jl")   JPD
#JPD
include("../src/isf_benchmark_values.jl")          # tool used in other benchmarking functions   
#JPD

# defines the groups of instances per types
# list_rand = ["rand_5_10", "rand_4_11", "rand_6_12", "rand_5_13", "rand_7_14", "rand_7_15", "rand_8_16", "rand_9_17"] # "rand_2_8", "rand_4_8", "rand_4_9", 
# list_2d = ["2d_5", "2d_6", "2d_7", "2d_8"] # "2d_4", 
# list_srand = ["srand_8_20_2", "srand_8_20_4", "srand_8_20_6"]
# list_perm = ["perm_5", "perm_6", "perm_7", "perm_8"]
# list_ratio = ["ratio_5_20_07", "ratio_5_20_09", "ratio_6_20_07", "ratio_6_20_09", "ratio_7_20_07", "ratio_7_20_09"] # "ratio_3_20_07", "ratio_3_20_09", "ratio_4_20_07", "ratio_4_20_09", 
# list_threshold = ["threshold_4", "threshold_5", "threshold_6"]
# list_resonance = ["resonance_4", "resonance_5", gen_resonance_6]
# list_demicube = ["demicube_5", "demicube_6", "demicube_7"]
# list_crosspoly = ["crosspoly_6", "crosspoly_7", "crosspoly_8", "crosspoly_9", "crosspoly_10", "crosspoly_11", "crosspoly_12", "crosspoly_13"]

# # list with all the instances
# list_total = [list_rand ; list_2d ; list_srand ; list_perm ; list_ratio ; list_threshold ; list_resonance ; list_demicube ; list_crosspoly]
# # list with less of the instances
# list_short = ["rand_5_10", "rand_6_12", "2d_5", "2d_6", "srand_8_20_2", "perm_5", "perm_6", "ratio_5_20_07", "ratio_5_20_09", "threshold_4", "resonance_4", "demicube_5", "crosspoly_6", "crosspoly_7"]


function tables_affine(list_instances)

    ### For the algorithms 6 and 7, where all the stem vectors are computed, the code "chooses" what seems best for the svsprod and wechelon options, based on heuristics.

    options_0 = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)    # RC
    # options_1 = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_2 = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_3 = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)     # P
    # options_4 = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_5 = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)     # PD
    # options_6 = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_7 = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false)   # D

    options_8 = Options(8, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)    # RC-C
    # options_9 = Options(9, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_10 = Options(10, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_11 = Options(11, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)   # P-C
    # options_12 = Options(12, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_13 = Options(13, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)   # PD-C
    # options_14 = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_15 = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false) # D-C

    options_0s = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)    # RC-S
    # options_1s = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    # options_2s = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_3s = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)     # P-S
    # options_4s = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_5s = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, true, 100000*eps(), 1000*eps(), false, false, true)     # PD-S
    # options_6s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_7s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, false, false)   # D-S


    DataFrames_times = DataFrame(Problems = Any[], RC = Float64[], RC_comp = Float64[], RC_ratio_comp = Float64[],
                                    P = Float64[], P_ratio = Float64[], P_comp = Float64[], P_ratio_comp = Float64[], 
                                    PD = Float64[], PD_ratio = Float64[], PD_comp = Float64[], PD_ratio_comp = Float64[], 
                                    D = Float64[], D_ratio = Float64[], D_comp = Float64[], D_ratio_comp = Float64[])

    DataFrames_times_sym = DataFrame(Problems = Any[], RC = Float64[], RC_comp = Float64[], RC_ratio_comp = Float64[],
                                    P = Float64[], P_ratio = Float64[], P_comp = Float64[], P_ratio_comp = Float64[], 
                                    PD = Float64[], PD_ratio = Float64[], PD_comp = Float64[], PD_ratio_comp = Float64[], 
                                    D = Float64[], D_ratio = Float64[], D_comp = Float64[], D_ratio_comp = Float64[],
                                    RCs = Float64[], RCs_ratio = Float64[], Ps = Float64[], Ps_ratio = Float64[], 
                                    PDs = Float64[], PDs_ratio = Float64[], Ds = Float64[], Ds_ratio = Float64[])

    DataFrames_tableA2 = DataFrame(Problems = Any[], RC_C = Float64[], ratio_RC_RCC = Float64[], ratio_RC_RCC_ = Float64[], P_C = Float64[], ratio_P_PC = Float64[], ratio_RC_PC = Float64[], 
                                            PD_C = Float64[], ratio_PD_PDC = Float64[], ratio_RC_PDC = Float64[], D_C = Float64[], ratio_D_DC = Float64[], ratio_RC_DC = Float64[])

    DataFrames_table72 = DataFrame(Problems = Any[], n = Int[], p = Int[], Max_Circuits = Int[], Sym_stems = Int[], Asym_stems = Int[], Aug_stems = Int[], Aug_stems_bound = Int[], Chambers = Int[], Chambers_bound = Int[])

    DataFrames_tableA3 = DataFrame(Problems = Any[], RC = Float64[], RC_C = Float64[], ratio_RC_RCC = Float64[], RC_S = Float64[], ratio_RC_RCS = Float64[], 
                                    Primal = Float64[], ratio_RC_P = Float64[], Primal_C = Float64[], ratio_RC_P_C = Float64[], Primal_S = Float64[], ratio_RC_P_S = Float64[],
                                    PD = Float64[], ratio_RC_PD = Float64[], PD_C = Float64[], ratio_RC_PD_C = Float64[], PD_S = Float64[], ratio_RC_PD_S = Float64[],
                                    Dual = Float64[], ratio_RC_D = Float64[], Dual_C = Float64[], ratio_RC_D_C = Float64[], Dual_S = Float64[], ratio_RC_D_S = Float64[])

    for DF_index in eachindex(list_instances)

        name = list_instances[DF_index]

        if name in [list_threshold ; list_resonance ; list_crosspoly ; list_demicube]
            (matrice, sym) = name(false)    # retrieve data, by adding zeros for the experiments below
            sym = true
        else
            (matrice, sym) = name()         # retrieve data
        end

        # if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros"
        #     fullname = "linear_data/" * name * ".txt"
        #     matrice = readdlm(fullname)
        #     n, p = size(matrice)
        #     matrice = [matrice ; zeros(Int,1,p)]
        # else
        #     fullname = "affine_data/" * name * ".txt"
        #     matrice = readdlm(fullname)
        # end

        # if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros" || name[1:4] == "perm"
        #     sym = true
        # else
        #     sym = false
        # end

        println(name,"\n") # prints a function name 

        nplus, p = size(matrice) # here all the matrices are with a right-hand side \tau

        if true
            info_0 = isf(matrice, options_0)
            count_0 = 1

            time_0, FLP_0, ILP_0, time_LP_0, avg_LP_0, prop_LP_0, syms_0, asyms_0, dupli_syms_0, dupli_asyms_0, time_sv_0, prop_sv_0, checks_0, detect_ratio_0, time_cover_0, avg_cover_0, prop_cover_0 = isf_benchmark_values(info_0, options_0)

            while count_0 < 3
                count_0 += 1
                info_0 = isf(matrice, options_0)
                time_0       += info_0.cput_total
                time_LP_0    += info_0.cput_lp
                prop_LP_0    += (info_0.cput_lp / info_0.cput_total)
                time_sv_0    += info_0.cput_sv
                prop_sv_0    += (info_0.cput_sv / info_0.cput_total)
                time_cover_0 += info_0.cput_cover
                prop_cover_0 += (info_0.cput_cover / info_0.cput_total)
            end
            time_0       /= count_0
            time_LP_0    /= count_0
            prop_LP_0    /= count_0
            time_sv_0    /= count_0
            prop_sv_0    /= count_0
            time_cover_0 /= count_0
            prop_cover_0 /= count_0
        end

        ##################################################
        ##################################################

        if true
            info_3 = isf(matrice, options_3)
            count_3 = 1
            
            time_3, FLP_3, ILP_3, time_LP_3, avg_LP_3, prop_LP_3, syms_3, asyms_3, dupli_syms_3, dupli_asyms_3, time_sv_3, prop_sv_3, checks_3, detect_ratio_3, time_cover_3, avg_cover_3, prop_cover_3 = isf_benchmark_values(info_3, options_3)
            
            while count_3 < 3
                count_3 += 1
                info_3 = isf(matrice, options_3)
                time_3       += info_3.cput_total
                time_LP_3    += info_3.cput_lp
                prop_LP_3    += (info_3.cput_lp / info_3.cput_total)
                time_sv_3    += info_3.cput_sv
                prop_sv_3    += (info_3.cput_sv / info_3.cput_total)
                time_cover_3 += info_3.cput_cover
                prop_cover_3 += (info_3.cput_cover / info_3.cput_total)
            end
            time_3       /= count_3
            time_LP_3    /= count_3
            prop_LP_3    /= count_3
            time_sv_3    /= count_3
            prop_sv_3    /= count_3
            time_cover_3 /= count_3
            prop_cover_3 /= count_3
        end
        
        ##################################################
        ##################################################

        if true
            info_5 = isf(matrice, options_5)
            count_5 = 1
            
            time_5, FLP_5, ILP_5, time_LP_5, avg_LP_5, prop_LP_5, syms_5, asyms_5, dupli_syms_5, dupli_asyms_5, time_sv_5, prop_sv_5, checks_5, detect_ratio_5, time_cover_5, avg_cover_5, prop_cover_5 = isf_benchmark_values(info_5, options_5)
            
            while count_5 < 3
                count_5 += 1
                info_5 = isf(matrice, options_5)
                time_5       += info_5.cput_total
                time_LP_5    += info_5.cput_lp
                prop_LP_5    += (info_5.cput_lp / info_5.cput_total)
                time_sv_5    += info_5.cput_sv
                prop_sv_5    += (info_5.cput_sv / info_5.cput_total)
                time_cover_5 += info_5.cput_cover
                prop_cover_5 += (info_5.cput_cover / info_5.cput_total)
            end
            time_5       /= count_5
            time_LP_5    /= count_5
            prop_LP_5    /= count_5
            time_sv_5    /= count_5
            prop_sv_5    /= count_5
            time_cover_5 /= count_5
            prop_cover_5 /= count_5
        end

        ##################################################
        ##################################################

        if true
            if name == gen_resonance_6
                count_7 = 50
                info_7 = isf_get_info()
                time_7, FLP_7, ILP_7, time_LP_7, avg_LP_7, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, avg_cover_7, prop_cover_7 = 10^8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            else
                info_7 = isf(matrice, options_7)
                count_7 = 1
            
                time_7, FLP_7, ILP_7, time_LP_7, avg_LP_7, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, avg_cover_7, prop_cover_7 = isf_benchmark_values(info_7, options_7)
            end
            
            while count_7 < 3
                count_7 += 1
                info_7 = isf(matrice, options_7)
                time_7       += info_7.cput_total
                time_LP_7    += info_7.cput_lp
                prop_LP_7    += (info_7.cput_lp / info_7.cput_total)
                time_sv_7    += info_7.cput_sv
                prop_sv_7    += (info_7.cput_sv / info_7.cput_total)
                time_cover_7 += info_7.cput_cover
                prop_cover_7 += (info_7.cput_cover / info_7.cput_total)
            end
            time_7       /= count_7
            time_LP_7    /= count_7
            prop_LP_7    /= count_7
            time_sv_7    /= count_7
            prop_sv_7    /= count_7
            time_cover_7 /= count_7
            prop_cover_7 /= count_7
        end

        ##################################################
        ##################################################

        if true
            info_8 = isf(matrice, options_8)
            count_8 = 1
            
            time_8, FLP_8, ILP_8, time_LP_8, avg_LP_8, prop_LP_8, syms_8, asyms_8, dupli_syms_8, dupli_asyms_8, time_sv_8, prop_sv_8, checks_8, detect_ratio_8, time_cover_8, avg_cover_8, prop_cover_8 = isf_benchmark_values(info_8, options_8)
            
            while count_8 < 3
                count_8 += 1
                info_8 = isf(matrice, options_8)
                time_8       += info_8.cput_total
                time_LP_8    += info_8.cput_lp
                prop_LP_8    += (info_8.cput_lp / info_8.cput_total)
                time_sv_8    += info_8.cput_sv
                prop_sv_8    += (info_8.cput_sv / info_8.cput_total)
                time_cover_8 += info_8.cput_cover
                prop_cover_8 += (info_8.cput_cover / info_8.cput_total)
            end
            time_8       /= count_8
            time_LP_8    /= count_8
            prop_LP_8    /= count_8
            time_sv_8    /= count_8
            prop_sv_8    /= count_8
            time_cover_8 /= count_8
            prop_cover_8 /= count_8
        end

        ##################################################
        ##################################################

        if true
            info_11 = isf(matrice, options_11)
            count_11 = 1
            
            time_11, FLP_11, ILP_11, time_LP_11, avg_LP_11, prop_LP_11, syms_11, asyms_11, dupli_syms_11, dupli_asyms_11, time_sv_11, prop_sv_11, checks_11, detect_ratio_11, time_cover_11, avg_cover_11, prop_cover_11 = isf_benchmark_values(info_11, options_11)
            
            while count_11 < 3
                count_11 += 1
                info_11 = isf(matrice, options_11)
                time_11       += info_11.cput_total
                time_LP_11    += info_11.cput_lp
                prop_LP_11    += (info_11.cput_lp / info_11.cput_total)
                time_sv_11    += info_11.cput_sv
                prop_sv_11    += (info_11.cput_sv / info_11.cput_total)
                time_cover_11 += info_11.cput_cover
                prop_cover_11 += (info_11.cput_cover / info_11.cput_total)
            end
            time_11       /= count_11
            time_LP_11    /= count_11
            prop_LP_11    /= count_11
            time_sv_11    /= count_11
            prop_sv_11    /= count_11
            time_cover_11 /= count_11
            prop_cover_11 /= count_11
        end

        ##################################################
        ##################################################

        if true
            info_13 = isf(matrice, options_13)
            count_13 = 1
            
            time_13, FLP_13, ILP_13, time_LP_13, avg_LP_13, prop_LP_13, syms_13, asyms_13, dupli_syms_13, dupli_asyms_13, time_sv_13, prop_sv_13, checks_13, detect_ratio_13, time_cover_13, avg_cover_13, prop_cover_13 = isf_benchmark_values(info_13, options_13)
            
            while count_13 < 3
                count_13 += 1
                info_13 = isf(matrice, options_13)
                time_13       += info_13.cput_total
                time_LP_13    += info_13.cput_lp
                prop_LP_13    += (info_13.cput_lp / info_13.cput_total)
                time_sv_13    += info_13.cput_sv
                prop_sv_13    += (info_13.cput_sv / info_13.cput_total)
                time_cover_13 += info_13.cput_cover
                prop_cover_13 += (info_13.cput_cover / info_13.cput_total)
            end
            time_13       /= count_13
            time_LP_13    /= count_13
            prop_LP_13    /= count_13
            time_sv_13    /= count_13
            prop_sv_13    /= count_13
            time_cover_13 /= count_13
            prop_cover_13 /= count_13
        end

        ##################################################
        ##################################################

        if true
            if name == gen_resonance_6
                count_15 = 50
                info_15 = isf_get_info()
                time_15, FLP_15, ILP_15, time_LP_15, avg_LP_15, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, avg_cover_15, prop_cover_15 = 10^8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            else
                info_15 = isf(matrice, options_15)
                count_15 = 1
                
                time_15, FLP_15, ILP_15, time_LP_15, avg_LP_15, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, avg_cover_15, prop_cover_15 = isf_benchmark_values(info_15, options_15)
            end
            
            while count_15 < 3
                count_15 += 1
                info_15 = isf(matrice, options_15)
                time_15       += info_15.cput_total
                time_LP_15    += info_15.cput_lp
                prop_LP_15    += (info_15.cput_lp / info_15.cput_total)
                time_sv_15    += info_15.cput_sv
                prop_sv_15    += (info_15.cput_sv / info_15.cput_total)
                time_cover_15 += info_15.cput_cover
                prop_cover_15 += (info_15.cput_cover / info_15.cput_total)
            end
            time_15       /= count_15
            time_LP_15    /= count_15
            prop_LP_15    /= count_15
            time_sv_15    /= count_15
            prop_sv_15    /= count_15
            time_cover_15 /= count_15
            prop_cover_15 /= count_15
        end

        if sym
            # to get the correct data
            matrice = matrice[1:(nplus-1),:] # remove the last line that is useless (centering for the perm and line of zeros for the others)

            if true
                info_0s = isf(matrice, options_0s)
                count_0s = 1

                time_0s, FLP_0s, ILP_0s, time_LP_0s, avg_LP_0s, prop_LP_0s, syms_0s, asyms_0s, dupli_syms_0s, dupli_asyms_0s, time_sv_0s, prop_sv_0s, checks_0s, detect_ratio_0s, time_cover_0s, avg_cover_0s, prop_cover_0s = isf_benchmark_values(info_0s, options_0s)

                while count_0s < 3
                    count_0s += 1
                    info_0s = isf(matrice, options_0s)
                    time_0s       += info_0s.cput_total
                    time_LP_0s    += info_0s.cput_lp
                    prop_LP_0s    += (info_0s.cput_lp / info_0s.cput_total)
                    time_sv_0s    += info_0s.cput_sv
                    prop_sv_0s    += (info_0s.cput_sv / info_0s.cput_total)
                    time_cover_0s += info_0s.cput_cover
                    prop_cover_0s += (info_0s.cput_cover / info_0s.cput_total)
                end
                time_0s       /= count_0s
                time_LP_0s    /= count_0s
                prop_LP_0s    /= count_0s
                time_sv_0s    /= count_0s
                prop_sv_0s    /= count_0s
                time_cover_0s /= count_0s
                prop_cover_0s /= count_0s
            end

            ##################################################
            ##################################################

            if true
                info_3s = isf(matrice, options_3s)
                count_3s = 1
                
                time_3s, FLP_3s, ILP_3s, time_LP_3s, avg_LP_3s, prop_LP_3s, syms_3s, asyms_3s, dupli_syms_3s, dupli_asyms_3s, time_sv_3s, prop_sv_3s, checks_3s, detect_ratio_3s, time_cover_3s, avg_cover_3s, prop_cover_3s = isf_benchmark_values(info_3s, options_3s)
                
                while count_3s < 3
                    count_3s += 1
                    info_3s = isf(matrice, options_3s)
                    time_3s       += info_3s.cput_total
                    time_LP_3s    += info_3s.cput_lp
                    prop_LP_3s    += (info_3s.cput_lp / info_3s.cput_total)
                    time_sv_3s    += info_3s.cput_sv
                    prop_sv_3s    += (info_3s.cput_sv / info_3s.cput_total)
                    time_cover_3s += info_3s.cput_cover
                    prop_cover_3s += (info_3s.cput_cover / info_3s.cput_total)
                end
                time_3s       /= count_3s
                time_LP_3s    /= count_3s
                prop_LP_3s    /= count_3s
                time_sv_3s    /= count_3s
                prop_sv_3s    /= count_3s
                time_cover_3s /= count_3s
                prop_cover_3s /= count_3s
            end
            
            ##################################################
            ##################################################

            if true
                info_5s = isf(matrice, options_5s)
                count_5s = 1
                
                time_5s, FLP_5s, ILP_5s, time_LP_5s, avg_LP_5s, prop_LP_5s, syms_5s, asyms_5s, dupli_syms_5s, dupli_asyms_5s, time_sv_5s, prop_sv_5s, checks_5s, detect_ratio_5s, time_cover_5s, avg_cover_5s, prop_cover_5s = isf_benchmark_values(info_5s, options_5s)
                
                while count_5s < 3
                    count_5s += 1
                    info_5s = isf(matrice, options_5s)
                    time_5s       += info_5s.cput_total
                    time_LP_5s    += info_5s.cput_lp
                    prop_LP_5s    += (info_5s.cput_lp / info_5s.cput_total)
                    time_sv_5s    += info_5s.cput_sv
                    prop_sv_5s    += (info_5s.cput_sv / info_5s.cput_total)
                    time_cover_5s += info_5s.cput_cover
                    prop_cover_5s += (info_5s.cput_cover / info_5s.cput_total)
                end
                time_5s       /= count_5s
                time_LP_5s    /= count_5s
                prop_LP_5s    /= count_5s
                time_sv_5s    /= count_5s
                prop_sv_5s    /= count_5s
                time_cover_5s /= count_5s
                prop_cover_5s /= count_5s
            end

            ##################################################
            ##################################################

            if true
                if name == gen_resonance_6
                    count_7s = 50
                    info_7s = isf_get_info()
                    time_7s, FLP_7s, ILP_7s, time_LP_7s, avg_LP_7s, prop_LP_7s, syms_7s, asyms_7s, dupli_syms_7s, dupli_asyms_7s, time_sv_7s, prop_sv_7s, checks_7s, detect_ratio_7s, time_cover_7s, avg_cover_7s, prop_cover_7s = 10^8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                else
                    info_7s = isf(matrice, options_7s)
                    count_7s = 1
                    
                    time_7s, FLP_7s, ILP_7s, time_LP_7s, avg_LP_7s, prop_LP_7s, syms_7s, asyms_7s, dupli_syms_7s, dupli_asyms_7s, time_sv_7s, prop_sv_7s, checks_7s, detect_ratio_7s, time_cover_7s, avg_cover_7s, prop_cover_7s = isf_benchmark_values(info_7s, options_7s)
                end
                
                while count_7s < 3
                    count_7s += 1
                    info_7s = isf(matrice, options_7s)
                    time_7s       += info_7s.cput_total
                    time_LP_7s    += info_7s.cput_lp
                    prop_LP_7s    += (info_7s.cput_lp / info_7s.cput_total)
                    time_sv_7s    += info_7s.cput_sv
                    prop_sv_7s    += (info_7s.cput_sv / info_7s.cput_total)
                    time_cover_7s += info_7s.cput_cover
                    prop_cover_7s += (info_7s.cput_cover / info_7s.cput_total)
                end
                time_7s       /= count_7s
                time_LP_7s    /= count_7s
                prop_LP_7s    /= count_7s
                time_sv_7s    /= count_7s
                prop_sv_7s    /= count_7s
                time_cover_7s /= count_7s
                prop_cover_7s /= count_7s
            end
        end

        ##################################################

        if sym == false

            push!(DataFrames_times, [name, time_0, time_8, time_0/time_8, time_3, time_0/time_3, time_11, time_0/time_11, 
                                    time_5, time_0/time_5, time_13, time_0/time_13, time_7, time_0/time_7, time_15, time_0/time_15])

            push!(DataFrames_tableA2, [name, time_8, time_0/time_8, time_0/time_8, time_11, time_3/time_11, time_0/time_11, time_13, time_5/time_13, time_0/time_13, time_15, time_7/time_15, time_0/time_15])
            push!(DataFrames_table72, [name, nplus-1, p, binomial(p, nplus), info_7.nb_stems_sym, info_7.nb_stems_asym, info_15.nb_stems_asym, binomial(p, nplus+1), info_0.ns, isf_noncentral_max(p, nplus-1)])

        elseif sym == true
            push!(DataFrames_times_sym, [name, time_0, time_8, time_0/time_8, time_3, time_0/time_3, time_11, time_0/time_11, 
                                        time_5, time_0/time_5, time_13, time_0/time_13, time_7, time_0/time_7, time_15, time_0/time_15,
                                        time_0s, time_0/time_0s, time_3s, time_0/time_3s, time_5s, time_0/time_5s, time_7s, time_0/time_7s])
            
            push!(DataFrames_tableA2, [name, time_8, time_0/time_8, time_0/time_8, time_11, time_3/time_11, time_0/time_11, time_13, time_5/time_13, time_0/time_13, time_15, time_7/time_15, time_0/time_15])
            push!(DataFrames_table72, [name, nplus-1, p, binomial(p, nplus), info_7.nb_stems_sym, info_7.nb_stems_asym, info_15.nb_stems_asym, binomial(p, nplus), info_0.ns, isf_central_max(p, nplus-1)])
            push!(DataFrames_tableA3, [name, time_0, time_8, time_0/time_8, time_0s, time_0/time_0s, time_3, time_0/time_3, time_11, time_0/time_11, time_3s, time_0/time_3s, 
                                        time_5, time_0/time_5, time_13, time_0/time_13, time_5s, time_0/time_5s, time_7, time_0/time_7, time_15, time_0/time_15, time_7s, time_0/time_7s])

        end

    end

    DataFrames_tableA1 = [DataFrames_times[:,[1,2,5,6,9,10,13,14]] ; DataFrames_times_sym[:,[1,2,5,6,9,10,13,14]]]

    meanP = mean(DataFrames_tableA1[:,4])
    medianP = median(DataFrames_tableA1[:,4])
    meanPD = mean(DataFrames_tableA1[:,6])
    medianPD = median(DataFrames_tableA1[:,6])
    meanD = mean(DataFrames_tableA1[:,8])
    medianD = median(DataFrames_tableA1[:,8])

    push!(DataFrames_tableA1, ["Mean", "", "", conv_nb_fig(meanP), "", conv_nb_fig(meanPD), "", conv_nb_fig(meanD)], promote=true)
    push!(DataFrames_tableA1, ["Median", "", "", conv_nb_fig(medianP), "", conv_nb_fig(medianPD), "", conv_nb_fig(medianD)], promote=true)

    for i in 1:length(list_instances) # to put 3 numbers in the values instead of Float
        for j in 2:8
            DataFrames_tableA1[i,j] = conv_nb_fig(DataFrames_tableA1[i,j])
        end
    end

    meanRCC1 = mean(DataFrames_tableA2[:,3])
    medianRCC1 = median(DataFrames_tableA2[:,3])
    meanRCC2 = mean(DataFrames_tableA2[:,4])
    medianRCC2 = median(DataFrames_tableA2[:,4])

    meanPC1 = mean(DataFrames_tableA2[:,6])
    medianPC1 = median(DataFrames_tableA2[:,6])
    meanPC2 = mean(DataFrames_tableA2[:,7])
    medianPC2 = median(DataFrames_tableA2[:,7])

    meanPDC1 = mean(DataFrames_tableA2[:,9])
    medianPDC1 = median(DataFrames_tableA2[:,9])
    meanPDC2 = mean(DataFrames_tableA2[:,10])
    medianPDC2 = median(DataFrames_tableA2[:,10])

    meanDC1 = mean(DataFrames_tableA2[:,12])
    medianDC1 = median(DataFrames_tableA2[:,12])
    meanDC2 = mean(DataFrames_tableA2[:,13])
    medianDC2 = median(DataFrames_tableA2[:,13])

    push!(DataFrames_tableA2, ["Mean", "", conv_nb_fig(meanRCC1), conv_nb_fig(meanRCC2), "", conv_nb_fig(meanPC1), conv_nb_fig(meanPC2), "", conv_nb_fig(meanPDC1), conv_nb_fig(meanPDC2), "", conv_nb_fig(meanDC1), conv_nb_fig(meanDC2)], promote=true)
    push!(DataFrames_tableA2, ["Median", "", conv_nb_fig(medianRCC1), conv_nb_fig(medianRCC2), "", conv_nb_fig(medianPC1), conv_nb_fig(medianPC2), "", conv_nb_fig(medianPDC1), conv_nb_fig(medianPDC2), "", conv_nb_fig(medianDC1), conv_nb_fig(medianDC2)], promote=true)

    for i in 1:length(list_instances)
        for j in 2:13
            DataFrames_tableA2[i,j] = conv_nb_fig(DataFrames_tableA2[i,j])
        end
    end

    meanRCC3 = mean(DataFrames_tableA3[:,4])
    medianRCC3 = median(DataFrames_tableA3[:,4])
    meanRCS = mean(DataFrames_tableA3[:,6])
    medianRCS = median(DataFrames_tableA3[:,6])

    meanRCP = mean(DataFrames_tableA3[:,8])
    medianRCP = median(DataFrames_tableA3[:,8])
    meanRCPC = mean(DataFrames_tableA3[:,10])
    medianRCPC = median(DataFrames_tableA3[:,10])
    meanRCPS = mean(DataFrames_tableA3[:,12])
    medianRCPS = median(DataFrames_tableA3[:,12])

    meanRCPD = mean(DataFrames_tableA3[:,14])
    medianRCPD = median(DataFrames_tableA3[:,14])
    meanRCPDC = mean(DataFrames_tableA3[:,16])
    medianRCPDC = median(DataFrames_tableA3[:,16])
    meanRCPDS = mean(DataFrames_tableA3[:,18])
    medianRCPDS = median(DataFrames_tableA3[:,18])

    meanRCD = mean(DataFrames_tableA3[:,20])
    medianRCD = median(DataFrames_tableA3[:,20])
    meanRCDC = mean(DataFrames_tableA3[:,22])
    medianRCDC = median(DataFrames_tableA3[:,22])
    meanRCDS = mean(DataFrames_tableA3[:,24])
    medianRCDS = median(DataFrames_tableA3[:,24])

    push!(DataFrames_tableA3, ["Mean", "", "", conv_nb_fig(meanRCC3), "", conv_nb_fig(meanRCS), "", conv_nb_fig(meanRCP), "", conv_nb_fig(meanRCPC), "", conv_nb_fig(meanRCPS), "", conv_nb_fig(meanRCPD), "", conv_nb_fig(meanRCPDC), "", conv_nb_fig(meanRCPDS), "", conv_nb_fig(meanRCD), "", conv_nb_fig(meanRCDC), "", conv_nb_fig(meanRCDS)], promote=true)
    push!(DataFrames_tableA3, ["Median", "", "", conv_nb_fig(medianRCC3), "", conv_nb_fig(medianRCS), "", conv_nb_fig(medianRCP), "", conv_nb_fig(medianRCPC), "", conv_nb_fig(medianRCPS), "", conv_nb_fig(medianRCPD), "", conv_nb_fig(medianRCPDC), "", conv_nb_fig(medianRCPDS), "", conv_nb_fig(medianRCD), "", conv_nb_fig(medianRCDC), "", conv_nb_fig(medianRCDS)], promote=true)

    for i in 1:(size(DataFrames_tableA3,1)-2)
        for j in 2:24
            DataFrames_tableA3[i,j] = conv_nb_fig(DataFrames_tableA3[i,j])
        end
    end

    # for more detailed results
    # return DataFrames_vector, DataFrames_times, DataFrames_times_sym, DataFrames_table72, DataFrames_tableA1, DataFrames_tableA2, DataFrames_tableA3
    return DataFrames_table72, DataFrames_tableA1, DataFrames_tableA2, DataFrames_tableA3
end
