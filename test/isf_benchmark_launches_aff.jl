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

using DelimitedFiles
using Printf

list_rand = ["rand_2_8", "rand_4_8", "rand_4_9", "rand_5_10", "rand_4_11", "rand_6_12", "rand_5_13", "rand_7_14", "rand_7_15", "rand_8_16", "rand_9_17"]
list_2d = ["2d_4", "2d_5", "2d_6", "2d_7", "2d_8"]
list_srand = ["srand_8_20_2", "srand_8_20_4", "srand_8_20_6"]
list_perm = ["perm_5", "perm_6", "perm_7", "perm_8"]
list_ratio = ["ratio_3_20_07", "ratio_3_20_09", "ratio_4_20_07", "ratio_4_20_09", "ratio_5_20_07", "ratio_5_20_09", "ratio_6_20_07", "ratio_6_20_09", "ratio_7_20_07", "ratio_7_20_09"]

list_threshold = ["threshold_4", "threshold_5", "threshold_6"]
list_resonance = ["resonance_4", "resonance_5", "resonance_6"]
list_demicube = ["demicube_5", "demicube_6", "demicube_7"]
list_crosspoly = ["crosspoly_6", "crosspoly_7", "crosspoly_8", "crosspoly_9", "crosspoly_10", "crosspoly_11", "crosspoly_12", "crosspoly_13"]

list_total = [list_rand ; list_2d ; list_srand ; list_perm ; list_ratio ; list_threshold ; list_resonance ; list_demicube ; list_crosspoly]

function benchmark_affine(list_instances)

    ### For the algorithms 6 and 7, where all the stem vectors are computed, the code "chooses" what seems best for the svsprod and wechelon options, based on heuristics.

    options_0 = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_1 = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_2 = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_3 = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_4 = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_4r = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
    options_5 = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_5r = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
    options_6 = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
    # options_6r = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
    # options_6w = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, true)
    # options_6wr = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, true)
    options_7 = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
    # options_7r = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
    # options_7w = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, false)
    # options_7wr = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, false)

    options_8 = Options(8, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_9 = Options(9, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_10 = Options(10, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_11 = Options(11, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
    options_12 = Options(12, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_12r = Options(12, 3, true, HiGHS.Optimizer, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, false, true)  # unused, no recursive covering when very few stem vectors 
    options_13 = Options(13, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)
    # options_13r = Options(13, 3, true, HiGHS.Optimizer, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, false, true)  # unused, no recursive covering when very few stem vectors 
    options_14 = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)  # w often faster
    # options_14r = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, true)  # w often faster
    # options_14w = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, true)
    # options_14wr = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, true)
    options_15 = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false)# w often faster
    # options_15r = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, false)# w often faster
    # options_15w = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, false)
    # options_15wr = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, false)

    options_0_s = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_1_s = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_2_s = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_3_s = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, false, true)
    options_4_s = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, true, 100000*eps(), 1000*eps(), false, false, true)
    # options_4r_s = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, true, true, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
    options_5_s = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, true, 100000*eps(), 1000*eps(), false, false, true)
    # options_5r_s = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, true, true, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
    options_6_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
    # options_6r_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
    # options_6w_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, false, true)
    # options_6wr_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, false, true)
    options_7_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
    # options_7r_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
    # options_7w_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, false, false)
    # options_7wr_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, false, false)

    options_l = Options(7, 0, false, HiGHS.Optimizer, false, true, 4, true, false, 100000*eps(), 1000*eps(), true, false, false)    # no tree structure

    algos = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] #,"l"
    algos_s = [0,1,2,3,4,5,6,7]

    DataFrame_instance = DataFrame(total_time = Float64[], FeasLP = Int64[], InfeasLP = Int64[], ratio_LP = Float64[], LP_time = Float64[], avg_LP = Float64[], prop_time_LP = Float64[],
                                    syms = Int64[], asyms = Int64[], dupli_syms = Int64[], dupli_asyms = Int64[], sv_time = Float64[], prop_time_sv = Float64[],
                                    checks = Int64[], detect_ratio = Float64[], cover_time = Float64[], avg_cover = Float64[], prop_time_cover = Float64[])

    DataFrames_vector = [copy(DataFrame_instance) for i in 1:length(list_instances)]
    # DataFrames_vector_comp = [DataFrame_instance for i in 1:length(list_instances)]

    DataFrames_times = DataFrame(Problems = String[], RC = Float64[], RC_comp = Float64[], RC_ratio_comp = Float64[],
                                    A = Float64[], A_ratio = Float64[], A_comp = Float64[], A_ratio_comp = Float64[], 
                                    B = Float64[], B_ratio = Float64[], B_comp = Float64[], B_ratio_comp = Float64[], 
                                    P = Float64[], P_ratio = Float64[], P_comp = Float64[], P_ratio_comp = Float64[], 
                                    D1 = Float64[], D1_ratio = Float64[], D1_comp = Float64[], D1_ratio_comp = Float64[], 
                                    PD = Float64[], PD_ratio = Float64[], PD_comp = Float64[], PD_ratio_comp = Float64[], 
                                    D3 = Float64[], D3_ratio = Float64[], D3_comp = Float64[], D3_ratio_comp = Float64[], 
                                    D = Float64[], D_ratio = Float64[], D_comp = Float64[], D_ratio_comp = Float64[])

    DataFrames_times_sym = DataFrame(Problems = String[], RC = Float64[], RC_comp = Float64[], RC_ratio_comp = Float64[],
                                    A = Float64[], A_ratio = Float64[], A_comp = Float64[], A_ratio_comp = Float64[], 
                                    B = Float64[], B_ratio = Float64[], B_comp = Float64[], B_ratio_comp = Float64[], 
                                    P = Float64[], P_ratio = Float64[], P_comp = Float64[], P_ratio_comp = Float64[], 
                                    D1 = Float64[], D1_ratio = Float64[], D1_comp = Float64[], D1_ratio_comp = Float64[], 
                                    PD = Float64[], PD_ratio = Float64[], PD_comp = Float64[], PD_ratio_comp = Float64[], 
                                    D3 = Float64[], D3_ratio = Float64[], D3_comp = Float64[], D3_ratio_comp = Float64[], 
                                    D = Float64[], D_ratio = Float64[], D_comp = Float64[], D_ratio_comp = Float64[],
                                    RCs = Float64[], RCs_ratio = Float64[], As = Float64[], As_ratio = Float64[], 
                                    Bs = Float64[], Bs_ratio = Float64[], Ps = Float64[], Ps_ratio = Float64[], 
                                    D1s = Float64[], D1s_ratio = Float64[], PDs = Float64[], PDs_ratio = Float64[],
                                    D3s = Float64[], D3s_ratio = Float64[], Ds = Float64[], Ds_ratio = Float64[])

    DataFrames_tableA2 = DataFrame(Problems = String[], RC_C = Float64[], ratio_RC_RCC = Float64[], ratio_RC_RCC_ = Float64[], P_C = Float64[], ratio_P_PC = Float64[], ratio_RC_PC = Float64[], 
                                            PD_C = Float64[], ratio_PD_PDC = Float64[], ratio_RC_PDC = Float64[], D_C = Float64[], ratio_D_DC = Float64[], ratio_RC_DC = Float64[])

    DataFrames_table72 = DataFrame(Problems = String[], n = Int[], p = Int[], Max_Circuits = Int[], Sym_stems = Int[], Asym_stems = Int[], Aug_stems = Int[], Aug_stems_bound = Int[], Chambers = Int[], Chambers_bound = Int[])

    DataFrames_tableA3 = DataFrame(Problems = String[], RC = Float64[], RC_C = Float64[], ratio_RC_RCC = Float64[], RC_S = Float64[], ratio_RC_RCS = Float64[], 
                                    Primal = Float64[], ratio_RC_P = Float64[], Primal_C = Float64[], ratio_RC_P_C = Float64[], Primal_S = Float64[], ratio_RC_P_S = Float64[],
                                    PD = Float64[], ratio_RC_PD = Float64[], PD_C = Float64[], ratio_RC_PD_C = Float64[], PD_S = Float64[], ratio_RC_PD_S = Float64[],
                                    Dual = Float64[], ratio_RC_D = Float64[], Dual_C = Float64[], ratio_RC_D_C = Float64[], Dual_S = Float64[], ratio_RC_D_S = Float64[])

    for DF_index in eachindex(list_instances)

        name = list_instances[DF_index]
        
        if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros"
            fullname = "linear_data/" * name * ".txt"
            matrice = readdlm(fullname)
            n, p = size(matrice)
            matrice = [matrice ; zeros(Int,1,p)]
        else
            fullname = "affine_data/" * name * ".txt"
            matrice = readdlm(fullname)
        end

        if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros" || name[1:4] == "perm"
            sym = true
        else
            sym = false
        end

        println(name,"\n")

        nplus, p = size(matrice) # here all the matrices are with a right-hand side \tau

        if 0 in algos
            info_0 = isf(matrice, options_0)
            count_0 = 1

            time_0, FLP_0, ILP_0, time_LP_0, avg_LP_0, prop_LP_0, syms_0, asyms_0, dupli_syms_0, dupli_asyms_0, time_sv_0, prop_sv_0, checks_0, detect_ratio_0, time_cover_0, avg_cover_0, prop_cover_0 = isf_benchmark_values(info_0, options_0)

            while count_0 < 3 || (time_0 < 50 && count_0 < 50)
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

            push!(DataFrames_vector[DF_index], [time_0, FLP_0, ILP_0, 1.0, time_LP_0, time_LP_0 / (FLP_0 + ILP_0), prop_LP_0, syms_0, asyms_0, dupli_syms_0, dupli_asyms_0, time_sv_0, prop_sv_0, checks_0, detect_ratio_0, time_cover_0, 0, prop_cover_0])
        end

        ##################################################
        ##################################################

        if 1 in algos
            info_1 = isf(matrice, options_1)
            count_1 = 1

            time_1, FLP_1, ILP_1, time_LP_1, avg_LP_1, prop_LP_1, syms_1, asyms_1, dupli_syms_1, dupli_asyms_1, time_sv_1, prop_sv_1, checks_1, detect_ratio_1, time_cover_1, avg_cover_1, prop_cover_1 = isf_benchmark_values(info_1, options_1)

            while count_1 < 3 || (time_1 < 50 && count_1 < 50)
                count_1 += 1
                info_1 = isf(matrice, options_1)
                time_1       += info_1.cput_total
                time_LP_1    += info_1.cput_lp
                prop_LP_1    += (info_1.cput_lp / info_1.cput_total)
                time_sv_1    += info_1.cput_sv
                prop_sv_1    += (info_1.cput_sv / info_1.cput_total)
                time_cover_1 += info_1.cput_cover
                prop_cover_1 += (info_1.cput_cover / info_1.cput_total)
            end
            time_1       /= count_1
            time_LP_1    /= count_1
            prop_LP_1    /= count_1
            time_sv_1    /= count_1
            prop_sv_1    /= count_1
            time_cover_1 /= count_1
            prop_cover_1 /= count_1

            push!(DataFrames_vector[DF_index], [time_1, FLP_1, ILP_1, (FLP_0+ILP_0)/(FLP_1+ILP_1) , time_LP_1, time_LP_1 / (FLP_1 + ILP_1), prop_LP_1, syms_1, asyms_1, dupli_syms_1, dupli_asyms_1, time_sv_1, prop_sv_1, checks_1, detect_ratio_1, time_cover_1, 0, prop_cover_1])
        end

        ##################################################
        ##################################################

        if 2 in algos
            info_2 = isf(matrice, options_2)
            count_2 = 1
            
            time_2, FLP_2, ILP_2, time_LP_2, avg_LP_2, prop_LP_2, syms_2, asyms_2, dupli_syms_2, dupli_asyms_2, time_sv_2, prop_sv_2, checks_2, detect_ratio_2, time_cover_2, avg_cover_2, prop_cover_2 = isf_benchmark_values(info_2, options_2)
            
            while count_2 < 3 || (time_2 < 50 && count_2 < 50)
                count_2 += 1
                info_2 = isf(matrice, options_2)
                time_2       += info_2.cput_total
                time_LP_2    += info_2.cput_lp
                prop_LP_2    += (info_2.cput_lp / info_2.cput_total)
                time_sv_2    += info_2.cput_sv
                prop_sv_2    += (info_2.cput_sv / info_2.cput_total)
                time_cover_2 += info_2.cput_cover
                prop_cover_2 += (info_2.cput_cover / info_2.cput_total)
            end
            time_2       /= count_2
            time_LP_2    /= count_2
            prop_LP_2    /= count_2
            time_sv_2    /= count_2
            prop_sv_2    /= count_2
            time_cover_2 /= count_2
            prop_cover_2 /= count_2
            
            push!(DataFrames_vector[DF_index], [time_2, FLP_2, ILP_2, (FLP_0+ILP_0)/(FLP_2+ILP_2) , time_LP_2, time_LP_2 / (FLP_2 + ILP_2), prop_LP_2, syms_2, asyms_2, dupli_syms_2, dupli_asyms_2, time_sv_2, prop_sv_2, checks_2, detect_ratio_2, time_cover_2, 0, prop_cover_2])
        end

        ##################################################
        ##################################################

        if 3 in algos
            info_3 = isf(matrice, options_3)
            count_3 = 1
            
            time_3, FLP_3, ILP_3, time_LP_3, avg_LP_3, prop_LP_3, syms_3, asyms_3, dupli_syms_3, dupli_asyms_3, time_sv_3, prop_sv_3, checks_3, detect_ratio_3, time_cover_3, avg_cover_3, prop_cover_3 = isf_benchmark_values(info_3, options_3)
            
            while count_3 < 3 || (time_3 < 50 && count_3 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_3, FLP_3, ILP_3, (FLP_0+ILP_0)/(FLP_3+ILP_3) , time_LP_3, time_LP_3 / (FLP_3 + ILP_3), prop_LP_3, syms_3, asyms_3, dupli_syms_3, dupli_asyms_3, time_sv_3, prop_sv_3, checks_3, detect_ratio_3, time_cover_3, 0, prop_cover_3])
        end

        ##################################################
        ##################################################

        if 4 in algos
            info_4 = isf(matrice, options_4)
            count_4 = 1
            
            time_4, FLP_4, ILP_4, time_LP_4, avg_LP_4, prop_LP_4, syms_4, asyms_4, dupli_syms_4, dupli_asyms_4, time_sv_4, prop_sv_4, checks_4, detect_ratio_4, time_cover_4, avg_cover_4, prop_cover_4 = isf_benchmark_values(info_4, options_4)
            
            while count_4 < 3 || (time_4 < 50 && count_4 < 50)
                count_4 += 1
                info_4 = isf(matrice, options_4)
                time_4       += info_4.cput_total
                time_LP_4    += info_4.cput_lp
                prop_LP_4    += (info_4.cput_lp / info_4.cput_total)
                time_sv_4    += info_4.cput_sv
                prop_sv_4    += (info_4.cput_sv / info_4.cput_total)
                time_cover_4 += info_4.cput_cover
                prop_cover_4 += (info_4.cput_cover / info_4.cput_total)
            end
            time_4       /= count_4
            time_LP_4    /= count_4
            prop_LP_4    /= count_4
            time_sv_4    /= count_4
            prop_sv_4    /= count_4
            time_cover_4 /= count_4
            prop_cover_4 /= count_4
            
            push!(DataFrames_vector[DF_index], [time_4, FLP_4, ILP_4, (FLP_0+ILP_0)/(FLP_4+ILP_4) , time_LP_4, time_LP_4 / (FLP_4 + ILP_4), prop_LP_4, syms_4, asyms_4, dupli_syms_4, dupli_asyms_4, time_sv_4, prop_sv_4, checks_4, detect_ratio_4, time_cover_4, time_cover_4 / checks_4, prop_cover_4])
        end
        
        ##################################################
        ##################################################

        if 5 in algos
            info_5 = isf(matrice, options_5)
            count_5 = 1
            
            time_5, FLP_5, ILP_5, time_LP_5, avg_LP_5, prop_LP_5, syms_5, asyms_5, dupli_syms_5, dupli_asyms_5, time_sv_5, prop_sv_5, checks_5, detect_ratio_5, time_cover_5, avg_cover_5, prop_cover_5 = isf_benchmark_values(info_5, options_5)
            
            while count_5 < 3 || (time_5 < 50 && count_5 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_5, FLP_5, ILP_5, (FLP_0+ILP_0)/(FLP_5+ILP_5) , time_LP_5, time_LP_5 / (FLP_5 + ILP_5), prop_LP_5, syms_5, asyms_5, dupli_syms_5, dupli_asyms_5, time_sv_5, prop_sv_5, checks_5, detect_ratio_5, time_cover_5, time_cover_5 / checks_5, prop_cover_5])
        end

        ##################################################
        ##################################################

        if 6 in algos
            if name == "resonance_6" # computing all stem vectors on this instance is not feasible
                count_6 = 10
                time_6, FLP_6, ILP_6, time_LP_6, avg_LP_6, prop_LP_6, syms_6, asyms_6, dupli_syms_6, dupli_asyms_6, time_sv_6, prop_sv_6, checks_6, detect_ratio_6, time_cover_6, avg_cover_6, prop_cover_6 = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            else
                info_6 = isf(matrice, options_6)
                count_6 = 1
                
                time_6, FLP_6, ILP_6, time_LP_6, avg_LP_6, prop_LP_6, syms_6, asyms_6, dupli_syms_6, dupli_asyms_6, time_sv_6, prop_sv_6, checks_6, detect_ratio_6, time_cover_6, avg_cover_6, prop_cover_6 = isf_benchmark_values(info_6, options_6)
            end

            while count_6 < 3 || (time_6 < 50 && count_6 < 50)
                count_6 += 1
                info_6 = isf(matrice, options_6)
                time_6       += info_6.cput_total
                time_LP_6    += info_6.cput_lp
                prop_LP_6    += (info_6.cput_lp / info_6.cput_total)
                time_sv_6    += info_6.cput_sv
                prop_sv_6    += (info_6.cput_sv / info_6.cput_total)
                time_cover_6 += info_6.cput_cover
                prop_cover_6 += (info_6.cput_cover / info_6.cput_total)
            end
            time_6       /= count_6
            time_LP_6    /= count_6
            prop_LP_6    /= count_6
            time_sv_6    /= count_6
            prop_sv_6    /= count_6
            time_cover_6 /= count_6
            prop_cover_6 /= count_6
            
            push!(DataFrames_vector[DF_index], [time_6, FLP_6, ILP_6, (FLP_0+ILP_0)/(FLP_6+ILP_6) , time_LP_6, time_LP_6 / (FLP_6 + ILP_6), prop_LP_6, syms_6, asyms_6, dupli_syms_6, dupli_asyms_6, time_sv_6, prop_sv_6, checks_6, detect_ratio_6, time_cover_6, time_cover_6 / checks_6, prop_cover_6])
        end

        ##################################################
        ##################################################

        if 7 in algos
            if name == "resonance_6"
                count_7 = 10
                time_7, FLP_7, ILP_7, time_LP_7, avg_LP_7, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, avg_cover_7, prop_cover_7 = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            else
                info_7 = isf(matrice, options_7)
                count_7 = 1
            
                time_7, FLP_7, ILP_7, time_LP_7, avg_LP_7, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, avg_cover_7, prop_cover_7 = isf_benchmark_values(info_7, options_7)
            end

            while count_7 < 3 || (time_7 < 50 && count_7 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_7, FLP_7, ILP_7, 0, time_LP_7, 0, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, time_cover_7 / checks_7, prop_cover_7])
        end

        ##################################################
        ##################################################

        if 8 in algos
            info_8 = isf(matrice, options_8)
            count_8 = 1
            
            time_8, FLP_8, ILP_8, time_LP_8, avg_LP_8, prop_LP_8, syms_8, asyms_8, dupli_syms_8, dupli_asyms_8, time_sv_8, prop_sv_8, checks_8, detect_ratio_8, time_cover_8, avg_cover_8, prop_cover_8 = isf_benchmark_values(info_8, options_8)
            
            while count_8 < 3 || (time_8 < 50 && count_8 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_8, FLP_8, ILP_8, (FLP_0+ILP_0)/(FLP_8+ILP_8) , time_LP_8, time_LP_8 / (FLP_8 + ILP_8), prop_LP_8, syms_8, asyms_8, dupli_syms_8, dupli_asyms_8, time_sv_8, prop_sv_8, checks_8, detect_ratio_8, time_cover_8, 0, prop_cover_8])
        end

        ##################################################
        ##################################################

        if 9 in algos
            info_9 = isf(matrice, options_9)
            count_9 = 1
            
            time_9, FLP_9, ILP_9, time_LP_9, avg_LP_9, prop_LP_9, syms_9, asyms_9, dupli_syms_9, dupli_asyms_9, time_sv_9, prop_sv_9, checks_9, detect_ratio_9, time_cover_9, avg_cover_9, prop_cover_9 = isf_benchmark_values(info_9, options_9)
            
            while count_9 < 3 || (time_9 < 50 && count_9 < 50)
                count_9 += 1
                info_9 = isf(matrice, options_9)
                time_9       += info_9.cput_total
                time_LP_9    += info_9.cput_lp
                prop_LP_9    += (info_9.cput_lp / info_9.cput_total)
                time_sv_9    += info_9.cput_sv
                prop_sv_9    += (info_9.cput_sv / info_9.cput_total)
                time_cover_9 += info_9.cput_cover
                prop_cover_9 += (info_9.cput_cover / info_9.cput_total)
            end
            time_9       /= count_9
            time_LP_9    /= count_9
            prop_LP_9    /= count_9
            time_sv_9    /= count_9
            prop_sv_9    /= count_9
            time_cover_9 /= count_9
            prop_cover_9 /= count_9
            
            push!(DataFrames_vector[DF_index], [time_9, FLP_9, ILP_9, (FLP_0+ILP_0)/(FLP_9+ILP_9) , time_LP_9, time_LP_9 / (FLP_9 + ILP_9), prop_LP_9, syms_9, asyms_9, dupli_syms_9, dupli_asyms_9, time_sv_9, prop_sv_9, checks_9, detect_ratio_9, time_cover_9, 0, prop_cover_9])
        end

        ##################################################
        ##################################################

        if 10 in algos
            info_10 = isf(matrice, options_10)
            count_10 = 1
            
            time_10, FLP_10, ILP_10, time_LP_10, avg_LP_10, prop_LP_10, syms_10, asyms_10, dupli_syms_10, dupli_asyms_10, time_sv_10, prop_sv_10, checks_10, detect_ratio_10, time_cover_10, avg_cover_10, prop_cover_10 = isf_benchmark_values(info_10, options_10)
            
            while count_10 < 3 || (time_10 < 50 && count_10 < 50)
                count_10 += 1
                info_10 = isf(matrice, options_10)
                time_10       += info_10.cput_total
                time_LP_10    += info_10.cput_lp
                prop_LP_10    += (info_10.cput_lp / info_10.cput_total)
                time_sv_10    += info_10.cput_sv
                prop_sv_10    += (info_10.cput_sv / info_10.cput_total)
                time_cover_10 += info_10.cput_cover
                prop_cover_10 += (info_10.cput_cover / info_10.cput_total)
            end
            time_10       /= count_10
            time_LP_10    /= count_10
            prop_LP_10    /= count_10
            time_sv_10    /= count_10
            prop_sv_10    /= count_10
            time_cover_10 /= count_10
            prop_cover_10 /= count_10
            
            push!(DataFrames_vector[DF_index], [time_10, FLP_10, ILP_10, (FLP_0+ILP_0)/(FLP_10+ILP_10) , time_LP_10, time_LP_10 / (FLP_10 + ILP_10), prop_LP_10, syms_10, asyms_10, dupli_syms_10, dupli_asyms_10, time_sv_10, prop_sv_10, checks_10, detect_ratio_10, time_cover_10, 0, prop_cover_10])
        end

        ##################################################
        ##################################################

        if 11 in algos
            info_11 = isf(matrice, options_11)
            count_11 = 1
            
            time_11, FLP_11, ILP_11, time_LP_11, avg_LP_11, prop_LP_11, syms_11, asyms_11, dupli_syms_11, dupli_asyms_11, time_sv_11, prop_sv_11, checks_11, detect_ratio_11, time_cover_11, avg_cover_11, prop_cover_11 = isf_benchmark_values(info_11, options_11)
            
            while count_11 < 3 || (time_11 < 50 && count_11 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_11, FLP_11, ILP_11, (FLP_0+ILP_0)/(FLP_11+ILP_11) , time_LP_11, time_LP_11 / (FLP_11 + ILP_11), prop_LP_11, syms_11, asyms_11, dupli_syms_11, dupli_asyms_11, time_sv_11, prop_sv_11, checks_11, detect_ratio_11, time_cover_11, 0, prop_cover_11])
        end

        ##################################################
        ##################################################

        if 12 in algos
            info_12 = isf(matrice, options_12)
            count_12 = 1
            
            time_12, FLP_12, ILP_12, time_LP_12, avg_LP_12, prop_LP_12, syms_12, asyms_12, dupli_syms_12, dupli_asyms_12, time_sv_12, prop_sv_12, checks_12, detect_ratio_12, time_cover_12, avg_cover_12, prop_cover_12 = isf_benchmark_values(info_12, options_12)
            
            while count_12 < 3 || (time_12 < 50 && count_12 < 50)
                count_12 += 1
                info_12 = isf(matrice, options_12)
                time_12       += info_12.cput_total
                time_LP_12    += info_12.cput_lp
                prop_LP_12    += (info_12.cput_lp / info_12.cput_total)
                time_sv_12    += info_12.cput_sv
                prop_sv_12    += (info_12.cput_sv / info_12.cput_total)
                time_cover_12 += info_12.cput_cover
                prop_cover_12 += (info_12.cput_cover / info_12.cput_total)
            end
            time_12       /= count_12
            time_LP_12    /= count_12
            prop_LP_12    /= count_12
            time_sv_12    /= count_12
            prop_sv_12    /= count_12
            time_cover_12 /= count_12
            prop_cover_12 /= count_12
            
            push!(DataFrames_vector[DF_index], [time_12, FLP_12, ILP_12, (FLP_0+ILP_0)/(FLP_12+ILP_12) , time_LP_12, time_LP_12 / (FLP_12 + ILP_12), prop_LP_12, syms_12, asyms_12, dupli_syms_12, dupli_asyms_12, time_sv_12, prop_sv_12, checks_12, detect_ratio_12, time_cover_12, time_cover_12 / checks_12, prop_cover_12])
        end

        ##################################################
        ##################################################

        if 13 in algos
            info_13 = isf(matrice, options_13)
            count_13 = 1
            
            time_13, FLP_13, ILP_13, time_LP_13, avg_LP_13, prop_LP_13, syms_13, asyms_13, dupli_syms_13, dupli_asyms_13, time_sv_13, prop_sv_13, checks_13, detect_ratio_13, time_cover_13, avg_cover_13, prop_cover_13 = isf_benchmark_values(info_13, options_13)
            
            while count_13 < 3 || (time_13 < 50 && count_13 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_13, FLP_13, ILP_13, (FLP_0+ILP_0)/(FLP_13+ILP_13) , time_LP_13, time_LP_13 / (FLP_13 + ILP_13), prop_LP_13, syms_13, asyms_13, dupli_syms_13, dupli_asyms_13, time_sv_13, prop_sv_13, checks_13, detect_ratio_13, time_cover_13, time_cover_13 / checks_13, prop_cover_13])
        end

        ##################################################
        ##################################################

        if 14 in algos
            if name == "resonance_6"
                count_14 = 10
                time_14, FLP_14, ILP_14, time_LP_14, avg_LP_14, prop_LP_14, syms_14, asyms_14, dupli_syms_14, dupli_asyms_14, time_sv_14, prop_sv_14, checks_14, detect_ratio_14, time_cover_14, avg_cover_14, prop_cover_14 = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0                
            else
                info_14 = isf(matrice, options_14)
                count_14 = 1
                
                time_14, FLP_14, ILP_14, time_LP_14, avg_LP_14, prop_LP_14, syms_14, asyms_14, dupli_syms_14, dupli_asyms_14, time_sv_14, prop_sv_14, checks_14, detect_ratio_14, time_cover_14, avg_cover_14, prop_cover_14 = isf_benchmark_values(info_14, options_14)
            end

            while count_14 < 3 || (time_14 < 50 && count_14 < 50)
                count_14 += 1
                info_14 = isf(matrice, options_14)
                time_14       += info_14.cput_total
                time_LP_14    += info_14.cput_lp
                prop_LP_14    += (info_14.cput_lp / info_14.cput_total)
                time_sv_14    += info_14.cput_sv
                prop_sv_14    += (info_14.cput_sv / info_14.cput_total)
                time_cover_14 += info_14.cput_cover
                prop_cover_14 += (info_14.cput_cover / info_14.cput_total)
            end
            time_14       /= count_14
            time_LP_14    /= count_14
            prop_LP_14    /= count_14
            time_sv_14    /= count_14
            prop_sv_14    /= count_14
            time_cover_14 /= count_14
            prop_cover_14 /= count_14
            
            push!(DataFrames_vector[DF_index], [time_14, FLP_14, ILP_14, (FLP_0+ILP_0)/(FLP_14+ILP_14) , time_LP_14, time_LP_14 / (FLP_14 + ILP_14), prop_LP_14, syms_14, asyms_14, dupli_syms_14, dupli_asyms_14, time_sv_14, prop_sv_14, checks_14, detect_ratio_14, time_cover_14, time_cover_14 / checks_14, prop_cover_14])
        end

        ##################################################
        ##################################################

        if 15 in algos
            if name == "resonance_6"
                count_15 = 10
                time_15, FLP_15, ILP_15, time_LP_15, avg_LP_15, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, avg_cover_15, prop_cover_15 = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            else
                info_15 = isf(matrice, options_15)
                count_15 = 1
                
                time_15, FLP_15, ILP_15, time_LP_15, avg_LP_15, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, avg_cover_15, prop_cover_15 = isf_benchmark_values(info_15, options_15)
            end

            while count_15 < 3 || (time_15 < 50 && count_15 < 50)
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
            
            push!(DataFrames_vector[DF_index], [time_15, FLP_15, ILP_15, 0, time_LP_15, 0, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, time_cover_15 / checks_15, prop_cover_15])
        end

        if sym
            # to get the correct data
            if name[1:4] == "perm" 
                fullname = "linear_data/" * "data_" * name * ".txt"
                matrice = readdlm(fullname)
            else
                matrice = matrice[1:n,:]
            end

            if 0 in algos_s
                info_0s = isf(matrice, options_0s)
                count_0s = 1

                time_0s, FLP_0s, ILP_0s, time_LP_0s, avg_LP_0s, prop_LP_0s, syms_0s, asyms_0s, dupli_syms_0s, dupli_asyms_0s, time_sv_0s, prop_sv_0s, checks_0s, detect_ratio_0s, time_cover_0s, avg_cover_0s, prop_cover_0s = isf_benchmark_values(info_0s, options_0s)

                while count_0s < 3 || (time_0s < 50 && count_0s < 50)
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

                push!(DataFrames_vector[DF_index], [time_0s, FLP_0s, ILP_0s, 1.0, time_LP_0s, time_LP_0s / (FLP_0s + ILP_0s), prop_LP_0s, syms_0s, asyms_0s, dupli_syms_0s, dupli_asyms_0s, time_sv_0s, prop_sv_0s, checks_0s, detect_ratio_0s, time_cover_0s, 0, prop_cover_0s])
            end

            ##################################################
            ##################################################

            if 1 in algos_s
                info_1s = isf(matrice, options_1s)
                count_1s = 1

                time_1s, FLP_1s, ILP_1s, time_LP_1s, avg_LP_1s, prop_LP_1s, syms_1s, asyms_1s, dupli_syms_1s, dupli_asyms_1s, time_sv_1s, prop_sv_1s, checks_1s, detect_ratio_1s, time_cover_1s, avg_cover_1s, prop_cover_1s = isf_benchmark_values(info_1s, options_1s)

                while count_1s < 3 || (time_1s < 50 && count_1s < 50)
                    count_1s += 1
                    info_1s = isf(matrice, options_1s)
                    time_1s       += info_1s.cput_total
                    time_LP_1s    += info_1s.cput_lp
                    prop_LP_1s    += (info_1s.cput_lp / info_1s.cput_total)
                    time_sv_1s    += info_1s.cput_sv
                    prop_sv_1s    += (info_1s.cput_sv / info_1s.cput_total)
                    time_cover_1s += info_1s.cput_cover
                    prop_cover_1s += (info_1s.cput_cover / info_1s.cput_total)
                end
                time_1s       /= count_1s
                time_LP_1s    /= count_1s
                prop_LP_1s    /= count_1s
                time_sv_1s    /= count_1s
                prop_sv_1s    /= count_1s
                time_cover_1s /= count_1s
                prop_cover_1s /= count_1s

                push!(DataFrames_vector[DF_index], [time_1s, FLP_1s, ILP_1s, (FLP_0s+ILP_0s)/(FLP_1s+ILP_1s) , time_LP_1s, time_LP_1s / (FLP_1s + ILP_1s), prop_LP_1s, syms_1s, asyms_1s, dupli_syms_1s, dupli_asyms_1s, time_sv_1s, prop_sv_1s, checks_1s, detect_ratio_1s, time_cover_1s, 0, prop_cover_1s])
            end

            ##################################################
            ##################################################

            if 2 in algos_s
                info_2s = isf(matrice, options_2s)
                count_2s = 1
                
                time_2s, FLP_2s, ILP_2s, time_LP_2s, avg_LP_2s, prop_LP_2s, syms_2s, asyms_2s, dupli_syms_2s, dupli_asyms_2s, time_sv_2s, prop_sv_2s, checks_2s, detect_ratio_2s, time_cover_2s, avg_cover_2s, prop_cover_2s = isf_benchmark_values(info_2s, options_2s)
                
                while count_2s < 3 || (time_2s < 50 && count_2s < 50)
                    count_2s += 1
                    info_2s = isf(matrice, options_2s)
                    time_2s       += info_2s.cput_total
                    time_LP_2s    += info_2s.cput_lp
                    prop_LP_2s    += (info_2s.cput_lp / info_2s.cput_total)
                    time_sv_2s    += info_2s.cput_sv
                    prop_sv_2s    += (info_2s.cput_sv / info_2s.cput_total)
                    time_cover_2s += info_2s.cput_cover
                    prop_cover_2s += (info_2s.cput_cover / info_2s.cput_total)
                end
                time_2s       /= count_2s
                time_LP_2s    /= count_2s
                prop_LP_2s    /= count_2s
                time_sv_2s    /= count_2s
                prop_sv_2s    /= count_2s
                time_cover_2s /= count_2s
                prop_cover_2s /= count_2s
                
                push!(DataFrames_vector[DF_index], [time_2s, FLP_2s, ILP_2s, (FLP_0s+ILP_0s)/(FLP_2s+ILP_2s) , time_LP_2s, time_LP_2s / (FLP_2s + ILP_2s), prop_LP_2s, syms_2s, asyms_2s, dupli_syms_2s, dupli_asyms_2s, time_sv_2s, prop_sv_2s, checks_2s, detect_ratio_2s, time_cover_2s, 0, prop_cover_2s])
            end

            ##################################################
            ##################################################

            if 3 in algos_s
                info_3s = isf(matrice, options_3s)
                count_3s = 1
                
                time_3s, FLP_3s, ILP_3s, time_LP_3s, avg_LP_3s, prop_LP_3s, syms_3s, asyms_3s, dupli_syms_3s, dupli_asyms_3s, time_sv_3s, prop_sv_3s, checks_3s, detect_ratio_3s, time_cover_3s, avg_cover_3s, prop_cover_3s = isf_benchmark_values(info_3s, options_3s)
                
                while count_3s < 3 || (time_3s < 50 && count_3s < 50)
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
                
                push!(DataFrames_vector[DF_index], [time_3s, FLP_3s, ILP_3s, (FLP_0s+ILP_0s)/(FLP_3s+ILP_3s) , time_LP_3s, time_LP_3s / (FLP_3s + ILP_3s), prop_LP_3s, syms_3s, asyms_3s, dupli_syms_3s, dupli_asyms_3s, time_sv_3s, prop_sv_3s, checks_3s, detect_ratio_3s, time_cover_3s, 0, prop_cover_3s])
            end

            ##################################################
            ##################################################

            if 4 in algos_s
                info_4s = isf(matrice, options_4s)
                count_4s = 1
                
                time_4s, FLP_4s, ILP_4s, time_LP_4s, avg_LP_4s, prop_LP_4s, syms_4s, asyms_4s, dupli_syms_4s, dupli_asyms_4s, time_sv_4s, prop_sv_4s, checks_4s, detect_ratio_4s, time_cover_4s, avg_cover_4s, prop_cover_4s = isf_benchmark_values(info_4s, options_4s)
                
                while count_4s < 3 || (time_4s < 50 && count_4s < 50)
                    count_4s += 1
                    info_4s = isf(matrice, options_4s)
                    time_4s       += info_4s.cput_total
                    time_LP_4s    += info_4s.cput_lp
                    prop_LP_4s    += (info_4s.cput_lp / info_4s.cput_total)
                    time_sv_4s    += info_4s.cput_sv
                    prop_sv_4s    += (info_4s.cput_sv / info_4s.cput_total)
                    time_cover_4s += info_4s.cput_cover
                    prop_cover_4s += (info_4s.cput_cover / info_4s.cput_total)
                end
                time_4s       /= count_4s
                time_LP_4s    /= count_4s
                prop_LP_4s    /= count_4s
                time_sv_4s    /= count_4s
                prop_sv_4s    /= count_4s
                time_cover_4s /= count_4s
                prop_cover_4s /= count_4s
                
                push!(DataFrames_vector[DF_index], [time_4s, FLP_4s, ILP_4s, (FLP_0s+ILP_0s)/(FLP_4s+ILP_4s) , time_LP_4s, time_LP_4s / (FLP_4s + ILP_4s), prop_LP_4s, syms_4s, asyms_4s, dupli_syms_4s, dupli_asyms_4s, time_sv_4s, prop_sv_4s, checks_4s, detect_ratio_4s, time_cover_4s, time_cover_4s / checks_4s, prop_cover_4s])
            end
            
            ##################################################
            ##################################################

            if 5 in algos_s
                info_5s = isf(matrice, options_5s)
                count_5s = 1
                
                time_5s, FLP_5s, ILP_5s, time_LP_5s, avg_LP_5s, prop_LP_5s, syms_5s, asyms_5s, dupli_syms_5s, dupli_asyms_5s, time_sv_5s, prop_sv_5s, checks_5s, detect_ratio_5s, time_cover_5s, avg_cover_5s, prop_cover_5s = isf_benchmark_values(info_5s, options_5s)
                
                while count_5s < 3 || (time_5s < 50 && count_5s < 50)
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
                
                push!(DataFrames_vector[DF_index], [time_5s, FLP_5s, ILP_5s, (FLP_0s+ILP_0s)/(FLP_5s+ILP_5s) , time_LP_5s, time_LP_5s / (FLP_5s + ILP_5s), prop_LP_5s, syms_5s, asyms_5s, dupli_syms_5s, dupli_asyms_5s, time_sv_5s, prop_sv_5s, checks_5s, detect_ratio_5s, time_cover_5s, time_cover_5s / checks_5s, prop_cover_5s])
            end

            ##################################################
            ##################################################

            if 6 in algos_s
                if name == "resonance_6"
                    count_6s = 10
                    time_6s, FLP_6s, ILP_6s, time_LP_6s, avg_LP_6s, prop_LP_6s, syms_6s, asyms_6s, dupli_syms_6s, dupli_asyms_6s, time_sv_6s, prop_sv_6s, checks_6s, detect_ratio_6s, time_cover_6s, avg_cover_6s, prop_cover_6s = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                else
                    info_6s = isf(matrice, options_6s)
                    count_6s = 1
                    
                    time_6s, FLP_6s, ILP_6s, time_LP_6s, avg_LP_6s, prop_LP_6s, syms_6s, asyms_6s, dupli_syms_6s, dupli_asyms_6s, time_sv_6s, prop_sv_6s, checks_6s, detect_ratio_6s, time_cover_6s, avg_cover_6s, prop_cover_6s = isf_benchmark_values(info_6s, options_6s)
                end
                
                while count_6s < 3 || (time_6s < 50 && count_6s < 50)
                    count_6s += 1
                    info_6s = isf(matrice, options_6s)
                    time_6s       += info_6s.cput_total
                    time_LP_6s    += info_6s.cput_lp
                    prop_LP_6s    += (info_6s.cput_lp / info_6s.cput_total)
                    time_sv_6s    += info_6s.cput_sv
                    prop_sv_6s    += (info_6s.cput_sv / info_6s.cput_total)
                    time_cover_6s += info_6s.cput_cover
                    prop_cover_6s += (info_6s.cput_cover / info_6s.cput_total)
                end
                time_6s       /= count_6s
                time_LP_6s    /= count_6s
                prop_LP_6s    /= count_6s
                time_sv_6s    /= count_6s
                prop_sv_6s    /= count_6s
                time_cover_6s /= count_6s
                prop_cover_6s /= count_6s
                
                push!(DataFrames_vector[DF_index], [time_6s, FLP_6s, ILP_6s, (FLP_0s+ILP_0s)/(FLP_6s+ILP_6s) , time_LP_6s, time_LP_6s / (FLP_6s + ILP_6s), prop_LP_6s, syms_6s, asyms_6s, dupli_syms_6s, dupli_asyms_6s, time_sv_6s, prop_sv_6s, checks_6s, detect_ratio_6s, time_cover_6s, time_cover_6s / checks_6s, prop_cover_6s])
            end

            ##################################################
            ##################################################

            if 7 in algos_s
                if name == "resonance_6"
                    count_7s = 10
                    time_7s, FLP_7s, ILP_7s, time_LP_7s, avg_LP_7s, prop_LP_7s, syms_7s, asyms_7s, dupli_syms_7s, dupli_asyms_7s, time_sv_7s, prop_sv_7s, checks_7s, detect_ratio_7s, time_cover_7s, avg_cover_7s, prop_cover_7s = 10^10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                else
                    info_7s = isf(matrice, options_7s)
                    count_7s = 1
                    
                    time_7s, FLP_7s, ILP_7s, time_LP_7s, avg_LP_7s, prop_LP_7s, syms_7s, asyms_7s, dupli_syms_7s, dupli_asyms_7s, time_sv_7s, prop_sv_7s, checks_7s, detect_ratio_7s, time_cover_7s, avg_cover_7s, prop_cover_7s = isf_benchmark_values(info_7s, options_7s)
                end
                
                while count_7s < 3 || (time_7s < 50 && count_7s < 50)
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
                
                push!(DataFrames_vector[DF_index], [time_7s, FLP_7s, ILP_7s, 0, time_LP_7s, 0, prop_LP_7s, syms_7s, asyms_7s, dupli_syms_7s, dupli_asyms_7s, time_sv_7s, prop_sv_7s, checks_7s, detect_ratio_7s, time_cover_7s, time_cover_7s / checks_7s, prop_cover_7s])
            end

        end

        ##################################################

        if algos == [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] && sym == false

            push!(DataFrames_times, [name, time_0, time_8, time_0/time_8, time_1, time_0/time_1, time_9, time_0/time_9, time_2, time_0/time_2, time_10, time_0/time_10, time_3, time_0/time_3, time_11, time_0/time_11, 
                                    time_4, time_0/time_4, time_12, time_0/time_12, time_5, time_0/time_5, time_13, time_0/time_13, time_6, time_0/time_6, time_14, time_0/time_14, time_7, time_0/time_7, time_15, time_0/time_15])
            ### for the complete set of results
            # println(" & \\multicolumn{4}{c|}{RC} & \\multicolumn{4}{c|}{A} & \\multicolumn{4}{c|}{AB} & \\multicolumn{4}{c|}{ABC}    \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp \\\\ \\hline ")
            # println("time & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_8 * "\$} & " * string_ratio_8 * " & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_9 * "\$} & " * string_ratio_9 * " & \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_10 * "\$} & " * string_ratio_10 * " & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_11 * "\$} & " * string_ratio_11 * "\\\\ \\hline")
            # println("LP & \\multicolumn{1}{c|}{\$" * string_avg_LP_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_8 * "\$} & " * string_ratio_LP_8 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_9 * "\$} & " * string_ratio_LP_9 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_10 * "\$} & " * string_ratio_LP_10 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_11 * "\$} & " * string_ratio_LP_11 * "\\\\ \\hline")
            # println("\\%LP & \\multicolumn{1}{c|}{\$" * string_ppt_LP_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_8 * "\$} & " * string_ratio_ppt_LP_8 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_9 * "\$} & " * string_ratio_ppt_LP_9 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_10 * "\$} & " * string_ratio_ppt_LP_10 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_11 * "\$} & " * string_ratio_ppt_LP_11 * "\\\\ \\hline")
            # println("sv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline ")
            # println("\\%sv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline")
            # println("cv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline")
            # println("\\%cv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\\\ \\hline ")
            # println("\n")
            # println(" & \\multicolumn{4}{c|}{ABCD1} & \\multicolumn{4}{c|}{ABCD2} & \\multicolumn{4}{c|}{ABCD3}& \\multicolumn{4}{c|}{AD4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp \\\\ \\hline ")
            # println("time & \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_12 * "\$} & " * string_ratio_12 * " & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_13 * "\$} & " * string_ratio_13 * " & \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_14 * "\$} & " * string_ratio_14 * " & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_15 * "\$} & " * string_ratio_15 * "\\\\ \\hline")
            # println("LP & \\multicolumn{1}{c|}{\$" * string_avg_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_12 * "\$} & " * string_ratio_LP_12 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_13 * "\$} & " * string_ratio_LP_13 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_14 * "\$} & " * string_ratio_LP_14 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_15 * "\$} & " * string_ratio_LP_15 * "\\\\ \\hline")
            # println("\\%LP & \\multicolumn{1}{c|}{\$" * string_ppt_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_12 * "\$} & " * string_ratio_ppt_LP_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_13 * "\$} & " * string_ratio_ppt_LP_13 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_14 * "\$} & " * string_ratio_ppt_LP_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_15 * "\$} & " * string_ratio_ppt_LP_15 * "\\\\ \\hline")
            # println("sv & \\multicolumn{1}{c|}{\$" * string_sv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_12 * "\$} & " * string_ratio_sv_12 * " & \\multicolumn{1}{c|}{\$" * string_sv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_13 * "\$} & " * string_ratio_sv_13 * "& \\multicolumn{1}{c|}{\$" * string_sv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_14 * "\$} & " * string_ratio_sv_14 * " & \\multicolumn{1}{c|}{\$" * string_sv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_15 * "\$} & " * string_ratio_sv_15 * "\\\\ \\hline")
            # println("\\%sv & \\multicolumn{1}{c|}{\$" * string_ppt_sv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_12 * "\$} & " * string_ratio_ppt_sv_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_sv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_13 * "\$} & " * string_ratio_ppt_sv_13 * "& \\multicolumn{1}{c|}{\$" * string_ppt_sv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_14 * "\$} & " * string_ratio_ppt_sv_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_sv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_15 * "\$} & " * string_ratio_ppt_sv_15 * "\\\\ \\hline")
            # println("cv & \\multicolumn{1}{c|}{\$" * string_cv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_12 * "\$} & " * string_ratio_cv_12 * " & \\multicolumn{1}{c|}{\$" * string_cv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_13 * "\$} & " * string_ratio_cv_13 * "& \\multicolumn{1}{c|}{\$" * string_cv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_14 * "\$} & " * string_ratio_cv_14 * " & \\multicolumn{1}{c|}{\$" * string_cv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_15 * "\$} & " * string_ratio_cv_15 * "\\\\ \\hline")
            # println("\\%cv & \\multicolumn{1}{c|}{\$" * string_ppt_cv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_12 * "\$} & " * string_ratio_ppt_cv_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_cv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_13 * "\$} & " * string_ratio_ppt_cv_13 * "& \\multicolumn{1}{c|}{\$" * string_ppt_cv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_14 * "\$} & " * string_ratio_ppt_cv_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_cv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_15 * "\$} & " * string_ratio_ppt_cv_15 * "\\\\ \\hline")
            ### end of the complete results
            # println("\n")

            # println("Problem & \\multicolumn{2}{c|}{RC} & \\multicolumn{2}{c|}{A} & \\multicolumn{2}{c|}{AB} & \\multicolumn{2}{c|}{ABC} & \\multicolumn{2}{c|}{D1} & \\multicolumn{2}{c|}{D2} & \\multicolumn{2}{c|}{D3} & \\multicolumn{2}{c|}{D4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & " * " " * " & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & {\\color{blue}{" * string_ratio_1* "}}& \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & {\\color{blue}{" * string_ratio_2* "}} & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & {\\color{blue}{" * string_ratio_3* "}}& \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & {\\color{blue}{" * string_ratio_4* "}} & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & {\\color{blue}{" * string_ratio_5* "}}& \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & {\\color{blue}{" * string_ratio_6* "}} & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & {\\color{blue}{" * string_ratio_7* "}}\\\\ \\hline")
            # println("\n")
            # println("Problem & \\multicolumn{2}{c|}{RC} & \\multicolumn{2}{c|}{A} & \\multicolumn{2}{c|}{AB} & \\multicolumn{2}{c|}{ABC} & \\multicolumn{2}{c|}{D1} & \\multicolumn{2}{c|}{D2} & \\multicolumn{2}{c|}{D3} & \\multicolumn{2}{c|}{D4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_8 * "\$} & {\\color{blue}{" * string_ratio_8 * "}} & \\multicolumn{1}{c|}{\$" * string_avg_9 * "\$} & {\\color{blue}{" * string_ratio_9* "}}& \\multicolumn{1}{c|}{\$" * string_avg_10 * "\$} & {\\color{blue}{" * string_ratio_10* "}} & \\multicolumn{1}{c|}{\$" * string_avg_11 * "\$} & {\\color{blue}{" * string_ratio_11* "}}& \\multicolumn{1}{c|}{\$" * string_avg_12 * "\$} & {\\color{blue}{" * string_ratio_12* "}} & \\multicolumn{1}{c|}{\$" * string_avg_13 * "\$} & {\\color{blue}{" * string_ratio_13* "}}& \\multicolumn{1}{c|}{\$" * string_avg_14 * "\$} & {\\color{blue}{" * string_ratio_14* "}} & \\multicolumn{1}{c|}{\$" * string_avg_15 * "\$} & {\\color{blue}{" * string_ratio_15* "}}\\\\ \\hline")

            # println("\n\n")

            push!(DataFrames_tableA2, [name, time_8, time_0/time_8, time_0/time_8, time_11, time_3/time_11, time_0/time_11, time_13, time_5/time_13, time_0/time_13, time_15, time_7/time_15, time_0/time_15])
            push!(DataFrames_table72, [name, nplus-1, p, binomial(p, nplus), info_7.nb_stems_sym, info_7.nb_stems_asym, info_15.nb_stems_asym, binomial(p, nplus+1), info_0.ns, isf_noncentral_max(p, nplus-1)])

        elseif algos == [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] && sym == true && algos_s == [0,1,2,3,4,5,6,7]
            push!(DataFrames_times_sym, [name, time_0, time_8, time_0/time_8, time_1, time_0/time_1, time_9, time_0/time_9, time_2, time_0/time_2, time_10, time_0/time_10, time_3, time_0/time_3, time_11, time_0/time_11, 
                                        time_4, time_0/time_4, time_12, time_0/time_12, time_5, time_0/time_5, time_13, time_0/time_13, time_6, time_0/time_6, time_14, time_0/time_14, time_7, time_0/time_7, time_15, time_0/time_15,
                                        time_0s, time_0/time_0s, time_1s, time_0/time_1s, time_2s, time_2/time_2s, time_3s, time_0/time_3s, time_4s, time_0/time_4s, time_5s, time_0/time_5s, time_6s, time_0/time_6s, time_7s, time_0/time_7s])
            
            push!(DataFrames_tableA2, [name, time_8, time_0/time_8, time_0/time_8, time_11, time_3/time_11, time_0/time_11, time_13, time_5/time_13, time_0/time_13, time_15, time_7/time_15, time_0/time_15])
            push!(DataFrames_table72, [name, nplus-1, p, binomial(p, nplus), info_7.nb_stems_sym, info_7.nb_stems_asym, info_15.nb_stems_asym, binomial(p, nplus+1), info_0.ns, isf_central_max(p, nplus-1)])
            push!(DataFrames_tableA3, [name, time_0, time_8, time_0/time_8, time_0s, time_0/time_0s, time_3, time_0/time_3, time_11, time_0/time_11, time_3s, time_0/time_3s, 
                                        time_5, time_0/time_5, time_13, time_0/time_13, time_5s, time_0/time_5s, time_7, time_0/time_7, time_15, time_0/time_15, time_7s, time_0/time_7s])

        end

    end

    DataFrames_tableA1 = [DataFrames_times[:,[1,2,13,14,21,22,29,30]] ; DataFrames_times_sym[:,[1,2,13,14,21,22,29,30]]]

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

    return DataFrames_vector, DataFrames_times, DataFrames_times_sym, DataFrames_table72, DataFrames_tableA1, DataFrames_tableA2, DataFrames_tableA3
end
