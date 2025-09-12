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

# types of instances
instances_rand = ["rand_2_8", "rand_4_8", "rand_4_9", "rand_5_10", "rand_4_11", "rand_6_12", "rand_5_13", "rand_7_14", "rand_7_15", "rand_8_16", "rand_9_17"]
instances_2d = ["2d_4", "2d_5", "2d_6", "2d_7", "2d_8"]
instances_srand = ["srand_8_20_2", "srand_8_20_4", "srand_8_20_6"]
instances_perm = ["perm_5", "perm_6", "perm_7", "perm_8"]
instances_ratio = ["ratio_3_20_07", "ratio_3_20_09", "ratio_4_20_07", "ratio_4_20_09", "ratio_5_20_07", "ratio_5_20_09", "ratio_6_20_07", "ratio_6_20_09", "ratio_7_20_07", "ratio_7_20_09"]
instances_threshold = ["threshold_4", "threshold_5", "threshold_6"]
instances_resonance = ["resonance_4", "resonance_5"] # , "resonance_6" # resonance_6 is excluded "a priori" since the computation of its stem vectors takes too long. 
instances_demicube = ["demicube_5", "demicube_6", "demicube_7"]
instances_crosspoly = ["crosspoly_6", "crosspoly_7", "crosspoly_8", "crosspoly_9", "crosspoly_10", "crosspoly_11", "crosspoly_12", "crosspoly_13"]

# example with all instances
# list_instances = [instances_rand ; instances_2d ; instances_srand ; instances_perm ; instances_ratio ; instances_threshold ; instances_resonance ; instances_demicube ; instances_crosspoly]
list_instances = [instances_perm ; instances_ratio ; instances_threshold ; instances_resonance ; instances_demicube ; instances_crosspoly]
# list_instances = ["R_8_2"]
# list_instances = ["threshold_3"]

# explicit possible options
options_0 = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_1 = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_2 = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_3 = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_4 = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_4r = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
options_5 = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_5r = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, false, true)    # unused, no recursive covering when very few stem vectors 
options_6 = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
options_6r = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, true)    # w often faster
options_6w = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, true)
options_6wr = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, true)
options_7 = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
options_7r = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, false)  # w often faster
options_7w = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, false)
options_7wr = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, false)

options_8 = Options(8, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_9 = Options(9, 0, false, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_10 = Options(10, 0, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_11 = Options(11, 3, true, HiGHS.Optimizer, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_12 = Options(12, 3, true, HiGHS.Optimizer, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_12r = Options(12, 3, true, HiGHS.Optimizer, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, false, true)  # unused, no recursive covering when very few stem vectors 
options_13 = Options(13, 3, true, HiGHS.Optimizer, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, false, true)
options_13r = Options(13, 3, true, HiGHS.Optimizer, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, false, true)  # unused, no recursive covering when very few stem vectors 
options_14 = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, true)  # w often faster
options_14r = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, true)  # w often faster
options_14w = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, true)
options_14wr = Options(14, 3, true, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, true)
options_15 = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false, false)# w often faster
options_15r = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false, false)# w often faster
options_15w = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false, false)
options_15wr = Options(15, 0, false, HiGHS.Optimizer, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false, false)

options_0_s = Options(0, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_1_s = Options(1, 0, false, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_2_s = Options(2, 0, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_3_s = Options(3, 3, true, HiGHS.Optimizer, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_4_s = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_4r_s = Options(4, 3, true, HiGHS.Optimizer, false, true, 1, true, true, 100000*eps(), 1000*eps(), false, true, true)    # unused, no recursive covering when very few stem vectors 
options_5_s = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, false, true, 100000*eps(), 1000*eps(), false, true, true)
options_5r_s = Options(5, 3, true, HiGHS.Optimizer, false, true, 2, true, true, 100000*eps(), 1000*eps(), false, true, true)    # unused, no recursive covering when very few stem vectors 
options_6_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, true, true)    # w often faster
options_6r_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, true, true)    # w often faster
options_6w_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, true, true)
options_6wr_s = Options(6, 3, true, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, true, true)
options_7_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, true, false)  # w often faster
options_7r_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, true, false)  # w often faster
options_7w_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, true, false)
options_7wr_s = Options(7, 0, false, HiGHS.Optimizer, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, true, false)

options_l = Options(7, 0, false, HiGHS.Optimizer, false, true, 4, true, false, 100000*eps(), 1000*eps(), true, true, false)    # no tree structure

# example with all the tree algorithms
algos=["0","1","2","3","4","4r","5","5r","6","6r","6wr","6w","7","7r","7wr","7w","8","9","10","11","12","12r","13","13r","14","14r","14wr","14w","15","15r","15wr","15w"] #,"l"
algos_s = ["0s", "1s", "2s", "3s", "4s", "4rs", "5s", "5rs", "6s", "6rs", "6ws", "6wrs", "7s", "7rs", "7ws", "7wrs"]
# Depending on how many algorithms and how many instances, the rest of the code may be very long. 

for name in list_instances

    if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros"
        fullname = "linear_data/" * name * ".txt"
        matrice = readdlm(fullname)
        n, p = size(matrice)
        matrice = [matrice ; zeros(Int,1,p)]
    else
        fullname = "affine_data/data_aff_" * name * ".txt"
        matrice = readdlm(fullname)
    end

    if name[1:4] == "thre" || name[1:4] == "reso" || name[1:4] == "demi" || name[1:4] == "cros" || name[1:4] == "perm"
        sym = true
    else
        sym = false
    end

    println(name,"\n")

    # This part computes the different information for each algorithm on a given instance. 
    if "0" in algos
        info_0 = isf(matrice, options_0)
        println("fin 0,", info_0.cput_total)
    end
    
    if "1" in algos
        info_1 = isf(matrice, options_1)
        println("fin 1,", info_1.cput_total)
    end
    
    if "2" in algos
        info_2 = isf(matrice, options_2)
        println("fin 2,", info_2.cput_total)
    end
    
    if "3" in algos
        info_3 = isf(matrice, options_3)
        println("fin 3,", info_3.cput_total)
    end
    
    if "4" in algos
        info_4 = isf(matrice, options_4)
        println("fin 4,", info_4.cput_total)
    end
    
    if "4r" in algos
        info_4r = isf(matrice, options_4r)
        println("fin 4r,", info_4r.cput_total)
    end
    
    if "5" in algos
        info_5 = isf(matrice, options_5)
        println("fin 5,", info_5.cput_total)
    end
    
    if "5r" in algos
        info_5r = isf(matrice, options_5r)
        println("fin 5r,", info_5r.cput_total)
    end
    
    if "6" in algos
        info_6 = isf(matrice, options_6)
        println("fin 6,", info_6.cput_total)
    end
    
    if "6r" in algos
        info_6r = isf(matrice, options_6r)
        println("fin 6r,", info_6r.cput_total)
    end
    
    if "6wr" in algos
        info_6wr = isf(matrice, options_6wr)
        println("fin 6wr,", info_6wr.cput_total)
    end
    
    if "6w" in algos
        info_6w = isf(matrice, options_6w)
        println("fin 6w,", info_6w.cput_total)
    end
    
    if "7" in algos
        info_7 = isf(matrice, options_7)
        println("fin 7,", info_7.cput_total)
    end
    
    if "7r" in algos
        info_7r = isf(matrice, options_7r)
        println("fin 7r,", info_7r.cput_total)
    end
    
    if "7wr" in algos
        info_7wr = isf(matrice, options_7wr)
        println("fin 7wr,", info_7wr.cput_total)
    end
    
    if "7w" in algos
        info_7w = isf(matrice, options_7w)
        println("fin 7w,", info_7w.cput_total)
    end

    if "8" in algos
        info_8 = isf(matrice, options_8)
        println("fin 8,", info_8.cput_total)
    end
    
    if "9" in algos
        info_9 = isf(matrice, options_9)
        println("fin 9,", info_9.cput_total)
    end
    
    if "10" in algos
        info_10 = isf(matrice, options_10)
        println("fin 10,", info_10.cput_total)
    end
    
    if "11" in algos
        info_11 = isf(matrice, options_11)
        println("fin 11,", info_11.cput_total)
    end
    
    if "12" in algos
        info_12 = isf(matrice, options_12)
        println("fin 12,", info_12.cput_total)
    end
    
    if "12r" in algos
        info_12r = isf(matrice, options_12r)
        println("fin 12r,", info_12r.cput_total)
    end
    
    if "13" in algos
        info_13 = isf(matrice, options_13)
        println("fin 13,", info_13.cput_total)
    end
    
    if "13r" in algos
        info_13r = isf(matrice, options_13r)
        println("fin 13r,", info_13r.cput_total)
    end
    
    if "14" in algos
        info_14 = isf(matrice, options_14)
        println("fin 14,", info_14.cput_total)
    end
    
    if "14r" in algos
        info_14r = isf(matrice, options_14r)
        println("fin 14r,", info_14r.cput_total)
    end
    
    if "14wr" in algos
        info_14wr = isf(matrice, options_14wr)
        println("fin 14wr,", info_14wr.cput_total)
    end
    
    if "14w" in algos
        info_14w = isf(matrice, options_14w)
        println("fin 14w,", info_14w.cput_total)
    end
    
    if "15" in algos
        info_15 = isf(matrice, options_15)
        println("fin 15,", info_15.cput_total)
    end
    
    if "15r" in algos
        info_15r = isf(matrice, options_15r)
        println("fin 15r,", info_15r.cput_total)
    end
    
    if "15wr" in algos
        info_15wr = isf(matrice, options_15wr)
        println("fin 15wr,", info_15wr.cput_total)
    end
    
    if "15w" in algos
        info_15w = isf(matrice, options_15w)
        println("fin 15w,", info_15w.cput_total)
    end

    if "l" in algos
        info_l = isf(matrice, options_l)
        println("fin l,", info_l.cput_total)
    end

    if sym
        if name[1:4] == "perm"
            fullname = "linear_data/" * "data_" * name * ".txt"
            matrice = readdlm(fullname)
        else
            matrice = matrice[1:n,:]
        end

        if "0s" in algos_s
            info_0_s = isf(matrice, options_0_s)
            println("fin 0_s,", info_0_s.cput_total)
        end
        
        if "1s" in algos_s
            info_1_s = isf(matrice, options_1_s)
            println("fin 1_s,", info_1_s.cput_total)
        end

        if "2s" in algos_s
            info_2_s = isf(matrice, options_2_s)
            println("fin 2_s,", info_2_s.cput_total)
        end

        if "3s" in algos_s
            info_3_s = isf(matrice, options_3_s)
            println("fin 3_s,", info_3_s.cput_total)
        end

        if "4s" in algos_s
            info_4_s = isf(matrice, options_4_s)
            println("fin 4_r,", info_4_s.cput_total)
        end

        if "4rs" in algos_s
            info_4r_s = isf(matrice, options_4r_s)
            println("fin 4r_s,", info_4r_s.cput_total)
        end

        if "5s" in algos_s
            info_5_s = isf(matrice, options_5_s)
            println("fin 5_s,", info_5_s.cput_total)
        end

        if "5rs" in algos_s
            info_5r_s = isf(matrice, options_5r_s)
            println("fin 5r_s,", info_5r_s.cput_total)
        end

        if "6s" in algos_s
            info_6_s = isf(matrice, options_6_s)
            println("fin 6_s,", info_6_s.cput_total)
        end

        if "6rs" in algos_s
            info_6r_s = isf(matrice, options_6r_s)
            println("fin 6r_s,", info_6r_s.cput_total)
        end

        if "6ws" in algos_s
            info_6w_s = isf(matrice, options_6w_s)
            println("fin 6w_s,", info_6w_s.cput_total)
        end

        if "6wrs" in algos_s
            info_6wr_s = isf(matrice, options_6wr_s)
            println("fin 6wr_s,", info_6wr_s.cput_total)
        end

        if "7s" in algos_s
            info_7_s = isf(matrice, options_7_s)
            println("fin 7_s,", info_7_s.cput_total)
        end

        if "7rs" in algos_s
            info_7r_s = isf(matrice, options_7r_s)
            println("fin 7r_s,", info_7r_s.cput_total)
        end

        if "7ws" in algos_s
            info_7w_s = isf(matrice, options_7w_s)
            println("fin 7w_s,", info_7w_s.cput_total)
        end

        if "7wrs" in algos_s
            info_7wr_s = isf(matrice, options_7wr_s)
            println("fin 7wr_s,", info_7wr_s.cput_total)
        end
    end
    ##################################################
    ##################################################
    ##################################################
    ##################################################

    # quick verification of the number of chambers
    println(info_0.ns, " ", info_1.ns, " ", info_2.ns, " ", info_3.ns, " ", info_4.ns, " ", info_4r.ns, " ", info_5.ns, " ", info_5r.ns, " ", info_6.ns, " ", info_6r.ns, " ", info_6wr.ns, " ", info_6w.ns, " ", info_7.ns, " ", info_7r.ns, " ", info_7wr.ns, " ", info_7w.ns, " ", info_8.ns, " ", info_9.ns, " ", info_10.ns, " ", info_11.ns, " ", info_12.ns, " ", info_12r.ns, " ", info_13.ns, " ", info_13r.ns, " ", info_14.ns, " ", info_14r.ns, " ", info_14wr.ns, " ", info_14w.ns, " ", info_15.ns, " ", info_15r.ns, " ", info_15wr.ns, " ", info_15w.ns)

    if sym
        println(info_0_s.ns, " ", info_1_s.ns, " ", info_2_s.ns, " ", info_3_s.ns, " ", info_4_s.ns, " ", info_4r_s.ns, " ", info_5_s.ns, " ", info_5r_s.ns, " ", info_6_s.ns, " ", info_6r_s.ns, " ", info_6wr_s.ns, " ", info_6w_s.ns, " ", info_7_s.ns, " ", info_7r_s.ns, " ", info_7wr_s.ns, " ", info_7w_s.ns)
    end

    # comparisons ; the order is to try comparing the most similar sets (which should accelerate the comparison)
    # to adjust if not all algorithms are to be compared

    if info_0.ns != info_1.ns
        # println(length(setdiff(info_0.s, info_1.s)))
        # println(length(setdiff(info_1.s, info_0.s)))
        println("$(name) has an issue between 0 and 1")
    else
        sd = setdiff(info_0.s, info_1.s)
        if length(sd) == 0
            println("there is equality between 0 and 1")
        else
            println(length(sd))
            println("$(name) has an issue between 0 and 1")
        end
    end

    if info_1.ns != info_2.ns
        # println(length(setdiff(info_1.s, info_2.s)))
        # println(length(setdiff(info_2.s, info_1.s)))
        println("$(name) has an issue between 1 and 2")
    else
        sd = setdiff(info_1.s, info_2.s)
        if length(sd) == 0
            println("there is equality between 1 and 2")
        else
            println(length(sd))
            println("$(name) has an issue between 1 and 2")
        end
    end

    if info_2.ns != info_3.ns
        # println(length(setdiff(info_2.s, info_3.s)))
        # println(length(setdiff(info_3.s, info_2.s)))
        println("$(name) has an issue between 2 and 3")
    else
        sd = setdiff(info_2.s, info_3.s)
        if length(sd) == 0
            println("there is equality between 2 and 3")
        else
            println(length(sd))
            println("$(name) has an issue between 2 and 3")
        end
    end

    if info_3.ns != info_4.ns
        # println(length(setdiff(info_3.s, info_4.s)))
        # println(length(setdiff(info_4.s, info_3.s)))
        println("$(name) has an issue between 3 and 4")
    else
        sd = setdiff(info_3.s, info_4.s)
        if length(sd) == 0
            println("there is equality between 3 and 4")
        else
            println(length(sd))
            println("$(name) has an issue between 3 and 4")
        end
    end

    if info_4.s != info_4r.s
        println("$(name) has an issue between algorithms 4?")
    else
        println("there is equality between 4 and 4r")
    end

    if info_4.ns != info_5.ns
        # println(length(setdiff(info_4.s, info_5.s)))
        # println(length(setdiff(info_5.s, info_4.s)))
        println("$(name) has an issue between 4 and 5")
    else
        sd = setdiff(info_4.s, info_5.s)
        if length(sd) == 0
            println("there is equality between 4 and 5")
        else
            println(length(sd))
            println("$(name) has an issue between 4 and 5")
        end
    end

    if info_5.s != info_5r.s
        println("$(name) has an issue between algorithms 5?")
    else
        println("there is equality between 5 and 5r")
    end

    if info_5.ns != info_6.ns
        # println(length(setdiff(info_5.s, info_6.s)))
        # println(length(setdiff(info_6.s, info_5.s)))
        println("$(name) has an issue between 5 and 6")
    else
        sd = setdiff(info_5.s, info_6.s)
        if length(sd) == 0
            println("there is equality between 5 and 6")
        else
            println(length(sd))
            println("$(name) has an issue between 5 and 6")
        end
    end

    if (info_6.s != info_6r.s) || (info_6.s != info_6wr.s) || (info_6.s != info_6w.s) 
        println("$(name) has an issue between algorithms 6?")
    else
        println("there is equality between 6, 6r, 6rw, 6w")
    end

    if info_6.ns != info_7.ns
        # println(length(setdiff(info_6.s, info_7.s)))
        # println(length(setdiff(info_7.s, info_6.s)))
        println("$(name) has an issue between 6 and 7")
    else
        sd = setdiff(info_6.s, info_7.s)
        if length(sd) == 0
            println("there is equality between 6 and 7")
        else
            println(length(sd))
            println("$(name) has an issue between 6 and 7")
        end
    end

    if (info_7.s != info_7r.s) || (info_7.s != info_7wr.s) || (info_7.s != info_7w.s) 
        println("$(name) has an issue between algorithms 7?")
    else
        println("there is equality between 7, 7r, 7rw, 7w")
    end

    if info_8.ns != info_9.ns
        # println(length(setdiff(info_8.s, info_9.s)))
        # println(length(setdiff(info_9.s, info_8.s)))
        println("$(name) has an issue between 8 and 9")
    else
        sd = setdiff(info_8.s, info_9.s)
        if length(sd) == 0
            println("there is equality between 8 and 9")
        else
            println(length(sd))
            println("$(name) has an issue between 8 and 9")
        end
    end

    if info_9.ns != info_10.ns
        # println(length(setdiff(info_9.s, info_10.s)))
        # println(length(setdiff(info_10.s, info_9.s)))
        println("$(name) has an issue between 9 and 10")
    else
        sd = setdiff(info_9.s, info_10.s)
        if length(sd) == 0
            println("there is equality between 9 and 10")
        else
            println(length(sd))
            println("$(name) has an issue between 9 and 10")
        end
    end

    if info_10.ns != info_11.ns
        # println(length(setdiff(info_10.s, info_11.s)))
        # println(length(setdiff(info_11.s, info_10.s)))
        println("$(name) has an issue between 10 and 11")
    else
        sd = setdiff(info_10.s, info_11.s)
        if length(sd) == 0
            println("there is equality between 10 and 11")
        else
            println(length(sd))
            println("$(name) has an issue between 10 and 11")
        end
    end

    if info_11.ns != info_12.ns
        # println(length(setdiff(info_11.s, info_12.s)))
        # println(length(setdiff(info_12.s, info_11.s)))
        println("$(name) has an issue between 11 and 12")
    else
        sd = setdiff(info_11.s, info_12.s)
        if length(sd) == 0
            println("there is equality between 11 and 12")
        else
            println(length(sd))
            println("$(name) has an issue between 11 and 12")
        end
    end

    if info_12.s != info_12r.s
        println("$(name) has an issue between algorithms 12?")
    else
        println("there is equality between 12 and 12r")
    end

    if info_12.ns != info_13.ns
        # println(length(setdiff(info_12.s, info_13.s)))
        # println(length(setdiff(info_13.s, info_12.s)))
        println("$(name) has an issue between 12 and 13")
    else
        sd = setdiff(info_12.s, info_13.s)
        if length(sd) == 0
            println("there is equality between 12 and 13")
        else
            println(length(sd))
            println("$(name) has an issue between 12 and 13")
        end
    end

    if info_13.s != info_13r.s
        println("$(name) has an issue between algorithms 13?")
    else
        println("there is equality between 13 and 13r")
    end

    if info_13.ns != info_14.ns
        # println(length(setdiff(info_13.s, info_14.s)))
        # println(length(setdiff(info_14.s, info_13.s)))
        println("$(name) has an issue between 13 and 14")
    else
        sd = setdiff(info_13.s, info_14.s)
        if length(sd) == 0
            println("there is equality between 13 and 14")
        else
            println(length(sd))
            println("$(name) has an issue between 13 and 14")
        end
    end

    if (info_14.s != info_14r.s) || (info_14.s != info_14wr.s) || (info_14.s != info_14w.s) 
        println("$(name) has an issue between algorithms 14?")
    else
        println("there is equality between 14, 14r, 14rw, 14w")
    end

    if info_14.ns != info_15.ns
        # println(length(setdiff(info_14.s, info_15.s)))
        # println(length(setdiff(info_15.s, info_14.s)))
        println("$(name) has an issue between 14 and 15")
    else
        sd = setdiff(info_14.s, info_15.s)
        if length(sd) == 0
            println("there is equality between 14 and 15")
        else
            println(length(sd))
            println("$(name) has an issue between 14 and 15")
        end
    end

    if (info_15.s != info_15r.s) || (info_15.s != info_15wr.s) || (info_15.s != info_15w.s) 
        println("$(name) has an issue between algorithms 15?")
    else
        println("there is equality between 15, 15r, 15rw, 15w")
    end

    if info_7.ns != info_15.ns
        # println(length(setdiff(info_7.s, info_15.s)))
        # println(length(setdiff(info_15.s, info_7.s)))
        println("$(name) has an issue between 7 and 15")
    else
        sd = setdiff(info_7.s, info_15.s)
        if length(sd) == 0
            println("there is equality between 7 and 15")
        else
            println(length(sd))
            println("$(name) has an issue between 7 and 15")
        end
    end

    if sym
        if info_0.ns != info_0_s.ns
            # println(length(setdiff(info_0.s, info_0_s.s)))
            # println(length(setdiff(info_0_s.s, info_0.s)))
            println("$(name) has an issue between 0 and 0s")
        else
            sd = setdiff(info_0_s.s, info_0.s)
            if length(sd) == 0
                println("there is equality between 0 and 0_s")
            else
                println(length(sd))
                println("$(name) has an issue between 0 and 0_s")
            end
        end

        if info_0_s.ns != info_1_s.ns
            # println(length(setdiff(info_0_s.s, info_1_s.s)))
            # println(length(setdiff(info_1_s.s, info_0_s.s)))
            println("$(name) has an issue between 0 and 1")
        else
            sd = setdiff(info_0_s.s, info_1_s.s)
            if length(sd) == 0
                println("there is equality between 0 and 1")
            else
                println(length(sd))
                println("$(name) has an issue between 0 and 1")
            end
        end

        if info_1_s.ns != info_2_s.ns
            # println(length(setdiff(info_1_s.s, info_2_s.s)))
            # println(length(setdiff(info_2_s.s, info_1_s.s)))
            println("$(name) has an issue between 1 and 2")
        else
            sd = setdiff(info_1_s.s, info_2_s.s)
            if length(sd) == 0
                println("there is equality between 1 and 2")
            else
                println(length(sd))
                println("$(name) has an issue between 1 and 2")
            end
        end

        if info_2_s.ns != info_3_s.ns
            # println(length(setdiff(info_2_s.s, info_3_s.s)))
            # println(length(setdiff(info_3_s.s, info_2_s.s)))
            println("$(name) has an issue between 2 and 3")
        else
            sd = setdiff(info_2_s.s, info_3_s.s)
            if length(sd) == 0
                println("there is equality between 2 and 3")
            else
                println(length(sd))
                println("$(name) has an issue between 2 and 3")
            end
        end

        if info_3_s.ns != info_4_s.ns
            # println(length(setdiff(info_3_s.s, info_4_s.s)))
            # println(length(setdiff(info_4_s.s, info_3_s.s)))
            println("$(name) has an issue between 3 and 4")
        else
            sd = setdiff(info_3_s.s, info_4_s.s)
            if length(sd) == 0
                println("there is equality between 3 and 4")
            else
                println(length(sd))
                println("$(name) has an issue between 3 and 4")
            end
        end

        if info_4_s.s != info_4r_s.s
            println("$(name) has an issue between algorithms 4?")
        else
            println("there is equality between 4 and 4r")
        end

        if info_4_s.ns != info_5_s.ns
            # println(length(setdiff(info_4_s.s, info_5_s.s)))
            # println(length(setdiff(info_5_s.s, info_4_s.s)))
            println("$(name) has an issue between 4 and 5")
        else
            sd = setdiff(info_4_s.s, info_5_s.s)
            if length(sd) == 0
                println("there is equality between 4 and 5")
            else
                println(length(sd))
                println("$(name) has an issue between 4 and 5")
            end
        end

        if info_5_s.s != info_5r_s.s
            println("$(name) has an issue between algorithms 5?")
        else
            println("there is equality between 5 and 5r")
        end

        if info_5_s.ns != info_6_s.ns
            # println(length(setdiff(info_5_s.s, info_6_s.s)))
            # println(length(setdiff(info_6_s.s, info_5_s.s)))
            println("$(name) has an issue between 5 and 6")
        else
            sd = setdiff(info_5_s.s, info_6_s.s)
            if length(sd) == 0
                println("there is equality between 5 and 6")
            else
                println(length(sd))
                println("$(name) has an issue between 5 and 6")
            end
        end

        if (info_6_s.s != info_6r_s.s) || (info_6_s.s != info_6wr_s.s) || (info_6_s.s != info_6w_s.s) 
            println("$(name) has an issue between algorithms 6?")
        else
            println("there is equality between 6, 6r, 6rw, 6w")
        end

        if info_6_s.ns != info_7_s.ns
            # println(length(setdiff(info_6_s.s, info_7_s.s)))
            # println(length(setdiff(info_7_s.s, info_6_s.s)))
            println("$(name) has an issue between 6 and 7")
        else
            sd = setdiff(info_6_s.s, info_7_s.s)
            if length(sd) == 0
                println("there is equality between 6 and 7")
            else
                println(length(sd))
                println("$(name) has an issue between 6 and 7")
            end
        end

        if (info_7_s.s != info_7r_s.s) || (info_7_s.s != info_7wr_s.s) || (info_7_s.s != info_7w_s.s) 
            println("$(name) has an issue between algorithms 7?")
        else
            println("there is equality between 7, 7r, 7rw, 7w")
        end

    end

    println("Done")
end