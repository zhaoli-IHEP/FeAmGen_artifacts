import  Pkg
Pkg.activate(".")

using FeAmGen
using SymEngine
using Test

@funs   SP
@vars   shat
half    =   Basic(1) / Basic(2)

nn              =   50
mom_list_str    =   join(["K$ii" for ii ∈ 1:nn], ", ")
mass2_list_str  =   join(["m$ii" for ii ∈ 1:nn], ", ")
ver_list_str    =   join(["ver$ii" for ii ∈ 1:(3 * nn - 10)], ", ")
(eval ∘ Meta.parse)(
    "@vars $mom_list_str, $mass2_list_str, $ver_list_str"
)
mom_list        =   [Basic("K$ii") for ii ∈ 1:nn]
mass2_list      =   [Basic("m$ii") for ii ∈ 1:nn]

@testset "Test generate_kin_relation_v2" begin

    # 1 -> 1
    n_inc   =   1
    n_out   =   1
    n_tot   =   n_inc + n_out

    kin_relation    =   FeAmGen.generate_kin_relation_v2(
        n_inc, n_out,
        mom_list[1:n_tot],
        mass2_list[1:n_tot]
    )
    kin_relation_bench  =   Dict{Basic, Basic}(
        SP(K1, K2)  =>  m1
    )
    for ii ∈ 1:n_tot
        kin_relation_bench[Basic("SP(K$ii, K$ii)")] =   Basic("m$ii")
    end
    @test kin_relation == kin_relation_bench
    # 1 -> 1

    # 1 -> 2
    n_inc   =   1
    n_out   =   2
    n_tot   =   n_inc + n_out

    kin_relation    =   FeAmGen.generate_kin_relation_v2(
        n_inc, n_out,
        mom_list[1:n_tot],
        mass2_list[1:n_tot]
    )
    kin_relation_bench  =   Dict{Basic, Basic}(
        # Basic("SP(K1, K2)") =>  Basic("(1/2)*m1 + (1/2)*m2 + (-1/2)*m3"),
        # Basic("SP(K1, K3)") =>  Basic("(1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
        # Basic("SP(K2, K3)") =>  Basic("(1/2)*m1 + (-1/2)*m2 + (-1/2)*m3"),
        SP(K1, K2)  =>  half * m1 + half * m2 - half * m3,
        SP(K1, K3)  =>  half * m1 - half * m2 + half * m3,
        SP(K2, K3)  =>  half * m1 - half * m2 - half * m3
    )
    for ii ∈ 1:n_tot
        kin_relation_bench[Basic("SP(K$ii, K$ii)")] =   Basic("m$ii")
    end
    @test kin_relation == kin_relation_bench
    # 1 -> 2

    # 2 -> 1
    n_inc   =   2
    n_out   =   1
    n_tot   =   n_inc + n_out

    kin_relation    =   FeAmGen.generate_kin_relation_v2(
        n_inc, n_out,
        mom_list[1:n_tot],
        mass2_list[1:n_tot]
    )
    kin_relation_bench  =   Dict{Basic, Basic}(
        # Basic("SP(K1, K2)") =>  Basic("(-1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
        # Basic("SP(K1, K3)") =>  Basic("(1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
        # Basic("SP(K2, K3)") =>  Basic("(-1/2)*m1 + (1/2)*m2 + (1/2)*m3")
        SP(K1, K2)  =>  - half * m1 - half * m2 + half * m3,
        SP(K1, K3)  =>  half * m1 - half * m2 + half * m3,
        SP(K2, K3)  =>  - half * m1 + half * m2 + half * m3
    )
    for ii ∈ 1:n_tot
        kin_relation_bench[Basic("SP(K$ii, K$ii)")] =   Basic("m$ii")
    end
    @test kin_relation == kin_relation_bench
    # 2 -> 1

    # 2 -> 2
    n_inc   =   2
    n_out   =   2
    n_tot   =   n_inc + n_out

    kin_relation    =   FeAmGen.generate_kin_relation_v2(
        n_inc, n_out,
        mom_list[1:n_tot],
        mass2_list[1:n_tot]
    )
    kin_relation_bench  =   Dict{Basic, Basic}(
        # Basic("SP(K1, K2)") =>  Basic("(1/2)*(-m1 - m2 + shat)"),
        # Basic("SP(K1, K3)") =>  Basic("(-1/2)*(-m1 - m3 + ver1)"),
        # Basic("SP(K1, K4)") =>  Basic("(-1/2)*m2 + (-1/2)*m3 + (1/2)*shat + (1/2)*ver1"),
        # Basic("SP(K2, K3)") =>  Basic("(-1/2)*m1 + (-1/2)*m4 + (1/2)*shat + (1/2)*ver1"),
        # Basic("SP(K2, K4)") =>  Basic("(1/2)*m2 + (1/2)*m4 + (-1/2)*ver1"),
        # Basic("SP(K3, K4)") =>  Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*shat")
        SP(K1, K2)  =>  half * (- m1 - m2 + shat),
        SP(K1, K3)  =>  - half * (- m1 - m3 + ver1),
        SP(K1, K4)  =>  - half * m2 - half * m3 + half * shat + half * ver1,
        SP(K2, K3)  =>  - half * m1 - half * m4 + half * shat + half * ver1,
        SP(K2, K4)  =>  half * m2 + half * m4 - half * ver1,
        SP(K3, K4)  =>  - half * m3 - half * m4 + half * shat
    )
    for ii ∈ 1:n_tot
        kin_relation_bench[Basic("SP(K$ii, K$ii)")] =   Basic("m$ii")
    end
    @test kin_relation == kin_relation_bench
    # 2 -> 2

    # 2 -> 3
    n_inc   =   2
    n_out   =   3
    n_tot   =   n_inc + n_out

    kin_relation    =   FeAmGen.generate_kin_relation_v2(
        n_inc, n_out,
        mom_list[1:n_tot],
        mass2_list[1:n_tot]
    )
    kin_relation_bench  =   Dict{Basic, Basic}(
        # # Basic("SP(K1, K2)") =>  Basic("(1/2)*(-m1 - m2 + shat)"),
        # Basic("SP(K1, K3)") =>  Basic("(-1/2)*(-m1 - m3 + ver1)"),
        # Basic("SP(K1, K4)") =>  Basic("(-1/2)*(-m1 - m4 + ver2)"),
        # Basic("SP(K1, K5)") =>  Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*ver1 + (1/2)*ver2 + SP(K1, K2)"),
        # Basic("SP(K2, K3)") =>  Basic("(-1/2)*(-m2 - m3 + ver3)"),
        # Basic("SP(K2, K4)") =>  Basic("(-1/2)*(-m2 - m4 + ver4)"),
        # Basic("SP(K2, K5)") =>  Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*ver3 + (1/2)*ver4 + SP(K1, K2)"),
        # Basic("SP(K3, K4)") =>  Basic("(1/2)*m1 + (1/2)*m2 + (1/2)*m3 + (1/2)*m4 + (1/2)*m5 + (-1/2)*ver1 + (-1/2)*ver2 + (-1/2)*ver3 + (-1/2)*ver4 - SP(K1, K2)"),
        # Basic("SP(K3, K5)") =>  Basic("(-1/2)*m3 + (-1/2)*m4 + (-1/2)*m5 + (1/2)*ver2 + (1/2)*ver4 + SP(K1, K2)"),
        # Basic("SP(K4, K5)") =>  Basic("(-1/2)*m3 + (-1/2)*m4 + (-1/2)*m5 + (1/2)*ver1 + (1/2)*ver3 + SP(K1, K2)")
        SP(K1, K3)  =>  - half * (- m1 - m3 + ver1),
        SP(K1, K4)  =>  - half * (- m1 - m4 + ver2),
        SP(K1, K5)  =>  - half * m3 - half * m4 + half * ver1 + half * ver2 + SP(K1, K2),
        SP(K2, K3)  =>  - half * (- m2 - m3 + ver3),
        SP(K2, K4)  =>  - half * (- m2 - m4 + ver4),
        SP(K2, K5)  =>  - half * m3 - half * m4 + half * ver3 + half * ver4 + SP(K1, K2),
        SP(K3, K4)  =>  half * m1 + half * m2 + half * m3 + half * m4 + half * m5 - half * ver1 - half * ver2 - half * ver3 - half * ver4 - SP(K1, K2),
        SP(K3, K5)  =>  - half * m3 - half * m4 - half * m5 + half * ver2 + half * ver4 + SP(K1, K2),
        SP(K4, K5)  =>  - half * m3 - half * m4 - half * m5 + half * ver1 + half * ver3 + SP(K1, K2)
    )
    for ii ∈ 1:n_tot
        kin_relation_bench[Basic("SP(K$ii, K$ii)")] =   Basic("m$ii")
    end
    @test kin_relation == kin_relation_bench
    # 2 -> 3

    # 2 -> (n - 2)
    n_inc   =   2
    for n_tot ∈ 6:nn   
        n_out   =   n_tot - n_inc

        sign_list   =   [ii ≤ n_inc ? (- 1) : 1 for ii ∈ 1:n_tot]

        kin_relation_bench  =   Dict{Basic, Basic}(
            [
                Basic("SP(K$ii, K$ii)") =>  Basic("m$ii")
                for ii ∈ 1:n_tot
            ]
        )
        for ii ∈ 3:(n_tot - 1)
            kin_relation_bench[SP(K1, Basic("K$ii"))]   =   first(sign_list) * sign_list[ii] * half * Basic("- m1 - m$ii + ver$(ii - 2)")
        end
        ver_index   =   n_tot - 2
        for ii ∈ 2:(n_tot - 3)
            for jj ∈ (ii + 1):(n_tot - 1)
                kin_relation_bench[Basic("SP(K$ii, K$jj)")] =   sign_list[ii] * sign_list[jj] * half * Basic("- m$ii - m$jj + ver$(ver_index)")
                ver_index   +=  1
            end
        end
        @assert ver_index   ==  n_tot * (n_tot - 3) / 2

        for ii ∈ 1:(n_tot - 3)
            tmp =   - Basic("m$ii") - sum(
                [
                    sign_list[ii] * sign_list[jj] * (
                        (ii == jj) ? 0 : (
                            (ii < jj) ? Basic("SP(K$ii, K$jj)") : Basic("SP(K$jj, K$ii)")
                        )
                    )
                    for jj ∈ 1:(n_tot - 1)
                ]
            )
            kin_relation_bench[Basic("SP(K$ii, K$n_tot)")]  =   sign_list[ii] * sign_list[n_tot] * (
                expand ∘ subs
            )(
                tmp, kin_relation_bench
            )
        end

        tmp =   inv(
            Basic[
                1 1 0
                1 0 1
                0 1 1
            ]
        ) * [
            - Basic("m$ii") - sign_list[ii] * sum(
                [
                    sign_list[jj] * Basic("SP(K$jj, K$ii)")
                    for jj ∈ 1:(n_tot - 3)
                ]
            )
            for ii ∈ (n_tot - 2):n_tot
        ]
        tmp =   (expand ∘ subs).(tmp, Ref(kin_relation_bench))
        push!(
            kin_relation_bench,
            Pair.(
                [
                    Basic("SP(K$(n_tot - 2), K$(n_tot - 1))"),
                    Basic("SP(K$(n_tot - 2), K$n_tot)"),
                    Basic("SP(K$(n_tot - 1), K$n_tot)")
                ],
                tmp
            )...
        )
        
        kin_relation    =   FeAmGen.generate_kin_relation_v2(
            n_inc, n_out,
            mom_list[1:n_tot],
            mass2_list[1:n_tot]
        )
        
        @test   length(kin_relation) == length(kin_relation_bench)
        key_pair_list   =   union(
            [
                [
                    (
                        Basic("SP(K$ii, K$jj)"),
                        Basic("SP(K$jj, K$ii)")
                    )
                    for jj ∈ (ii + 1):n_tot
                ]
                for ii ∈ 1:n_tot
            ]...
        )
        deleteat!(key_pair_list, 1)
        for key_pair ∈ key_pair_list
            @test   xor(
                first(key_pair) ∈ keys(kin_relation),
                last(key_pair) ∈ keys(kin_relation)
            )
            @test   xor(
                first(key_pair) ∈ keys(kin_relation_bench),
                last(key_pair) ∈ keys(kin_relation_bench)
            )
        end
        @test   prod(
            [
                try
                    (iszero ∘ expand)(
                        kin_relation[key] - kin_relation_bench[key]
                    )
                catch KeyError
                    key_index   =   findfirst(
                        key_pair -> first(key_pair) == key,
                        key_pair_list
                    )
                    new_key     =   last(key_pair_list[key_index])

                    (iszero ∘ expand)(
                        kin_relation[new_key] - kin_relation_bench[key]
                    )
                end
                for key ∈ keys(kin_relation_bench)
            ]
        )
    end
end # @testset

#mom_list = [ K1, K2, K3, K4 ]
#mass2_list = [ 0, 0, mt2, mt2 ]

#mom_list = [ K1, K2 ]
#mass2_list = [ mt2, mt2 ]
