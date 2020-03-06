#=
TODO: expand test coverage
=#

using Oxygen

println("Testing Periodic Table")

@testset "Periodic Table" begin
    @test Oxygen.periodic_table["H"] == 1
    @test Oxygen.periodic_table["C"] == 6
    @test Oxygen.periodic_table["Si"] == 14
    @test Oxygen.periodic_table["Cl"] == 17
    @test Oxygen.periodic_table["Ti"] == 22
    @test Oxygen.periodic_table["Fe"] == 26
    @test Oxygen.periodic_table["Ge"] == 32
    @test Oxygen.periodic_table["Kr"] == 36
    @test Oxygen.periodic_table["Nb"] == 41
    @test Oxygen.periodic_table["Pd"] == 46
    @test Oxygen.periodic_table["Sb"] == 51
    @test Oxygen.periodic_table["Ba"] == 56
    @test Oxygen.periodic_table["Pm"] == 61
    @test Oxygen.periodic_table["Dy"] == 66
    @test Oxygen.periodic_table["Lu"] == 71
    @test Oxygen.periodic_table["Os"] == 76
    @test Oxygen.periodic_table["Tl"] == 81
    @test Oxygen.periodic_table["Rn"] == 86
    @test Oxygen.periodic_table["Pa"] == 91
    @test Oxygen.periodic_table["Cm"] == 96
    @test Oxygen.periodic_table["Md"] == 101
    @test Oxygen.periodic_table["Sg"] == 106
    @test_throws KeyError Oxygen.periodic_table["LL"]
end

@testset "Reverse Periodic Table" begin
    @test Oxygen.periodic_table_reverse[1] == "H"
    @test Oxygen.periodic_table_reverse[6] == "C"
    @test Oxygen.periodic_table_reverse[14] == "Si"
    @test Oxygen.periodic_table_reverse[17] == "Cl"
    @test Oxygen.periodic_table_reverse[22] == "Ti"
    @test Oxygen.periodic_table_reverse[26] == "Fe"
    @test Oxygen.periodic_table_reverse[32] == "Ge"
    @test Oxygen.periodic_table_reverse[36] == "Kr"
    @test Oxygen.periodic_table_reverse[41] == "Nb"
    @test Oxygen.periodic_table_reverse[46] == "Pd"
    @test Oxygen.periodic_table_reverse[51] == "Sb"
    @test Oxygen.periodic_table_reverse[56] == "Ba"
    @test Oxygen.periodic_table_reverse[61] == "Pm"
    @test Oxygen.periodic_table_reverse[66] == "Dy"
    @test Oxygen.periodic_table_reverse[71] == "Lu"
    @test Oxygen.periodic_table_reverse[76] == "Os"
    @test Oxygen.periodic_table_reverse[81] == "Tl"
    @test Oxygen.periodic_table_reverse[86] == "Rn"
    @test Oxygen.periodic_table_reverse[91] == "Pa"
    @test Oxygen.periodic_table_reverse[96] == "Cm"
    @test Oxygen.periodic_table_reverse[101] == "Md"
    @test Oxygen.periodic_table_reverse[106] == "Sg"
    @test_throws KeyError Oxygen.periodic_table_reverse[109]
end

@testset "Periodic Table Mass" begin
    @test Oxygen.periodic_table_mass[1] == 1.0080
    @test Oxygen.periodic_table_mass[6] == 12.011
    @test Oxygen.periodic_table_mass[8] == 15.999
    @test Oxygen.periodic_table_mass[17] == 35.45
    @test Oxygen.periodic_table_mass[22] == 47.87
    @test Oxygen.periodic_table_mass[35] == 79.9
end

@testset "Periodic Table Electronegativity" begin
    @test Oxygen.periodic_table_electronegativity[1] == 2.2
    @test Oxygen.periodic_table_electronegativity[6] == 2.55
    @test Oxygen.periodic_table_electronegativity[8] == 3.44
    @test Oxygen.periodic_table_electronegativity[17] == 3.16
    @test Oxygen.periodic_table_electronegativity[22] == 1.54
    @test Oxygen.periodic_table_electronegativity[35] == 2.96
    @test ismissing(Oxygen.periodic_table_electronegativity[2])
end
