using PDEBase
using PDEBase: SupportedConditions, supported_conditions
using SciMLBase
using Test

struct UndeclaredDiscretization <: SciMLBase.AbstractDiscretization end
struct DeclaredDiscretization <: SciMLBase.AbstractDiscretization end

PDEBase.supported_conditions(::DeclaredDiscretization) = SupportedConditions(
    geometries = (:interval, :box), dirichlet = true, neumann = true, robin = true,
    periodic = true, interface = false, notes = "Rectangular grids only."
)

@testset "supported_conditions" begin
    @testset "Default declares nothing" begin
        sc = supported_conditions(UndeclaredDiscretization())
        @test sc isa SupportedConditions
        @test isempty(sc)
        @test sc.geometries == ()
        for kind in (:dirichlet, :neumann, :robin, :periodic, :interface, :time_dependent)
            @test getfield(sc, kind) === nothing
        end
        @test sc.notes == ""
        @test sc == SupportedConditions()
    end

    @testset "Override" begin
        sc = supported_conditions(DeclaredDiscretization())
        @test !isempty(sc)
        @test sc.geometries == (:interval, :box)
        @test sc.dirichlet === true
        @test sc.interface === false
        @test sc.time_dependent === nothing
        @test sc.notes == "Rectangular grids only."
        @test !isempty(SupportedConditions(notes = "Anything."))
        @test SupportedConditions(geometries = [:mesh]).geometries == (:mesh,)
    end

    @testset "Markdown table" begin
        md = repr(MIME"text/markdown"(), supported_conditions(DeclaredDiscretization()))
        @test occursin("| Geometries | interval, box |", md)
        @test occursin("| Dirichlet | yes |", md)
        @test occursin("| Interface | no |", md)
        @test occursin("| Time-dependent | unknown |", md)
        @test occursin("Rectangular grids only.", md)
        empty_md = repr(MIME"text/markdown"(), SupportedConditions())
        @test occursin("| Geometries | unknown |", empty_md)
    end
end
