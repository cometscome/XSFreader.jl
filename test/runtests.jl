using XSFreader
using Test

@testset "XSFreader.jl" begin
    filename = "structure0001.xsf"
    xsf = XSFdata(filename)
    #println(xsf)
    energy = get_energy(xsf)
    println(energy)
    display(xsf)
    # Write your tests here.
end
