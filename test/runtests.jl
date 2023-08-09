using XSFreader
using Test

@testset "XSFreader.jl" begin
    filename = "structure0001.xsf"
    xsf = XSFdata(filename)
    #println(xsf)
    energy = get_energy(xsf)
    println(energy)
    display(xsf)

    Rmax = 8
    ith_atom = 1
    println(xsf.haslocalinfo)
    @time R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j = get_atoms_inside_the_sphere(xsf,ith_atom,Rmax)
    xsf.haslocalinfo = false
    println(xsf.haslocalinfo)
    @time R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j = get_atoms_inside_the_sphere(xsf,ith_atom,Rmax)
    println(xsf.haslocalinfo)
    @time R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j = get_atoms_inside_the_sphere(xsf,ith_atom,Rmax)
    display(R_js)

    R_j = make_Rmatrix(R_js,200)
    println(R_j)
    # Write your tests here.
end
