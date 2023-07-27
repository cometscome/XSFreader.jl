module XSFreader
using LinearAlgebra
export XSFdata,get_energy,get_atoms_inside_the_sphere
# Write your package code here.
struct XSFdata{numatoms}
    R::Matrix{Float64} #3 x numatoms
    F::Matrix{Float64} #3 x numatoms
    comments::String
    filename::String
    cell::Matrix{Float64}
    kinds::Vector{String}
end

function get_energy(xsf::XSFdata)
    return parse(Float64,split(xsf.comments)[5])
end

function Base.display(xsf::XSFdata{numatoms}) where numatoms
    println(xsf.comments)
    println("unit cell:")
    display(xsf.cell)
    println("\t")
    println("----------------------------------------------------------------")
    println("num. of atoms: $(numatoms)")
    println("kind \t | X \t | Y \t | Z \t | F_x \t | F_y \t | F_z \t |")
    for i=1:numatoms
        R = xsf.R[:,i]
        F = xsf.F[:,i]
        println("$(xsf.kinds[i]) \t $(R[1]) \t $(R[2]) \t $(R[3]) \t $(F[1]) \t $(F[2]) \t $(F[3])")
    end
end

function XSFdata(filename)
    data = readlines(filename)
    comments = data[1]
    @assert data[3] == "CRYSTAL" "Only CRYSTAL format is supported. Now $(data[3])"
    idata = 4
    @assert data[idata] == "PRIMVEC" "Only PRIMVEC format is supported. Now $(data[4])"
    cell = zeros(3,3)    
    for k=1:3
        cell[k,:] = parse.(Float64,split(data[idata+k]))
    end
    idata = 8
    @assert data[idata] == "PRIMCOORD" "Only PRIMCORD format is supported. Now $(data[idata])"
    idata = 9
    numatoms = parse(Int64,split(data[idata])[1])
    R = zeros(3,numatoms)
    F = zeros(3,numatoms)
    kinds = Vector{String}(undef,numatoms)
    for i=1:numatoms
        idata = idata + 1
        u = split(data[idata])
        #println(u)
        R[:,i] = parse.(Float64,u[2:4])
        F[:,i] = parse.(Float64,u[5:end])
        kinds[i] = u[1]
    end
    return XSFdata{numatoms}(R,F,comments,filename,cell,kinds)
end

function cross_product!(a,b, c)
    c[1] = a[2] * b[3] - a[3] * b[2]
    c[2] = a[3] * b[1] - a[1] * b[3]
    c[3] = a[1] * b[2] - a[2] * b[1]
end

function calculate_pbcbox(xsf::XSFdata,Rmax)
    lattice_vector_a = xsf.cell[1,:]
    lattice_vector_b = xsf.cell[2,:]
    lattice_vector_c = xsf.cell[3,:]
    axb = zero(lattice_vector_a)
    bxc = zero(lattice_vector_a)
    axc = zero(lattice_vector_a)

    cross_product!(lattice_vector_a, lattice_vector_b, axb)
    axb ./= norm(axb)
    cross_product!(lattice_vector_a, lattice_vector_c, axc)
    axc ./= norm(axc)
    cross_product!(lattice_vector_b, lattice_vector_c, bxc)
    bxc ./= norm(bxc)

    project_to_a = abs(dot(lattice_vector_a, bxc))
    project_to_b = abs(dot(lattice_vector_b, axc))
    project_to_c = abs(dot(lattice_vector_c, axb))

    pbcbox_x = 0
    pbcbox_y = 0
    pbcbox_z = 0

    while pbcbox_x * project_to_a <= Rmax
        pbcbox_x += 1 
    end

    while pbcbox_y * project_to_b <= Rmax
        pbcbox_y += 1 
    end

    while pbcbox_z * project_to_c <= Rmax
        pbcbox_z += 1 
    end

    return pbcbox_x ,pbcbox_y ,pbcbox_z 

end

function get_atoms_inside_the_sphere(xsf::XSFdata{numatoms},ith_atom,Rmax,
    ) where {numatoms}
    R_i = zeros(Float64,3)
    for p=1:3
        R_i[p] = xsf.R[p,ith_atom]
    end
    atomkind_i = xsf.kinds[ith_atom]
    index_i = ith_atom
    Rmax2 = Rmax^2

    pbcbox = calculate_pbcbox(xsf,Rmax)

    count = 0
    box_min1 = -pbcbox[1]
    box_max1 = pbcbox[1]
    box_min2 = -pbcbox[2]
    box_max2 = pbcbox[2] 
    box_min3 = -pbcbox[3]
    box_max3 = pbcbox[3]

    lattice_vector_a = xsf.cell[1,:]
    lattice_vector_b = xsf.cell[2,:]
    lattice_vector_c = xsf.cell[3,:]
    R_js = Vector{Float64}[]
    atomkinds_j = String[]
    indices_j = Int64[]

    eps = 1e-4

    for j=1:numatoms
        for box_1= box_min1:box_max1
            for box_2= box_min2:box_max2
                for box_3= box_min3:box_max3
                    if (j == ith_atom && box_1 == 0 && box_2 == 0 && box_3 == 0) == false
                        R_j = xsf.R[:,j]
                        for k=1:3
                            R_j[k] += -box_1 * lattice_vector_a[k] - box_2 * lattice_vector_b[k] - box_3 * lattice_vector_c[k]
                        end

                        Rij2 = 0.0
                        for k=1:3
                            Rij2 += (R_i[k]-R_j[k])^2
                        end

                        if Rij2 < Rmax2 && Rij2 > eps
                            push!(indices_j,j)
                            push!(atomkinds_j,xsf.kinds[j])
                            push!(R_js,R_j)
                            count += 1
                        end

                    end
                end
            end
        end
    end

    return R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j


end

end
