module XSFreader
using LinearAlgebra
export XSFdata,get_energy,get_atoms_inside_the_sphere,make_Rmatrix,get_number,get_localRvectors
# Write your package code here.

struct LocalRdata
    R_i::Vector{Float64}
    atomkind_i::String
    index_i::Int64
    R_js::Vector{Vector{Float64}}
    atomkinds_j::Vector{String}
    indices_j::Vector{Int64}
end

mutable struct XSFdata{numatoms}
    const R::Matrix{Float64} #3 x numatoms
    const F::Matrix{Float64} #3 x numatoms
    const comments::String
    const filename::String
    const cell::Matrix{Float64}
    const kinds::Vector{String}
    haslocalinfo::Bool
    localinfo::Union{Nothing,LocalRdata}
    haslocalRvectors::Bool
    localRvectors::Union{Nothing,Matrix{Float64}}
end

function get_energy(xsf::XSFdata)
    return parse(Float64,split(xsf.comments)[5])
end

function get_number(xsf::XSFdata{num}) where num
    return num
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

function get_position(string,data)
    idata = 1
    for i =1:length(data)
        u = split(data[i])
        if length(u) != 0
            if u[1] == string
                idata = i 
                break
            end
        end
        if i == length(data)
            error("$string  is not found")
        end
    end
    return idata
end



function XSFdata(filename)
    data = readlines(filename)
    comments = data[1]
    idata = 1
    idata = get_position("CRYSTAL",data)
    #@assert data[idata] == "CRYSTAL" "Only CRYSTAL format is supported. Now $(data[3])"
    idata = get_position("PRIMVEC",data)
    #@assert data[idata] == "PRIMVEC" "Only PRIMVEC format is supported. Now $(data[4])"
    cell = zeros(3,3)    
    for k=1:3
        cell[k,:] = parse.(Float64,split(data[idata+k]))
    end
    idata = get_position("PRIMCOORD",data)
    #@assert data[idata] == "PRIMCOORD" "Only PRIMCORD format is supported. Now $(data[idata])"
    #idata = 9
    idata += 1
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
    return XSFdata{numatoms}(R,F,comments,filename,cell,kinds,false,nothing,false,nothing)
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

function make_Rmatrix(R_js::Vector{Vector{T}},natoms) where {T<:Real}
    @assert length(R_js) <= natoms "length(R_js) should be smaller than natoms"
    R_j = zeros(T,3,natoms)
    for i=1:length(R_js) 
        R_j[:,i] .= R_js[i]
    end
    return R_j
end

function make_Rmatrix(R_js::Vector{Vector{T}},atomkinds_j,natoms,kinds) where {T<:Real}
    @assert length(R_js) <= natoms "length(R_js) should be smaller than natoms"
    R_j = zeros(T,4,natoms)
    for i=1:length(R_js)
        for k=1:3   
            R_j[k,i] = R_js[i][k]
        end
        R_j[4,i] = findfirst(x -> x == atomkinds_j[i],kinds)
        #R_j[4,i] = ifelse(atomkinds_j[i]=="H",1,-1)
    end
    return R_j
end

function get_localRvectors(xsf::XSFdata{numatoms},ith_atom,Rmax,natoms,haskinds=false)
    if xsf.haslocalRvectors
    else
        if xsf.haslocalRvectors
        else
            get_atoms_inside_the_sphere(xsf,ith_atom,Rmax)
            xsf.haslocalRvectors = true
        end
        if haskinds
            R_j = make_Rmatrix(R_js,xsf.localinfo.atomkinds_j,natoms,xsf.kinds)
        else
            R_j = make_Rmatrix(R_js,natoms) 
        end
        xsf.localRvectors = R_j
        xsf.haslocalRvectors = true
    end

    return xsf.localRvectors
end

function get_atoms_inside_the_sphere(xsf::XSFdata{numatoms},ith_atom,Rmax,
    ) where {numatoms}

    if xsf.haslocalinfo
        #return xsf.localinfo.R_i,xsf.localinfo.atomkind_i,xsf.localinfo.index_i,xsf.localinfo.R_js, xsf.localinfo.atomkinds_j, xsf.localinfo.indices_j
    else
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
                                push!(R_js,R_j-R_i)
                                count += 1
                            end
    
                        end
                    end
                end
            end
        end
        xsf.localinfo = LocalRdata(R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j)
        xsf.haslocalinfo = true
        #return R_i,atomkind_i,index_i,R_js, atomkinds_j, indices_j
    end

    return xsf.localinfo.R_i,xsf.localinfo.atomkind_i,xsf.localinfo.index_i,xsf.localinfo.R_js, xsf.localinfo.atomkinds_j, xsf.localinfo.indices_j

    


end

end
