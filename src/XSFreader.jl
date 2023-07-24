module XSFreader
export XSFdata,get_energy
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


end
