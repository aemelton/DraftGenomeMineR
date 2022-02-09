using Base: String
using CSV, DataFrames, LinearAlgebra, Plots, Statistics, Distances, Clustering,Tables 

begin
    file = "Data/Sagebrush_Promoter/promoter_catCount_by_Scaff.csv"
    df = DataFrame(CSV.File(file; transpose=true, normalizenames=true))
end

rename!(df, :water_stress => :σ_water, 
    :temperature_stress => :σ_temp, :stress_hormone => :σ_hormone,:other_stress => :σ_other);

begin
    split_names = split.(df[:,1],"_");
    number_of_scaffolds = length(split_names)

    scaffold_names = []
    gene_names = []
    gene_type = []
    gene_category = []

    for i ∈ 1:number_of_scaffolds
        values = split_names[i,:][1]
        push!(scaffold_names, values[1])
        push!(gene_names, values[2])
        push!(gene_category, values[2][1:3])
        if length(values[2]) > 3
            push!(gene_type, values[2][1:4])
        else
            push!(gene_type, values[2][1:3])
        end
    end
end

begin
    df.scaffolds = scaffold_names
    df.genes = gene_names
    df.gene_type = gene_type
    df.gene_category = gene_category
end;


# PCA

begin
    matrix = Matrix(df[:,2:8]);
    
    F₁ = svd(matrix);
    F₁.U;
    F₁.S;
    percent_explained₁ = []
    total = sum(F₁.S)

    for i ∈ 1:length(F₁.S)
        push!(percent_explained₁,F₁.S[i]/total*100)
    end
    
    𝑉 = F₁.V

end

gr()

#K-means
begin
    k_cluster = kmeans(transpose(matrix),3;maxiter=100, display=:iter);
    grp = assignments(k_cluster)
    c = counts(k_cluster)
    M = k_cluster.centers 
end

# Plot
begin
    families = []
    for i ∈ eachindex(df[:,10])
        push!(families, df[i,10][1:3])
    end

    gene_code = Vector{Int64}()
    for i ∈ eachindex(gene_category)
        if gene_category[i] == "NIP"
            push!(gene_code,1::Int64) 
        elseif gene_category[i] == "PIP" 
            push!(gene_code,2::Int64) 
        else
            push!(gene_code,3::Int)
        end
    end

    x_val = Matrix(df[:,2:8]) * 𝑉[:,1]
    y_val =  Matrix(df[:,2:8]) * 𝑉[:,2]
    pointx = transpose(M)*𝑉[:,1]
    pointy = transpose(M)*𝑉[:,2]


    k_plot₃ = scatter(x_val,
        y_val,
        #xlims = (-30,-10),
        marker_z = gene_code,
        color =:rainbow,
        size=(800,800),
        markersize=8,
        #markershape = gene_code,
        legend=:none,
        title ="PC 1 and PC 2 grouped by k-means clusters, k=3",
        xlabel="PC 1", 
        ylabel="PC 2"
        )
    
    annotate!(x_val.+.075,(y_val.+0.20),gene_names,10)

    plot!([-25,-17],[3,-.75], linecolor = :gray)
    plot!([-17,-17],[-.75,-4], linecolor = :gray)
    plot!([-17,-14.5],[-.75,4], linecolor = :gray)
    
    #png("PCA_Aquaporins1")

end

kmeansdata = [x_val y_val grp]

#CSV.write("kmeans_data", Tables.table(kmeansdata))

clusterdf = DataFrame(kmeansdata[1:end,:],:auto)
rename!(clusterdf, :x1 => :xpos, :x2 => :ypos, :x3 => :grp)

using LazySets
begin
    hull1 = Matrix(filter(row -> row.grp == 1,clusterdf)[:,1:2])
    hull2 = Matrix(filter(row -> row.grp == 2,clusterdf)[:,1:2])
    hull3 = Matrix(filter(row -> row.grp == 3,clusterdf)[:,1:2])


    h1 = Vector{Vector{Float64}}()
    h2 = Vector{Vector{Float64}}()
    h3 = Vector{Vector{Float64}}()

    for i ∈ eachindex(hull1[:,1])    
        vec = [hull1[i,1], hull1[i,2]]
        push!(h1, vec)
    end
    
    for i ∈ eachindex(hull2[:,1])    
        vec = [hull2[i,1], hull2[i,2]]
        push!(h2, vec)
    end

    for i ∈ eachindex(hull3[:,1])    
        vec = [hull3[i,1], hull3[i,2]]
        push!(h3, vec)
    end

    h1_1 = convex_hull(h1)
    h2_1 = convex_hull(h2)
    h3_1 = convex_hull(h3)
end

begin
    k_plot = scatter(x_val,
        y_val,
        marker_z = gene_code,
        color =:rainbow,
        size=(800,800),
        markersize=8,
        #markershape = gene_code,
        legend=:none,
        title ="PC 1 and PC 2 grouped by k-means clusters, k=3",
        xlabel="PC 1", 
        ylabel="PC 2"
        )

    annotate!(x_val.+.075,(y_val.+0.20),gene_names,10)
    
    p1 = plot!([Singleton(pt1) for pt1 in h1_1], alpha=0)
    plot!(p1, VPolygon(h1_1), alpha=0.2)
    
    p2 = plot!([Singleton(pt2) for pt2 in h2_1], alpha=0)
    plot!(p2, VPolygon(h2_1), alpha=0.2)
    
    p3 = plot!([Singleton(pt3) for pt3 in h3_1], alpha=0)
    plot!(p3, VPolygon(h3_1), alpha=0.2)


end
png("PCA_Aquaporins2")