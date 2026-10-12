module moireband_particle_hole_sym
    using LinearAlgebra
    using Plots
    using ProgressMeter
    using SpecialFunctions
    using Polynomials
    using HypergeometricFunctions
    using BlockArrays
    using Base.Threads
    using Primes
    using JLD2
    using Combinatorics
    using SparseArrays
    using Memoization
    using Dates
    using Arpack
    using Rotations
    include("parameters.jl")
    using .parameters  # 加载同目录下的 parameters.jl

    @memoize function effH1(k,τ)
        result = (hbar^2*((k[1]+τ*kappa1[1])^2+(k[2]+τ*kappa1[2])^2))/(2m*am^2)/e*1e+3
        return result
    end

    @memoize function effH2(k,τ)
        result = (hbar^2*((k[1]+τ*kappa2[1])^2+(k[2]+τ*kappa2[2])^2))/(2m*am^2)/e*1e+3
        return result
    end

    function b(n)
        b11 = b1
        result = Vector(RotMatrix(pi/3*(n-1))*b11)
        return result
    end

    function delta(x)
        # 使用绝对值比较，确保正负小数值都能被检测到
        return abs(x) < 1e-6 ? 1 : 0
    end

    function Δ1(Gp,G)
        V1 = V*exp(im*ψ)
        V2 = V*exp(-im*ψ)
        Delta1 = V1*sum([delta(norm(G+b(n)-Gp)) for n = 1:2:5])
        Delta2 = V2*sum([delta(norm(G+b(n)-Gp)) for n = 2:2:6])
        result = Delta1 + Delta2
        return result
    end

    function Δ2(Gp,G)
        V1 = V*exp(-im*ψ)
        V2 = V*exp(im*ψ)
        Delta1 = V1*sum([delta(norm(G+b(n)-Gp)) for n = 1:2:5])
        Delta2 = V2*sum([delta(norm(G+b(n)-Gp)) for n = 2:2:6])
        result = Delta1 + Delta2
        return result
    end

    function ΔT(Gp,G,τ)
        Delta = ω*delta(norm(Gp-G)) + ω*sum([delta(norm(G-τ*b(n)-Gp)) for n = 2:3])
        return Delta
    end

    function ΔTconj(Gp,G,τ)
        Delta = ω*delta(norm(Gp-G)) + ω*sum([delta(norm(G+τ*b(n)-Gp)) for n = 2:3])
        return Delta
    end

    ````
    创建G截断
    ```
    const range_val = 4  # 假设 range 为 2（根据实际需求修改）
    const G1 = b(1) # 等价于 (4π)/√3 {1, 0}
    const G2 = b(2)  # 等价于 (4π)/√3 {1/2, √3/2}

    # 生成所有整数组合 {m1, m2}，其中 m1, m2 ∈ [-range_val, range_val]
    combinations = []
    for m1 in -range_val:range_val
        for m2 in -range_val:range_val
            push!(combinations, (m1, m2))  # 存储为元组
        end
    end


    # 筛选满足六重对称性条件的组合：-range_val ≤ m1 + m2 ≤ range_val
    filteredCombinations = [
        (m1, m2) for (m1, m2) in combinations 
        if -range_val ≤ m1 + m2 ≤ range_val
    ]

    # 生成倒格点：Gall = m1*G1 + m2*G2（删除重复点）
    const Gshell = [m1*G1 + m2*G2 for (m1, m2) in filteredCombinations]
    const normG = [norm(m1*G1 + m2*G2) for (m1, m2) in filteredCombinations]
    sort!(normG)

    unique!(Gshell)  # 等价于 Mathematica 的 DeleteDuplicates

    # 按范数（Norm）升序排序
    sort!(Gshell, by = x -> norm(x))
    cutoff = length(Gshell)


    ````
    根据G截断构建矩阵
    ```
    function Htotal(k)
        dim = length(Gshell)
        cutoff = length(Gshell)
        bigH1 = BlockArray(zeros(ComplexF64, 2*dim, 2*dim), 
                       fill(2, dim),fill(2, dim))
        bigH2= BlockArray(zeros(ComplexF64, 2*dim, 2*dim), 
                        fill(2, dim),fill(2, dim))
        for np = 1:cutoff
            for n = 1:cutoff
                dynamic_up1 = effH1(k+Gshell[n],1)* Int(np == n)
                dynamic_up2 = effH2(k+Gshell[n],1)* Int(np == n)
                dynamic_down1 = effH1(k+Gshell[n],-1)* Int(np == n)
                dynamic_down2 = effH2(k+Gshell[n],-1)* Int(np == n)
                inter_upcj = ΔTconj(Gshell[np],Gshell[n],1)
                inter_up = ΔT(Gshell[np],Gshell[n],1)
                inter_downcj = ΔTconj(Gshell[np],Gshell[n],-1)
                inter_down = ΔT(Gshell[np],Gshell[n],-1)
                potential1 = Δ1(Gshell[np],Gshell[n])
                potential2 = Δ2(Gshell[np],Gshell[n])
                H_up = [dynamic_up1 - potential1 inter_upcj ; inter_up dynamic_up2 - potential2]
                H_down = [dynamic_down1 - potential1 inter_downcj ; inter_down dynamic_down2 - potential2]
                # Hblock = [Htop 0;0 Hbottom]
                bigH1[Block(np,n)] = H_up
                bigH2[Block(np,n)] = H_down
            end
        end
        return bigH1,bigH2
    end

    const gg =  [0.0,0.0]
    const kp = (b(2)+b(3))/3
    const kk = (b(1)+b(2))/3
    const kpp = -kk
    const mm = (kk + kp)/2
    # GKP = range(kpp,gg,length = 60)
    # KPM = range(kp,mm,length = 30)
    # MKK = range(mm,kk,length = 30)
    # KKG = range(kk,GG,length = 60)
    # kpath = vcat(GKP, KPM ,MKK,KKG)
    kppg = range(kpp,gg,length=40)
    gk = range(gg,kk,length=40)
    kkp = range(kk,kp,length = 40)
    kpkpp = range(kp,kpp,length=70)
    kpath = vcat(kppg, gk ,kkp, kpkpp)
    const kpoints = length(kpath)
    const highsymPoints = [1,40,80,120,190]

    @memoize function evaluate_eigen(valley)
        if valley == 1
            allM1 = [Htotal(kpath[i])[1] for i = 1:kpoints]
        else
            allM1 = [Htotal(kpath[i])[2] for i = 1:kpoints]
        end
        dimMatrix = size(allM1[1])[1]
        energies1 = zeros(kpoints, dimMatrix)
        eigvectors1 = zeros(ComplexF64, kpoints, dimMatrix, dimMatrix);
        @threads for i = 1:kpoints
            eig1 = eigen(allM1[i])
            # 对本征值排序并获取排序索引
            perm1 = sortperm(real(eig1.values))
            energies1[i,:]=sort(real(eig1.values))
            eigvectors1[i,:,:]=eig1.vectors[:,perm1]
        end
        return energies1,eigvectors1
    end

    function plotSpectrum()  # 默认绘制两者
        # 计算能带数据（保持原逻辑）
        energies1,_ = evaluate_eigen(1)
        energies2,_ = evaluate_eigen(2)
        
        # 假设已定义的变量（确保实际代码中有这些变量的定义）
        # highsymPoints、highsymLabels、kpoints 需要在函数内或外部定义
        n_k, n_bands = size(energies1)
        highsymLabels = ["κ'", "γ", "κ''", "κ", "κ'"]
        ylims = (-7, 27)
        
        # 设置绘图参数
        gr(dpi=300, size=(1000, 600))
        
        # 创建空图
        p = plot(
            xlabel="k-path",
            ylabel="Energy (meV)",
            title="WSe2",  # 标题显示模式
            ylims=ylims,
            xticks=(highsymPoints, highsymLabels),
            linewidth=1,
            legend=false
        )
        
        # 根据模式绘制能带
        n_bands1 = size(energies1, 2)
        for band in 1:n_bands1
            plot!(p, 1:kpoints, -energies1[:, band], color=:blue, alpha=0.8)
        end
        for band in 1:n_bands1
            plot!(p, 1:kpoints, -energies2[:, band], color=:red, alpha=0.8)
        end

        
        # 标记高对称点垂直线
        vline!(p, highsymPoints, color=:black, linestyle=:dash, linewidth=1)
        
        # 显示并保存图形（注意：num_layer需要定义，否则会报错）
        display(p)
        savefig(p, "energyspectrum_WSE2_VD_tb.pdf")
    end

    

    export evaluate_eigen,plotSpectrum,Htotal
end