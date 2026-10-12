module moireband
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
    function effH(k)
        result = -(hbar^2*(k[1]^2+k[2]^2))/(2m*am^2)/e*1e+3
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

    function Δ(Gp,G)
        V1 = V*exp(im*ψ)
        V2 = V*exp(-im*ψ)
        Delta1 = V1*sum([delta(norm(G+b(n)-Gp)) for n = 1:2:5])
        Delta2 = V2*sum([delta(norm(G+b(n)-Gp)) for n = 2:2:6])
        result = Delta1 + Delta2
        return result
    end

    
    function Htotal(k)
        dim = length(Gshell)
        cutoff = length(Gshell)
        bigH1 = Matrix{ComplexF64}(undef,dim,dim)
        for np = 1:cutoff
            for n = 1:cutoff
                bigH1[np,n] = effH(k+Gshell[n])*delta(norm(Gshell[np]-Gshell[n])) + Δ(Gshell[np],Gshell[n])
            end
        end
        return bigH1
    end

    const GG =  [0.0,0.0]
    const kp = (b(1)+b(2))/3
    const kk = (b(2)+b(3))/3
    const mm = (kk + kp)/2
    GKP = range(GG,kp,length = 50)
    KPM = range(kp,mm,length =25)
    MKK = range(mm,kk,length =25)
    KKG = range(kk,GG,length =50)
    kpath = vcat(GKP, KPM ,MKK,KKG)
    const kpoints = length(kpath)
    const highsymPoints = [1,50,75,100,150]

    @memoize function evaluate_eigen()
        allM = [Htotal(kpath[i]) for i = 1:kpoints]
        dimMatrix = size(allM[1])[1]
        energies = zeros(kpoints, dimMatrix);
        @threads for i = 1:kpoints
            eig_vals = eigvals(allM[i])  # 只返回本征值向量
            # 排序本征值（无需单独perm，直接sort!原地排序）
            sort!(eig_vals, by=real)  # 按实部排序
            energies[i, :] = real(eig_vals)  # 存储实部（若需复数可直接存eig_vals）
        end
        return energies
    end

    function plotSpectrum()  # 默认绘制两者
        # 计算能带数据（保持原逻辑）
        energies1 = evaluate_eigen()
        
        # 假设已定义的变量（确保实际代码中有这些变量的定义）
        # highsymPoints、highsymLabels、kpoints 需要在函数内或外部定义
        n_k, n_bands = size(energies1)
        highsymLabels = ["γ", "κ'", "m", "κ","γ"]
        ylims = (-80, 10)
        
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
            plot!(p, 1:kpoints, energies1[:, band], color=:blue, alpha=0.8)
        end
        
        # 标记高对称点垂直线
        vline!(p, highsymPoints, color=:black, linestyle=:dash, linewidth=1)
        
        # 显示并保存图形（注意：num_layer需要定义，否则会报错）
        display(p)
        savefig(p, "energyspectrum_WSE2.pdf")
    end

    

    export evaluate_eigen,plotSpectrum,Htotal
end