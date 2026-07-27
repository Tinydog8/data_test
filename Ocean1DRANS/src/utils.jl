# 共享插值 / LES CSV 读取（供 LESNutClosure 与 validation 使用）

function _interp_linear(x::AbstractVector, f::AbstractVector, xq::Real)
    if xq >= maximum(x)
        # y/H 序列通常为 0 → -1（递减）；若 xq 在表面之上
        return x[1] >= x[end] ? f[1] : f[end]
    elseif xq <= minimum(x)
        return x[1] >= x[end] ? f[end] : f[1]
    end
    for i in 1:length(x)-1
        x1, x2 = x[i], x[i + 1]
        if (x1 >= xq >= x2) || (x1 <= xq <= x2)
            t = (xq - x1) / (x2 - x1)
            return f[i] + t * (f[i + 1] - f[i])
        end
    end
    return f[end]
end

function load_les_nut_csv(path::AbstractString)
    y = Float64[]
    n02 = Float64[]
    n03 = Float64[]
    open(path, "r") do io
        readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(strip(line), ',')
            push!(y, parse(Float64, parts[1]))
            push!(n02, parse(Float64, parts[2]))
            push!(n03, parse(Float64, parts[3]))
        end
    end
    return y, n02, n03
end
