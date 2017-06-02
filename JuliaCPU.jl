using StaticArrays

e = 0.6
q10 = 1 - e
q20 = 0.0
p10 = 0.0
p20 = sqrt((1 + e) / (1 - e))

h = 0.01

x1 = SVector{2,Float64}(q10, q20)
x2 = SVector{2,Float64}(p10, p20)
x3 = SVector{2,SVector{2,Float64}}(x1,x2)

@inline function oneStep(h, prev)
    h2 = h / 2
    @inbounds qsPrev = prev[1]
    @inbounds psPrev = prev[2]
    @inline function nablaQQ(qs)
        @inbounds q1 = qs[1]
        @inbounds q2 = qs[2]
        r = (q1^2 + q2^2) ^ (3/2)
        return SVector{2,Float64}(q1 / r, q2 / r)
    end
    @inline function nablaPP(ps)
        return ps
    end
    p2 = psPrev - h2 * nablaQQ(qsPrev)
    qNew = qsPrev + h * nablaPP(p2)
    pNew = p2 - h2 * nablaQQ(qNew)
    return SVector{2,SVector{2,Float64}}(qNew, pNew)
end

function manyStepsFinal(n,h,prev)
    for i in 1:n
        prev = oneStep(h,prev)
    end
    return prev
end

final = manyStepsFinal(100000000,h,x3)
print(final)
