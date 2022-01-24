export normalize

const mmgs = RefUnits(
    1u"mm",
    1u"g",
    1u"s",
    1u"A",
    1u"K",
    1u"c",
    1u"mol"
)

################

ref(::typeof(Unitful.𝐋),refU) = refU.ref_length
ref(::typeof(Unitful.𝐌),refU) = refU.ref_mass
ref(::typeof(Unitful.𝐓),refU) = refU.ref_time
ref(::typeof(Unitful.𝐈),refU) = refU.ref_current
ref(::typeof(Unitful.𝚯),refU) = refU.ref_temp
ref(::typeof(Unitful.𝐉),refU) = refU.ref_lum
ref(::typeof(Unitful.𝐍),refU) = refU.ref_mat

function ref(::Quantity{T,D,U}, refU::RefUnits) where {U,D,T}
    N = length(typeof(D).parameters[1])
    normalization = 1
    for i in 1:N
        typ = Unitful.Dimensions{typeof(D).parameters[1][i:i]}()
        power = typeof(D).parameters[1][i:i][1].power
        normalization *= ref(typ^(1/power),refU)^power
    end
    return normalization
end

function inrefunit(x::Quantity{T,D,U}) where {U,D,T}
    global mmgs
    inrefunit(x,mmgs)
end

function inrefunit(x::Quantity{T,D,U}, refU::RefUnits) where {U,D,T}
    return x/ref(x, refU) |> NoUnits
end


function inunit(x::Quantity{T,D,U}, ) where {U,D,T}
    return x*ref(x)
end
