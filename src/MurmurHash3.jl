"""
MurmurHash3 was written by Austin Appleby, and is placed in the public
domain. The author hereby disclaims copyright to this source code.

This version was translated into Julia by Scott P. Jones
It is licensed under the MIT license

Note - The x86 and x64 versions do _not_ produce the same results, as the
algorithms are optimized for their respective platforms. You can still
compile and run any of them on any platform, but your performance with the
non-native version will be less than optimal.
"""
module MurmurHash3
export mmhash128_a, mmhash128_u, mmhash128_c, mmhash32

u8(val)   = val%UInt8
u32(val)  = val%UInt32
u64(val)  = val%UInt64
u128(val) = val%UInt128

@inline rotl(x::Unsigned, r) = (x << r) | (x >>> (sizeof(typeof(x))*8 - r))

@inline xor33(k::UInt64) = xor(k, k >>> 33)

@inline rotl27(k) = rotl(k, 27)
@inline rotl31(k) = rotl(k, 31)
@inline rotl33(k) = rotl(k, 33)

#-----------------------------------------------------------------------------
# Finalization mix - force all bits of a hash block to avalanche

export dbf
const dbf = Ref(false)

@inline fmix(k::UInt64) = xor33(xor33(xor33(k) * 0xff51afd7ed558ccd) * 0xc4ceb9fe1a85ec53)

const c1 = 0x87c37b91114253d5
const c2 = 0x4cf5ad432745937f

@inline mhtail1(h1, k1) = xor(h1, rotl31(k1 * c1) * c2)
@inline mhtail2(h2, k2) = xor(h2, rotl33(k2 * c2) * c1)

@inline function mhblock(h1, h2, k1, k2)
    # dbf[] && print("mhblock($(repr(h1)), $(repr(h2)), $(repr(k1)), $(repr(k2))) => ")
    h1 = (rotl27(mhtail1(h1, k1)) + h2) * 5 + 0x52dce729
    h2 = (rotl31(mhtail2(h2, k2)) + h1) * 5 + 0x38495ab5
    # dbf[] && println(repr(h1), ", ", repr(h2))
    h1, h2
end

@inline mhblock(h1, h2, k1) = mhblock(h1, h2, u64(k1), u64(k1 >>> 64))

@inline function mhbody(nblocks, pnt, h1, h2)
    for i = 1:nblocks
        h1, h2 = mhblock(h1, h2, unsafe_load(pnt), unsafe_load(pnt + 8))
        pnt += 16
    end
    pnt, h1, h2
end

@inline function mhfin(len, h1, h2)
    h1 = xor(h1, u64(len))
    h2 = xor(h2, u64(len))

    h1 += h2
    h2 += h1

    h1 = fmix(h1)
    h2 = fmix(h2)

    h1 += h2
    h1, h1 + h2
end

#---------------------------------------------------------------------------

# mmhash128_8 implementation for strings stored little-ending packed in an Unsigned value

# simplified case for when length == 0
@inline function mhfin(h1)
    h2 = h1
    h1 += h2
    h2 += h1
    h1 = fmix(h1)
    h2 = fmix(h2)
    h1 += h2
    h1, h1 + h2
end

mask_val(val, left) = u64(val) & ((UInt64(1) << ((left & 7) << 3)) - 0x1)

# For val that is stored in up to 128 bits, we handle it without any loops
# len can only be 1-7 in this case
_mmhash128(len::Integer, val::Union{UInt32,UInt64}, seed::UInt32) =
    mhfin(len, mhtail1(u64(seed), mask_val(val, len)), u64(seed))

# For val that is stored in up to 128 bits, we handle it without any loops
# len can only be 1-15 in this case
_mmhash128(len::Integer, val::UInt128, seed::UInt32) =
    mhfin(len, mhtail1(u64(seed), len < 8 ? mask_val(val, len) : u64(val)),
          len > 8 ? mhtail2(u64(seed), mask_val(val >> 64, len)) : u64(seed))

@inline function mhbody(nblocks, val::Unsigned, h1, h2)
    for i = 1:nblocks
        h1, h2 = mhblock(h1, h2, u64(val), u64(val>>>64))
        val >>>= 128
    end
    val, h1, h2
end

# Handle values that are more than 128 bits long
function _mmhash128(len::Integer, val::Unsigned, seed::UInt32)
    val, h1, h2 = mhbody(len >>> 4, val, u64(seed), u64(seed))
    if (left = len & 15) > 0
        h1 = mhtail1(h1, left < 8 ? mask_val(val, left) : u64(val))
        left > 8 && (h2 = mhtail2(h2, mask_val(val >> 64, left)))
    end
    mhfin(len, h1, h2)
end

mmhash128_8_a(len::Integer, val::Unsigned, seed::UInt32) =
    len === 0 ? mhfin(seed) : _mmhash128(len, val, seed)

mmhash128_8_c(len::Integer, val::Unsigned, seed::UInt32) = mmhash128_8_a(len, val, seed)

#---------------------------------------------------------------------------

up8(val)  = u32(val) << 8
up16(val) = u32(val) << 16
up24(val) = u32(val) << 24
up32(val) = u64(val) << 32
up40(val) = u64(val) << 40
up48(val) = u64(val) << 48
up56(val) = u64(val) << 56

up13b(val) = u128(val) << 104
up14b(val) = u128(val) << 112
up15b(val) = u128(val) << 120

dn6(val) = u8(val >>> 6)
dn12(val) = u8(val >>> 12)
dn18(val) = u8(val >>> 18)

msk6(val) = u8(val & 0x3f)

# Support functions for UTF-8 handling
@inline get_utf8_2(ch) = 0x000080c0 | dn6(ch) | up8(msk6(ch))
@inline get_utf8_3(ch) = 0x008080e0 | dn12(ch) | up8(msk6(dn6(ch))) | up16(msk6(ch))
@inline get_utf8_4(ch) =
    0x808080f0 | dn18(ch) | up8(msk6(dn12(ch))) | up16(msk6(dn6(ch))) | up24(msk6(ch))

# Optimized in-place conversion to UTF-8 for hashing compatibly with isequal / String
@inline shift_n(v, n) = u128(v) << (n<<3)

# if cnt == 0 - 4, bytes must fit in k1
# cnt between 5 - 8, may overflow into k2
# if l == 8 - 12,  bytes must fit in k2
# cnt between 12 - 15, may overflow into k3

@inline function add_utf8(cnt, ch::Char, k1::UInt128)
    v = bswap(reinterpret(UInt32, ch))
    (cnt + ifelse(v < 0x80, 1, ifelse((v & 0xff) < 0xe0, 2, 3 + ((v & 0xff) >= 0xf0))),
     k1 | shift_n(v, cnt))
end

@inline function add_utf8_split(cnt, ch::Char, k1::UInt128)
    v = bswap(reinterpret(UInt32, ch))
    v < 0x80 && return cnt + 1, k1 | shift_n(ch, cnt), u64(0)
    if (v & 0xff) < 0xe0
        nc = cnt + 2
        v1, v2 = cnt == 15 ? (up15b(v), u64(v) >>> 8) : (shift_n(v, cnt), u64(0))
    else 
        nc = cnt + 3 + ((v & 0xff) >= 0xf0)
        v1, v2 = cnt == 13 ? (up13b(v), u64(v) >>> 24) :
            cnt == 14 ? (up14b(v), u64(v) >>> 16) :
            (up15b(v), u64(v) >>> 8)
    end
    return (nc, k1 | v1, v2)
end

@inline function add_utf8(cnt, chr, k1::UInt128)
    ch = u32(chr)
    if ch <= 0x7f
        cnt + 1, k1 | shift_n(ch, cnt)
    elseif ch <= 0x7ff
        cnt + 2, k1 | shift_n(get_utf8_2(ch), cnt)
    elseif ch <= 0xffff
        cnt + 3, k1 | shift_n(get_utf8_3(ch), cnt)
    else
        cnt + 4, k1 | shift_n(get_utf8_4(ch), cnt)
    end
end

@inline function add_utf8_split(cnt, chr, k1::UInt128)
    ch = u32(chr)
    ch <= 0x7f && return (cnt + 1, k1 | shift_n(ch, cnt), u64(0))
    # dbf[] && print("add_utf_split($cnt, $(repr(ch)), $(repr(k1))")
    if ch <= 0x7ff
        nc = cnt + 2
        v = get_utf8_2(ch)
        v1, v2 = cnt == 15 ? (up15b(v), u64(v) >>> 8) : (shift_n(v, cnt), u64(0))
    elseif ch <= 0xffff
        nc = cnt + 3
        v = get_utf8_3(ch)
        v1, v2 = cnt == 13 ? (up13b(v), u64(0)) :
                 cnt == 14 ? (up14b(v), u64(v >>> 16)) :
                 (up15b(v), u64(v) >>> 8)
    else
        # This will always go over, may be 1, 2, 3 bytes in second word
        nc = cnt + 4
        v = get_utf8_4(ch)
        # dbf[] && println(" : cnt=$cnt, v=$(repr(v))")
        v1, v2 = cnt == 13 ? (up13b(v), u64(v) >>> 24) :
                 cnt == 14 ? (up14b(v), u64(v) >>> 16) :
                 (up15b(v), u64(v) >>> 8)
    end
    # dbf[] && println(" -> ($nc, $(repr(v1)) => $(repr(k1|v1)), $(repr(v2)))")
    return (nc, k1 | v1, v2)
end

#-----------------------------------------------------------------------------

# AbstractString MurmurHash3, converts to UTF-8 on the fly
function mmhash128_8_c(str::AbstractString, seed::UInt32)
    k1 = UInt128(0)
    h1 = h2 = u64(seed)
    cnt = len = 0
    @inbounds for ch in str
        if cnt < 13
            cnt, k1 = add_utf8(cnt, ch, k1)
            if cnt == 16
                h1, h2 = mhblock(h1, h2, k1)
                k1 = 0%UInt128
                len += 16
                cnt = 0
            end
        else
            cnt, k1, k2 = add_utf8_split(cnt, ch, k1)
            # When k1 and k2 are full, then hash another block
            if cnt > 15
                h1, h2 = mhblock(h1, h2, k1)
                k1 = k2%UInt128
                len += 16
                cnt &= 15
            end
        end
    end
    # We should now have characters in k1 and k2, and total length in len
    if cnt != 0
        h1 = mhtail1(h1, u64(k1))
        cnt > 8 && (h2 = mhtail2(h2, u64(k1>>>64)))
    end
    mhfin(len + cnt, h1, h2)
end

#----------------------------------------------------------------------------

# Note: this is designed to work on the Str/String types, where I know in advance that
# the start of the strings are 8-byte aligned, and it is safe to access a full
# 8-byte chunk always at the end (simply masking off the remaining 1-7 bytes)

@inline mask_load(pnt, left) = unsafe_load(pnt) & ((UInt64(1) << ((left & 7) << 3)) - 0x1)

function mmhash128_8_a(len::Integer, pnt::Ptr, seed::UInt32)
    pnt8, h1, h2 = mhbody(len >>> 4, reinterpret(Ptr{UInt64}, pnt), u64(seed), u64(seed))
    if (left = len & 15) > 0
        h1 = mhtail1(h1, left < 8 ? mask_load(pnt8, left) : unsafe_load(pnt8))
        left > 8 && (h2 = mhtail2(h2, mask_load(pnt8 + 8, left)))
    end
    mhfin(len, h1, h2)
end

function mmhash128_8_a(seed::Integer)
    h1 = fmix(2 * u64(seed))
    h2 = fmix(3 * u64(seed))
    h1 + h2, h1 + 2 * h2
end

#----------------------------------------------------------------------------

# Combine bits from k1, k2, k3
@inline shift_mix(shft, k1::T, k2::T, k3::T) where {T<:Union{UInt32,UInt64}} =
    k1 >>> shft | (k2 << (sizeof(T)*8 - shft)),
    k2 >>> shft | (k3 << (sizeof(T)*8 - shft)),
    k3 >>> shft

#----------------------------------------------------------------------------

function mmhash128_8_u(len::Integer, unaligned_pnt::Ptr, seed::UInt32)
    # Should optimize handling of short (< 16 byte) unaligned strings
    ulp = reinterpret(UInt, unaligned_pnt)
    pnt = reinterpret(Ptr{UInt64}, ulp & ~UInt(7))
    fin = reinterpret(Ptr{UInt64}, (ulp + len + 0x7) & ~UInt(7)) - 8
    shft = (ulp & UInt(7))<<3
    h1 = h2 = u64(seed)
    k1 = unsafe_load(pnt) # Pick up first 1-7 bytes
    k2 = u64(0)
    while pnt < fin
        k1, k2, k3 = shift_mix(shft, k1, unsafe_load(pnt += 8), unsafe_load(pnt += 8))
        h1, h2 = mhblock(h1, h2, k1, k2)
        k1 = k3
    end
    # We should now have characters in k1 and k2, and total length in len
    if (len & 15) != 0
        h1 = mhtail1(h1, k1)
        (len & 15) > 8 && (h2 = mhtail2(h2, k2))
    end
    mhfin(len, h1, h2)
end

#----------------------------------------------------------------------------

# 32-bit MurmurHash3 (see MurmurHash3_x86_32)

@inline xor16(k::UInt32) = xor(k, k >>> 16)
@inline xor13(k::UInt32) = xor(k, k >>> 13)

# I don't think these help the generated code anymore (but they make the code easier to read)
@inline rotl13(k) = rotl(k, 13)
@inline rotl15(k) = rotl(k, 15)
@inline rotl16(k) = rotl(k, 16)
@inline rotl17(k) = rotl(k, 17)
@inline rotl18(k) = rotl(k, 18)
@inline rotl19(k) = rotl(k, 19)

# Constants for mmhash_32
const d1 = 0xcc9e2d51
const d2 = 0x1b873593

@inline fmix(h::UInt32) = xor16(xor13(xor16(h) * 0x85ebca6b) * 0xc2b2ae35)

@inline mhblock(h1, k1) = u32(rotl13(xor(h1, rotl15(k1 * d1) * d2))*0x00005) + 0xe6546b64

@inline function mhbody(nblocks, pnt, h1)
    for i = 1:nblocks
        h1 = mhblock(h1, unsafe_load(pnt))
        pnt += 4
    end
    pnt, h1
end

@inline mhtail32(h, v) = xor(h, rotl15(v * d1) * d2)
@inline mask32(v, res) = v & ifelse(res==1, 0x000ff, ifelse(res==2, 0x0ffff, 0xffffff))

@inline function calc32(len, pnt::Ptr, seed)
    res = len & 3
    res != 0 && (seed = mhtail32(seed, mask32(unsafe_load(pnt), res)))
    fmix(xor(seed, u32(len)))
end

@inline function calc32(len, val::UInt32, seed)
    res = len & 3
    res != 0 && (seed = mhtail32(seed, mask32(val, res)))
    fmix(xor(seed, u32(len)))
end
    

mmhash32(len, pnt::Ptr, seed::UInt32) =
    calc32(len, mhbody(len >>> 2, reinterpret(Ptr{UInt32}, pnt), seed)...)

# length must be 0-3
mmhash32(len, val::UInt32, seed::UInt32) = calc32(len, val, seed)

# length must be 0-7
mmhash32(len, val::UInt64, seed::UInt32) =
    (len > 3
     ? calc32(len, u32(val>>>32), mhblock(seed, u32(val)))
     : calc32(len, u32(val), seed))

function mmhash32(len, val::Unsigned, seed::UInt32)
    for i = 1:(len>>>2)
        seed = mhblock(seed, u32(val))
        val >>>= 32
    end
    calc32(len, u32(val), seed)
end

@inline function mhfin(len, h1, h2, h3, h4)
    h1 = xor(h1, u32(len))
    h2 = xor(h2, u32(len))
    h3 = xor(h3, u32(len))
    h4 = xor(h4, u32(len))

    h1 += h2; h1 += h3; h1 += h4; h2 += h1; h3 += h1; h4 += h1

    h1 = fmix(h1)
    h2 = fmix(h2)
    h3 = fmix(h3)
    h4 = fmix(h4)

    h1 += h2; h1 += h3; h1 += h4; h2 += h1; h3 += h1; h4 += h1

    up32(h2) | h1, up32(h4) | h3
end

#-----------------------------------------------------------------------------

# Calculate MurmurHash for 32-bit platforms

# Constants for mmhash128_4
const e1 = 0x239b961b
const e2 = 0xab0e9789
const e3 = 0x38b34ae5
const e4 = 0xa1e38b93

mhtail4_1(h, v) = xor(h, rotl15(v * e1) * e2)
mhtail4_2(h, v) = xor(h, rotl16(v * e2) * e3)
mhtail4_3(h, v) = xor(h, rotl17(v * e3) * e4)
mhtail4_4(h, v) = xor(h, rotl18(v * e4) * e1)

@inline mask_v32(val, left) = u32(val) & ((UInt32(1) << ((left & 3) << 3)) - 0x1)

@inline function mhblock(h1, h2, h3, h4, k1, k2, k3, k4)
    h1 = (rotl19(mhtail4_1(h1, k1)) + h2)*5 + 0x561ccd1b
    h2 = (rotl17(mhtail4_2(h2, k2)) + h3)*5 + 0x0bcaa747
    h3 = (rotl15(mhtail4_3(h3, k3)) + h4)*5 + 0x96cd1c35
    h4 = (rotl13(mhtail4_4(h4, k4)) + h1)*5 + 0x32ac3b17
    h1, h2, h3, h4
end

@inline function mhbody(nblocks, pnt, h1, h2, h3, h4)
    for i = 1:nblocks
        h1, h2, h3, h4 =
            mhblock(h1, h2, h3, h4,
                    unsafe_load(pnt), unsafe_load(pnt+4), unsafe_load(pnt+8), unsafe_load(pnt+12))
        pnt += 16
    end
    pnt, h1, h2, h3, h4
end

# degenerate case, hash for 0 length strings, based entirely on seed
function mmhash128_4(seed::UInt32)
    h = fmix(5*seed)*5
    up32(h) | fmix(4*seed)*4, up32(h) | h
end

# For val that is stored in up to 32 bits, we handle it without any loops
# len can only be 1-3 in this case
function mmhash128_4(len::Integer, val::UInt32, seed::UInt32)
    len == 0 && return mmhash128_4(seed)
     mhfin(len, mhtail4_1(seed, mask_v32(val, len)), seed, seed, seed)
end

# For val that is stored in up to 64 bits, we handle it without any loops
# len can only be 1-7 in this case
function mmhash128_4(len::Integer, val::UInt64, seed::UInt32)
    len == 0 && return mmhash128_4(seed)
    mhfin(len, mhtail4_1(seed, len < 4 ? mask_v32(val, len) : u32(val)),
          len > 4 ? mhtail4_2(seed, mask_v32(val>>>32, len)) : seed,
          seed, seed)
end

# For val that is stored in up to 128 bits, we handle it without any loops
# len can only be 1-15 in this case
function mmhash128_4(len::Integer, val::UInt128, seed::UInt32)
    len == 0 && return mmhash128_4(seed)
    h2 = h3 = h4 = seed
    if len > 4
        val >>>= 32
        h2  = mhtail4_2(h2, len < 8 ? mask_v32(val, len) : u32(val))
        if len > 8
            val >>>= 32
            h3  = mhtail4_3(h3, len < 12 ? mask_v32(val, len) : u32(val))
            len > 12 && (h4 = mhtail4_4(h4, mask_v32(val>>>32, len)))
        end
    end
    mhfin(len, mhtail4_1(seed, len < 4 ? mask_v32(val, len) : u32(val)), h2, h3, h4)
end

function mmhash128_4(len::Integer, val::Unsigned, seed::UInt32)
    h1 = h2 = h3 = h4 = seed
    for i = 1:(len>>>4)
        h1, h2, h3, h4 =
            mhblock(h1, h2, h3, h4, u32(val), u32(val>>>32), u32(val>>>64), u32(val>>>96))
        val >>>= 128
    end
    if (left = len & 15) != 0
        # Pick up 32-bit
        h1  = mhtail4_1(h1, left < 4 ? mask_v32(val, left) : u32(val))
        if left > 4
            val >>>= 32
            h2  = mhtail4_2(h2, left < 8 ? mask_v32(val, left) : u32(val))
            if left > 8
                val >>>= 32
                h3  = mhtail4_3(h3, left < 12 ? mask_v32(val, left) : u32(val))
                left > 12 && (h4 = mhtail4_4(h4, mask_v32(val>>>32, left)))
            end
        end
    end
    mhfin(len, h1, h2, h3, h4)
end

function mmhash128_4(len::Integer, pnt::Ptr{UInt32}, seed::UInt32)
    pnt, h1, h2, h3, h4 = mhbody(len >>> 4, pnt, seed, seed, seed, seed)
    if (left = len & 15) != 0
        h1  = mhtail4_1(h1, unsafe_load(pnt))
        if left > 4
            h2  = mhtail4_2(h2, unsafe_load(pnt+4))
            if left > 8
                h3  = mhtail4_3(h3, unsafe_load(pnt+8))
                left > 12 && (h4  = mhtail4_4(h4, unsafe_load(pnt+12)))
            end
        end
    end
    mhfin(len, h1, h2, h3, h4)
end

mmhash128_4(len::Integer, pnt::Ptr, seed::UInt32) =
    mmhash128_4(len, reinterpret(Ptr{UInt32}, pnt), seed)

# Handle value stored in an unsigned value
import Base.GC: @preserve

# AbstractString MurmurHash3, converts to UTF-8 on the fly (not optimized yet!)
function mmhash128_4(s::AbstractString, seed::UInt32)
    str = string(s)
    @preserve str mmhash128_4(sizeof(str), pointer(str), seed)
end

@inline function get_utf8(cnt, ch)
    if ch <= 0x7f
        cnt + 1, u32(ch)
    elseif ch <= 0x7ff
        cnt + 2, get_utf8_2(ch)
    elseif ch <= 0xffff
        cnt + 3, get_utf8_3(ch)
    else
        cnt + 4, get_utf8_4(ch)
    end
end

const mmhash128_c = @static sizeof(Int) == 8 ? mmhash128_8_c : mmhash128_4
const mmhash128_a = @static sizeof(Int) == 8 ? mmhash128_8_a : mmhash128_4
const mmhash128_u = @static sizeof(Int) == 8 ? mmhash128_8_c : mmhash128_4 # Todo: fix unaligned

end # module MurmurHash3
